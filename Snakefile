import os
import pandas as pd

CONDAENV= "../envs"
TMPDIR = '/tmp'

if not 'sra_table' in config:
    from glob import glob
    potential_tables = glob('*_info.tab.txt')
    if len(potential_tables)==1:
        config['sra_table']= potential_tables[0]
    else:
        raise IOException("Define 'sra_table' in config file")
config['simplejob_threads'] =4
config['runtime']= {"simple_job": 1 }

include: "rules/sra.smk"


def load_and_validate_sra_table(path):

    RunTable= pd.read_csv(path, sep='\t', index_col=0)

    # validate sra table
    Expected_headers= ['LibraryLayout','LibrarySource','LibrarySelection','LibraryStrategy','BioSample']
    for header in Expected_headers:
        if not header in RunTable.columns:

            logger.error(f"Didn't found expected header {header} in sra table {path}"
                            )
            exit(1)


    assert all(RunTable.index.str[1:3]=="RR"), "Expect runs as index"

    Expected_library_values= {#"LibrarySelection": "RANDOM",
                              "LibraryStrategy": "WGS" ,
                              "LibrarySource": "METAGENOMIC"
                             }

    for key in Expected_library_values:
        if not all( RunTable[key] == Expected_library_values[key] ):
            logger.critical(f"'{key}' should be '{Expected_library_values[key]}' for all reads")

            exit(1)


    assert RunTable.LibraryLayout.isin(['PAIRED','SINGLE']).all() , "LibraryLayout should be paired or single"

    assert RunTable.BioSample.str.startswith('SAM').all(), "BioSample should start with 'SAMN'"


    Nlanes= pd.crosstab(RunTable.BioSample,RunTable.LibraryLayout)
    assert ~ Nlanes.isnull().any().any()
    assert (Nlanes.diff(axis=1).iloc[:,1]==0).all(), f"Not all samples have equal lanes {Nlanes}"

    return RunTable


def get_sra_fractions():

    LibraryLayouts= list( RunTable.LibraryLayout.unique())

    if 'PAIRED' in LibraryLayouts:
        Fractions = ['1','2']
    else:
        Fractions = []
    if 'SINGLE' in LibraryLayouts:
        Fractions += ['se']

    assert len(Fractions)>0

    return Fractions





RunTable = load_and_validate_sra_table(config['sra_table'])
BioSamples = list(RunTable.BioSample.unique())



rule download_all_reads:
    input:
        expand("SRAreads/{sample}_{fraction}.fastq.gz",
               sample=BioSamples,
               fraction=get_sra_fractions()
               )


#calculate list of fractions for biosample
# Biosample_Fraction= RunTable.groupby(['BioSample','LibraryLayout']).size().index.to_frame()
# Biosample_Fraction['Fractions']= Biosample_Fraction.LibraryLayout.map({'SINGLE':['_se'],'PAIRED':['_1','_2']})
# Biosample_Fraction= Biosample_Fraction.groupby(level=0).Fractions.sum()


def get_input_for_merging(wildcards):

    if (wildcards.fraction=='1') or (wildcards.fraction=='2'):
        LibraryLayout= 'PAIRED'
    elif wildcards.fraction=='se':
        LibraryLayout= 'SINGLE'
    else:
        Exception(f"Cannot map fraction {wildcards.fraction}")

    logger.debug(f"Sample = {wildcards.sample}")
    logger.debug(f"LibraryLayout = {LibraryLayout}")


    #get_runs_for_biosample
    runs= RunTable.query(f" BioSample =='{wildcards.sample}' & LibraryLayout=='{LibraryLayout}'"
                   ).index

    assert len(runs) >0, "Found no runs for ths sample"

    logger.debug(f"N runs = {len(runs)}")

    if LibraryLayout=='SINGLE':
        return expand("SRAreads/{sra_run}.fastq.gz",
                            sra_run= runs,
                            )
    else:
        return expand("SRAreads/{sra_run}_{fraction}.fastq.gz",
                        sra_run= runs,
                        fraction= ['1','2']
                        )


localrules: merge_runs_to_sample
rule merge_runs_to_sample:
    input:
        get_input_for_merging
    output:
        "SRAreads/{sample}_{fraction}.fastq.gz"
    shell:
        "cat {input} > {output}"
