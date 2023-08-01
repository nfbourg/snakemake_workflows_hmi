from collections import defaultdict
import os
import shutil
import pandas as pd
from scripts.utils import *

configfile: "config.yaml"

config = init_default_config(config)

workdir: config['workdir']

from snakemake.utils import min_version
min_version("6.0")

start = config['start_format']

sample_sheets = config['illumina_sample_sheets'].values()

#------------------ Get Sample info --------------------#
df_list = []
for sample_sheet in sample_sheets:
    df_t = read_sample_sheet(sample_sheet).set_index("Sample_Name", drop=False)
    df_list.append(df_t)

sample_ref = pd.concat(df_list)
samps = sample_ref.index.values

rule all:
    input:        # 
        expand("data/aligned/{sample}.bam",sample=samps)

if start == 'raw':
        
    module miseq:
        snakefile: "../miseq/Snakefile.smk"
        config: config

    use rule * from miseq as miseq_*

else:

    rule pseudo_fastq:
        input:
            bam = '{sample}.fastq.gz'
        output:
            pseudo_bam = 'data/aligned/{sample}.bam'
        shell:
            "ln -rs {input.bam} {output.pseudo_bam}"

include: 'rules/align.smk'




# samples = pd.read_csv(config["data_table"]).set_index("Sample_Name", drop=False)
# dictionary to correspond each of the directories with their matching sample sheet
# samp_sheet_dict, samp_dict, samps = match_sample_sheets(miseq_dirs,sample_sheets,config['bad_run'])


