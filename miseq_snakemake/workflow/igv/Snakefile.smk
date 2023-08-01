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
refs = sample_ref['Reference'].unique()
print(refs)
# refs = pull_refs(config['reference'], config['workdir'])
#------------------ Get Sample info --------------------#

rule all:
    input:        # 
        expand("data/igv/{sample}.bam",sample=samps),
        expand("{igv_loc}/summary/{project}_summary.json",project=config['project'],igv_loc=config['igv_output']),
        "qc/multiqc.html",




include: 'rules/igv.smk'

if start == 'raw':

    module miseq:
        snakefile: "../preprocessing/miseq/Snakefile.smk"
        config: config

    use rule * from miseq as miseq_*

else:
    rule pseudo_bam:
        input:
            bam = '{sample}.bam'
        output:
            pseudo_bam = 'data/aligned/{sample}.bam'
        shell:
            "ln -rs {input.bam} {output.pseudo_bam}"
