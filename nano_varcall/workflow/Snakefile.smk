from collections import defaultdict
import shutil
import pandas as pd
from scripts.utils import *

configfile: "config.yaml"

config = init_default_config(config)

workdir: config['workdir']

####--------------- Get Sample info -----------------####
# samples = pd.read_csv(os.path.join(config['workdir'],config["data_table"])).set_index("Sample_Name", drop=False)
# in_dirs = list(config['input_data_dirs'].values())
# path=os.path.join(config['pacbio_output'],in_dirs[0])

sample_sheets = config['sample_sheets'].values()

df_list = []
# barcodes = []
for sample_sheet in sample_sheets:
    df_t = pd.read_csv(sample_sheet).set_index("Sample_Name", drop=False)
    df_list.append(df_t)
sample_df = pd.concat(df_list)
samps = sample_df['Sample_Name'].tolist()
sample_df = sample_df.set_index('Sample_Name')

# barcodes = 
# sample_ref = pd.concat(df_list)
# samps = sample_ref.index.values
# refs = sample_ref['Reference'].unique()

####--------------- Get Sample info -----------------####
rule all:
    input:
        expand('data/aligned/{sample}.bam',sample=samps),
        expand("{igv_loc}/summary/{project}_summary.json",project=config['project'],igv_loc=config['igv_output']),
        expand("results/vcf/{sample}",sample=samps),
        
include: 'rules/align.smk'   
include: 'rules/clair3.smk'   

module ref:
    snakefile: "../../reference_preparation/workflow/Snakefile.smk"
    config: config

use rule * from ref as ref_*

module igv:
    snakefile: "../../igv/workflow/Snakefile_nano.smk"
    config: config

use rule * from igv as igv_*

# print(config['reference'])
# refs = pull_refs(config['reference'],config['workdir']) # copy the references to the wordking directory

