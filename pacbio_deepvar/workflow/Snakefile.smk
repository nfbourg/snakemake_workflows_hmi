from collections import defaultdict
import shutil
import pandas as pd
from scripts.utils import *

configfile: "config.yaml"

config = init_default_config(config)

workdir: config['workdir']

####--------------- Get Sample info -----------------####
# samples = pd.read_csv(os.path.join(config['workdir'],config["data_table"])).set_index("Sample_Name", drop=False)
in_dirs = list(config['input_data_dirs'].values())

sample_sheets = config['sample_sheets'].values()

# wildcard_constraints:
    # requried for agg step, miseq file format of {6#'s}_M{5#'s}_{4#'s}_{9zeros}-{5char}
    # miseq_dir="\d{6}_M\d{5}_\d{4}_0{9}-.{5}" 

# dictionary to correspond each of the directories with their matching sample sheet
# samp_sheet_dict, samp_dict, samps = match_sample_sheets(miseq_dirs,sample_sheets,config['bad_run'])

df_list = []
barcodes = []
path=os.path.join(config['pacbio_output'],in_dirs[0])
for sample_sheet in sample_sheets:
    df_t = pd.read_csv(sample_sheet).set_index("Sample_Name", drop=False)
    df_list.append(df_t)
sample_df = pd.concat(df_list)
samps = sample_df['Sample_Name'].tolist()
sample_df = sample_df.set_index('Sample_Name')
barcode_dict = sample_df['Barcode'].to_dict()
# barcodes = 
# sample_ref = pd.concat(df_list)
# samps = sample_ref.index.values
# refs = sample_ref['Reference'].unique()
####--------------- Get Sample info -----------------####
rule all:
    input:
        expand('data/aligned/{sample}.bam',sample=samps),
        expand("results/vcf/{sample}",sample=samps),
        
include: 'rules/align-pbmm2.smk'   
include: 'rules/deepvar.smk'   

module ref:
    snakefile: "../../reference_preparation/workflow/Snakefile.smk"
    config: config

use rule * from ref as ref_*

# module igv:
#     snakefile: "../../igv/workflow/Snakefile.smk"
#     config: config

# use rule * from igv as igv_*

# Run Align


# print(config['reference'])
# refs = pull_refs(config['reference'],config['workdir']) # copy the references to the wordking directory

