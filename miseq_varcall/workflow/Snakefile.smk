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
for sample_sheet in sample_sheets:
    df_t = read_sample_sheet(sample_sheet).set_index("Sample_Name", drop=False)
    df_list.append(df_t)

sample_df = pd.concat(df_list)
samps = sample_df.index.values
refs = sample_df['Reference'].unique()

# barcodes = 
# sample_ref = pd.concat(df_list)
# samps = sample_ref.index.values
# refs = sample_ref['Reference'].unique()

####--------------- Get Sample info -----------------####
rule all:
    input:
        expand("results/vcf/{sample}",sample=samps),
        expand("{igv_loc}/summary/{project}_summary.json",project=config['project'],igv_loc=config['igv_output']),
        "qc/multiqc.html",

include: 'rules/clair3.smk'   


module miseq:
    snakefile: "../../miseq_preprocess/workflow/Snakefile.smk"
    config: config

use rule * from miseq as miseq_*

module igv:
    snakefile: "../../igv/workflow/Snakefile.smk"
    config: config

use rule * from igv as igv_*

# module igv:
#     snakefile: "../../igv/workflow/Snakefile.smk"
#     config: config

# use rule * from igv as igv_*

# Run Align


# print(config['reference'])
# refs = pull_refs(config['reference'],config['workdir']) # copy the references to the wordking directory

