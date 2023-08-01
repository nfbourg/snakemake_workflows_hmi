from collections import defaultdict
import shutil
import pandas as pd
from scripts.utils import *

configfile: "config.yaml"

config = init_default_config(config)

workdir: config['workdir']

#------------------ Get Sample info --------------------#
# samples = pd.read_csv(os.path.join(config['workdir'],config["data_table"])).set_index("Sample_Name", drop=False)
miseq_dirs = config['input_data_dirs'].values()

if config['use_default_sheets']:
    sample_sheets = []
    for miseq_dir in miseq_dirs:
        miseq_home = config["miseq_output"]
        sample_sheet = os.path.join(miseq_home,miseq_dir,'SampleSheet.csv')
    sample_sheets.append(sample_sheet)
else:
    sample_sheets = config['sample_sheets'].values()

wildcard_constraints:
    miseq_dir="\d{6}_M\d{5}_\d{4}_0{9}-.{5}" #requried for agg step

# dictionary to correspond each of the directories with their matching sample sheet
samp_sheet_dict, samp_dict, samps = match_sample_sheets(miseq_dirs,sample_sheets,config['bad_run'])

df_list = []
for sample_sheet in sample_sheets:
    df_t = read_sample_sheet(sample_sheet).set_index("Sample_Name", drop=False)
    df_list.append(df_t)

sample_ref = pd.concat(df_list)
samps = sample_ref.index.values
#------------------ Get Sample info --------------------#


rule all:
    input:
        # expand("data/raw/fastq/{sample}_R{read}.fastq.gz",sample=samps,read=[1,2]),
        # "qc/multiqc.html",
        expand("data/aligned/{sample}.bam",sample=samps),
        # 'qc/all_samples.csv'

# Run Demux
# include: 'rules/fastq.smk'
# include: 'rules/fastqc.smk'

# Run Align
include: 'rules/mini_align.smk'   
print(config['reference'])
refs = pull_refs(config['reference'],config['workdir'])
include: '../../refernce_format/rules/fmt_ref.smk'   
