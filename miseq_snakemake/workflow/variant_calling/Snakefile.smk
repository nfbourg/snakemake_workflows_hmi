from collections import defaultdict
import shutil
import pandas as pd
from scripts.utils import *

configfile: "config.yaml"

config = init_default_config(config)

workdir: config['workdir']

#------------------ Get Sample info --------------------#
miseq_dirs = config['input_data_dirs'].values()

if config['use_default_sheets']:
    sample_sheets = []
    for miseq_dir in miseq_dirs:
        miseq_home = config["miseq_output"]
        sample_sheet = os.path.join(miseq_home,miseq_dir,'SampleSheet.csv')
    sample_sheets.append(sample_sheet)
else:
    sample_sheets = config['illumina_sample_sheets'].values()

wildcard_constraints:
    # requried for agg step, miseq file format of {6#'s}_M{5#'s}_{4#'s}_{9zeros}-{5char}
    miseq_dir="\d{6}_M\d{5}_\d{4}_0{9}-.{5}" ,
    # sample="[^.*\/.*]"

df_list = []
for sample_sheet in sample_sheets:
    df_t = read_sample_sheet(sample_sheet).set_index("Sample_Name", drop=False)
    df_list.append(df_t)

sample_ref = pd.concat(df_list)
samps = sample_ref.index.values
refs = sample_ref['Reference'].unique()
################



rule all:
    input:        # 
        # expand("results/{sample}{ext}",sample=samps,
        #          ext=[".featureCounts",
        #          ".featureCounts.summary",
        #          ".featureCounts.jcounts"]),
        expand("results/{sample}.genes.results",sample=samps),
        "genes.results.tsv",
        "isoforms.results.tsv",
        "qc/multiqc.html",

include: 'rules/truseq.smk'

module miseq:
    snakefile: "../preprocessing/miseq/Snakefile.smk"
    config: config

use rule * from miseq as miseq_*


