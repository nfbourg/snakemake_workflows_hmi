from collections import defaultdict
import os
import shutil
import pandas as pd
from scripts.utils import *

configfile: "config.yaml"
workdir: config['workdir']

from snakemake.utils import min_version
min_version("6.0")

start = config['start_format']

#------------------ Get Sample info --------------------#
sample_ref = pd.read_csv(os.path.join(config['workdir'],config["data_table"])).set_index("Sample_Name", drop=False)
samps = sample_ref.index.values
refs = pull_refs(config['reference'])

rule all:
    input:        # 
        'results/counts.tsv',
        'qc/metrics.tsv'

include: 'rules/ampliseq.smk'

module igv:
    snakefile: "../igv/Snakefile.smk"
    config: config

use rule * from igv as igv_*


