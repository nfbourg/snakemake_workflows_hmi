from collections import defaultdict
import pysam
import os
import shutil
import pandas as pd
from scripts.utils import *

configfile: "config.yaml"

miseq_dirs = config['input_data_dirs'].values()
if config['use_default_sheets']:
    sample_sheets = []
    for miseq_dir in miseq_dirs:
        sample_sheet = os.path.join(config["miseq_output"],miseq_dir,'SampleSheet.csv')
    sample_sheets.append(sample_sheet)
else:
    sample_sheets = config['sample_sheets'].values()

wildcard_constraints:
    miseq_dir="\d{6}_M\d{5}_\d{4}_0{9}-.{5}" #requried for agg step


# dictionary to correspond each of the directories with their matching sample sheet
samp_sheet_dict, samp_dict, samps = match_sample_sheets(miseq_dirs,sample_sheets,config['bad_run'])

ref = check_local_ref(config['reference'], multiext(config['reference'], ".amb", ".ann", ".bwt", ".pac", ".sa"))
# gen_coords_file(config['genes'])

def pick_run():
    inputs = []
    if config['ampliseq']:
        inputs.append('qc/metrics.tsv')
        inputs.append('results/counts.tsv')
    if config['snp_calling']:
        inputs.append(expand("results/{sample}.counts.tsv", sample=samps))
        inputs.append(expand("results/{sample}.countsfreq.tsv", sample=samps))
        inputs.append(expand("results/{sample}.fig.pdf", sample=samps))
        inputs.append(expand("results/{sample}.m.pdf", sample=samps))
        #inputs.append(expand("results/{sample}.allele.txt",sample=samps))
    if config['only_fastq']:
        print(miseq_dirs)

        return inputs
    if True:
        inputs.extend(expand('data/aligned/{sample}.{reads}.bam.bai', sample=samps, reads=[1,2]))
    else:
        inputs.extend(expand("data/filtered/{sample}.bam",sample=samps))
    return inputs



include: 'rules/fastq.smk'
include: 'rules/qc.smk'

if not config['only_fastq']:
    include: 'rules/align.smk'

    if config['ampliseq']:
        include: 'rules/ampliseq.smk'

    if config['snp_calling']:
        include: 'rules/snp.smk'

include: 'rules/genotype.smk'


rule all:
    input:
        "qc/multiqc.html",
        pick_run()