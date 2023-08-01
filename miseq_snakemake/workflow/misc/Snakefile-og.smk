from collections import defaultdict
import pysam
import os
import shutil
import pandas as pd
from scripts.utils import *

configfile: "config.yaml"
workdir: config['workdir']

samples = pd.read_csv(config["data_table"]).set_index("Sample_Name", drop=False)

miseq_dirs = config['input_data_dirs'].values()

if config['use_default_sheets']:
    sample_sheets = []
    for miseq_dir in miseq_dirs:
        sample_sheet = os.path.join(config["miseq_output"],miseq_dir,'SampleSheet.csv')
    sample_sheets.append(sample_sheet)
else:
    sample_sheets = config['illumina_sample_sheets'].values()

wildcard_constraints:
    miseq_dir="\d{6}_M\d{5}_\d{4}_0{9}-.{5}" #requried for agg step


# dictionary to correspond each of the directories with their matching sample sheet
samp_sheet_dict, samp_dict, samps = match_sample_sheets(miseq_dirs,sample_sheets,config['bad_run'])

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
    if config['igv']:
        inputs.extend(expand("{igv_loc}/summary/{project}_summary.json",project=config['project'],igv_loc=config['igv_output']))
        # inputs.extend(expand("{igv_loc}/{sample}_track.json",sample=samps,igv_loc=config['igv_output']))
        # inputs.extend(expand("{igv_loc}/}.fai",sample=samps,project=config['project'],
        #                                                         igv_loc=config['igv_output']))
    if True: #genotyping
        # inputs.extend(expand('data/aligned/{sample}.{reads}.bam.bai', sample=samps, reads=[1,2]))
        inputs.extend(
                      expand('analysis/plots/{sample}_{plot}_plot.png',
                             sample=samps,plot=['pos','delsize'])
                    )
        inputs.extend(expand("data/filtered/{sample}.bam",sample=samps))
    else:
        inputs.extend(expand("data/aligned/{sample}.bam",sample=samps))
    return inputs

include: 'rules/fastq.smk'
include: 'rules/fastqc.smk'

if not config['only_fastq']:
    refs = pull_refs(config['reference'])

    include: 'rules/align.smk'
    include: 'rules/genotype.smk'
    include: 'rules/igv.smk'
    
    if config['ampliseq']:
        include: 'rules/ampliseq.smk'

    if config['snp_calling']:
        include: 'rules/snp.smk'

        include: 'rules/genotype.smk'

rule all:
    input:        # 
        "qc/multiqc.html",
        pick_run(),
        'qc/all_samples.csv'
