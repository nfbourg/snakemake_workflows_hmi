
def match_barcode(wildcards):
    """creates the input for the ampliseq_quant function while also removing the 
    ill_id from the file names"""
    match_barcode = sample_df.loc[wildcards.sample,'Barcode']

    # samp_sheet = samp_sheet_dict[miseq_dir]
    # if config['bad_run']:
    #     return 'data/aligned/Undetermined_S0.bam'

    # else:
    #     sub_samps = read_sample_sheet(samp_sheet)['Sample_ID'].tolist()
    #     ill_id = sub_samps.index(wildcards.sample) + 1
    #     return expand('data/raw/fastq/{{sample}}_S{ill_id}_L001_R{{read}}_001.fastq.gz',ill_id=ill_id)[0]

    final_path = os.path.join(path, f'demultiplex.{match_barcode}.bam')
    return final_path


def match_ref(wildcards):
    sample = wildcards.sample
    match_ref = sample_df.loc[sample,'Reference']
    fasta = f'ref/{match_ref}.fa'
    return fasta

rule ln_demux:
    input:
        demux_in = match_barcode
    output:
        bam='data/demux/{sample}.bam',
        bam_idx='data/demux/{sample}.bam.pbi'
        # bam='data/demux/demultiplex.{barcodes}.bam',
        # bam_idx='data/demux/demultiplex.{barcodes}.bam.pbi'
    params:
        demux_loc = path
    shell:
        # for path in config['pacbio_output']:
        #     full_path = os.path.join('')
        """
        ln -s {input.demux_in} {output.bam};
        ln -s {input.demux_in}.pbi {output.bam_idx};"""
        # """
        # ln -s {params.demux_loc}/demultiplex.{wildcards.barcodes}.bam {output.bam};
        # ln -s {params.demux_loc}/demultiplex.{wildcards.barcodes}.bam.pbi {output.bam_idx};"""

rule pbmm2_align:
    input:
        reference=match_ref, # can be either genome index or genome fasta
        query="data/demux/{sample}.bam", # can be either unaligned bam, fastq, or fasta
    output:
        bam="data/aligned/{sample}.bam",
        index="data/aligned/{sample}.bam.bai",
    log:
        "logs/pbmm2_align/{sample}.log",
    params:
        preset="CCS", # SUBREAD, CCS, HIFI, ISOSEQ, UNROLLED
        sample="{sample}", # sample name for @RG header
        extra="--sort", # optional additional args
        loglevel="INFO",
    threads: 12
    wrapper:
        "v1.23.4/bio/pbmm2/align"