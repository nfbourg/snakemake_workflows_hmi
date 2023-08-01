

def match_ref(wildcards):
    sample = wildcards.sample
    match_ref = sample_df.loc[sample,'Reference']
    fasta = f'ref/{match_ref}.fa'
    return fasta

rule minimap2_bam_sorted:
    input:
        target=match_ref, # can be either genome index or genome fasta
        query="data/fastq/{sample}_L800_Q80.fastq.gz",
    output:
        "data/aligned/{sample}.bam",
    log:
        "logs/minimap2/{sample}.log",
    params:
        extra="-x map-ont -A 1 -B 2 -m 2 --eqx",  # optional
        sorting="coordinate",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra="",  # optional: extra arguments for samtools/picard
    threads: 8
    wrapper:
        "v1.25.0/bio/minimap2/aligner"