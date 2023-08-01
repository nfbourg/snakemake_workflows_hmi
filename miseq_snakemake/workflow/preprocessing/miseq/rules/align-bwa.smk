localrules: bwa_mem

def match_bwa_ref(wildcards):
    sample = wildcards.sample
    match_ref = sample_ref.loc[sample,'Reference']
    fasta = f'ref/{match_ref}.fa'
    return multiext(fasta, ".amb", ".ann", ".bwt", ".pac", ".sa")
        
rule match_bwa_ref:
    input:
        reads = expand('data/raw/fastq/{{sample}}_R{read}.fastq.gz',read=[1,2]),
        idx=match_bwa_ref,
    output:
        bam_file = 'data/aligned/{sample}.bam'
    log:
        "logs/bwa_mem/{sample}.log",
    params:
        extra="",
        sorting="samtools",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        # sort_extra="",  # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "v1.5.0/bio/bwa/mem"

