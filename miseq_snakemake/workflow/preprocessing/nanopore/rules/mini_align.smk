localrules: mini_align

def match_fa_ref(wildcards):
    sample = wildcards.sample
    match_ref = sample_ref.loc[sample,'Reference']
    fasta = f'ref/{match_ref}.fa'
    return fasta
        
rule mini_align:
    input:
        sample = 'data/raw/fastq/{sample}.fastq.gz',
        ref=match_fa_ref,
    output:
        bam_file = 'data/aligned/{sample}.bam'
    log:
        "logs/mini_align/{sample}.log",
    conda:
        'envs/pomoxis.yaml'
    params:
        m="2"
    threads: 8
    shell:
        "mini_align -t {threads} -r {input.ref} -i {input.sample}  -p data/aligned/{wildcards.sample} -m {params.m} > {log}"