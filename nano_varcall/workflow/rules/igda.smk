def match_ref(wildcards):
    sample = wildcards.sample
    match_ref = sample_df.loc[sample,'Reference']
    fasta = f'ref/{match_ref}.fa'
    return fasta
    
    
rule realign_bams:
    input:
        bam='data/aligned/{sample}.bam',
        ref= match_ref
    output:
        bam='data/aligned/{sample}.neg.bam'
    threads: 4
    conda:
        'envs/igda.yml'
    shell:
        'igda_align_ont {input.bam} {input.ref} {output.bam} {threads}'

rule detect_snvs:
    input:
        bam='data/aligned/{sample}.sorted.neg.bam',
        ref=match_ref
    output:
        bam='data/aligned/{sample}.neg.bam'
    threads: 4
    conda:
        'envs/igda.yml'
    shell:
        'igda_pipe_detect -m ont {input.bam} {input.ref} {output.bam} {threads}'

       reffile contextmodel outdir