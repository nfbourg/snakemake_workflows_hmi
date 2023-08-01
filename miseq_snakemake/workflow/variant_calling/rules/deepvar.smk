rule deepvariant:
    input:
        bam="mapped/{sample}.bam",
        ref="genome/genome.fasta"
    output:
        vcf="results/calls/{sample}.vcf.gz"
    params:
        model="pacbio",   # {wgs, wes, pacbio, hybrid}
        sample_name=lambda w: w.sample, # optional
        extra=""
    threads: 2
    log:
        "logs/deepvariant/{sample}/stdout.log"
    wrapper:
        "v1.23.4/bio/deepvariant"


rule deepvariant_gvcf:
    input:
        bam="mapped/{sample}.bam",
        ref="ref/genome.fasta"
    output:
        vcf="results/gvcf_calls/{sample}.vcf.gz",
        gvcf="results/gvcf_calls/{sample}.g.vcf.gz"
    params:
        model="pacbio",   # {wgs, wes, pacbio, hybrid}
        extra=""
    threads: 2
    log:
        "logs/deepvariant/{sample}/stdout.log"
    wrapper:
        "v1.23.4/bio/deepvariant"