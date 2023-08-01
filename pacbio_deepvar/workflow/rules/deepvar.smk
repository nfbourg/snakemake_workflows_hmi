def match_ref(wildcards):
    sample = wildcards.sample
    match_ref = sample_df.loc[sample,'Reference']
    fasta = f'ref/{match_ref}.fa'
    return fasta

def match_ref_idx(wildcards):
    sample = wildcards.sample
    match_ref = sample_df.loc[sample,'Reference']
    fai = f'ref/{match_ref}.fai'
    return fai
    
def get_model_dir(platform):
    path = workflow.source_path(f"resources/models/{platform}/full_alignment.index")
    print(path)
    path = '/'.join(path.split('/')[:-1])
    print(path)
    return(path)

rule clair3:
    input:
        bam="data/aligned/{sample}.bam",
        ref=match_ref,
        idx=match_ref_idx
    output:
        vcf=directory("results/vcf/{sample}")
    params:
        platform="hifi",   # {wgs, wes, pacbio, hybrid}
        contig=config['contig'],
        extra=''
    threads: 12
    conda:
        'envs/clair3.yml'
    log:
        "logs/deepvariant/{sample}/stdout.log"
    shell:
        """run_clair3.sh \
        --bam_fn={input.bam} \
        --ref_fn={input.ref} \
        --threads={threads} \
        --platform={params.platform} \
        --model_path=/grid/home/nbourgeois/snakemake_workflows/pacbio_deepvar/workflow/rules/resources/models/hifi \
        --output={output.vcf} \
        --ctg_name={params.contig}"""
# rule deepvariant:
#     input:
#         bam="data/aligned/{sample}.bam",
#         ref=match_ref
#     output:
#         vcf="results/calls/{sample}.vcf.gz"
#     params:
#         model="pacbio",   # {wgs, wes, pacbio, hybrid}
#         sample_name=lambda w: w.sample, # optional
#         extra=""
#     threads: 2
#     conda:
#         'envs/bcl2fastq2.yaml'
#     log:
#         "logs/deepvariant/{sample}/stdout.log"
#     wrapper:
#         "v1.23.0/bio/deepvariant"


# rule deepvariant_gvcf:
#     input:
#         bam="data/aligned/{sample}.bam",
#         ref=match_ref
#     output:
#         vcf="results/gvcf_calls/{sample}.vcf.gz",
#         gvcf="results/gvcf_calls/{sample}.g.vcf.gz"
#     params:
#         model="pacbio",   # {wgs, wes, pacbio, hybrid}
#         extra=""
#     conda:
#         'envs/bcl2fastq2.yaml'
#     threads: 2
#     log:
#         "logs/deepvariant/{sample}/stdout.log"
#     wrapper:
#         "v1.22.0/bio/deepvariant"