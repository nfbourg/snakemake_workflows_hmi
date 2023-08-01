def match_ref_fa(wildcards):
    sample = wildcards.sample
    match_ref = sample_df.loc[sample,'Reference']
    fasta = f'ref/{match_ref}.fa'
    return fasta

def match_ref_idx(wildcards):
    sample = wildcards.sample
    match_ref = sample_df.loc[sample,'Reference']
    fai = f'ref/{match_ref}.fai'
    return fai

def match_ref_contig(wildcards):
    sample = wildcards.sample
    match_ref = sample_df.loc[sample,'Reference']
    contig = config['reference'][match_ref]['contig']
    return contig
    
def get_model_dir(platform):
    path = workflow.source_path(f"resources/models/{platform}/full_alignment.index")
    print(path)
    path = '/'.join(path.split('/')[:-1])
    print(path)
    return(path)

def clair3_path(ext):
    path = os.getcwd()
    full_path = os.path.join(path,'data','aligned','{sample}' + ext)
    return full_path
    
def clair3_bai():
    path = os.getcwd()
    full_path = os.path.join(path,'data','aligned','{sample}')
    return full_path

rule clair3:
    input:
        bam="data/aligned/{sample}.bam",
        bamidx="data/aligned/{sample}.bai",
        # bam=clair3_path('.bam'),
        # bamidx=clair3_path('.bai'),
        ref=match_ref_fa,
        idx=match_ref_idx
    output:
        vcf=directory("results/vcf/{sample}")
    params:
        mnt_path=os.getcwd(),
        platform="ilmn",  
        contig=match_ref_contig,
        extra=''
    threads: 12
    # conda:
    #     'envs/clair3.yml'
    log:
        "logs/clair3/{sample}/stdout.log"
    shell:
        """singularity exec --bind {params.mnt_path}:/mnt /software/compbio/singularity/clair3_latest.sif \
        /opt/bin/run_clair3.sh \
        --bam_fn=/mnt/{input.bam} \
        --ref_fn=/mnt/{input.ref} \
        --threads={threads} \
        --platform={params.platform} \
        --model_path=/opt/models/{params.platform} \
        --output=/mnt/{output.vcf} \
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