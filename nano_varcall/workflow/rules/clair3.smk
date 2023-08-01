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
        platform="ont",   # {wgs, wes, pacbio, hybrid}
        mnt_path=os.getcwd(),
        contig=config['contig'],
        extra=''
    threads: 12
    log:
        "logs/clair3/{sample}/stdout.log"
    shell:
        """singularity exec --bind {params.mnt_path}:/mnt /software/compbio/singularity/clair3_latest.sif \
        /opt/bin/run_clair3.sh \
        --bam_fn=/mnt/{input.bam} \
        --ref_fn=/mnt/{input.ref} \
        --threads={threads} \
        --platform={params.platform} \
        --model_path=/opt/models/r941_prom_sup_g5014 \
        --output=/mnt/{output.vcf} \
        --ctg_name={params.contig}"""