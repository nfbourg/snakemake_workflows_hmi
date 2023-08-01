localrules: snap_to_fa, snap_to_coords, bwa_index, bwa_mem, remove_ill_id

def match_ref(wildcards):
    sample = wildcards.sample
    match_ref = sample_ref.loc[sample,'Reference']
    fasta = f'ref/{match_ref}.fa'
    return multiext(fasta, ".amb", ".ann", ".bwt", ".pac", ".sa")

# rule ln_ref_loc:
#     input:
#         in_ref = 'ref/{ref}.bed',
#     output:
#         ln_bed = "{igv_loc}/{project}/{ref}.bed",
#     shell:
#         """
#         cwd=$(pwd);
#         ln -s $cwd/{input.in_bed} {output.ln_bed};
#         """

rule snap_to_fa:
    input:
        in_dna='ref/{ref}.dna'
    output:
        out_fa='ref/{ref}.fa',
    log:
        "logs/snap_to_fa/{ref}.log"
    params:
        fmt='fasta'
    conda:
        "envs/snap.yaml"
    script:
        "../scripts/snap_to_fa.py" 
        
rule snap_to_coords:
    input:
        in_dna='ref/{ref}.dna'
    output:
        out_coords='ref/{ref}.coords.txt',
    log:
        "logs/snap_to_coords/{ref}.log"
    conda:
        "envs/snap.yaml"
    params:
        feature = lambda wildcards: config['reference'][wildcards.ref]['feature_label']
    script:
        "../scripts/snap_to_coords.py" 

rule snap_to_bed:
    input:
        in_dna='ref/{ref}.dna'
    output:
        out_bed='ref/{ref}.bed',
    log:
        "logs/snap_to_bed/{ref}.log"
    params:
        fmt='bed'
    conda:
        "envs/snap.yaml"
    script:
        "../scripts/snap_to_fa.py" 

rule bwa_index:
    input:
        'ref/{ref}.fa',
    output:
        idx=multiext('ref/{ref}.fa', ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index/{ref}.bwa.log",
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.5.0/bio/bwa/index"
        
rule bwa_mem:
    input:
        reads = expand('data/raw/fastq/{{sample}}_R{read}.fastq.gz',read=[1,2]),
        idx=match_ref,
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

