localrules: snap_to_fa, snap_to_coords, bwa_index

# rule source_ref:
#     input:
#         in_dna=config['reference']
#     output:
#         out_fa='source/source.txt',
#     log:
#         "logs/snap_to_fa/{ref}.log"
#     params:
#         fmt='fasta'
#     conda:
#         "envs/snap.yaml"
#     run:
#         pull_refs(config['reference'],config['workdir']) 
        
# rule ln_source:

rule snap_to_fa:
    input:
        in_dna='ref/{ref}.dna'
    output:
        out_fa='ref/{ref}.fa',
    log:
        "logs/snap_to_fa/{ref}.log"
    params:
        fmt='fasta',
        infmt='snapgene'
    conda:
        "envs/snap.yaml"
    script:
        "../scripts/snap_to_fa.py" 
        
rule gb_to_fa:
    input:
        in_dna='ref/{ref}.gb'
    output:
        out_fa='ref/{ref}.fa',
    log:
        "logs/snap_to_fa/{ref}.log"
    params:
        fmt='fasta',
        infmt='gb'
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
        fmt='bed',
        infmt='snapgene'
    conda:
        "envs/snap.yaml"
    script:
        "../scripts/snap_to_fa.py" 

rule gb_to_bed:
    input:
        in_dna='ref/{ref}.gb'
    output:
        out_bed='ref/{ref}.bed',
    log:
        "logs/snap_to_bed/{ref}.log"
    params:
        fmt='bed',
        infmt='gb'
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


rule fa_idx:
    input:
        'ref/{ref}.fa',
    output:
        'ref/{ref}.fai',
    log:
        "logs/faidx/{ref}.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v1.23.4/bio/samtools/faidx"