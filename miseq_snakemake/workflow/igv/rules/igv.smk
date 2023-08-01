def check_bed(wildcards):
    ref = wildcards.ref
    if config['reference'][ref]['bed_data']:
        return expand("{{igv_loc}}/{project}/{{ref}}.bed",project=config['project'])
    else:
        return []

rule samtools_index_fasta:
    input:
        'ref/{ref}.fa',
    output:
        'ref/{ref}.fa.fai',
    log:
        "logs/samtools/{ref}.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v1.17.2/bio/samtools/faidx"

rule subsample_bam:
    input:
        'data/aligned/{sample}.bam',
    output:
        "data/igv/{sample}.bam",
    log:
        "logs/subsample/{sample}.log",
    params:
        extra="-s 0.25",  # optional params string
    threads: 2
    wrapper:
        "v1.17.3/bio/samtools/view"

rule samtools_index_sub:
    input:
        "data/igv/{sample}.bam",
    output:
        "data/igv/{sample}.bai",
    log:
        "logs/samtools_subindex/{sample}.log",
    params:
        extra="",  # optional params string
    threads: 2  # This value - 1 will be sent to -@
    wrapper:
        "v1.5.0/bio/samtools/index"

# rule genomecov_bam:
#     input:
#         'data/aligned/{sample}.bam',
#     output:
#         "data/igv/{sample}.genomecov",
#     log:
#         "logs/genomecov_bam/{sample}.log"
#     params:
#         "-bg"  # optional parameters
#     wrapper:
#         "v1.20.0/bio/bedtools/genomecov"

rule snapbed_link_to_igv:
    input:
        in_bed = 'ref/{ref}.bed',
    output:
        ln_bed = "{igv_loc}/{project}/{ref}.bed",
    shell:
        """
        cwd=$(pwd);
        ln -s $cwd/{input.in_bed} {output.ln_bed};
        """


# rule genomecov_bed:
#     input:
#         # for genome file format please see:
#         # https://bedtools.readthedocs.io/en/latest/content/general-usage.html#genome-file-format
#         bed="bed_input/{sample}.sorted.bed",
#         ref="bed_input/genome_file"
#     output:
#         "genomecov_bed/{sample}.genomecov"
#     log:
#         "logs/genomecov_bed/{sample}.log"
#     params:
#         "-bg"  # optional parameters
#     wrapper:
#         "v1.20.0/bio/bedtools/genomecov"


# rule gc_link_to_igv:
#     input:
#         in_gc = 'data/igv/{sample}.genomecov',
#     output:
#         ln_gc = "{igv_loc}/{project}/samples/{sample}.genomecov",
#     shell:
#         """
#         cwd=$(pwd);
#         ln -s $cwd/{input.in_gc} {output.ln_gc};
#         """

rule ref_link_to_igv:
    input:
        in_fa = 'ref/{ref}.fa',
        in_fai = 'ref/{ref}.fa.fai',
    output:
        ln_fa = "{igv_loc}/{project}/{ref}.fa",
        ln_fai = "{igv_loc}/{project}/{ref}.fai",
    shell:
        """
        cwd=$(pwd);
        ln -s $cwd/{input.in_fa} {output.ln_fa};
        ln -s $cwd/{input.in_fai} {output.ln_fai};
        """

rule sample_link_to_igv:
    input:
        in_bam= "data/igv/{sample}.bam",
        in_bai = "data/igv/{sample}.bai",
    output:
        ln_bam = "{igv_loc}/{project}/samples/{sample}.bam",
        ln_bai = "{igv_loc}/{project}/samples/{sample}.bai",
    shell:
        """
        cwd=$(pwd);
        ln -s $cwd/{input.in_bam} {output.ln_bam};
        ln -s $cwd/{input.in_bai} {output.ln_bai};
        """

def check_bed(wildcards):

    ref = wildcards.ref
    bool_bed = config['reference'][ref]['bed_data']
    if bool_bed:
        bed = expand("{{igv_loc}}/{project}/{{ref}}.bed",project=config['project'])
    else:
        bed = []
    return(bed)

rule ref_track:
    input:
        ln_fa = expand("{{igv_loc}}/{project}/{{ref}}.fa",project=config['project']),
        ln_fai = expand("{{igv_loc}}/{project}/{{ref}}.fai",project=config['project']),
        ln_bed = check_bed,

    output:
        out_loc = "{igv_loc}/{project}/{ref}_track.json",
    params:
        chrom="{ref}",
        type='ref'
    script:
        "../scripts/igv.py"      

rule bam_track:
    input:
        ln_bam = "{igv_loc}/{project}/samples/{sample}.bam",
        ln_bai = "{igv_loc}/{project}/samples/{sample}.bai",
        # ln_gc = "{igv_loc}/{project}/samples/{sample}_genomecov_track.json"
    output:
        out_loc = "{igv_loc}/{project}/{sample}_track.json"
    params:
        type='bam'
    script:
        "../scripts/igv.py" 

# rule bed_track:
#     input:
#         ln_bed = "{igv_loc}/{project}/samples/{sample}.bam",
#         ln_bai = "{igv_loc}/{project}/samples/{sample}.bai"
#     output:
#         out_loc = "{igv_loc}/{project}/{sample}_genomecov_track.json"
#     params:
#         type='bed'
#     script:
#         "../scripts/igv.py" 

# rule gc_track:
#     input:
#         ln_bed = "{igv_loc}/{project}/samples/{sample}.genomecov",
#     output:
#         out_loc = "{igv_loc}/{project}/samples/{sample}_track_genomecov.json"
#     params:
#         type='bed'
#     script:
#         "../scripts/igv.py" 

# rule snap_track:
#     input:
#         ln_bed = "{igv_loc}/{project}/{ref}.bed",
#     output:
#         out_loc = "{igv_loc}/{project}/{ref}_bed_track.json"
#     params:
#         type='bed'
#     script:
#         "../scripts/igv.py" 

rule summary_json:
    input:
        refs = expand("{{igv_loc}}/{{project}}/{ref}_track.json",ref=refs),
        samps = expand("{{igv_loc}}/{{project}}/{sample}_track.json",sample=samps),
        # gcs = expand("{{igv_loc}}/{{project}}/samples/{sample}_track_genomecov.json",sample=samps),
    output:
        out_loc = "{igv_loc}/summary/{project}_summary.json",
    params:
        type='summary'
    script:
        "../scripts/igv.py" 


