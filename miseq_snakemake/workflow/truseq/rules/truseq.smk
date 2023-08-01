localrules: feature_counts

rule calculate_expression:
    input:
        # input.bam or input.fq_one must be specified (and if input.fq_one, optionally input.fq_two if paired-end)
        # an aligned to transcriptome BAM
        bam="data/aligned/{sample}.bam", # list of sam or bam files
        # Index files created by rsem-prepare-reference
        reference=multiext(config['reference']['hg38']['rsem'], ".grp", ".ti", ".transcripts.fa", ".seq", ".idx.fa", ".n2g.idx.fa"),
        # reference_bowtie: Additionally needed for FASTQ input; Index files created (by bowtie-build) from the reference transcriptome
        # reference_bowtie=multiext("index/reference", ".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"),
    output:
        # genes_results must end in .genes.results; this suffix is stripped and passed to rsem as an output name prefix
        # this file contains per-gene quantification data for the sample
        genes_results="results/{sample}.genes.results",
        # isoforms_results must end in .isoforms.results and otherwise have the same prefix as genes_results
        # this file contains per-transcript quantification data for the sample
        isoforms_results="results/{sample}.isoforms.results",
    params:
        # optional, specify if sequencing is paired-end
        paired_end=False,
        # additional optional parameters to pass to rsem, for example,
        extra="--seed 42",
    log:
        "logs/rsem/calculate_expression/{sample}.log",
    threads: 8
    wrapper:
        "v1.22.0/bio/rsem/calculate-expression"

rule rsem_generate_gene_matrix:
    input:
        # one or more expression files created by rsem-calculate-expression
        expand("results/{sample}.genes.results",sample=samps),

    output:
        # a tsv containing each sample in the input as a column
        "genes.results.tsv",
    params:
        # optional additional parameters
        extra="",
    log:
        "logs/rsem/generate_data_matrix.log",
    wrapper:
        "v1.22.0/bio/rsem/generate-data-matrix"

rule rsem_generate_isoform_matrix:
    input:
        # one or more expression files created by rsem-calculate-expression
        expand("results/{sample}.isoforms.results",sample=samps),

    output:
        # a tsv containing each sample in the input as a column
        "isoforms.results.tsv",
    params:
        # optional additional parameters
        extra="",
    log:
        "logs/rsem/generate_data_matrix.log",
    wrapper:
        "v1.22.0/bio/rsem/generate-data-matrix"


# rule feature_counts:
#     input:
#         # list of sam or bam files
#         sam="data/aligned/{sample}.bam", # list of sam or bam files
#         annotation=config['reference']['hg38']['gtf']
#         # optional input
#         #chr_names="",           # implicitly sets the -A flag
#         #fasta="genome.fasta"    # implicitly sets the -G flag
#     output:
#         multiext(
#             "results/{sample}",
#             ".featureCounts",
#             ".featureCounts.summary",
#             ".featureCounts.jcounts",
#         ),
#     threads: 2
#     params:
#         strand=1,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
#         r_path="",  # implicitly sets the --Rpath flag
#         extra="-O --fracOverlap 0.2 -J",
#     log:
#         "logs/{sample}.log",
#     wrapper:
#         "0.72.0/bio/subread/featurecounts"

