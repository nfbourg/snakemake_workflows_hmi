localrules: bwa_mem


# rule star_pe_multi:
#     input:
#         # use a list for multiple fastq files for one sample
#         # usually technical replicates across lanes/flowcells
#         fq1=expand('data/raw/fastq/{{sample}}_R{read}.fastq.gz',read=[1,2]),
#         # paired end reads needs to be ordered so each item in the two lists match
#         fq2=["reads/{sample}_R2.1.fastq", "reads/{sample}_R2.2.fastq"],  #optional
#         # path to STAR reference genome index
#         idx="index",
#     output:
#         # see STAR manual for additional output files
#         aln="star/pe/{sample}/pe_aligned.sam",
#         log="logs/pe/{sample}/Log.out",
#         sj="star/pe/{sample}/SJ.out.tab",
#     log:
#         "logs/pe/{sample}.log",
#     params:
#         # optional parameters
#         extra="",
#     threads: 8
#     wrapper:
#         "v1.22.0/bio/star/align"
rule trimmomatic:
    input:
        "data/raw/fastq/{sample}_R1.fastq.gz"  # input and output can be uncompressed or compressed
    output:
        "data/trimmed/fastq/{sample}_R1.fastq.gz"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["TRAILING:3"],
        # optional parameters
        extra="",
        # optional compression levels from -0 to -9 and -11
        compression_level="-9"
    threads:
        8
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=2048
    wrapper:
        "v1.22.0/bio/trimmomatic/se"

rule star_se:
    input:
        fq1=expand('data/trimmed/fastq/{{sample}}_R{read}.fastq.gz',read=reads),
        # path to STAR reference genome index
        idx=config['reference']['hg38']['file_loc'],
    output:
        # see STAR manual for additional output files
        aln='data/aligned/{sample}Aligned.sortedByCoord.out.bam',
        transcriptome='data/aligned/{sample}Aligned.toTranscriptome.out.bam'
    log:
        "logs/star/{sample}.log",
    conda:
        'envs/star.yaml'
    params:
        # optional parameters
        prefix='data/aligned/{sample}',
        tmpdir=temp('data/tmp/{sample}'),
        extra="--outSAMtype BAM SortedByCoordinate --alignSoftClipAtReferenceEnds No --quantMode TranscriptomeSAM",
    threads: 32
    resources:
         mem_mb=40000
    shell:
        """STAR \
        --runThreadN {threads} \
        --genomeDir {input.idx} \
        --readFilesIn {input.fq1} --readFilesCommand  gunzip -c \
        {params.extra} \
        --outFileNamePrefix {params.prefix}
        """

rule rename_star_transcriptome:
    input:
        transcriptome='data/aligned/{sample}Aligned.toTranscriptome.out.bam',
        aln='data/aligned/{sample}Aligned.sortedByCoord.out.bam',

    output:
        transcriptome_new='data/aligned/{sample}.bam',
        aln_new='data/aligned/{sample}.sorted.bam',
    shell:
        """
        mv {input.transcriptome} {output.transcriptome_new};
        mv {input.aln} {output.aln_new};
        """

# rule star_se:
#     input:
#         fq1="reads/{sample}_R1.1.fastq",
#         # path to STAR reference genome index
#         idx="index",
#     output:
#         # see STAR manual for additional output files
#         aln="star/se/{sample}/se_aligned.bam",
#         log="logs/se/{sample}/Log.out",
#         log_final="logs/se/{sample}/Log.final.out",
#     log:
#         "logs/se/{sample}.log",
#     params:
#         # optional parameters
#         extra="--outSAMtype BAM Unsorted",
#     threads: 8
#     wrapper:
#         "v1.22.0/bio/star/align"