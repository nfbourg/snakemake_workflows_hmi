localrules: samtools_index, pileup

NUM = 1
#splits the aligned reads by read group 1
rule r1_stview:
    input:
        bam_file = 'data/aligned/{sample}.bam'
    output:
        R1='data/aligned/{sample}.1.bam',
    params:
        extra="-bf 64"
    wrapper:
        "v1.5.0/bio/samtools/view"

#splits the aligned reads by read group 2
rule r2_stview:
    input:
        bam_file = 'data/aligned/{sample}.bam'
    output:
        R2='data/aligned/{sample}.2.bam'
    params:
        extra="-bf 128"
    wrapper:
        "v1.5.0/bio/samtools/view"

# rule dummy_rule:
#     input:
#         bam_file_R2='data/aligned/{sample}.2.bam',
#         bam_file_R1='data/aligned/{sample}.1.bam'
#     output:
#         "data/aligned/{sample}.{reads}.bam.bai",

#creates indices for the reads
rule samtools_index:
    input:
        "data/aligned/{sample}.{reads}.bam",
    output:
        "data/aligned/{sample}.{reads}.bam.bai",
    log:
        "logs/samtools_index/{sample}.{reads}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v1.5.0/bio/samtools/index"

#creates allele.txt output
rule pileup:
    input:
        bam='data/aligned/{sample}.{reads}.bam',
        bai='data/aligned/{sample}.{reads}.bam.bai',
        coords_file = config['coords_file']
    output:
        allele="results/{sample}.{reads}.allele.txt"
    params:
        script='/grid/data/software/gcerqueira_legacy/ngs/check_pileup.py'
    run:
        # loop through lines in coord file 
        for line in open(input.coords_file, 'r'): 
            chrom = line.split(" ")[0]
            coord = line.split(" ")[1] 
            coord = coord.replace("\n", "")
            out = output.allele
            shell(f"python {params.script} {input.bam} {chrom} {coord} >> {output.allele}")

#Creates a tsv of all the nucleotides at the coordinates of interest, compiled to include all samples
rule summarize:
    input:
        allele=expand('results/{sample}.{reads}.allele.txt', sample=samps, reads=[1,2])
    output:
        res = "results/nt_by_corrdinate.tsv"
    conda:
        "envs/tools.yaml"
    params:
        lowq = config['bad_run']
    script:
        "../scripts/summarize.py"

#Identifies any mismatches by comparing the bases with the reference bases at the coordinates of interest, R = expected
rule pivot:
    input:
        res = "results/nt_by_corrdinate.tsv",
        coords_file = config['coords_file'], 
        reference = ref  
    output:
        tab = "results/nt_ref_summary.tsv"
    conda:
        "envs/tools.yaml"
    params:
        lowq = config['bad_run']
    script:
        "../scripts/pivot.py"
    
#produces tsv of the proportion of each base at each coordinate, file per sample
rule freqs:
    input:
        df = "results/nt_ref_summary.tsv",
        coords_file = config['coords_file']
    output:
        output = expand("results/{sample}.counts.tsv", sample=samps)
    conda:
        "envs/tools.yaml"
    script:
        "../scripts/freqs.py"

#produces tsv of the proportion of each base at each coordinate, filtered by number of mutations, file per sample
#mutation number can be changed 
rule fcounts:
    input:
        df = "results/nt_ref_summary.tsv",
        coords_file = config['coords_file']
    output:
        output = expand("results/{sample}.countsfreq.tsv", sample=samps)
    params:
        m_num = NUM
    conda:
        "envs/tools.yaml"
    script:
        "../scripts/fcounts.py"

#produces pdf of stacked proportion bar plot of %nt per sample
rule plotnt:
    input:
        reference = ref,
        coords_file = config['coords_file'],
        files = expand("results/{sample}.counts.tsv", sample=samps)
    output:
        res = expand("results/{sample}.fig.pdf", sample=samps)
    conda:
        "envs/tools.yaml"
    params:
        lowq = config['bad_run']
    script:
        "../scripts/plotnt.py"

#produces histogram of mutation density per sample
rule plot_m:
    input:
        tab = "results/nt_ref_summary.tsv"
    output:
        out = expand("results/{sample}.m.pdf", sample=samps)
    conda:
        "envs/tools.yaml"
    script:
        "../scripts/plot_m.py"