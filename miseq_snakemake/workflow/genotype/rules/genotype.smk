def match_coords(wildcards):
    sample = wildcards.sample
    match_ref = samples.loc[sample,'Reference']
    coords_file = f'ref/{match_ref}.coords.txt'
    return coords_file

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

rule samtools_index_aligned:
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

rule filter_bam:
    input:
        bam_in = 'data/aligned/{sample}.{reads}.bam',
        bam_index = 'data/aligned/{sample}.{reads}.bam.bai',
        coords_file = match_coords
    output:
        filtered_bam = temp('data/filtered/temp_{sample}.{reads}.bam'),
        filtering_metric = temp('qc/filtering_{sample}.{reads}.csv')
    log:
        "logs/{sample}.{reads}"
    conda:
        "envs/genotype.yaml"
    script:
        "../scripts/genotype/filter_bam.py"

rule merge_filtered:
    input:
        expand('data/filtered/temp_{{sample}}.{reads}.bam',reads=[1,2])
    output:
        "data/filtered/{sample}.bam",
    log:
        "logs/samtools_merge/{sample}"
    threads: 8
    wrapper:
        "v1.17.2/bio/samtools/merge"

rule samtools_index_filtered:
    input:
        "data/filtered/{sample}.bam",
    output:
        "data/filtered/{sample}.bam.bai",
    log:
        "logs/samtools_index/{sample}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v1.5.0/bio/samtools/index" 

rule calculate_deletions:
    input:
        filtered_bam = "data/filtered/{sample}.bam",
        bam_index =  "data/filtered/{sample}.bam.bai",
        metric_df_locs = expand('qc/filtering_{{sample}}.{reads}.csv',reads=[1,2]),
        coords_file = match_coords
    output:
        del_csv = 'analysis/tables/{sample}_del.csv',
        pos_csv = 'analysis/tables/{sample}_pos.csv',
        freq_csv = 'analysis/tables/{sample}_freq.csv',
        delsize_csv = 'analysis/tables/{sample}_delsize.csv',
        metric_df_updated = 'qc/final_{sample}.csv'
    log:
        "logs/{sample}"
    conda:
        "envs/genotype.yaml"
    script:
        "../scripts/genotype/indels.py"

rule gen_plot:
    input:
        pos_csv = 'analysis/tables/{sample}_{plot}.csv'
    output:
        plot_loc = 'analysis/plots/{sample}_{plot}_plot.png'
    log:
        "logs/plots/{plot}/{sample}"
    conda:
        "envs/genotype.yaml"
    params:
        plot = '{plot}'
    script:
        "../scripts/genotype/gen_plot.py"

rule combine_metrics:
    input:
        paths = expand('qc/final_{sample}.csv',sample=samps)
    output:
        out = 'qc/all_samples.csv'
    run:
        import pandas as pd 

        dfs = []
        for path in input.paths:
            df_temp = pd.read_csv(path)
            dfs.append(df_temp)
        df = pd.concat(dfs)
        df = df.set_index('Sample')
        df.to_csv(output.out)

        

# rule plot_delsize:
#     input:
#         pos_csv = 'analysis/tables/{sample}_pos.csv'
#     output:
#         plot_loc = 'analysis/plots/{sample}_pos_plot.png'
#     log:
#         "logs/plots/pos/{sample}"
#     conda:
#         "envs/genotype.yaml"
#     script:
#         "../scripts/genotype/indels.py"
