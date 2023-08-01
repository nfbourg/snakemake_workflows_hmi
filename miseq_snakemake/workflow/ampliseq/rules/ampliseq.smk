localrules: ampliseq_quant, combine_counts, combine_counts

rule ampliseq_quant:
    input:
        # bam_file = ill_id_removal
        bam_file = "data/aligned/{sample}.bam"

    output:
        count_output = 'results/counts/{sample}.csv',
        ambig_output = 'results/counts/{sample}_ambig.csv',
        metric_output = 'qc/metrics/{sample}_metric.csv'
    conda:
        "envs/quant.yaml"
    params:
        ref=config['reference']
    script:
        "../scripts/quant.py" #file /grid/home/nbourgeois/nbourgeois_scripts/snakemake/miseq/workflow/python ./

rule combine_counts:
    input:
        counts=expand('results/counts/{sample}.csv', sample=samps)
    output:
        tsv_loc='results/counts.tsv'
    run:
        first = True
        for count_loc in input.counts:
            df_t = pd.read_csv(count_loc)
            df_t = df_t.set_index(['Gene','Transcript ID'])
            sample_name = '_'.join(count_loc.split('/')[-1].split('.')[:-1])
            df_t = df_t.rename(columns={'val':sample_name})
            
            if first:
                df = df_t.copy()
                first = False
            else:
                df = df.merge(df_t,left_index=True, right_index=True,how='outer')

        df.to_csv(output.tsv_loc,sep='\t',index=True)

rule combine_metrics:
    input:
        metrics=expand('qc/metrics/{sample}_metric.csv', sample=samps)
    output:
        tsv_loc='qc/metrics.tsv'
    run:
        df_list = []
        for metric_loc in input.metrics:
            print(metric_loc)
            sample = ''.join(metric_loc.split('/')[-1].split('.')[:-1])
            df_t = pd.read_csv(metric_loc)
            df_t.columns = ['Metric','Count']
            df_t['Sample'] = sample
            df_list.append(df_t)
        df = pd.concat(df_list)
        df.to_csv(output.tsv_loc,index=False)