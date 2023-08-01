import pandas as pd
import glob
print('\n')
counts = glob.glob('/grid/home/nbourgeois/nbourgeois_scripts/snakemake/miseq/tests/ampliseq_run_2201/results/counts/*.csv')
counts = [x for x in counts if 'ambig' not in x]
counts = [x for x in counts if '.csvcounted' not in x]
counts = [x for x in counts if '.csvunmapped' not in x]

tsv_loc='/grid/home/nbourgeois/nbourgeois_scripts/snakemake/miseq/tests/ampliseq_run_2201/results/counts.tsv'
first = True
for count_loc in counts:
    print(count_loc)
    df_t = pd.read_csv(count_loc)
    df_t = df_t.set_index(['Gene','Transcript ID'])
    sample_name = '_'.join(count_loc.split('/')[-1].split('.')[:-1])
    df_t = df_t.rename(columns={'val':sample_name})
    
    if first:
        df = df_t.copy()
        first = False
    else:
        df = df.merge(df_t,left_index=True, right_index=True,how='outer')
    print(len(df))
        
print(tsv_loc)
df.to_csv(tsv_loc,sep='\t',index=False)
