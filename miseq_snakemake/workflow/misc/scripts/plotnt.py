#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt

colors=[(106/200,118/200,132/200,1),
        (0/255,50/255,90/255,1),
        (246/255,126/255,94/255,1),
        (116/255,81/255,173/255,1),
        (20/255,103/255,172/255,1),
        (41/255,171/255,226/255,1),
        (51/255,51/255,51/255,1),
       (0/255,50/255,90/255,1)]

def coords_list(coords_file):
    file1 = open(coords_file, 'r') 
    cols = []

    # Using for loop 
    for line in file1: 
        coord = line.split(" ")[1] 
        coord = coord.replace("\n", "")
        cols.append(coord)
    return cols

def create_dict(reference):
    with open(reference) as f:
        line = f.readline()
        line = f.readline()
    d = {}
    for i in cols:
        d[i] = line[int(i)]

    return d

def plotnt(df,fig_title,out):

    # fig_title = 'results/' + f.split('.')[0][8:] + '.fig.pdf'
    xsize = round(df.shape[0]/5+5.5,1)

    # fig, ax = plt.subplots()
    ax = df.plot.bar(x='CoordinateBase', stacked=True, title=fig_title, color=colors, figsize=(xsize, 6))
    plt.tight_layout(pad=5.5)

    ax.set_ylabel("Percent")
    vals = ax.get_yticks()
    ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.figure.savefig(out)

def main(files,outputs,lowq):
    for f,out in zip(files,outputs):
        dat = pd.read_csv(f, sep="\t")
        df1 = dat.transpose()

        df1.columns = df1.iloc[0]
        df1 = df1[1:]
        df1.index.name = 'Coordinate'
        df1.reset_index(inplace=True)
        first_column = df1.pop('R')
        df1.insert(1, 'R', first_column)
        
        #added for bases in x ticks
        cols = coords_list(coord_file)
        ref_dict = create_dict(reference)
        ref_bases = []
        for i in dict1.values():
            ref_bases += i
        df1['RefBase'] = ref_bases
        df1['CoordinateBase'] = df1['Coordinate'] + ' (' + df1['RefBase'] + ')'

        #outputting the vertical df
        df1.to_csv(f.split('.')[0][8:] + '.nt_ref_vertical.tsv', sep="\t", index=False)

        fig_title = f.split('.')[0][8:] + ' %NT'
        if lowq:
            n_out = ''.join(out.split('.')[:-1]) + '_with-N' + out.split('.')[-1]
            plotnt(df1,fig_title,n_out)
            cols_no_n = [x for x in df1.columns if x != 'N']
            df1_no_n = df1[cols_no_n].fillna(0)
            sums = df1_no_n[cols_no_n[1:]].sum(axis=1)
            df1_no_n[cols_no_n[1:]] = (df1_no_n[cols_no_n[1:]].T/sums.T).T
            plotnt(df1_no_n,fig_title,out)
        else:
            plotnt(df1,fig_title,out)


if __name__ == "__main__":
    coord_file = smakemake.input.coords_file
    reference = snakemake.input.reference
    files = snakemake.input.files
    outputs = snakemake.output.res
    lowq = snakemake.params.lowq
    main(files,outputs,lowq)
