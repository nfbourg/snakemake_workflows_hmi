import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

colors=[(0/255,50/255,90/255,1),
        (246/255,126/255,94/255,1),
        (106/255,118/255,132/255,1),
        (116/255,81/255,173/255,1),
        (20/255,103/255,172/255,1),
        (41/255,171/255,226/255,1),
        (51/255,51/255,51/255,1),
       (0/255,50/255,90/255,1)]

pal =sns.color_palette(colors)
sns.set_palette(pal)


def plot_delsize(df, id_name, plot_out):
    df['Read'] = df['Read'].astype(int)
    del_perc = 100-df.set_index('Size').loc[0,'Percent']
    fig = plt.figure()
    ax = sns.barplot(data=df,x='Size',y='Percent')
    ax.set_title(label=f"{id_name} ({round(del_perc,1)}% Reads w/ Dels) ")
    fig.savefig(plot_out)


def plot_pos(df, id_name, plot_out):

    plt.figure(figsize=(10,5))
    fig = plt.figure()

    ax = sns.barplot(data=df,x='Pos',y='Read')
    ax.set_ylabel('Count')
    ax.set_title(label=f"{id_name} del positions")
    ax.set_xticklabels(df['Pos'].unique(), rotation = 90)
    fig.savefig(plot_out)

def main(csv_loc, plot_out, plot):
    id_name = os.path.basename(csv_loc).split('.csv')[0]
    df = pd.read_csv(csv_loc)
    if plot == 'pos':
        id_name = id_name.strip('_pos')
        plot_pos(df, id_name, plot_out)
    elif plot == 'delsize':
        id_name = id_name.strip('_delsize')
        plot_delsize(df, id_name, plot_out)
    else:
        for i in range(100):
            print('ERROR\n')

if __name__ == '__main__':
    pos_csv = snakemake.input.pos_csv
    plot_loc = snakemake.output.plot_loc
    plot = snakemake.params.plot

    main(pos_csv, plot_loc, plot)