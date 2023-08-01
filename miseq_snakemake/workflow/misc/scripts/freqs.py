#!/usr/bin/env python
import pandas as pd

if __name__ == "__main__":
    df_loc = snakemake.input.df
    coords = snakemake.input.coords_file
    output = snakemake.output.output

#creates a list of coordinates of interest
def coords_list(coords_file):
    file1 = open(coords_file, 'r') 
    cols = []
    for line in file1: 
        coord = line.split(" ")[1] 
        coord = coord.replace("\n", "")
        cols.append(coord)
    return cols

cols = coords_list(coords)


data = pd.read_csv(df_loc, sep="\t")
cols = [col for col in cols if col in data.columns]

for o in output:
    s = data[data["Sample"] == str(o.split(".")[0][8:])]
    #df of bases
    d = {'Base': ['A', 'C', 'G', 'T', 'R','N']}
    Base_df = pd.DataFrame(data=d)

    #loop creates df of prop of each base for single sample then appends it to the bases df
    for i in cols:
        x = s.groupby(i).count()
        if '-' in x.index:
            x = x.drop('-')
        y = x['index']
        y = y.reindex()  
        y = pd.DataFrame(y)
        y.rename(columns={int(i): "Base", "index": "Count"}, inplace = True)

        y[int(i)] = y.Count/y.Count.sum()
        y.drop(columns=['Count'], inplace = True)

        Base_df = pd.merge(Base_df, y, "left", left_on="Base", right_index=True)
    Base_df.to_csv(o, sep="\t", index=False)