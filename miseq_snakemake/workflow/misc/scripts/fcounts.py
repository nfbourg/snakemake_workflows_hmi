#!/usr/bin/env python
from collections import defaultdict
import pysam
import pandas as pd
import numpy as np
import sys

nums = 1

if __name__ == "__main__":
    df = snakemake.input.df
    coords = snakemake.input.coords_file
    output = snakemake.output.output

def coords_list(coords_file):
    file1 = open(coords_file, 'r') 
    cols = []
    for line in file1: 
        coord = line.split(" ")[1] 
        coord = coord.replace("\n", "")
        cols.append(coord)
    return cols

cols = coords_list(coords)

data = pd.read_csv(df, sep="\t")
cols = [col for col in cols if col in data.columns]

#filtering for mismatches
data = data[data["Mismatches"] > 0]

for o in output:
    s = data[data["Sample"] == str(o.split(".")[0][8:])]
    #filtering for mismatches
    s = s[s["Mismatches"] == nums]

    d = {'Base': ['A', 'C', 'G', 'T', 'R']}
    Base_df = pd.DataFrame(data=d)

    for i in cols:
        x = s.groupby(i).count()

        y = x['index']
        y = y.reindex()
        y = pd.DataFrame(y)
        y.rename(columns={int(i): "Base", "index": "Count"}, inplace = True)

        y[int(i)] = y.Count/y.Count.sum()
        y.drop(columns=['Count'], inplace = True)
    
        Base_df = pd.merge(Base_df, y, "left", left_on="Base", right_index=True)
    Base_df.to_csv(o, sep="\t", index=False)