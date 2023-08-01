#!/usr/bin/env python
from collections import defaultdict
import pysam
import pandas as pd
import numpy as np
import sys


if __name__ == "__main__":
    res = snakemake.input.res
    coords_file = snakemake.input.coords_file
    reference = snakemake.input.reference
    output = snakemake.output.tab
    lowq = snakemake.params.lowq


def label_SNP_speedup (item,col):
    #not sure what the d is referring to 
    if item.upper() == d.get(col).upper(): 
        return 'R'
    return item

def create_dict(reference):
    with open(reference) as f:
        line = f.readline()
        line = f.readline()
    d = {}
    for i in cols:
        d[i] = line[int(i)-1]
    return d


def coords_list(coords_file):
    file1 = open(coords_file, 'r') 
    cols = []
    for line in file1: 
        coord = line.split(" ")[1] 
        coord = coord.replace("\n", "")
        cols.append(coord)
    return cols

Pivot = pd.read_csv(res, sep="\t")
cols = coords_list(coords_file)
d = create_dict(reference)
cols = [col for col in cols if col in Pivot.columns]


for col in Pivot[cols]:
    Pivot[col] = [label_SNP_speedup(x,col) for x in Pivot[col]]

for col in cols:
    if col == cols[0]:
        column = Pivot[col].astype(str)
    else:
        column = column + Pivot[col].astype(str)
Pivot['stringAll'] = column

for col in cols[2:-2]:
    if col == cols[0]:
        column = Pivot[col].astype(str)
    else:
        column = column + Pivot[col].astype(str)
Pivot['stringBlock'] = column

data = pd.concat([Pivot,pd.DataFrame(columns=['Count_R', 'Count_N'])])
data = data.reset_index()
data['Count_R'] = (Pivot[cols]=='R').sum(axis=1)
data['Count_N'] = (Pivot[cols]=='N').sum(axis=1)

if not lowq: #keep N's for lowq runs
    data = data[data["Count_N"] == 0]

data["Mismatches"] = len(cols) - data.Count_R
data.to_csv(output, sep="\t", index=False)