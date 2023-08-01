#!/usr/bin/env python
from collections import defaultdict
import pysam
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

colors=[(0/255,50/255,90/255,1),
        (246/255,126/255,94/255,1),
        (106/255,118/255,132/255,1),
        (116/255,81/255,173/255,1),
        (20/255,103/255,172/255,1),
        (41/255,171/255,226/255,1),
        (51/255,51/255,51/255,1),
       (0/255,50/255,90/255,1)]

if __name__ == "__main__":
    tab = snakemake.input.tab
    output = snakemake.output.out

dat = pd.read_csv(tab, sep="\t")

i = 0

for o in output:
    fdat = dat[dat["Sample"] == o.split('.')[0][8:]]
    sns.set_context('talk')
    plt = sns.displot(data = fdat, x="Mismatches", discrete = True, stat = "density", common_norm = False, color=colors[i])
    i += 1
    if i == len(colors):
        i=0
    plt.set(title=o.split('.')[0][8:]+' Mismatches')
    plt.figure.savefig(o)