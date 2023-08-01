#!/usr/bin/env python
from Bio import SeqIO
import os

def get_guides(snap_file,guide_name):
    # assumes there is only 1 record in the snapgene file
    record = [x for x in SeqIO.parse(snap_file, "snapgene").records][0]
    guide = [x for x in record.features if x.qualifiers['label'][0] == guide_name][0]
    return(int(guide.location.start), int(guide.location.end))

def main(infile, outfile, feature):
    start, stop = get_guides(infile,feature)
    with open(outfile,'w') as fileo:
        fileo.write(f'Start Stop\n{start} {stop}')

if __name__ == '__main__':
    infile = snakemake.input.in_dna
    outfile = snakemake.output.out_coords
    feature = snakemake.params.feature

    main(infile, outfile, feature)
        