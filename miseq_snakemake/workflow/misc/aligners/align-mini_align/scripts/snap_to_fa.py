#!/usr/bin/env python

from Bio import SeqIO
import os

def convert_snapgene_fa(infile, outfile):

    parsed = SeqIO.read(infile, 'snapgene')

    basename = os.path.basename(infile)
    parsed = SeqIO.read(infile, 'snapgene')
    parsed.id = os.path.splitext(basename)[0]
    parsed.description = '' 
    SeqIO.write(parsed, outfile, 'fasta')

def convert_snapgene_bed(infile, outfile):

    parsed = SeqIO.parse(infile, 'snapgene')
  
    features = set()
    basename = ''.join(os.path.basename(outfile).split('.')[:-1])
    with open(outfile,'w') as fileo:
        for record in parsed:
            print(record)
            for feature in record.features:
                    # if 
                    start = feature.location.start.position
                    stop = feature.location.end.position

                    name = feature.qualifiers['label'][0]
                    if feature.strand < 0:
                        strand = "-"
                    else:
                        strand = "+"
                    bed_line = f"{basename}\t{start}\t{stop}\t{name}\t1000\t{strand}\t{start}\t{stop}\t65,105,225\n"
                    if bed_line in features:
                        continue
                    else:
                        features.add(bed_line)
                        fileo.write(bed_line)

def main(inputs, outputs, fmt):

    if fmt == 'fasta':
        infile = inputs.in_dna
        outfile = outputs.out_fa
        convert_snapgene_fa(infile, outfile)

    
    if fmt == 'bed':
        infile = inputs.in_dna
        outfile = outputs.out_bed
        convert_snapgene_bed(infile, outfile)


    

if __name__ == '__main__':
    inputs = snakemake.input
    outputs = snakemake.output
    fmt = snakemake.params.fmt

    main(inputs, outputs, fmt)
    
        
