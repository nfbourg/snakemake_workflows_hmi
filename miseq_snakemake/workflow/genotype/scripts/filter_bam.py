from distutils import extension
import numpy as np
import pysam
import os
import pandas as pd
from geno_utils import get_coords, get_contig


def coverage_array(guidestart,guidestop,ref_pos):
    pos =  [x-guidestart for x in ref_pos if x in range(guidestart,guidestop+1)]
    coverage = np.zeros(guidestop+1-guidestart)
    coverage[pos] = coverage[pos] + 1
    return (coverage)

def filter_bam(samfile, start, stop, chrom, buffer=1, q_filter = 5):
    '''Filter the bam file to look for reads meeting the fiollowing criteria:
    - start/end of aligned ref pos are before/after the guide region
    - none of the bases have a quality score below threshold   
    '''
    # init vars
    oor_count = 0
    filt_count = 0
    data = []
    
    start = start-buffer
    
    coverage_lowq = np.zeros(stop+1-start+buffer*2)

    filtered_reads = []
    for read in samfile.fetch(chrom,start=start):
        # only look at the reference region
        if read.reference_start > start-buffer or \
                read.reference_end < stop+buffer+1:
            oor_count +=1
            continue

        ref_pos = read.get_reference_positions()

        # mark poor quality 
        lowq_pos = [x for x,y in zip(ref_pos,read.query_alignment_qualities)
                    if x in range(start,stop+1) and y < q_filter]

        if len(lowq_pos) == 0:
            filt_count+=1
            # samfile_filt.write(read)
            filtered_reads.append(read)

        else:
            coverage_lowq = coverage_lowq + coverage_array(start-buffer,stop+buffer,lowq_pos)
        
    total_reads = samfile.count()
    
    
    metrics = {
        'Total Reads':[total_reads],
        'In Range Reads':[total_reads-oor_count],
        'Filtered Reads':[filt_count],
      #  'lowq_coverage':coverage_lowq
    }
    


    return filtered_reads, metrics

    # print(f'{oor_count} reads were Out Of Range')
    # print(f'{total_reads-filt_count} reads contained poor quality bases')

def main(coords_file,bam_in,bam_out,metric_out):

 
    # guidestart, guidestop = get_guides(snap_file,guide_name)
    guidestart, guidestop = get_coords(coords_file)
    contig = get_contig(coords_file)

    samfile = pysam.AlignmentFile(bam_in, "rb")
    # write an "inrange" file that does not filter poor quality bases

    # write a file that does  filter poor quality bases and is in range
    samfile_filt = pysam.AlignmentFile(bam_out, "wb", template=samfile)

    try:
        filtered_reads, metrics = filter_bam(samfile, guidestart, guidestop, contig)
        assert(len(filtered_reads) > 1000)
    except:
        filtered_reads, metrics = filter_bam(samfile, guidestart, guidestop, contig,q_filter=0)
        print('LOW QUALITY')

    for read in filtered_reads:
        samfile_filt.write(read)
    
    samfile_filt.close()
    samfile.close()
    
    pd.DataFrame(metrics).to_csv(metric_out)




if __name__ == "__main__":

    # snap_file = snakemake.input.snap_file
    coords_file = snakemake.input.coords_file
    bam_in = snakemake.input.bam_in

    bam_out = snakemake.output.filtered_bam
    metric_out = snakemake.output.filtering_metric

    main(coords_file,bam_in,bam_out,metric_out)


