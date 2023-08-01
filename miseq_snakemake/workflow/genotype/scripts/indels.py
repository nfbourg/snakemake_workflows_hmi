import pysam
import os
import pandas as pd
from geno_utils import get_coords, get_contig

def calculate_deletions(samfile, guidestart, guidestop, chrom):

    # id_name = os.path.splitext(os.path.basename(filt_bam))[0].split('_R1')[0]
    dels = []
    min_qual = 15
    for pileupcolumn in samfile.pileup(chrom,max_depth=1_000_000,min_base_quality =min_qual):
        if pileupcolumn.pos not in range(guidestart,guidestop+1):
            continue
        pileupcolumn.set_min_base_quality(min_qual)

        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del and not pileupread.is_refskip:
                dels.append(
                    {
                        'Read':pileupread.alignment.query_name,
                        'Pos':str(pileupcolumn.pos),
    #                     'Ref':pileupread.alignment.reference_name                
                    }
                )

    return(dels)  

def get_del_size(samfile):
    
    del_sizes = []
    for read in samfile.fetch():
        sizes = [size for block,size in read.cigartuples if block==2]
        for size in sizes:
            del_size = {
                'Read':read.qname,
                'Size':size
            }        
            del_sizes.append(del_size)

    samfile.close()
    return pd.DataFrame(del_sizes)


def write_freq_df(del_df, total_reads, freq_out):
    freq_df = pd.DataFrame()
    freq_df['# Dels'] = del_df.groupby('Read').count()
    freq_df = freq_df.reset_index()
    freq_df = freq_df.groupby('# Dels').count()
    freq_df.loc[0,'Read']=total_reads-freq_df['Read'].sum()
    freq_df['Percent'] = freq_df['Read'] / total_reads * 100
    freq_df = freq_df.reset_index()
    freq_df.to_csv(freq_out)

def write_delsize_df(samfile, total_reads, delsize_out):
    del_sizes = get_del_size(samfile)
    del_size_max = del_sizes.groupby('Read').max()
    del_size_max = pd.DataFrame(del_size_max.reset_index().groupby('Size').count())
    del_size_max = del_size_max.reset_index()
    total_dels = del_size_max['Read'].sum()
    del_size_max['Percent'] = del_size_max['Read']/total_reads * 100
    del_size_max = del_size_max.loc[del_size_max['Percent']>.1]
    del_0 = {'Size':0,
             'Read':total_reads-total_dels,
             'Percent': (total_reads-total_dels)/total_reads *100}
    del_size_max = del_size_max.append(del_0,ignore_index=True)
    del_size_max.to_csv(delsize_out)
    return total_dels

def combine_metric_df(metric_df_locs):
    '''improvement: switch cases, check for bad cases'''
    if len(metric_df_locs) == 1:
        # fill in pairwise:
        exit
    else:
        df = pd.read_csv(metric_df_locs[0]).add(pd.read_csv(metric_df_locs[1]))
    return df
        

def main(coords_file,bam_in,del_out,pos_out,freq_out,delsize_out,metric_df_locs,metric_df_updated):

    guidestart, guidestop = get_coords(coords_file)
    contig = get_contig(coords_file)
    
    samfile = pysam.AlignmentFile(bam_in, "rb")
    total_reads = samfile.count()

    dels = calculate_deletions(samfile, guidestart, guidestop, contig)

    del_df = pd.DataFrame(dels)
    del_df.to_csv(del_out)

    pos_df = del_df.groupby('Pos').count()
    pos_df.to_csv(pos_out)

    write_freq_df(del_df, total_reads, freq_out)

    total_dels = write_delsize_df(samfile, total_reads, delsize_out)


    metric_df = combine_metric_df(metric_df_locs)
    metric_df['Sample'] = os.path.basename(bam_in).split('.')[0]
    metric_df['Percent Reads with Deletion'] = round(total_dels/total_reads*100,1)
    metric_df.to_csv(metric_df_updated)

if __name__ == '__main__':

    coords_file = snakemake.input.coords_file
    bam_in = snakemake.input.filtered_bam
    metric_df_locs = snakemake.input.metric_df_locs
    

    del_out = snakemake.output.del_csv
    pos_out = snakemake.output.pos_csv
    freq_out = snakemake.output.freq_csv
    freq_out = snakemake.output.freq_csv
    delsize_out = snakemake.output.delsize_csv
    metric_df_updated = snakemake.output.metric_df_updated

    main(coords_file,bam_in,del_out,pos_out,freq_out,delsize_out,metric_df_locs,metric_df_updated)