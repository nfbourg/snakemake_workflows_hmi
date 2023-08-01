from collections import defaultdict
import pysam
import pandas as pd
import numpy as np

ampliseq_manifest_loc = '/data/reference_genome/illumina_ampliseq/amplicon_manifest.txt'
truseq_manifest_loc = '/data/reference_genome/illumina_truseq/amplicon_manifest.txt'

def sort_reads(df,verbose=False, mapq_threshold=20):

    # fitler reads
    df_prim = df.loc[~(df['Supp'] | df['Secondary'])]
    lowq_reads = df_prim.loc[df_prim['MapQ']<mapq_threshold]['Read'].unique()
    df_qual = df_prim.loc[~df_prim['Read'].isin(lowq_reads)]
    df_unmapped = df_qual.loc[df_qual['Ref'].isna()]
    df_mapped = df_qual.loc[~df_qual['Read'].isin(df_unmapped['Read'].unique())]

    # Collect proper pairs, these will be the counts
    df_prop = df_mapped.loc[df_mapped['Proper']]
    # collect improper pairs, these will be ambig counts
    df_improp = df_mapped.loc[~df_mapped['Proper']]
    # collect unmapped, lowq score + mapped to "None"
    df_all_unmapped = df_prim.loc[df_prim['Read'].isin(np.append(df_unmapped['Read'], lowq_reads)),]

    if verbose: # output the read dfs if verbose
        df_prop.to_csv(count_out+'counted_reads.csv')
        df_improp.to_csv(count_out+'ambig_reads.csv')
        df_all_unmapped.to_csv(count_out+'unmapped_reads.csv')

    # divide counts by 2 since the reads are paired end
    ser_counts = df_prop.groupby(['Ref']).count()['Read']/2
    ser_ambig = df_improp.groupby(['Ref'])['Read'].count()/2
    ser_metrics = pd.Series({'Counts':ser_counts.sum(),
                        'Ambiguous':ser_ambig.sum() ,
                        'Unmapped':len(df_all_unmapped['Read'].unique())})

    return(ser_counts,ser_ambig,ser_metrics)

def create_alignment_df(alignment_path):
    
    # read in samfile
    samfile = pysam.AlignmentFile(alignment_path, "rb")
    algnmnts = [ algnmnt for algnmnt in samfile.fetch(until_eof=True)]    

    df = pd.DataFrame()
    df['Read'] = ['-'.join(algnment.query_name.split('/')) for algnment in algnmnts]
    df['Secondary']= [algnment.is_secondary for algnment in algnmnts]
    df['MapQ']= [algnment.mapping_quality for algnment in algnmnts]
    df['Supp'] = [algnment.is_supplementary for algnment in algnmnts]
    df['Ref']= [algnment.reference_name for algnment in algnmnts]
    df['Proper'] = [algnment.is_proper_pair for algnment in algnmnts]
    
    return(df)

def srs_to_df(ser_vals,manifest):

    val_table = manifest[['Gene','Transcript ID']].copy()

    # algined to target ref
    if len([x for x in ser_vals.index if x in val_table.index]) > 0:
        val_table.loc[ser_vals.index,'val'] = ser_vals.values
        # removing transcripts without matching gene
        val_table = val_table.loc[~val_table['Gene'].isna()] 
        val_table = val_table.fillna(0) 
        return val_table

    #aligned to wg ref
    else:
        val_table = val_table.groupby('Transcript ID').agg(pd.Series.mode)
        val_table.loc[ser_vals.index,'val'] = ser_vals.values
        val_table = val_table.fillna(0) 
        val_table = val_table.reset_index()
        return val_table

def main(alignment_path,count_out,ambig_out,metric_out,library):

    # function if using additional transcripts
    def clean_transcipt_df(ampliseq_manifest):
        # get transcripts that appear multiple times in the manifest
        df_tran = ampliseq_manifest.reset_index().groupby('Transcript ID').count()
        multi_trans = df_tran.loc[df_tran['Target ID']>1].index.values

        for transcript in multi_trans:
            amp_t = ampliseq_manifest.loc[ampliseq_manifest['Transcript ID']==transcript]
            gene_list = amp_t['Gene']
            if any(~gene_list.isna()):
                gene = gene_list[~gene_list.isna()][0]
                ampliseq_manifest.loc[ampliseq_manifest['Transcript ID']==transcript,'Gene'] = gene
        return(ampliseq_manifest)

    print('a')
    if library=='ampliseq':
        ampliseq_manifest = pd.read_csv(ampliseq_manifest_loc,sep='\t',skiprows=9)
        manifest = ampliseq_manifest.set_index('Target ID')
        print('b')
    elif library=='truseq':
        truseq_manifest = pd.read_csv(truseq_manifest_loc,sep='\t',skiprows=9)
        manifest = truseq_manifest.set_index('Target ID')
    else:
        raise('ERROR, WRONG LIBRARY PARAMETER')
        
    # convert bam reads into a pandas df
    df = create_alignment_df(alignment_path)

    ser_counts,ser_ambig,ser_metrics = sort_reads(df)

    # Save data
    counts_table = srs_to_df(ser_counts,manifest)
    counts_table.to_csv(count_out,index=False)

    ambig_table = srs_to_df(ser_ambig,ampliseq_manifest)
    ambig_table.to_csv(ambig_out,index=False)

    ser_metrics.to_csv(metric_out)


if __name__ == "__main__":

    bam_file = snakemake.input.bam_file
    count_out = snakemake.output.count_output
    ambig_out = snakemake.output.ambig_output
    metric_out = snakemake.output.metric_output
    library = snakemake.params.library
    print('aa')
    main(bam_file,count_out,ambig_out,metric_out,library)

    # def gen_amp_reads_dict_orig(alignment_path,outpath):
    
#     def check_readnum(read):
#         if read.is_read1:
#             return('1')
#         else:
#             return('2')
        
#     def check_primary(read):
#         prim_sec = 'Primary'
#         if read.is_secondary:
#             prim_sec = 'Secondary'
#         if read.is_supplementary:
#             prim_sec = 'Supplementary'

#     """Generate df with info for qc, as well as info to qc"""
    
#     out_loc = alignment_path.split('.')[0] +'_rawmapping.tsv.gz'
 
#     samfile = pysam.AlignmentFile(alignment_path, "rb")
#     gene_dict = {}
#     bad_reads = []
#     read_dict = {}

#     for read in samfile.fetch(until_eof=True):
#         readname = '-'.join(read.query_name.split('/'))
#         if read.is_secondary or read.is_supplementary:
#             bad_reads.append(readname)
#         if read.is_proper_pair:
#             gene_dict[readname] = read.reference_name


#     for bad_read in bad_reads:
#         if bad_read in gene_dict.keys():
#             del gene_dict[bad_read] 
        
#     gene_counts = defaultdict(lambda: 0)
    
#     for read in gene_dict.keys():
#         gene = gene_dict[read]
#         gene_counts[gene] += 1
    
#     fname = outpath
#     pd.Series(data = gene_counts).to_csv(fname)

# def bam_to_algnment_dict(alignment_path):
#     samfile = pysam.AlignmentFile(alignment_path, "rb")

#     prim_dict = defaultdict(lambda : [None,None])
#     sec_dict = defaultdict(lambda : [set(),set()])
#     for algnmnt in samfile.fetch(until_eof=True):
#         read_name = '-'.join(algnmnt.query_name.split('/'))
#         if algnmnt.is_secondary:
#             if algnmnt.mapping_quality >= 20:
#                 sec_dict[read_name][int(algnmnt.is_read2)].add(algnmnt.reference_name)
                
#         else:
#             if algnmnt.mapping_quality < 20:
#                 prim_dict[read_name][int(algnmnt.is_read2)] = None
#             else:
#                 prim_dict[read_name][int(algnmnt.is_read2)] = algnmnt.reference_name
    
#     return(prim_dict,sec_dict)
   
# def bin_reads(prim_dict,sec_dict):
    
#     def check_R1R2_match(read):
        
        
#         if prim_dict[read][0] == prim_dict[read][1]:
#             if prim_dict[read][0] == None:
#                 return(None,'Unmapped')
#             else:
#                 return(prim_dict[read][0],'Mapped')

#         P1andS2 = False
#         P2andS1 = False

#         # check if read 1 has a matching secondary alignment
#         if prim_dict[read][0] in sec_dict[read][1]:
#             P1andS2=True

#         # check if read 2 has a matching secondary alignment
#         if prim_dict[read][1] in sec_dict[read][0]:
#             P2andS1=True

#         if P2andS1 ^ P1andS2:
#             if P1andS2:
#                 return (prim_dict[read][0],'Mapped')
#             else:
#                 return (prim_dict[read][1],'Mapped')

#         elif P2andS1 or P1andS2:
#             amb_reads = [prim_dict[read][0],prim_dict[read][1]]
#             return (amb_reads,'Ambiguous')

#         else:
#             return (None,'Unmapped')
    
#     counts_dict = defaultdict(lambda: 0)
#     amb_dict = defaultdict(lambda: 0)
#     unmapped = 0
#     counter = 0
    
#     for read in prim_dict.keys():
#         target, map_id = check_R1R2_match(read)
        
#         if map_id=='Mapped':
#             counts_dict[target] +=1

#         elif map_id=='Ambiguous':
#             for targ in target:
#                 amb_dict[targ] +=.5
            
#         else:
#             unmapped+=1
#         counter+=1   
#     return(counts_dict, amb_dict, unmapped)

# def quantify_ampliseq_old(alignment_path,count_out,ambig_out,metric_out):

#     def dict_to_df(val_dict):
#         ser_vals = pd.Series(val_dict)
#         val_table = ampliseq_manifest[['Gene','Transcript ID']].copy()
#         val_table.loc[ser_vals.index,'val'] = ser_vals.values
#         val_table = val_table.fillna(0) 
#         return val_table


#     ampliseq_manifest = pd.read_csv(manifest_loc,sep='\t',skiprows=9)
#     ampliseq_manifest = ampliseq_manifest.set_index('Target ID')

#     primary_dict, secondary_dict = bam_to_algnment_dict(alignment_path)
#     counts_dict, ambig_dict, unmapped = bin_reads(primary_dict, secondary_dict)

#     counts_table = dict_to_df(counts_dict)
#     counts_table.to_csv(count_out,sep='\t',index=False)

#     ambig_table = dict_to_df(ambig_dict)
#     ambig_table.to_csv(ambig_out,index=False)

#     metrics = pd.Series({'Counts':counts_table['val'].sum(),
#                          'Ambiguous':counts_table['val'].sum(),
#                          'Unmapped':unmapped})
#     metrics.to_csv(metric_out,sep='\t')