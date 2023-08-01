from re import L
import pandas as pd
import os
import shutil

def init_default_config(config):
    default_vals = {
        'workdir' : os.getcwd(),
        'bcl2fastq_flags': [],
        'bad_run': False,
        'legacy' : False,
        'use_default_sheets': False,
        'start_format': 'raw',
        'miseq_output' : '/data/Genomics/Upcoming_MiSeq_Runs',
        'igv_output' : '/grid/compbio/igv_data',
        'pacbio_output': '/grid/genomics/data/ccs_demux/results/'
    }

    for key,val in default_vals.items(): 
        if key not in config.keys(): config[key] = val    

    return config

def read_sample_sheet(loc):
    
    with open(loc) as fileo:
        for line in fileo:
    #         print(line)
            if line.startswith('[Data]'):
                return (pd.read_csv(fileo))
    raise('Error with Sample Sheet: No [Data] Header')

def match_sample_sheets(miseq_dirs,sample_sheets,bad_run=False):
    samp_sheet_dict = {}
    samp_dict = {}
    samps = []
    for indir, samp_sheet in zip(miseq_dirs,sample_sheets):
        samp_sheet_dict[indir] = samp_sheet
        if not bad_run:
            sub_samps =  read_sample_sheet(samp_sheet)['Sample_ID'].tolist()
        else:
            sub_samps = ['Undetermined']
        for samp in sub_samps:
            samp_dict[samp] = indir 
        samp_dict[indir] = sub_samps 
        samps.extend(sub_samps) 
    return samp_sheet_dict, samp_dict, samps


def pull_refs(ref_dict,wdir):
    
    if not os.path.exists('ref'):
        os.mkdir('ref')

    for refname in ref_dict:
        ref_loc = ref_dict[refname]['file_loc']
        ext = os.path.splitext(ref_loc)[1]
        target_path = f'{wdir}/ref/{refname}{ext}'
        shutil.copy2(ref_loc,target_path)

    return list(ref_dict.keys())

######## to readd ###########
# def check_local_ref(ref_loc,base):

#     # check if bwa index exists
#     ref_files = []
#     if os.path.splitext(ref_loc) == 'fa':
#         for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
#             ref_files.append(ref_loc+ext)

#     ref = 'ref/' + ref_loc.split('/')[-1]

#     if not os.path.exists('ref'):
#         os.mkdir('ref')
#     if not os.path.exists(ref):
#         shutil.copy2(ref_loc,ref)
#     for a_file in ref_files:
#         target = 'ref/' +a_file.split('/')[-1]
#         if os.path.exists(a_file) and not os.path.exists(target): 
#             age_afile = os.path.getmtime(a_file)
#             age_ref = os.path.getmtime(ref_loc)
#             if (age_afile > age_ref):
#                 shutil.copy2(a_file,target)

def get_flags(value):
    if value is None:
        return ''
    else:
        return value


def gen_coords_file(genes):
    with open('coords_file.txt', 'w') as f:
        for i in gs:
            name = i[0]
            for x in range(1,len(i)):
                nums = i[x].split('-')
                int_list = list(range(int(nums[0]), int(nums[-1])+1))
                for j in int_list:
                    f.write(name + ' ' + str(j))
                    f.write('\n')


    