import json
import os

def format_ref_track(fa_loc, fai_loc,ln_bed, chrom):
    ref = {
           'id': chrom,
           'name': chrom,
        #    'project': 'testP',
           'relativeLink': True,
           'fastaURL': '/'.join(str(fa_loc).split('/')[-3:]),
           'indexURL': '/'.join(str(fai_loc).split('/')[-3:]),
    }

    if len(ln_bed) > 0:
        ref['tracks'] = [
                {
                    'name': chrom,
                    'type': "annotation",
                    'format': "bed",
                    'url': '/'.join(str(ln_bed).split('/')[-3:]),
                    'displayMode': "EXPANDED"
                }
            ]

    # locus = '{chrom}'

    # json = {}
    # json['reference'] = ref
    # json['locus'] = locus
    # json['tracks'] = tracks
    
    return(ref)

def format_bam_track(bam_loc,bam_idx,gc_loc=None):
    id_name = os.path.basename(str(bam_loc))
    id_name = os.path.basename(str(bam_idx))
    if gc_loc is not None:
        data = {
                "name": id_name,
                'relativeLink': True,
                "url": '/'.join(str(bam_loc).split('/')[-3:]),
                "indexURL": '/'.join(str(bam_idx).split('/')[-3:]),
                "format": "bam",
                "type": "alignment",
                "annotation": '/'.join(str(gc_loc).split('/')[-3:]),
                }
    else:
        data = {
                "name": id_name,
                'relativeLink': True,
                "url": '/'.join(str(bam_loc).split('/')[-3:]),
                "indexURL": '/'.join(str(bam_idx).split('/')[-3:]),
                "format": "bam",
                "type": "alignment"
                }
    return(data)

def format_bed_track(bed_loc):
    id_name = os.path.basename(str(bed_loc))
    data = {
            "name": id_name,
            'relativeLink': True,
            "url": '/'.join(str(bed_loc).split('/')[-3:]),
            "format": "bed",
            "type": "annotation"
            }
    return(data)

def format_summary(refs, samples, gcs):
    refs = ['/'.join(str(ref).split('/')[-2:]) for ref in refs]
    samples = ['/'.join(str(samp).split('/')[-2:]) for samp in samples]
    if gcs is not None:
        gcs = ['/'.join(str(bed).split('/')[-2:]) for bed in gcs]

        data = {
            "references": refs,
            "samples": samples,
            "genome_coverages": gcs
            }
    else:
        data = {
            "references": refs,
            "samples": samples
            }
    return(data)

def main(inputs, outputs, params):
    
    type = params.type

    if type == 'ref':
        ln_path = inputs.ln_fa
        ln_bed = inputs.ln_bed
        ln_path_idx = inputs.ln_fai
        chrom = params.chrom
        jsonData = format_ref_track(ln_path, ln_path_idx, ln_bed, chrom)  

    elif type == 'bam':
        ln_path = inputs.ln_bam
        ln_path_idx = inputs.ln_bai
        try:
            ln_gc = inputs.ln_gc
        except:
            ln_gc = None
        jsonData = format_bam_track(ln_path, ln_path_idx, ln_gc)

    elif type == 'bed':
        ln_path = inputs.ln_bed
        # ln_path_idx = inputs.ln_bai
        jsonData = format_bed_track(ln_path)

    elif type == 'summary':
        refs = inputs.refs
        samples = inputs.samps
        # gcs = inputs.gcs
        jsonData = format_summary(refs, samples, gcs=None)

    refJson = open(outputs.out_loc,'w')
    refJson.write(json.dumps(jsonData))



if __name__ == '__main__':

    inputs = snakemake.input
    outputs = snakemake.output
    params = snakemake.params

    main(inputs, outputs, params)