from Bio import SeqIO

def get_guides(snap_file,guide_name):
    record = [x for x in SeqIO.parse(snap_file, "snapgene").records][0]
    guide = [x for x in record.features if x.qualifiers['label'][0] == guide_name][0]
    return(int(guide.location.start), int(guide.location.end))