import pandas as pd
import os

def get_coords(coords_file):
    df = pd.read_csv(coords_file, sep=' ')
    guidestart, guidestop= df.values[0]
    return(guidestart, guidestop)

def get_contig(coords_file):
    coords_file_base = os.path.basename(coords_file)
    contig = coords_file_base.split('.coords')[0]
    return contig