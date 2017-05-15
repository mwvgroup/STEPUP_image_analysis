import numpy as np

def create_mflat():
    flats = get_flats(dircalib)
    mbias = create_mbias()
    mdark = create_mdark()
    mflat = ((np.median(flats, 2) - mdark - mbias)/np.mean(flat, 2))
    return mflat
