'''INPUTS:
(dirdark) - module variable, directory in which all darks are stored:
/home/depot/STEPUP/raw/caliration/Dark/default
'''

'''OUTPUTS:
(darks) - 3D numpy array of all dark images located within dirdark
'''
from astropy.io import fits
import numpy as np
import os

def get_darks():
    files = [f for f in os.listdir(dirdark) if f.endswith('.fit')]
    hudlist = fits.open(files)
    darks = []
    for image in files:
        if hdulist[0].header['IMAGETYP'] == 'Dark Frame':
            darks.append(image)
    hdulist.close()
    darks = np.array(darks, dtype = float)
    return darks
