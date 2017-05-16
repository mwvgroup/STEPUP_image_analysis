from astropy.io import fits
import numpy as np
import os
import glob

def get_darks(dirdark):
    '''NAME: get_darks

    PURPOSE:
    Retrieve all dark images from dirdark.

    EXPLANATION:
    Searches in primary HDU headers of all files in dirdark for files with
    keyword "IMAGETYP" to point to "Dark Frame". It then puts all of these 
    FITS files into a 3D nump array and returns it to the caller.

    INPUTS:
    (dirdark) - module variable, directory in which all darks are stored:
    /home/depot/STEPUP/raw/caliration/Dark/default

    OUTPUTS:
    (darks) - 3D numpy array of all dark images located within dirdark

    RESTRICTIONS:
    Must accept one parameter, dirdark, which must be a string.
    '''
    
    files = glob.glob(os.path.join(dirdark, '*.fit'))
    darks = []
    for file in files:
        hdulist = fits.open(file)
        if hdulist[0].header['IMAGETYP'] == 'Dark Frame':
            image = fits.getdata(file)
            darks.append(image)
        hdulist.close()
    darks = np.array(darks, dtype = float)
    return darks
