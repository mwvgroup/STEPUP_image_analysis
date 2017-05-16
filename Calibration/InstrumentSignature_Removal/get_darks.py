from astropy.io import fits
import numpy as np
import os
import glob

def get_darks(dirdark):
    """Retrieves all dark images from dirdark.

    Searches in primary HDU headers of all files in dirdark for files with
    keyword "IMAGETYP" that point to "Dark Frame". It then puts all darks
    found into an array which is returned to the caller.

    Parameters
    ----------
    dirdark : str
        Module variable that is the directory in which all darks are
        stored. /home/depot/STEPUP/raw/caliration/Dark/default

    Returns
    -------
    darks : numpy array
        3D numpy array of all dark images located within dirdark
    """
    
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
