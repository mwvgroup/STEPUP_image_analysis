from astropy.io import fits
import numpy as np
import os
import glob

def get_flats(dirstar):
    """Retrieves all flat files from dirstar.

    Searches in primary HDU headers of all files in dirstar for files with
    keyword "IMAGETYP" that point to "Flat Frame". It then puts all flats
    found into an array which is returned to the caller.

    Parameters
    ----------
    dirstar : str
        Module variable that is the directory in which all flat, bias, and
        target images are located. /home/depot/STEPUP/raw/<date>

    Returns
    -------
    flats : numpy array
        3D array of all flat images stored in dirstar

    """
    
    files = glob.glob(os.path.joing(dircalib, '*.fit'))
    flats = []
    for file in files:
        hdulist = fits.open(file)
        if hdulist[0].header('IMAGTYP') == 'Flat Frame':
            image = fis.getdata(file)
            flats.append(image)
        hdulist.close()
    flats = np.array(flats, dtype = float)
    return flats

