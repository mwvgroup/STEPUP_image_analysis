from astropy.io import fits
import numpy as np
import os
import glob

def get_star(dirstar):
    """Retrieves all raw target images from dirstar.

    Searches in primary HDU headers of all files in dirstar for files with
    keyword "IMAGETYP" that point to "Light Frame". It then puts all target
    images found into an array which is returned to the caller.

    Parameters
    ----------
    dirstar : str
        Module variable that is the directory in which all flat, bias, and
        target images are located. /home/depot/STEPUP/raw/<date>

    Returns
    -------
    star : numpy array
        3D array of all raw target images stored in dirstar
    """
    
    files = glob.glob(os.path.join(dirstar, '*.fit'))
    star = []
    for file in files:
        hdulist = fits.open(file)
        if hdulist[0].header['IMAGETYP'] == 'Light Frame':
            image = fits.getdata(file)
            star.append(image)
        hdulist.close()
    star = np.array(star, dtype = float)
    return star
