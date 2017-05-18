from astropy.io import fits
import glob
import numpy as np
import os

def get_science_images(dirtarget):
    """Retrieves all raw target images from dirtarget.

    Searches in primary HDU headers of all files in dirtarget for files with
    keyword "IMAGETYP" that points to "Light Frame". It then puts all target
    images found into an array which is returned to the caller.

    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and target images.

    Returns
    -------
    target : numpy array
        3D array of all raw target images stored in dirstar.
    """
    
    files = glob.glob(os.path.join(dirtarget, '*.fit'))
    target = []
    for file in files:
        hdulist = fits.open(file)
        if hdulist[0].header['IMAGETYP'] == 'Light Frame':
            image = fits.getdata(file)
            target.append(image)
        hdulist.close()
    target = np.array(target, dtype = float)
    return target
