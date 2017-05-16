from astropy.io import fits
import numpy as np
import os
import glob

def get_biases(dirstar):
    """Retrieves all bias images from dirstar.

    Extended Summary
    ----------------
    Searches in primary HDU headers of all files in dirstar for files with
    keyword "IMAGETYP" that point to "Bias Frame". It then puts all biases
    found into an array which is returned to the caller.

    Parameters
    ----------
    dirstar : string
        module variable which is the directory that
        image, flat, and bias files are saved in
        /home/depot/STEPUP/raw/<date>

    Returns
    -------
    biases : numpy array
        3D array of all bias images stored in dirstar
    """
    
    files = glob.glob(os.path.join(dirstar, '*.fit'))
    biases = []
    for file in files:
        hdulist = fits.open(file)
        if hdulist[0].header['IMAGETYP'] == 'Bias Frame':
            image = fits.getdata(file)
            biases.append(image)
        hdulist.close()
    biases = np.array(biases, dtype = float)
    return biases

