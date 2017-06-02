import numpy as np
from astropy.io import fits
import os
import glob

def get_calibimages(dirtarget, dirdark):
    """Retrieves flats and biases from dirtarget and darks from dirdark.

    Searches for all files in dirtarget/dirdark with primary HDU header
    keyword, "IMAGETYP", to equal "Bias Frame", "Flat Field", or "Dark Frame".
    It then appends each file to an array, biases, flats, or darks. It also
    retrieves the exposure time for the dark images (dark_exptime). It then
    converts each array to a numpy array. These arrays are returned to the
    caller in the order biases, darks, flats.

    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and science images.
    dirdark : str
        Directory containing all dark images.

    Returns
    -------
    biases : numpy.ndarray
        3D array containing all bias images found in dirtarget.
    darks : numpy.ndarray
        3D array containing all dark images found in dirdark.
    flats : numpy.ndarray
        3D array containing all flat images found in dirtarget.
    dark_exptime : float
        Exposure time of dark images.

    """

    files = glob.glob(os.path.join(dirtarget, '*.fit'))
    d_files = glob.glob(os.path.join(dirdark, '*.fit'))

    biases = []
    darks = []
    flats = []
    dark_exptime = None
    
    for file in files:
        hdulist = fits.open(file)
        if hdulist[0].header['IMAGETYP'] == 'Bias Frame':
            image = fits.getdata(file)
            biases.append(image)
        if hdulist[0].header['IMAGETYP'] == 'Flat Field':
            image = fits.getdata(file)
            flats.append(image)
        hdulist.close()
        
    for file in d_files:
        hdulist = fits.open(file)
        if hdulist[0].header['IMAGETYP'] == 'Dark Frame':
            dark_exptime = hdulist[0].header['EXPTIME']
            image = fits.getdata(file)
            darks.append(image)

    biases = np.array(biases, dtype=float)
    darks = np.array(darks, dtype=float)
    flats = np.array(flats, dtype=float)

    return(biases, darks, flats, dark_exptime)
