import os
import numpy as np
from astropy.io import fits
import glob

def get_calibimages(dirtarget, dirdark):
    """Retrieves all bias, flat, and dark images from dirtarget/dirbias.

    Searches in primary HDU headers of all files in dirtarget/dirdark for files
    with keyword "IMAGETYP" that points to "Bias Frame", "Flat Field", or
    "Dark Frame". It then appends each image type to its respective array.

    Parameters
    ----------
    dirtarget : str
        Directory containing all flat, bias, and target images.
    dirdark : str
        Directory containing all dark images.

    Returns
    -------
    biases : numpy array
        3D array of all bias images stored in dirtarget.
    darks : numpy array
        3D numpy array of all dark images located within dirdark.
    flats : numpy array
        3D array of all flat images stored in dirtarget.
    """
    
    files = glob.glob(os.path.join(dirtarget, '*.fit'))
    d_files = glob.glob(os.path.join(dirdark, '*.fit'))
    biases = []
    flats = []
    darks = []

    for file in files:
        hdulist = fits.open(file)
        if hdulist[0].header['IMAGETYP'] == 'Bias Frame':
            image = fits.getdata(file)
            biases.append(image)
        if hdulist[0].header['IMAGETYP'] == 'Flat Field':
            image = fits.getdata(file)
            flats.append(image)
        hdulist.close()

    biases = np.array(biases, dtype=float)
    flats = np.array(flats, dtype=float)

    dark_exptime=None

    for file in d_files:
        hdulist = fits.open(file)
        if hdulist[0].header['IMAGETYP'] == 'Dark Frame':
            dark_exptime = hdulist[0].header['EXPTIME']
            image = fits.getdata(file)
            darks.append(image)
            
        hdulist.close()

    darks = np.array(darks, dtype = float)

    return (biases, darks, flats, dark_exptime)
