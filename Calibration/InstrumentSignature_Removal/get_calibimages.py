from astropy.io import fits
import numpy as np
import os
import glob

def get_biases(dirtarget):
    """Retrieves all bias images from dirtarget.

    Searches in primary HDU headers of all files in dirstar for files with
    keyword "IMAGETYP" that points to "Bias Frame". It then puts all biases
    found into an array which is returned to the caller.

    Parameters
    ----------
    dirtarget : string
        Directory containing all flat, bias, and target images.

    Returns
    -------
    biases : numpy array
        3D array of all bias images stored in dirtarget.
    """
    
    files = glob.glob(os.path.join(dirtarget, '*.fit'))
    biases = []
    for file in files:
        hdulist = fits.open(file, memmap=True)
        if hdulist[0].header['IMAGETYP'] == 'Bias Frame':
            image = fits.getdata(file)
            biases.append(image)
        hdulist.close()
    biases = np.array(biases, dtype = float)
    return biases

def get_darks(dirdark):
    """Retrieves all dark images from dirdark.

    Searches in primary HDU headers of all files in dirdark for files with
    keyword "IMAGETYP" that points to "Dark Frame". It then puts all darks
    found into an array which is returned to the caller.

    Parameters
    ----------
    dirdark : str
        Directory containing all dark images.

    Returns
    -------
    darks : numpy array
        3D numpy array of all dark images located within dirdark.
    """
    
    files = glob.glob(os.path.join(dirdark, '*.fit'))
    darks = []
    for file in files:
        hdulist = fits.open(file, memmap=True)
        if hdulist[0].header['IMAGETYP'] == 'Dark Frame':
            image = fits.getdata(file)
            darks.append(image)
        hdulist.close()
    darks = np.array(darks, dtype = float)
    return darks

def get_flats(dirtarget):
    """Retrieves all flat files from dirtarget.

    Searches in primary HDU headers of all files in dirtarget for files with
    keyword "IMAGETYP" that points to "Flat Field". It then puts all flats
    found into an array which is returned to the caller.

    Parameters
    ----------
    dirtarget : str
        Directory in which all flat, bias, and target images are stored.

    Returns
    -------
    flats : numpy array
        3D array of all flat images stored in dirtarget.

    """
    
    files = glob.glob(os.path.join(dirtarget, '*.fit'))
    flats = []
    for file in files:
        hdulist = fits.open(file, memmap=True)
        if hdulist[0].header['IMAGETYP'] == 'Flat Field':
            image = fits.getdata(file)
            flats.append(image)
        hdulist.close()
    flats = np.array(flats, dtype = float)
    return flats
