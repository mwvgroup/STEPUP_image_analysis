import numpy as np
from astropy.io import fits
import os
import glob

def get_calibimages(dirtarget, dirdark):
    """Retrieves flats and biases from dirtarget and darks from dirdark.

    Searches for all files in dirtarget/dirdark with primary HDR header
    keyword, "IMAGETYP", to equal "Bias Frame", "Flat Field", or "Dark Frame".
    It then appends each file to an array, biases, flats, or darks. It also
    retrieves the exposure times for the dark images (dark_exptime) and
    light frame images (exptime). It then gets the primary HDRs for each
    type of calibration image. It then converts each array to a numpy array.
    These arrays are returned to the caller in the order biases, darks, flats.

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
        Exposure time of dark frame in seconds.
    exptime : float
        Exposure time of light frame in seconds.
    bias_prihdr : astropy.io.fits.header.Header
        Primary HDR from HDU of bias images.
    dark_prihdr : astropy.io.fits.header.Header
        Primary HDR from HDU of dark images.
    flat_prihdr : astropy.io.fits.header.Header
        HDR from HDU of flat images.   
    """

    files = glob.glob(os.path.join(dirtarget, '*.fit'))
    d_files = glob.glob(os.path.join(dirdark, '*.fit'))

    biases = []
    darks = []
    r_flats = []
    b_flats = []
    v_flats = []

    dark_exptime = None
    exptime = None
    bias_prihdr = None
    dark_prihdr = None
    r_flat_prihdr = None
    b_flat_prihdr = None
    v_flat_prihdr = None
    r_filter = False
    b_filter = False
    v_filter = False
    
    for file in files:
        hdulist = fits.open(file)
        if hdulist[0].header['IMAGETYP'] == 'Bias Frame':
            bias_prihdr = hdulist[0].header
            image = fits.getdata(file)
            biases.append(image)
        if hdulist[0].header['IMAGETYP'] == 'Light Frame':
            exptime = hdulist[0].header['EXPTIME']
        if hdulist[0].header['IMAGETYP'] == 'Flat Field':
            if hdulist[0].header.['FILTER'] == 'r':
                r_filter = True
                r_flat_prihdr = hdulist[0].header
                image = fits.getdata(file)
                r_flats.append(image)
            if hdulist[0].header.['FILTER'] == 'b':
                b_filter = True
                b_flat_prihdr = hdulist[0].header
                image = fits.getdata(file)
                b_flats.append(image)
            if hdulist[0].header.['FILTER'] == 'v':
                v_filter = True
                b_flat_prihdr = hdulist[0].header
                image = fits.getdata(file)
                v_flats.append(image)
        hdulist.close()

    for file in d_files:
        hdulist = fits.open(file)
        if hdulist[0].header['IMAGETYP'] == 'Dark Frame':
            dark_prihdr = hdulist[0].header
            dark_exptime = hdulist[0].header['EXPTIME']
            dark_exptime = float(dark_exptime)
            image = fits.getdata(file)
            darks.append(image)

    biases = np.array(biases, dtype=float)
    darks = np.array(darks, dtype=float)
    r_flats = np.array(r_flats, dtype=float)
    b_flats = np.array(b_flats, dtype=float)
    v_flats = np.array(v_flats, dtype=float)

    return (biases, darks, r_flats, b_flats, v_flats, dark_exptime,
            exptime, bias_prihdr, dark_prihdr, r_flat_prihdr, b_flat_prihdr,
            v_flat_prihdr, r_filter, b_filter, v_filter)
