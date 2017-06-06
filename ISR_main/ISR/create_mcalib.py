import numpy as np
from astropy.io import fits
import os

def create_mbias(biases, bias_prihdr, dirtarget):
    """Creates master bias array.

    Takes median of biases along the first axis. This creates
    the master bias image that is used to reduce the raw
    target imgaes.

    Parameters
    ----------
    biases : numpy.ndarray
        3D array containing all bias images found in dirtarget.        

    Returns
    -------
    mbias : numpy.ndarray
        2D array containing master bias image.
    """
    mbias = np.median(biases, 0)

    os.mkdir(dirtarget + '/mcalib')

    hdu = fits.PrimaryHDU(mbias, header=bias_prihdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(dirtarget + '/mcalib/mbias.fits', overwrite=True)
        
    return mbias

def create_mdark(darks, mbias, dark_prihdr, dirtarget, dark_exptime, exptime):
    """Creates array of master dark images.
    
    Takes median of darks along first axis. Then subtracts mbias.
    This creates the master dark image that is used to reduce the
    raw target images.
    
    Parameters
    ----------
    darks : numpy.ndarray
        3D array containing all dark images found in dirdark.
    mbias : numpy.ndarray
        2D array containing master bias image.
    
    Returns
    -------
    mdark : numpy.ndarray
        2D array containing master dark image.
    """      
    bias_subtracted_darks = []
    for dark in darks:
        dark -= mbias
        dark /= dark_exptime
        dark *= exptime
        bias_subtracted_darks.append(dark)

    bias_subtracted_darks = np.array(bias_subtracted_darks, dtype=float)

    mdark = np.median(darks, 0) - mbias

    hdu = fits.PrimaryHDU(mdark, header=dark_prihdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(dirtarget + '/mcalib/mdark.fits', overwrite=True)

    return mdark

def create_mflat(flats, mbias, flat_prihdr, dirtarget):
    """Creates master flat array.
    
    Takes median of flats along first axis. Then subtracts mbias
    and mdark. Then divides by the median, taken along the first
    axis, of flats. This creates the master flat image that is used
    to reduce the raw target images.

    Parameters
    ----------
    flats : numpy.ndarray
        3D array containing all flat images found in dirtarget.
    mbias : numpy.ndarray
        2D array containing master bias image.
    mdark : numpy.ndarray
        2D array containing master dark image.

    Returns
    -------
    mflat : numpy.ndarray
        2D array containing master flat image.
    """

    bias_subtracted_flats = []
    for flat in flats:
        flat -= mbias
        bias_subtracted_flats.append(flat)
    bias_subtracted_flats = np.array(bias_subtracted_flats, dtype=float)
    
    mflat = np.median(bias_subtracted_flats, 0)/np.mean(bias_subtracted_flats, 0)

    hdu = fits.PrimaryHDU(mflat, header=flat_prihdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(dirtarget + '/mcalib/mflat.fits', overwrite=True)
                 
    return mflat

