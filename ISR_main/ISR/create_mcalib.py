import numpy as np
from astropy.io import fits

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

    hdu = fits.PrimaryHDU(mbias, header=bias_prihdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(dirtarget + '/mbias.fits', overwrite=True)
        
    return mbias

def create_mdark(darks, mbias, dark_prihdr, dirtarget):
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

    mdark = np.median(darks, 0) - mbias

    hdu = fits.PrimaryHDU(mdark, header=dark_prihdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(dirtarget + + '/mdark.fits', overwrite=True)

    return mdark

def create_mflat(flats, mbias, mdark, flat_prihdr, dirtarget):
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

    mflat = (np.median(flats, 0) - mbias - mdark)/np.mean(flats, 0)

    hdu = fits.PrimaryHDU(mflat, header=flat_prihdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(dirtarget + '/mflat.fits', overwrite=True)
                 
    return mflat

