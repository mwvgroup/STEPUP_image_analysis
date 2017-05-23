import numpy as np
from astropy.io import fits

def add_saturation(mbias, mdark, mflat, dark_exptime):
    """Adds all saturations of reduced images to their headers.

    Calls calibstar to get reduced array of target images. It then subtracts
    the median of the master bias and time-corrected master dark, and divides
    by the median of the master flat. This sequence *=0.97 yields the final
    saturation of each image. It then opens the primary HDR of each image
    and adds the saturation.

    Parameters
    ----------
    mbias : numpy array
        2D array containing master bias image.
    mdark : numpy array
        2D array containing master dark image.
    mflat : numpy array
        2D array containing master flat image.
    dark_exptime : float
        Module variable that is the exposure time of the dark images.

    Returns
    -------
    calibstar : numpy array
        3D array of reduced target images with saturation levels in headers.
    """
    
    calibstar = calibstar(mbias, mflat, mdark, dark_exptime)

    saturation = 65535
    saturation -= np.median(mbias)
    saturation -= dark_exptime*np.median(mdark)
    saturation /= np.median(mflat)
    saturation *= 0.97
    
    for file in calibstar:
        hdulist = fits.open(file)
        prihdr = hdulist[0].header
        prihdr['SATLEVEL'] = saturation

    return calibstar
