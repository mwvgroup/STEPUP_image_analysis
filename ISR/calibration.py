import numpy as np
from astropy.io import fits

def add_saturation(mbias, mdark, mflat, dark_exptime, isr_scimage):
    """Adds all saturations of reduced images to their headers.

    Subtracts the median of the master bias from isr_scimage and
    time-corrected master dark, and divides by the median of the master flat.
    This sequence *=0.97 yields the expected saturation of each image.
    It then opens the primary HDR of each image and adds the saturation.

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
    isr_scimage : numpy array
        3D array containing instrument signature removed science images.

    Returns
    -------
    isr_scimage : numpy array
        3D array containing instrument signature removed science images
        with saturation levels in headers.
    """


    saturation = 65535
    saturation -= np.median(mbias)
    saturation -= dark_exptime*np.median(mdark)
    saturation /= np.median(mflat)
    saturation *= 0.97
    
    for file in isr_scimage:
        hdulist = fits.open(file)
        prihdr = hdulist[0].header
        prihdr['SATLEVEL'] = saturation

    return isr_scimage
