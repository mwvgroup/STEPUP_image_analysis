from astropy.io import fits
import numpy as np

def isr_removal(dirtarget, mbias, mdark, mflat, dark_exptime):
    """Creates reduced array of target images.

    Loops through images in dirtarget

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
    isr_science_images : numpy array
        3D array containing instrument signature removed array of
        target images.
    """

    n = 0

    for image in dirtarget:
        n += 1
        hdulist = fits.open(image)
        if hdulist[0].header['IMAGETYP'] == 'Light Frame':
            image -= dark_exptime*mdark
            image -= mbias
            image /= mflat
            hdulist.writeto('scimage_isr{}.fits'.format(n))

    
             
            
