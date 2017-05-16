from astropy.io import fits
import numpy as np

def calibstar(mbias, mdark, mflat, dark_exptime):
    """Creates reduced array of target images.

    Calls get_star to retrieve all raw target images and subtracts the time-
    corrected master dark, the masterbias, and divides by the master flat.
    This removes the instrument signatures from the raw dataset.

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
    calibstar - numpy array
        3D array containing reduced array of target images.
    """
    
    star = get_star(dirstar)
    calibstar = []

    for image in star:
        image -= dark_exptime*mdark
        image -= mbias
        image /= mflat
        calibstar.append(image)

    calibstar = np.array(image)
    return calibstar
