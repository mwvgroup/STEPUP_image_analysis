from astropy.io import fits
import numpy as np

def calibtarget(target, mbias, mdark, mflat, dark_exptime):
    """Creates reduced array of target images.

    Calls get_target to retrieve all raw target images and subtracts the time-
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
    calibtarget - numpy array
        3D array containing reduced array of target images.
    """
    calibtarget = []

    for image in target:
        image -= dark_exptime*mdark
        image -= mbias
        image /= mflat
        calibtarget.append(image)

    calibtarget = np.array(calibtarget)
    return calibtarget
