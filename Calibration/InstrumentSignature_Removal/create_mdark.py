import numpy as np

def create_mdark(mbias):
    """Creates array of master dark images.

    Extended Summary
    ----------------
    Calls get_darks to retrieve all dark images and takes the median
    along the 3rd axis. Then subtracts the master bias image. This creates
    the master dark image that will later be used to reduce the raw
    target images.

    Parameters
    ----------
    mbias : numpy array
        2D array containing master bias image.

    Returns
    -------
    mdark : numpy array
        2D array containing master dark image.
    """
    
    darks = get_darks()
    mdark = np.median(darks, 2) - mbias
    return mdark
    
