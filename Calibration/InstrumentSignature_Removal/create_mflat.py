import numpy as np

def create_mflat(mbias, mdark):
    """Creates master flat array.

    Calls get_flats to retrieve all flat images and takes the median along
    the third axis. Then subtracts master bias and master dark images. Then
    divides by the median of the flat array. This creates the master flat
    image that will later be used to reduce the raw target images.

    Parameters
    ----------
    mbias : numpy array
        2D array containing master bias image.
    mdark : numpy array
        2D array containing master dark image.

    Returns
    -------
    mflat : numpy array
        2D array containing master flat image.
    """
    
    flats = get_flats()
    mflat = ((np.median(flats, 2) - mdark - mbias)/np.mean(flats, 2))
    return mflat

