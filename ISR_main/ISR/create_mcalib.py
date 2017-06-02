import numpy as np

def create_mbias(biases):
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
    return mbias

def create_mdark(darks, mbias):
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
    return mdark

def create_mflat(flats, mbias, mdark):
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
    
    mflat = ((np.median(flats, 0) - mdark - mbias)/np.mean(flats, 0))
    return mflat

