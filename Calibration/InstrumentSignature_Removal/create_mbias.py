import numpy as np

def create_mbias():
    """Creates master bias array.

    Calls get_biases to retrieve all bias iamges and takes the median
    along the 3rd axis. This creates the master bias image that will
    later be used to reduce the raw target images.

    Parameters
    ----------
    None

    Returns
    -------
    mbias : numpy array
        2D array containing master bias image.
    """
    
    biases = get_biases(dirstar)
    mbias = np.median(biases, 2)
    return mbias
