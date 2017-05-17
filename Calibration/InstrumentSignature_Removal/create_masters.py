import numpy as np

def create_mcalib(dirstar, dirdark):
    """Creates master flat, bias, and dark arrays.

    Calls get_images to get 3D arrays of biases, darks, flats and target
    images. Each array of calibration images is passed into their respective
    functions in order to create the master bias, dark, and flat. Each function
    is called in order (bias, dark, flat) as each preceding master image is
    used in the creation of the next type. Each array of each type of master
    calibration image is then returned to the user in a tuple ordered mbias,
    mdark, mflat, star.
    
    Parameters
    ----------
    dirstar : str
        Directory in which all bias, flat and target images are stored.
    dirdark : str
        Directory in which all dark images are stored.

    Returns
    -------
    mbias : numpy array
        2D numpy array containing the master bias image.
    mdark : numpy array
        2D numpy array containing the master dark image.
    mflat : numpy array
        2D numpy array containing the master flat image.
    star : numpy array
        3D numpy array containing all raw target images in dirstar.
    """
    
    def create_mbias(biases):
        """Creates master bias array.

        Takes median of biases along the first axis. This creates
        the master bias image that is used to reduce the raw
        target imgaes.

        Parameters
        ----------
        biases : numpy array
            3D array containing all bias images found in dirstar.        
    
        Returns
        -------
        mbias : numpy array
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
        darks : numpy array
            3D array containing all dark images found in dirdark.
        mbias : numpy array
            2D array containing master bias image.
    
        Returns
        -------
        mdark : numpy array
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
        flats : numpy array
            3D array containing all flat images found in dirstar.
        mbias : numpy array
            2D array containing master bias image.
        mdark : numpy array
            2D array containing master dark image.
    
        Returns
        -------
        mflat : numpy array
            2D array containing master flat image.
        """
        
        mflat = ((np.median(flats, 0) - mdark - mbias)/np.mean(flats, 0))
        return mflat

    (biases, darks, flats, star) = get_images(dirstar, dirdark)
    
    mbias = create_mbias(biases)
    mdark = create_mdark(mbias, darks)
    mflat = create_mflat(mbias, mdark, flats)
    
    return (mbias, mdark, mflat, star)

