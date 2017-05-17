import numpy as np

import sys
sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/Calibration/InstrumentSignature_Removal/')

from get_images import get_images

dirstar = dirdark = '/Users/helenarichie/HAT-P-3_Test_Images'

def create_mcalib(dirstar, dirdark, biases, flats, darks):
    
    def create_mbias(dirstar, biases):
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
        
        #biases = get_biases(dirstar)
        mbias = np.median(biases, 0)
        return mbias
    
    
    def create_mdark(dirdark, mbias, darks):
        """Creates array of master dark images.
    
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
        
        #darks = get_darks(dirdark)
        mdark = np.median(darks, 0) - mbias
        return mdark
    
    
    def create_mflat(dirstar, mbias, mdark, biases):
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
        
        #flats = get_flats(dirstar)
        mflat = ((np.median(flats, 0) - mdark - mbias)/np.mean(flats, 0))
        return mflat
    
    
    mbias = create_mbias(dirstar, biases)
    mdark = create_mdark(dirdark, mbias, darks)
    mflat = create_mflat(dirstar, mbias, mdark, flats)
    
    return (mbias, mdark, mflat)

(biases, darks, flats, star) = get_images(dirstar, dirdark)

(mbias, mdark, mflat) = create_mcalib(dirstar, dirdark, biases, flats, darks)
