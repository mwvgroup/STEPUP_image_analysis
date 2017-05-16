import numpy as np

def create_mdark(mbias):
    '''
    NAME: create_mdark

    PURPOSE:
    Create master dark array that will be used to subtract imperfections from 
    raw target image data.

    EXPLANATION:
    Calls get_darks in order to begin sequence. Then takes the
    median of the dark images along the 3rd axis and subtracts the master bias
    to remove those imperfections.

    INPUTS:
    (mbias) - this is the module variable equal to the result of calling
    create_mbias in the main function.

    OUTPUTS:
    (mdark) - 2D numpy array containing the master dark array that will be
    subtracted from the raw target dataset.

    RESTRICTIONS:
    Function must accept one parameter, mbias, which is a 2D numpy array.
    '''
    
    darks = get_darks()
    mdark = np.median(darks, 2) - mbias
    return mdark
    
