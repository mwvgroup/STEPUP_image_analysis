'''
NAME: create_mflat

PURPOSE:
Create master flat array that will be subtracted from each image in raw target
dataset to remove imperfections.

EXPLANATION:
Calls get_flats and removes imperfections by subtracting mdark and mbias. Also
takes the median and mean along 3rd axis. Returns array, mflat, to user.

INPUTS:
(mbias) - this is the module variable equal to the result of calling
create_mbias in the main function.
(mdark) - this is the module variable equal to the result of calling
create_mdark in the main function.

OUTPUTS:
(mflat) - a 2D numpy array containing the master flat image which will be
subtracted from each image in the raw target dataset.

RESTRICTIONS:
Function must accept two parameters; mdark and mbias, which are 2D numpy arrays.
'''

import numpy as np

def create_mflat(mbias, mdark):
    flats = get_flats()
    mflat = ((np.median(flats, 2) - mdark - mbias)/np.mean(flats, 2))
    return mflat
