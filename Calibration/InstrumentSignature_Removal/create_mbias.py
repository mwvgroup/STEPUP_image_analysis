'''
NAME: create_mbias

PURPOSE:
Create master bias array that will be used to subtract imperfections from raw
target image data.

EXPLANATION:
Calls get_biases to get the array of all bias images. It then takes the median
of this array along the third axis and returns the array to the caller.

INPUTS:
none

OUTPUTS:
(mbias) - 2D array of master bias image that will be subtracted from each image
in the raw target dataset.
'''

import numpy as np

def create_mbias():
    biases = get_biases()
    mbias = np.median(biases, 2)
    return mbias
