'''
NAME: get_biases

PURPOSE:
Retrieves all biase images from dirstar.

EXPLANATION:
Looks at primary HDU headers of files in dirstar for ones with keyword
"IMAGETYP" pointing to "Bias Frame". It then puts all of these files into a
3D numpy array and returns that array to the caller.

INPUTS:
(dirstar) - module variable which is thedirectory that
image, flat, and bias files are saved in: /home/depot/STEPUP/raw/<date>

OUTPUTS:
(biases) - 3D numpy array of all bias images within dirstar

RESTRICTIONS:
dirstar must be a string.
'''

from astropy.io import fits
import numpy as np
import os
import glob

def get_biases():
    files = glob.glob(os.path.join(dirstar, '*.fit'))
    biases = []
    for file in files:
        hdulist = fits.open(file)
        if hdulist[0].header['IMAGETYP'] == 'Bias Frame':
            image = fits.getdata(file)
            biases.append(image)
        hdulist.close()
    biases = np.array(biases, dtype = float)
    return biases
