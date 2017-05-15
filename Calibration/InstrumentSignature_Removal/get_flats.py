'''
NAME: get_flats

PURPOSE:
Get all flat files from dirstar to use for image reduction.

EXPLANATION:
Searches in the primary HDU headers of all files in dirstar to find ones with
keyword "IMAGETYP" to point to "Flat Frame". It then puts these FITS files into
a 3D numpy array and returns it to the caller.

INPUTS:
(dirstar) - module variable which is the directory that 
image, flat, and bias files are stored in: /home/depot/STEPUP/raw/<date>

OUTPUTS:
(flats) - 3D numpy array of all flat images stored in dirstar

RESTRICTIONS:
dirstar must be a string.
'''

from astropy.io import fits
import numpy as np
import os
import glob

def get_flats():
    files = glob.glob(os.path.joing(dircalib, '*.fit'))
    flats = []
    for file in files:
        hdulist = fits.open(file)
        if hdulist[0].header('IMAGTYP') == 'Flat Frame':
            image = fis.getdata(file)
            flats.append(image)
        hdulist.close()
    flats = np.array(flats, dtype = float)
    return flats

