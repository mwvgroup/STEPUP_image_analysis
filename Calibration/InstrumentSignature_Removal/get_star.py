from astropy.io import fits
import numpy as np
import os
import glob

def get_star(dirstar):
    '''NAME: get_star

    PURPOSE:
    Get all target images from dirstar.

    EXPLANATION:
    Searches in the primary HDU headers of all files in dirstar for files with
    keyword "IMAGETYP" to point to "Light Frame". It then puts these FITS files
    into a 3D numpy array and returns it to the caller.

    INPUTS:
    (dirstar) - module variable, directory in which all flat, bias, and
    target images are located: /home/depot/STEPUP/raw/<date>

    OUTPUTS:
    (star) - 3D numpy array of all target images located in dirstar

    RESTRICTIONS:
    Must accept on parameter, dirstar, which must be a string.
    '''
    
    files = glob.glob(os.path.join(dirstar, '*.fit'))
    star = []
    for file in files:
        hdulist = fits.open(file)
        if hdulist[0].header['IMAGETYP'] == 'Light Frame':
            image = fits.getdata(file)
            star.append(image)
        hdulist.close()
    star = np.array(star, dtype = float)
    return star
