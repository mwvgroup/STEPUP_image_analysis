'''INPUTS:
(dirstar) - module variable, directory in which all flat, bias, and
target images are located: /home/depot/STEPUP/raw/<date>
'''

'''OUTPUTS:
(star) - 3D numpy array of all target images located in dirstar
'''

from astropy.io import fits
import numpy as np
import os

def get_star():
    files = [f for f in os.listdir(dirstar) if f.endswith('.fit')]
    hdulist = fits.open(files)
    star = []
    for image in files:
        if hdulist[0].header['IMAGETYP'] == 'Light Frame':
            star.append(image)
    hdulist.close()
    star = np.array(star, dtype = float)
    return star
