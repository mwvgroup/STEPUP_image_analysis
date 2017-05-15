'''INPUTS:
(dirstar) - module variable which is thedirectory that
image, flat, and bias files are saved in: /home/depot/STEPUP/raw/<date>
'''

'''OUTPUTS:
(flats) - 3D numpy array of all flat images stored in dirstar
'''

from astropy.io import fits
import numpy as np
import os

def get_flats():
    files = [f for f in os.listdir(dirstar) if f.endswith('.fit')]
    hdulist = fits.open(files)
    flats = []
    for image in files:
        if hdulist[0].header('IMAGTYP') == 'Flat Frame':
            flats.append(image)
    hdulist.close()
    flats = np.array(flats, dtype = float)
    return flats

