'''INPUTS:
(dirstar) - module variable which is thedirectory that
image, flat, and bias files are saved in: /home/depot/STEPUP/raw/<date>
'''

'''OUTPUTS:
(biases) - 3D numpy array of all bias images within dirstar
'''

from astropy.io import fits
import numpy as np
import os

def get_biases():
    files = [f for f in os.listdir(dirstar) if f.endswith('.fit')]
    hdulist = fits.open(files)
    biases = []
    for image in files:
        if hdulist[0].header['IMAGETYP'] == 'Bias Frame':
            biases.append(image)
    hdulist.close()
    biases = np.array(biases, dtype = float)
    return biases
