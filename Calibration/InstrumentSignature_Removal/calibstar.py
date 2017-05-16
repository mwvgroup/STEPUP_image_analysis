'''
NAME: reduction

PURPOSE:
Remove all instrument-signatures that may appear on raw array of target images.

EXPLANATION:
Calls get_darks and get_star to get all needed arrays. Also gets EXPTIME for
dark images in order to time correct the master dark array. Then loops
through array of target images and subtracts the imperfections
(instrument-singatures) from the images by subtracting the
mdarks, mbias, and mflat from each image and appends it to a new
array, calibstar. It then returns this 3D numpy array to the caller.

INPUTS:
(mbias) - this is the module variable equal to the result of calling
create_mbias in the main function.
(mdark) - this is the module variable equal to the result of calling
create_mdark in the main function.
(mflat) - this is the module variable equal to the result of calling
create_mflat in the main function.

OUTPUTS:
(calibstar) - 3D numpy array of all reduced target images.

RESTRICTIONS:
Function must accept three parameters; mbias, mflat, and mdark, which are
2D numpy arrays.
'''

from astropy.io import fits
import numpy as np

def calibstar(mbias, mflat, mdark):
    star = get_star()
    darks = get_darks
    hdulist = fits.open(darks)
    exptime = hdulist.header['EXPTIME']
    calibstar = []

    for image in star:
        image -= exptime*mdark
        image -= mbias
        image /= mflat
        calibstar.append(image)

    calibstar = np.array(image)
    return calibstar
