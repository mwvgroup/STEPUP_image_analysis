from astropy.io import fits
import numpy as np

def calibstar():
    mbias = get_mbias()
    mflat = get_mflat()
    mdark = get_mdark()
    star = get_star()
    hdulist = fits.open(darks)
    exptime = hdulist,header['EXPTIME']
    calibstar = []

    for image in star:
        image -= exptime*mdark
        image -= mbias
        image /= mflat
        calibstar.append(image)

    calibstar = np.array(image)
    return calibstar
