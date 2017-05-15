import numpy as np
from astropy.io import fits

def add_saturation():
    mbias = create_mbias()
    mflat = create_mflat()
    mdark = create_mdark()
    darks = get_darks()
    star = get_star()
    hdulist = fits.open(darks)
    exptime = hdulist.header['EXPTIME']

    saturation = 65535
    saturation -= np.median(mbias)
    saturation -= exptime*np.median(mdark)
    saturation /= np.median(mflat)
    saturation *= 0.97

    prihdr = hdulist[0].header
    prihdr['SATLEVEL'] = saturation
