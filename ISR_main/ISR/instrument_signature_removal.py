import os
import glob
from astropy.io import fits
import numpy as np

def instrument_signature_removal(dirtarget, mbias, mdark, mflat, dark_exptime, target):
    """Removes instrument signatures from science images.

    Searches in dirtarget for all FITS files with "IMAGETYP" keyword that
    equals "Light Frame". It then subtracts the master calibration images
    from all "Light Frame" images and then appends them to an array. It then
    adds the saturation to the header of each image and saves the reduced
    image to a new directory.

    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    mbias : numpy.ndarray
        2D array containing master bias image.
    mdark : numpy.ndarray
        2D array containing master dark image.
    mfalt : numpy.ndarray
        2D array containing master flat image.

    Returns
    -------
    isr_scimages : numpy.ndarray
        3D array containing all science images with instrument signatures
        removed.
    """

    saturation = 65535
    saturation -= np.median(mbias)
    saturation -= np.median(mdark*dark_exptime)
    saturation /= np.median(mflat)
    saturation *= 0.97

    files = glob.glob(os.path.join(dirtarget, '*.fit'))
    scimages = []
    prihdr = None

    for file in files:
        hdulist = fits.open(file)
        if hdulist[0].header['IMAGETYP'] == 'Light Frame':
            hdulist[0].header['SATLEVEL'] = (saturation)
            exptime = hdulist[0].header['EXPTIME']
            prihdr = hdulist[0].header
            image = fits.getdata(file)
            scimages.append(image)
        hdulist.close()
        
    scimages = np.array(scimages, dtype=float)
    isr_scimages = []
    n = 0

    for file in scimages:
            file -= mbias
            file -= mdark
            file /= mflat
            isr_scimages.append(file)
            
    isr_scimages = np.array(isr_scimages, dtype = float)

    os.mkdir(dirtarget + '/ISR_Images')

    for i in isr_scimages:
        n+=1
        hdu = fits.PrimaryHDU(i,header=prihdr)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(dirtarget + '/ISR_Images/' + target + '_{}.fits'.format(n),
                        overwrite=True)
    return isr_scimages
