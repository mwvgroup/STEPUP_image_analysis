import os
import glob
from astropy.io import fits
import numpy as np

def get_unfiltered_calibimages(dirtarget, dirdark):
    """Creates and saves master bias and dark.
    
    Searches in dirtarget for all bias frames and creates master bias by taking
    the median of the array of those images. It then saves the master bias to 
    dirtarget/mcalib. It then searches in dirdark for all dark frames and 
    subtracts the mbias and time-corrects each image. It then creates the 
    master dark by taking the median of the array of those images. It then 
    saves the master dark to dirtarget/mcalib.

    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and science images.
    dirdark : str
        Directory containing all dark images.

    Returns
    -------
    n/a
    """
    t_files = glob.glob(os.path.join(dirtarget, '*.fit'))

    d_files = glob.glob(os.path.join(dirdark, '*.fit'))
    
    os.mkdir(dirtarget + '/mcalib')    
    biases = []
    bias_prihdr = None

    darks = []
    dark_prihdr = None
    dark_exptime = None
    exptime = None

    for file in t_files:

        hdulist = fits.open(file)

        if hdulist[0].header['IMAGETYP'] == 'Bias Frame':
            image = fits.getdata(file)
            biases.append(image)
            bias_prihdr = hdulist[0].header

    biases = np.array(biases, dtype=float)

    mbias = np.median(biases, 0)

    hdu = fits.PrimaryHDU(mbias, header=bias_prihdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(dirtarget + '/mcalib/mbias.fits')
            

    for file in d_files:

        hdulist = fits.open(file)

        if hdulist[0].header['IMAGETYP'] == 'Dark Frame':
            image = fits.getdata(file)
            darks.append(image)
            dark_exptime = hdulist[0].header['EXPTIME']
            dark_prihdr = hdulist[0].header                 

        if hdulist[0].header['IMAGETYP'] == 'Light Frame':
            exptime = hdulist[0].header['EXPTIME']

    darks = np.array(darks, dtype=float)

    for dark in darks:
        dark -= mbias
        dark /= dark_exptime
        dark *= exptime
        
    mdark = np.median(darks, 0)

    hdu = fits.PrimaryHDU(mdark, header=dark_prihdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(dirtarget + '/mcalib/mdark.fits', overwrite=True)       
