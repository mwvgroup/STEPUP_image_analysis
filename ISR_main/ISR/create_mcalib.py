import numpy as np
from astropy.io import fits
import os

def create_mbias(biases, bias_prihdr, dirtarget):
    """Creates and saves master bias.

    Takes median of biases along the first axis. This creates
    the master bias array that is used to reduce the science imgaes.
    It then saves that FITS file in a folder "mcalib" within dirtarget.

    Parameters
    ----------
    biases : numpy.ndarray
        3D array containing all bias images found in dirtarget.
    bias_prihdr : astropy.io.fits.header.Header
        Primary HDR from HDU of bias images.
    dirtarget : str
        Directory containing all bias, flat, and science images.

    Returns
    -------
    mbias : numpy.ndarray
        2D array containing master bias image.
    """
    mbias = np.median(biases, 0)

    os.mkdir(dirtarget + '/mcalib')

    hdu = fits.PrimaryHDU(mbias, header=bias_prihdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(dirtarget + '/mcalib/mbias.fits')
        
    return mbias

def create_mdark(darks, mbias, dark_prihdr, dirtarget, dark_exptime, exptime):
    """Creates and saves master dark.
    
    Subtracts mbias out of each individual dark. Then, takes median of
    darks of this array along the first axis. It then time-corrects each image
    to match the light frame exposure time. This creates the master dark
    array that is used to reduce the science images. It then saves that
    FITS file in a folder "mcalib" within dirtarget.
    
    Parameters
    ----------
    darks : numpy.ndarray
        3D array containing all dark images found in dirdark.
    mbias : numpy.ndarray
        2D array containing master bias image.
    dark_prihdr : astropy.io.fits.header.Header
        Primary HDR from HDU of dark images.
    dirtarget : str
        Directory containing all bias, flat and science images.
    dark_exptime : float
        Exposure time of dark frame in seconds.
    exptime : float
        Exposure time of light frame in seconds.
    
    Returns
    -------
    mdark : numpy.ndarray
        2D array containing master dark image.
    """      
    bias_subtracted_darks = []
    for dark in darks:
        dark -= mbias
        dark /= dark_exptime
        dark *= exptime
        bias_subtracted_darks.append(dark)

    bias_subtracted_darks = np.array(bias_subtracted_darks, dtype=float)

    mdark = np.median(darks, 0)

    hdu = fits.PrimaryHDU(mdark, header=dark_prihdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(dirtarget + '/mcalib/mdark.fits', overwrite=True)

    return mdark

def create_mflat(r_flats, b_flats, v_flats, mbias, r_flat_prihdr,
                 b_flat_prihdr, v_flat_prihdr, r_filter, b_filter,
                 v_filter, dirtarget):
    """Creates and saves master flat.
    
    Subtracts mbias out of and normalizes each individual flat.
    It then takes the average along the first axis. This creates the master
    flat array that is used to reduce the science images. It then saves
    that FITS file in a folder "mcalib" within dirtarget.

    Parameters
    ----------
    flats : numpy.ndarray
        3D array containing all flat images found in dirtarget.
    mbias : numpy.ndarray
        2D array containing master bias image.
    flat_prihdr : astropy.io.fits.header.Header
        HDR from HDU of flat images.
    dirtarget : str
        Directory containing all bias, flat, and science images.

    Returns
    -------
    mflat : numpy.ndarray
        2D array containing master flat image.
    """
    r_mflat = None
    b_mflat = None
    v_mflat = None
    if r_filter:
        normalized_r_flats = []
        for flat in r_flats:
            flat -= mbias
            flat /= np.average(r_flat)
            normalized_r_flats.append(flat)
        normalized_r_flats = np.array(normalized_r_flats, dtype=float)
    
        r_mflat = np.average(normalized_r_flats, 0)

        hdu = fits.PrimaryHDU(mflat, header=r_flat_prihdr)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(dirtarget + '/mcalib/r_mflat.fits', overwrite=True)

    if b_filter:
        normalized_b_flats = []
        for flat in b_flats:
            flat -= mbias
            flat /= np.average(b_flat)
            normalized_b_flats.append(flat)
        normalized_b_flats = np.array(normalized_b_flats, dtype=float)
    
        b_mflat = np.average(normalized_b_flats, 0)

        hdu = fits.PrimaryHDU(mflat, header=b_flat_prihdr)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(dirtarget + '/mcalib/b_mflat.fits', overwrite=True)

    if v_filter:
        normalized_v_flats = []
        for flat in v_flats:
            flat -= mbias
            flat /= np.average(b_flat)
            normalized_b_flats.append(flat)
        normalized_b_flats = np.array(normalized_b_flats, dtype=float)
    
        b_mflat = np.average(normalized_b_flats, 0)

        hdu = fits.PrimaryHDU(mflat, header=b_flat_prihdr)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(dirtarget + '/mcalib/b_mflat.fits', overwrite=True)
                 
    return r_mflat, b_mlat, v_mflat

