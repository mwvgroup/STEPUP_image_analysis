import os
import glob
from astropy.io import fits
import numpy as np

def instrument_signature_removal(dirtarget, target, mbias, mdark, mflat, dark_exptime, exptime):
    """Removes instrument signatures from science images.

    Searches in dirtarget for all FITS files with "IMAGETYP" keyword that
    equals "Light Frame". It then subtracts the master calibration images
    from all "Light Frame" images and then appends them to an array. It then
    adds the saturation to the header of each image and saves the reduced
    image to a folder "ISR_Images" in dirtarget.

    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    target : str
        Name of target object.
    mbias : numpy.ndarray
        2D array containing master bias image.
    mdark : numpy.ndarray
        2D array containing master dark image.
    mflat : numpy.ndarray
        2D array containing master flat image.
    dark_exptime : float
        Exposure time of dark frame in seconds.
    exptime : float
        Exposure time of light frame in seconds.

    Returns
    -------
    isr_scimages : numpy.ndarray
        3D array containing all science images with instrument signatures
        removed.
    """
    r_saturation = None
    b_saturation = None
    v_saturation = None
    r_scimages = []
    b_scimages = []
    v_scimages = []
    rim_prihdr = None
    bim_prihdr = None
    vim_prihdr = None
    r_isr_scimages = None
    b_isr_scimages = None
    v_isr_scimages = None

    os.mkdir(dirtarget + '/ISR_Images')
    
    files = glob.glob(os.path.join(dirtarget, '*.fit'))
    
    if r_filter:
        r_saturation = 65535
        r_saturation -= np.median(mbias)
        r_saturation -= np.median((mdark*exptime)/dark_exptime)
        r_saturation /= np.average(r_mflat)
        r_saturation *= 0.97

        for file in files:
            hdulist = fits.open(file)
            if hdulist[0].header['IMAGETYP'] == 'Light Frame':
                if hdulist[0].header['FILTER'] == 'r':
                    hdulist[0].header['SATLEVEL'] = r_saturation
                    exptime = hdulist[0].header['EXPTIME']
                    rim_prihdr = hdulist[0].header
                    image = fits.getdata(file)
                    r_scimages.append(image)
            hdulist.close()

        r_scimages = np.array(r_scimages, dtype=float)
        r_isr_scimages = []
        n = 0

        for file in scimages:
                file -= mbias
                file -= mdark
                file /= r_mflat
                r_isr_scimages.append(file)
                
        r_isr_scimages = np.array(r_isr_scimages, dtype = float)

        for i in r_isr_scimages:
            n+=1
            hdu = fits.PrimaryHDU(i,header=rim_prihdr)
            hdulist = fits.HDUList([hdu])
            hdulist.writeto(dirtarget + '/ISR_Images/' + target + '-R_{}.fits'.format(n),
                            overwrite=True)

    if b_filter:
        b_saturation = 65535
        b_saturation -= np.median(mbias)
        b_saturation -= np.median((mdark*exptime)/dark_exptime)
        b_saturation /= np.average(b_mflat)
        b_saturation *= 0.97

        for file in files:
            hdulist = fits.open(file)
            if hdulist[0].header['IMAGETYP'] == 'Light Frame':
                if hdulist[0].header['FILTER'] == 'b':
                    hdulist[0].header['SATLEVEL'] = b_saturation
                    exptime = hdulist[0].header['EXPTIME']
                    bim_prihdr = hdulist[0].header
                    image = fits.getdata(file)
                    b_scimages.append(image)
            hdulist.close()

        b_scimages = np.array(b_scimages, dtype=float)
        b_isr_scimages = []
        n = 0

        for file in scimages:
                file -= mbias
                file -= mdark
                file /= b_mflat
                b_isr_scimages.append(file)
                
        b_isr_scimages = np.array(b_isr_scimages, dtype = float)

        for i in b_isr_scimages:
            n+=1
            hdu = fits.PrimaryHDU(i,header=bim_prihdr)
            hdulist = fits.HDUList([hdu])
            hdulist.writeto(dirtarget + '/ISR_Images/' + target + '-B_{}.fits'.format(n),
                            overwrite=True)

    if v_filter:
        v_saturation = 65535
        v_saturation -= np.median(mbias)
        v_saturation -= np.median((mdark*exptime)/dark_exptime)
        v_saturation /= np.average(b_mflat)
        v_saturation *= 0.97

        for file in files:
            hdulist = fits.open(file)
            if hdulist[0].header['IMAGETYP'] == 'Light Frame':
                if hdulist[0].header['FILTER'] == 'v':
                    hdulist[0].header['SATLEVEL'] = v_saturation
                    exptime = hdulist[0].header['EXPTIME']
                    prihdr = hdulist[0].header
                    image = fits.getdata(file)
                    v_scimages.append(image)
            hdulist.close()
        
        v_scimages = np.array(v_scimages, dtype=float)
        v_isr_scimages = []
        n = 0

        for file in scimages:
                file -= mbias
                file -= mdark
                file /= v_mflat
                v_isr_scimages.append(file)
                
        v_isr_scimages = np.array(v_isr_scimages, dtype = float)

        for i in v_isr_scimages:
            n+=1
            hdu = fits.PrimaryHDU(i,header=vim_prihdr)
            hdulist = fits.HDUList([hdu])
            hdulist.writeto(dirtarget + '/ISR_Images/' + target + '-V_{}.fits'.format(n),
                            overwrite=True)

    return r_isr_scimages, b_isr_scimages, v_isr_scimages
