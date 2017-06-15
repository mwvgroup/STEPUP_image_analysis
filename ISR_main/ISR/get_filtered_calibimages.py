import glob
import os
from astropy.io import fits
import numpy as np

def get_filtered_calibimages(dirtarget):
    """Creates and saves master flat for each filter.

    Creates a list of all unique filter names on flat and light
    images. It then loops through that list looking for all flats
    with the same filter name for each iteration of the loop. It
    then subtracts the bias and normalizes each individual flat and
    averages that array. The master flat is then saved as
    "<filter-name>_mflat.fits" to dirtarget/mcalib.

    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and science images.

    Returns
    -------
    n/a
    """
    image_filters = []
    mbias = []

    files = glob.glob(os.path.join(dirtarget, '*.fit'))
    calib_files = glob.glob(os.path.join(dirtarget + '/mcalib/', '*.fits'))

    for file in calib_files:
        hdulist = fits.open(file)
        if hdulist[0].header['IMAGETYP'] == 'Bias Frame':
            image = fits.getdata(file)
            mbias.append(image)
            
    mbias = np.array(mbias, dtype=float)

    for image in files:
        
        hdulist = fits.open(image)
        
        if hdulist[0].header['IMAGETYP'] == 'Flat Field':
            image_filters.append(hdulist[0].header['FILTER'])
            
        elif hdulist[0].header['IMAGETYP'] == 'Light Frame': # Unsure of whether or not I really need this line
            image_filters.append(hdulist[0].header['FILTER'])

    image_filters = np.unique(image_filters)
    
    for i in image_filters:
        flats = []
        flat_prihdr = None
        
        for file in files:
        
            hdulist = fits.open(file)

            if hdulist[0].header['IMAGETYP'] == 'Flat Field':

                if hdulist[0].header['FILTER'] == i:
                    image = fits.getdata(file)
                    flats.append(image)
                    flat_prihdr = hdulist[0].header

            hdulist.close()
        
        flats = np.array(flats, dtype=float)

        for flat in flats:
            flat -= mbias[0]
            flat /= np.average(flats)

        mflat = np.average(flats, 0)

        hdu = fits.PrimaryHDU(mflat, header=flat_prihdr)
        hdulist= fits.HDUList([hdu])
        hdulist.writeto(dirtarget + '/mcalib/' + i + '_mflat.fits', overwrite=True)
    
            
                    

                    
