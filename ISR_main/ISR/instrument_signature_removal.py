import os
import glob
from astropy.io import fits
import numpy as np

def instrument_signature_removal(dirtarget, target, exptime, image_filters):
    """Removes instrument signatures from raw science images.
    """

    os.mkdir(dirtarget + '/ISR_Images')

    files = glob.glob(os.path.join(dirtarget, '*.fit'))
    calib_files = glob.glob(os.path.join(dirtarget + '/mcalib/', '*.fits'))

    for i in image_filters:

        science_images = []
        science_image_prihdr = None
        mbias = []
        mdark = []
        mflat = []
        dark_exptime = None

        for file in calib_files:

            hdulist = fits.open(file)

            if hdulist[0].header['IMAGETYP'] == 'Bias Frame':
                image = fits.getdata(file)
                mbias.append(image)

            if hdulist[0].header['IMAGETYP'] == 'Dark Frame':
                dark_exptime = hdulist[0].header['EXPTIME']
                image = fits.getdata(file)
                mdark.append(image)

            if hdulist[0].header['IMAGETYP'] == 'Flat Field':

                if hdulist[0].header['FILTER'] == i:
                    image = fits.getdata(file)
                    mflat.append(image)
            hdulist.close()

        mbias = np.array(mbias, dtype=float)
        mdark = np.array(mdark, dtype=float)
        mflat = np.array(mflat, dtype=float)
        
        mbias = mbias[0]
        mdark = mdark[0]
        mflat = mflat[0]

        saturation = 65535
        saturation -= np.median(mbias)
        saturation -= np.median((mdark*exptime)/dark_exptime)
        saturation /= np.average(mflat)
        saturation *= 0.97
        saturation = int(saturation)
        

        for file in files:

            hdulist = fits.open(file)

            if hdulist[0].header['IMAGETYP'] == 'Light Frame':

                if hdulist[0].header['FILTER'] == i:
                    hdulist[0].header['SATLEVEL'] = saturation
                    image = fits.getdata(file)
                    science_images.append(image)
                    science_image_prihdr = hdulist[0].header

            hdulist.close()

        science_images = np.array(science_images, dtype=float)

        for file in science_images:
                file -= mbias
                file -= mdark
                file /= mflat

        n = 0

        for j in science_images:
            n += 1
            hdu = fits.PrimaryHDU(j, header=science_image_prihdr)
            hdulist = fits.HDUList([hdu])
            hdulist.writeto(dirtarget + '/ISR_Images/' + target + '_' + i + '_{}.fits'.format(n),
                            overwrite=True)
