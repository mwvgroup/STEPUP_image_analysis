import subprocess
import glob
import os
from astropy.io import fits


def perform_astrometry(dirtarget, filters):
    """Performs astrometry on dataset.

    Finds coordinates of stars using imstar and writes information to a
    table called new-image.tab. It then adds WCS information to the headers
    of all the images in the dataset to get a preliminary estimate of the
    true WCS information. It writes these files and then uses imstar to
    correct the WCS information in each header. These files are written to
    the same directory.

    Parameters
    ----------
    filters : list
        List of strings containing all filters used for observation.

    Returns
    -------
    None
    """
    dirtarget = dirtarget + '/ISR_Images'
    # Retrieve new-image.fits using API and save it to ISR_Images.
    # Value for -i argument?
    os.chdir(dirtarget)
    # Creates new-images.tab file with positions of stars.
    subprocess.call(['imstar', '-vhi', '700', '-tw', 'new-image.fits'])

    # Gets header with WCS information to append to all images.
    # Shorten this line.
    wcsim_hdu = fits.open(dirtarget + '/new-image.fits')
    wcsim_header = wcsim_hdu[0].header

    # for i, fil enumerate
    for fil in filters:
        os.mkdir(dirtarget + '/{}/WCS'.format(fil))
        # Shorten this line.
        isr_images = glob.glob(os.path.join(dirtarget + '/{}'.format(fil),
                                            '*.fits'))
        n = 0
        for image in isr_images:
            n += 1
            print('\niteration....\n')
            other_hdu = fits.open(image)
            imagedata = other_hdu[0].data
            other_header = other_hdu[0].header
            # Finds all uncommon header keywords.
            diff = fits.HeaderDiff(wcsim_header, other_header).diff_keywords
            diff = diff[0]

            for i in diff:
                if i == 'COMMENT' or i == 'HISTORY':
                    print("skipping....")
                else:
                    # Adds uncommon keywords and their value to image.
                    other_header.set(i, wcsim_header[i])

            # Writes file.
            hdu = fits.PrimaryHDU(imagedata, header=other_header)
            hdulist = fits.HDUList([hdu])
            # Shorten this line.
            hdulist.writeto(dirtarget + '/{}/WCS/wcs{}.fits'.format(fil, n),
                            overwrite=True)

        # Corrects WCS information in image header using known star
        # coordinates in new-image.tab
        for i in range(1, n):
            os.chdir(dirtarget + '/{}/WCS'.format(fil))
            subprocess.call(['imwcs', '-wv', '-i', '100', '-c',
                             'new-image.tab', 'wcs{}.fits'.format(i)])

        os.mkdir(dirtarget + '/{}/WCS/accurate_WCS'.format(fil))
        subprocess.call(['mv', dirtarget + '/new-image.tab'.format(fil),
                         dirtarget + '/{}/WCS'.format(fil)])
        for path in os.listdir(dirtarget + "/{}/WCS/".format(fil)):
            if path.endswith("w.fits"):
                subprocess.call(['mv', dirtarget + '/{}/WCS/'.format(fil) +
                                 path, dirtarget + '/{}/WCS/accurate_WCS'
                                 .format(fil)])
