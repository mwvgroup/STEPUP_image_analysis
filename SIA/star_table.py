from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
import glob
import numpy as np
import os
from photutils import DAOStarFinder
import warnings


def star_table(image, FWHM):
    """Create table of positions and brightnesses of all sources in an image.

    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    fil : str
        Name of filter used for images which are currently being processed.

    Returns
    -------
    None
    """
    o_file = os.path.join(image)
    hdulist = fits.open(o_file)
    header = hdulist[0].header
    data = hdulist[0].data

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        w = WCS(header)

    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    daofind = DAOStarFinder(fwhm=FWHM, threshold=5.*std)
    sources = daofind(data - median)

    x_cent = sources['xcentroid']
    y_cent = sources['ycentroid']
    ra, dec = w.wcs_pix2world(x_cent, y_cent, 1, ra_dec_order=True)
    sources['ra'] = ra
    sources['dec'] = dec

    print('\n', sources)

    path_file = os.path.split(os.path.abspath(image))
    path = path_file[0]
    filename = path_file[1]
    filechars = np.array(list(filename))
    n = np.where(filechars == '.')[0][0]
    filename = filename[:n]

    out_path = os.path.join(path, '{}_startable.txt'.format(filename))
    sources.write(out_path, format='csv', overwrite=True)


def main():
    """Execute star_table.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    t_file = input('\nInput global filename: ')
    while not os.path.exists((os.path.join(t_file))):
        t_file = input('\nThis file does not seem to exist. Check to ensure' +
                       ' that you have entered the correct filename.\n\n' +
                       'Input global filename or enter "Q" to quit: ')
        if t_file.lower() == 'q':
            return None

    FWHM = input('\nInput image FWHM: ')

    def is_float(var):
        try:
            float(var)
            return True
        except ValueError:
            return False
    float_check = is_float(FWHM)
    while not float_check:
        FWHM = input('\nCould not convert input to float. Enter valid ' +
                     'value for FWHM or enter "Q" to quit: ')
        if str(FWHM).lower() == 'q':
            return None
        float_check = is_float(FWHM)

    FWHM = float(FWHM)

    star_table(t_file, FWHM)


main()
