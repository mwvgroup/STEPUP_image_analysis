from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
import glob
import numpy as np
import os
from photutils import DAOStarFinder
import warnings
from photutils import CircularAperture
from astropy.visualization.mpl_normalize import ImageNormalize
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.fontset': 'stixsans'})
plt.rcParams.update({'axes.linewidth': 1.5})
plt.rcParams.update({'xtick.direction': 'in'})
plt.rcParams.update({'xtick.major.size': 5})
plt.rcParams.update({'xtick.major.width': 1.25})
plt.rcParams.update({'xtick.minor.size': 2.5})
plt.rcParams.update({'xtick.minor.width': 1.25})
plt.rcParams.update({'ytick.direction': 'in'})
plt.rcParams.update({'ytick.major.size': 5})
plt.rcParams.update({'ytick.major.width': 1.25})
plt.rcParams.update({'ytick.minor.size': 2.5})
plt.rcParams.update({'ytick.minor.width': 1.25})


def star_table(image, FWHM):
    """Create table of positions and brightnesses of all sources in an image.

    Parameters
    ----------
    image : str
        Absolute filename of image whose sources are being extracted.
    FWHM : float
        FWHM of image's PSF.

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

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r=15.)

    print('\n', sources)

    path_file = os.path.split(os.path.abspath(image))
    path = path_file[0]
    filename = path_file[1]
    filechars = np.array(list(filename))
    n = np.where(filechars == '.')[0][0]
    filename = filename[:n]

    fig = plt.figure(figsize=(10, 8))
    norm = ImageNormalize(stretch=SqrtStretch())
    plt.imshow(data, cmap='Greys', origin='lower', norm=norm)
    apertures.plot(color='blue', lw=1.5, alpha=0.5)
    plt.savefig(os.path.join(path,'{}_startable.pdf'.format(filename)))

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

    t_file = os.path.join(t_file)

    # Option to set FWHM using user input (in pixels).
    """
    FWHM = input('\nInput image FWHM (in pixels): ')

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
    """
    # Set FWHM using default value (in pixels).
    FWHM = 6

    star_table(t_file, FWHM)


main()
