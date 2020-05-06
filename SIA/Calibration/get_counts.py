from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.time import Time
from astropy import units as u
from astropy.wcs import WCS
import glob
from matplotlib.animation import FuncAnimation
from matplotlib import gridspec
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
from matplotlib import use
import numpy as np
import os
from photutils import aperture_photometry
import photutils.centroids as c
from photutils import DAOStarFinder
from photutils import SkyCircularAperture, SkyCircularAnnulus
import warnings

use('agg')
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.fontset':'stixsans'})
plt.rcParams.update({'axes.linewidth': 1.5})
plt.rcParams.update({'xtick.direction':'in'})
plt.rcParams.update({'xtick.major.size': 5})
plt.rcParams.update({'xtick.major.width': 1.25 })
plt.rcParams.update({'xtick.minor.size': 2.5})
plt.rcParams.update({'xtick.minor.width': 1.25 })
plt.rcParams.update({'ytick.direction':'in'})
plt.rcParams.update({'ytick.major.size': 5})
plt.rcParams.update({'ytick.major.width': 1.25 })
plt.rcParams.update({'ytick.minor.size': 2.5})
plt.rcParams.update({'ytick.minor.width': 1.25 })

def get_counts(dirtarget, ra, dec, fil, aper_rad, ann_in_rad, ann_out_rad, name,
               set_rad=False, centroid_plot=False):
    """Generates background-substracted aperture sums for a source in an image.

    Reads in all files from dirtarget and initializes arrays for various outputs
    based on the number of images in dirtarget. Then, for each image, the source
    located at ra, dec is checked to ensure that more than 10 stars were matched
    in the WCS calculation, that it is entirely in the image, and that it is not
    saturated. A 2D Gaussian fit is performed to centroid the aperture on the
    source by searching in a 1600 square pixel region surrounding the pixel
    corresponding to the RA and Dec. of the source. Then, aperture_photometry
    sums the counts within this aperture. The background level per pixel is
    calculated and used to subtract the background level from the aperture sum.

    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    ra : list
        List of string(s) of right ascension(s) of object(s) to be processed.
    dec : list
        List of string(s) of declination(s) of object(s) to be processed.
    fil : str
        Name of filter used for images which are currently being processed.
    set_rad : Boolean
        Determine whether user would like to use default aperture/annulus radii
        or specify their own.
    aper_rad : float
        User-specified aperture radius in arcseconds.
    ann_in_rad : float
        User-specified annulus inner radius in arcseconds.
    ann_out_rad : float
        User-specified annulus outer radius in arcseconds.
    name : str
        Type of star that get_counts is being ran on (e.g. target, comp, check).
    centroid_plot : Boolean
        Whether or not to plot centroid shifts for object.

    Returns
    -------
    total_sum : list
        List of numpy arrays containing aperture sum of each RA and dec in
        rightascension and declination
    err : numpy.ndarray
        Array of error values for each aperture sum for the last RA and dec
        in rightascension and declination.
    date_obs : numpy.ndarray
        Array of times each image was taken in Julian Days.
    alts_out : numpy.ndarray
        Array of object altitudes for each image.
    saturated : list
        List of string of image paths that have objects in them that are at the
        saturation value for the ISR images.
    exptimes : list
        List of floats corresponding to exposure time of a given image.
    cent_coord_list : numpy.ndarray
    init_coord_list : numpy.ndarray
    image_num : list
    pix_radius :
    image_arr :
    """
    dirtarget_wcs = os.path.join(dirtarget, 'ISR_Images', fil, 'WCS')
    files = sorted(glob.glob(os.path.join(dirtarget_wcs, '*.fits')))
    size_files = len(files)
    size_sources = len(ra)

    # Initialize output arrays.
    aper_sum = np.empty([size_sources, size_files])
    aper_sum[:][:] = np.nan
    err = np.empty([size_sources, size_files])
    err[:][:] = np.nan
    date_obs = np.empty([size_sources, size_files])
    date_obs[:][:] = np.nan
    altitudes = np.empty([size_sources, size_files])
    altitudes[:][:] = np.nan
    saturated = []
    exptimes = np.empty([size_sources, size_files])
    exptimes[:][:] = np.nan
    init_coords = np.empty([size_sources, size_files], dtype=object)
    init_coords[:][:] = 'temp'
    cent_coords = np.empty([size_sources, size_files], dtype=object)
    cent_coords[:][:] = 'temp'
    image_num = np.empty([size_sources, size_files])
    image_num[:][:] = np.nan
    pix_radius = np.empty(size_sources)
    pix_radius[:] = np.nan
    image_arr = np.empty([size_sources, size_files])
    image_arr[:][:] = np.nan

    sat_qual = np.empty([size_sources, size_files])
    sat_qual[:][:] = 0
    cent_qual = np.empty([size_sources, size_files])
    cent_qual[:][:] = 0

    for i, (ra_i, dec_i) in enumerate(zip(ra, dec)):
        print('\nProcessing {} star ({})...'.format(name, fil))

        im = None
        p1 = None
        p2 = None
        circ = None
        fig = None
        ax = None
        pix_radius = None

        init_x = np.empty(size_files)
        init_x[:] = np.nan
        init_y = np.empty(size_files)
        init_y[:] = np.nan
        cent_x = np.empty(size_files)
        cent_x[:] = np.nan
        cent_y = np.empty(size_files)
        cent_y[:] = np.nan

        saturated_i = []

        # Find nine evenly-spaced indices to plot a centroiding summary.
        cent_ind = np.linspace(0, size_files-1, 9).astype(int)
        if centroid_plot:
            fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(10,8))
            ax = ax.flatten()
            for k in range(0, 9):
                ax[k].set_title('Image {}'.format(cent_ind[k]), size=10)
                for l in (ax[k].get_xticklabels() + ax[k].get_yticklabels()):
                    l.set_fontsize(8)

        try:
            os.mkdir(os.path.join(dirtarget_wcs, 'output', name))
        except FileExistsError:
            pass

        # fig_anim, ax_anim = plt.subplots(nrows=1, ncols=1, figsize=(10,8))
        # im_anim = plt.imshow()
        # def animate(i):
        for j, item in enumerate(files):
            o_file = os.path.join(dirtarget_wcs, item)
            hdulist = fits.open(o_file)
            header = hdulist[0].header
            data = hdulist[0].data

            # Read output information from header of image being processed.
            exptimes[i][j] = float(header['EXPTIME'])
            date = header['DATE-OBS']
            t = Time(date)
            time = t.jd
            date_obs[i][j] = float(time)
            altitudes[i][j] = float(header['OBJCTALT'])

            # Check if WCS solution was successful.
            if header['WCSMATCH'] < 10:
                print('\nLess than 10 stars matched in WCS calculation.')
                cent_qual[i][j] = 1
                continue

            w = WCS(header)
            coords = SkyCoord(ra_i, dec_i, unit=(u.hourangle, u.deg))
            px_dec, py_dec = w.wcs_world2pix(coords.ra.deg, coords.dec.deg, 1)
            px, py = int(px_dec), int(py_dec)
            init_coords_pix = (px, py)
            init_coords[i][j] = ','.join(map(str, init_coords_pix))

            star = data[(py - 19):(py + 21), (px - 19):(px + 21)]

            # Check that region lies entirely within image.
            if ((py - 19) < 0) or ((py + 21) > 2084):
                print('\n{} star not entirely in the image.'.format(name))
                cent_qual[i][j] = 1
                continue
            if ((px - 19) < 0) or ((px + 21) > 3072):
                print('\n{} star not entirely in the image.'.format(name))
                cent_qual[i][j] = 1
                continue

            star_flat = star.reshape(1, len(star) * len(star[0]))

            # Check that source does not exceed the expected saturation level.
            max_pix = int(np.amax(star_flat))
            if max_pix >= header['SATLEVEL']:
                print('\nSaturation value: {}'.format(header['SATLEVEL']))
                print('\nMax aperture value: {}'.format(max_pix))
                saturated_i.append(item)
                sat_qual[i][j] = 1
                continue

            image_num[i][j] = j

            FWHM = 4
            try:
                mean, median, std = sigma_clipped_stats(star, sigma=3.0)
                daofind = DAOStarFinder(fwhm=FWHM, threshold=5.*std)
                sources = daofind(star - median)
                px_cent = np.average(sources['xcentroid'])
                py_cent = np.average(sources['ycentroid'])
            except TypeError:
                print('\nAperture centroiding failed for image {}.'.format(j))
                cent_qual[i][j] = 1
                continue

            x_cent = int(px - 19 + px_cent)
            y_cent = int(py - 19 + py_cent)

            cent_coords_pix = (x_cent, y_cent)
            cent_coords[i][j] = ','.join(map(str, cent_coords_pix))

            cent_equa = SkyCoord.from_pixel(x_cent, y_cent, w)

            # Define aperture and annulus radii.
            radius = None
            r_in = None
            r_out = None
            if set_rad:
                radius = float(aper_rad) * u.arcsec
                r_in = float(ann_in_rad) * u.arcsec
                r_out = float(ann_out_rad) * u.arcsec

            else:
                radius = 4 * u.arcsec
                r_in = 25 * u.arcsec
                r_out = 27 * u.arcsec

            # Create SkyCircularAperture and SkyCircularAnnulus objects
            # centered at the position of the star whose counts are being
            # summed.
            aperture = SkyCircularAperture(cent_equa, radius)
            annulus = SkyCircularAnnulus(cent_equa, r_in=r_in, r_out=r_out)

            apers = (aperture, annulus)

            secpix1 = abs(hdulist[0].header['SECPIX1'])

            # Determine the area of the aperture and annulus using the
            # arcseconds per pixel in the horizontal dimension header keyword.
            aper_area = np.pi * (radius / secpix1) ** 2
            area_out = np.pi * (r_out / secpix1) ** 2
            area_in = np.pi * (r_in / secpix1) ** 2
            annulus_area = area_out - area_in

            pix_radius = radius.value / secpix1

            # Call aperture_photometry function in order to sum all of the
            # counts in both the aperture and annulus for item.
            phot_table = aperture_photometry(hdulist, apers)

            # Remove the background level from the aperture sum.
            bkg_mean = phot_table['aperture_sum_1'] / annulus_area
            bkg_sum = bkg_mean * aper_area
            final_sum = phot_table['aperture_sum_0'] - bkg_sum
            phot_table['residual_aperture_sum'] = final_sum

            # Determine the error in the aperture sum and background level.
            # source_err = np.sqrt(phot_table['residual_aperture_sum'])
            source = phot_table['residual_aperture_sum']

            bkg_err = np.sqrt(bkg_sum)

            aper_sum[i][j] = phot_table['residual_aperture_sum'][0]
            err[i][j] = np.sqrt(source + bkg_sum)

            hdulist.close()
            del data

        saturated.append(saturated_i)

    return aper_sum, err, date_obs, altitudes, saturated, exptimes, \
           init_coords, cent_coords, image_num, sat_qual, cent_qual
