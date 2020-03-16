import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy import units as u
from photutils import SkyCircularAperture, SkyCircularAnnulus
import photutils.centroids as c
from photutils import aperture_photometry
from astropy.time import Time
from matplotlib import gridspec
from matplotlib.patches import Circle
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


def perform_photometry(target, dirtarget, filters, date, coords, comp_ra,
                       comp_dec, comp_mags, clabel, cra, cdec, set_rad,
                       aper_rad, ann_in_rad, ann_out_rad):
    """Runs photometry part of image analysis routine.

    Calls photometry, which calls get_counts to get the aperture sums for the
    target, comparison stars, and check star as well as the error on
    the target aperture sum. get_counts also returns a list of all of the times
    and altitudes for each image. Then these aperture sums are converted in
    counts_to_mag from count values to magnitude values using the aperture sums
    of the comparison stars and their mangitudes. This process is repeated for
    the target and check aperture sums as well as the target error.
    Then, mag_plot is called, which plots the target magnitudes with error over
    time and check magnitudes over time. Then write_file is called to write an
    output text file that includes a summary of the analysis results.

    Parameters
    ----------
    target : str
        Name of target.
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    filters : list
        List containing string of each filter keyword found in header of flat
        field and light frame images.
    date : str
        Date of observation.
    coords : list
        List of list of string of target right ascension and declination.
    comp_ra : list
        List of strings of comparison stars' right ascension.
    comp_dec : list
        List of strings of comparison stars' declination.
    comp_mags : list
        List of floats representing the magnitudes of the comparison stars.
    clabel : str
        Name of check star.
    cra : list
        List of string of right ascension of check star.
    cdec : list
        List of string of delination of check star.
    set_rad : Boolean
        Determine whether user would like to use default aperture/annulus radii
        or specify their own.
    aper_rad : float
        User-specified aperture radius in arcseconds.
    ann_in_rad : float
        User-specified annulus inner radius in arcseconds.
    ann_out_rad : float
        User-specified annulus outer radius in arcseconds.

    Returns
    -------
    None
    """
    for fil in filters:
        os.chdir(os.path.join(dirtarget, 'ISR_Images', fil, 'WCS'))

        try:
            os.mkdir('output')
        except FileExistsError:
            pass

        aper_sum, comp_aper_sums, check_aper_sum, err, date_obs, altitudes, final_comp_mags, saturated, exposure_times = photometry(dirtarget, fil, coords, comp_ra, comp_dec, cra, cdec, comp_mags, set_rad, aper_rad, ann_in_rad, ann_out_rad)

        write_net_counts(dirtarget, fil, date, comp_aper_sums, aper_sum, check_aper_sum, err, date_obs, altitudes, target, clabel)

        target_mags, target_err, check_mags = counts_to_mag(aper_sum, comp_aper_sums, err, final_comp_mags, check_aper_sum, fil, date_obs)

        mag_plot(target_mags, target_err, date_obs, target, date, fil, dirtarget, check_mags)

        write_file(target_mags, target_err, date_obs, target, dirtarget, fil, altitudes, clabel, check_mags, date)

        path = os.path.join(dirtarget, 'ISR_Images', fil, 'WCS',
                            'output/saturated_{}.txt'.format(fil))

        with open(path, 'w+') as f:
            f.write('#Files in which the target star met or exceeded the saturation level:\n')
            for o_file in saturated:
                f.write(str(o_file) + '\n')
            f.close()

        print('Saturated images ({}): {}'.format(fil, saturated))
        print('Exposure times ({}): {}'.format(fil, exposure_times))

    multi_filter_analysis(dirtarget, date, target, filters)


def photometry(dirtarget, fil, coords, comp_ra, comp_dec, cra, cdec,
               comp_mags, set_rad, aper_rad, ann_in_rad, ann_out_rad):
    """Gets aperture sums for target, comparison, and check stars.

    Calls get_counts for the list of right ascension(s) and declination(s) for
    the target star, each comparison star, and the check star.
    The err, date_obs, and altitudes are defined when get_counts is
    called for the target star. Then, if any of the comparison stars are not
    image (i.e. they have either nan or negative values in their aper_sum
    in the arrays), then they are removed from the comp_aper_sum array and the
    corresponding magnitude is removed from comp_mags. The same process is
    repreated for the check and star, except there is no corresponding
    magnitude to delete.

    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    fil : list
        List containing string of each filter keyword found in header of flat
        field and light frame images.
    coords : list
        List of list of string of target right ascension and declination.
    comp_ra : list
        List of strings of comparison stars' right ascension.
    comp_dec : list
        List of strings of comparison stars' right ascension.
    cra : list
        List of string of right ascension of check star.
    cdec : list
        List of string of delination of check star.
    comp_mags : list
        List of floats representing the magnitudes of the comparison stars.
    set_rad : Boolean
        Determine whether user would like to use default aperture/annulus radii
        or specify their own.
    aper_rad : float
        User-specified aperture radius in arcseconds.
    ann_in_rad : float
        User-specified annulus inner radius in arcseconds.
    ann_out_rad : float
        User-specified annulus outer radius in arcseconds.

    Returns
    -------
    aper_sum : numpy.ndarray
        Array of aperture sums for target star in counts.
    comp_apers : numpy.ndarray
        Filtered array of comparison star magnitudes.
    check_apers : numpy.ndarray
        Filtered array of aperture sum for check star in counts.
    err : numpy.ndarray
        Array of error values for each aperture sum of the target star in
        counts.
    date_obs : np.ndarray
        Array of times each image was taken in Julian Days.
    altitudes : np.ndarray
        Array of object altitudes for each image.
    comp_mags : numpy.ndarray
        Filtered array of comparison star magnitudes.
    saturated : list
        List of string of image paths that have objects in them that are at the
        saturation value for the ISR images.
    exposure_times : list
        List of floats corresponding to exposure time of a given image.
    """
    # Get aperture sum, error of aperture sum, times of data collection,
    # and altitudes for target.
    aper_sum, err, date_obs, altitudes, saturated, exposure_times = get_counts(dirtarget, coords[0], coords[1], fil, set_rad, aper_rad, ann_in_rad, ann_out_rad, "target", True)
    aper_sum = aper_sum[0]
    aper_sum = np.array(aper_sum, dtype=float)

    # Get aperture sums for each somparison star.
    comp_apers = (get_counts(dirtarget, comp_ra, comp_dec, fil, set_rad, aper_rad, ann_in_rad, ann_out_rad, "comp", False))[0]
    comp_apers = np.array(comp_apers, dtype=float)

    # Get aperture sum of the check star.
    check_aper_sum = (get_counts(dirtarget, cra, cdec, fil, set_rad, aper_rad, ann_in_rad, ann_out_rad, "check", True))[0]
    check_aper_sum = np.array(check_aper_sum, dtype=float)
    check_aper_sum = check_aper_sum[0]

    # To determine if any of the comparison stars are not in the image, the
    # comp_apers arrays are individually checked for NaN and negative values.
    # If an array contains a negative value, the the total number of images
    # containing negative numbers is considered before ruling out the star
    # completely.
    good_comp = []
    for i, o_aper in enumerate(comp_apers):
        no_nan = (len(np.argwhere(np.isnan(o_aper))) == 0)
        no_neg = (len(np.where(o_aper <= 0)[0]) == 0)
        if not no_neg:
            where_neg = np.where(o_aper < 0)[0]
            if len(where_neg) < (0.25 * len(o_aper)):
                no_neg = True
        if no_nan and no_neg:
            good_comp.append(i)
    comp_apers = np.array([comp_apers[i] for i in good_comp])
    comp_mags = np.array([comp_mags[i] for i in good_comp])

    # Remove any bad images from analysis.
    good_im = []
    good_im.extend(np.where(aper_sum != np.nan)[0])
    good_im.extend(np.where(aper_sum > 0)[0])
    good_im = np.unique(good_im)
    aper_sum = np.array([aper_sum[i] for i in good_im])
    err = np.array([err[i] for i in good_im])
    date_obs = np.array([date_obs[i] for i in good_im])
    altitudes = np.array([altitudes[i] for i in good_im])
    exposure_times = np.array([exposure_times[i] for i in good_im])
    comp_apers_n = []
    for obj in comp_apers:
        comp_apers_n.append([obj[i] for i in good_im])
    comp_apers = np.array(comp_apers_n)
    check_apers = np.array([check_aper_sum[i] for i in good_im])

    return aper_sum, comp_apers, check_apers, err, date_obs, altitudes, comp_mags, saturated, exposure_times


def write_net_counts(dirtarget, fil, date, comp_aper_sums, aper_sum,
                     check_aper_sum, t_err, date_obs, altitudes, target,
                     clabel):
    """*Incomplete* function to save file with net count values before
    converting to magnitudes.

    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    fil : str
        Name of filter used for images which are currently being processed.
    date : str
        Date of observation.
    comp_aper_sums : numpy.ndarray
        Filtered array of comparison star magnitudes.
    aper_sum : numpy.ndarray
        Array of aperture sums for target star in counts.
    check_aper_sum : numpy.ndarray
        Filtered array of aperture sum for check star in counts.
    t_err : numpy.ndarray
        Array of error values for each aperture sum of the target star in
        counts.
    date_obs : np.ndarray
        Array of times each image was taken in Julian Days.
    altitudes : np.ndarray
        Array of object altitudes for each image.
    target : str
        Name of target.
    clabel : str
        Name of check star.

    Returns
    -------
    None
    """
    path = os.path.join(dirtarget, 'ISR_Images', fil, 'WCS',
                        'output/net_counts_{}.txt'.format(date))

    with open(path, 'w+') as f:
        f.write('#SOFTWARE=STEPUP Image Analysis\n#DELIM=,\n#DATE=JD\n' +
                '#OBSTYPE=CCD\n')
        comp_n = len(comp_aper_sums)
        header_str = str('#TARGET NAME,DATE,TARGET COUNTS,ERR,FILTER,(C1,...,C{})'+
                         ' COUNTS,CHECK LABEL,CHECK COUNTS,' +
                         'AIRMASS\n'.format(comp_n))
        f.write(header_str)
        comp_sums = list(zip(*comp_aper_sums))
        for n, (date_i, tsum, err, csum, alt) in enumerate(zip(date_obs,
                                                               aper_sum,
                                                               t_err,
                                                               check_aper_sum,
                                                               altitudes)):
            csums = comp_sums[n]
            zenith = np.deg2rad(90 - alt)
            airmass = 1 / np.cos(zenith)
            input_list = [target, date_i, tsum, err, fil, csums, clabel, csum,
                          airmass]
            input_string = ",".join(map(str, input_list))
            f.write(input_string + '\n')
        f.close()


def counts_to_mag(aper_sum, comp_aper_sums, err, comp_mags, check_aper_sum,
                  fil, date_obs):
    """Scales the aperture sums from counts to magnitudes.

    Parameters
    ----------
    aper_sum : numpy.ndarray
        Array of aperture sums for target star in counts.
    comp_aper_sums : numpy.ndarray
        Filtered array of comparison star magnitudes.
    err : numpy.ndarray
        Array of error values for each aperture sum of the target star in
        counts.
    comp_mags : list
        List of floats representing the magnitudes of the comparison stars.
    check_aper_sum : numpy.ndarray
        Filtered array of aperture sum for check star in counts.
    fil : str
        Name of filter used for images which are currently being processed.
    date_obs : numpy.ndarray
        Array of times each image was taken in Julian Days.

    Returns
    -------
    target_mags : numpy.ndarray
        Array of magnitude values for the target.
    target_err : numpy.ndarray
        Array of error values for target magnitudes.
    check_mags_f : numpy.ndarray
        Array of magnitude values for the check star.
    """
    # Initialize arrays for target magnitudes, check magnitudes,
    # and errors in target magnitudes.
    scaled_mags = np.empty(comp_aper_sums.shape)
    scaled_mags[:] = np.nan
    check_mags = np.empty(comp_aper_sums.shape)
    check_mags[:] = np.nan

    for i, (obj, mag) in enumerate(zip(comp_aper_sums, comp_mags)):
        figure = plt.figure(figsize=(10,8))
        plt.plot(date_obs, 2.5 * np.log10(obj), "o", c='cadetblue', label='{}'.format(fil))
        plt.ylabel("Instrumental Magnitude")
        plt.xlabel("Time [JD]")
        plt.gca().invert_yaxis()
        plt.title("Comp {}, {} = {} mag".format(i + 1, fil, mag))
        plt.savefig("output/comp_{}_{}.pdf".format(i + 1, fil))
        plt.legend()

    for i, (mag, obj) in enumerate(zip(comp_mags, comp_aper_sums)):
        # Using magnitude value of comparison star (mag) and aperture sum
        # of comparison star (obj), each image's target count value
        # (aper_sum) is determined.
        scaled_mags[i] = mag - 2.5 * np.log10(aper_sum / obj)

        # If the check star is in the image, the magnitudes of each star are
        # determined for each image.
        if np.all(check_aper_sum != None):
            check_mags[i] = mag - 2.5 * np.log10(check_aper_sum / obj)

    target_err = (2.5 * np.log10((aper_sum + err) / aper_sum))

    # For each image, the scaled magnitude value for each comparison star is
    # averaged.
    target_mags = np.average(scaled_mags, axis=0)

    # If the check star is in the image, the calculated magnitudes are averaged
    # for each image.
    check_mags_f = None
    if np.all(check_aper_sum != None):
        check_mags = np.array(check_mags)
        check_mags_f = np.average(check_mags, axis=0)
    else:
        check_mags_f = 0

    print('\nTarget mags: ', target_mags)
    print('\nTarget error: ', target_err)
    print('\nCheck mags: ', check_mags_f)

    return target_mags, target_err, check_mags_f


def mag_plot(target_mags, target_err, date_obs, target, date, fil, dirtarget,
             check_mags):
    """Plots magnitudes of target star with error and check star over time.

    Creates a subplot with the target magntidues and their error over time in
    Julian Days on top and the check star magnitudes on the bottom with the
    same x-axis as the target star. A PDF of the lightcurve is written to the
    location of the ISR, plate-solved dataset.

    Parameters
    ----------
    target_mags : numpy.ndarray
        Array of magnitude values for the target.
    target_err : numpy.ndarray
        Array of error values for target magnitudes.
    date_obs : numpy.ndarray
        Array of times each image was taken in Julian Days.
    target : numpy.ndarray
        Name of target.
    date : str
        Date of observation.
    fil : str
        Name of filter used for images which are currently being processed.
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    check_mags : numpy.ndarray
        Array of magnitude values for the check star.

    Returns
    -------
    None
    """
    f, axarr = plt.subplots(2, sharex=True,
                            gridspec_kw={'height_ratios': [3, 1]},
                            figsize=(10, 8))
    axarr[0].errorbar(date_obs, target_mags, yerr=target_err, fmt='o',
                      c='cadetblue', label='{}'.format(fil))
    axarr[0].set_title('Light Curve of {}, {}'.format(target, date))
    axarr[0].set_ylabel('Magnitude')
    axarr[0].invert_yaxis
    axarr[0].legend()
    axarr[0].set_ylim(axarr[0].get_ylim()[::-1])
    axarr[1].scatter(date_obs, check_mags, c='cadetblue',
                     label='{}'.format(fil))
    axarr[1].set_ylim(axarr[1].get_ylim()[::-1])
    axarr[1].invert_yaxis
    axarr[1].legend()
    axarr[1].set_title('Check Star')
    axarr[1].set_ylabel('Magnitude')
    axarr[1].set_xlabel('Time [JD]')
    f.savefig(os.path.join(dirtarget, 'ISR_Images', fil, 'WCS',
                           'output/lightcurve_{}_{}.pdf'.format(date, fil)))


def write_file(target_mags, target_err, date_obs, target, dirtarget, fil,
               altitudes, clabel, cmags, date):
    """Writes file of observational data to submit to AAVSO.

    Formats data to be compatible for submission for the AAVSO Extended File
    Format, which is saved to the directory containing FITS files with accurate
    WCS information.

    Parameters
    ----------
    target_mags : numpy.ndarray
        Array of magnitude values of target star for each image.
    target_err : numpy.ndarray
        Array of error values of magntidue of target star for each image.
    date_obs : numpy.ndarray
        Array of time of observation in Julian Days.
    target : str
        Name of target.
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    fil : str
        Name of filter used for images which are currently being processed.
    altitudes : numpy.ndarray
        Array of altitudes of each image.
    clabel : str
        Name of check star.
    cmags : numpy.ndarray
        Array of count values of aperture sum of comparison star of closest
        location to the target star.
    date : str
        Date of observation.

    Returns
    -------
    None
    """
    path = os.path.join(dirtarget, 'ISR_Images', fil, 'WCS',
                        'output/output_{}_{}.txt'.format(date, fil))

    with open(path, 'w+') as f:
        f.write('#SOFTWARE=STEPUP Image Analysis\n#DELIM=,\n#DATE=JD\n' +
                '#OBSTYPE=CCD\n')
        f.write('TARGET,DATE,TARGET MAG,ERROR,FILTER,CHECK LABEL,CHECK MAG,' +
                'AIRMASS\n')
        for date_i, mag, err, cmag, alt in zip(date_obs, target_mags,
                                               target_err, cmags, altitudes):
            zenith = np.deg2rad(90 - alt)
            airmass = 1 / np.cos(zenith)
            input_list = [target, date_i, mag, err, fil, clabel, cmag, airmass]
            input_string = ",".join(map(str, input_list))
            f.write(input_string + '\n')
        f.close()


def get_counts(dirtarget, rightascension, declination, fil, set_rad, aper_rad,
               ann_in_rad, ann_out_rad, name, centroid_plt):
    """Determines count values for star(s).

    Reads in all files from dirtarget and initializes arrays for error,
    observation times, and altitudes corresponding to the number of images
    that will be processed as well as a list for the arrays of aperture sums
    that will be generated for each object get_counts has been called for.
    Then, for each object corresponding to each ra and dec pair in
    rightascension and declination, if WCS Tools has matches at least 20 stars
    in the image, a SkyCoord and SkyCircularAnnulus object are craeted, both
    centered at the object's location. The area of the aperture and annulus
    are calculated and phot_table is called in order to get an aperture and
    annulus sum for the object. The source error and background error are
    calculated in order to determine to aperture sum error and the aperture sum
    and its error are added to their corresponding arrays, as well as the time
    and altitude. The aperture is appended to total_sum and this process is
    repeated for each object that get_counts has been called for.

    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    rightascension : list
        List of string(s) of right ascension(s) of object(s) to be processed.
    declination : list
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
    centroid_plt : Boolean
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
    altitudes : numpy.ndarray
        Array of object altitudes for each image.
    saturated : list
        List of string of image paths that have objects in them that are at the
        saturation value for the ISR images.
    exposure_times : list
        List of floats corresponding to exposure time of a given image.
    """
    dirtarget_wcs = os.path.join(dirtarget, 'ISR_Images', fil, 'WCS')
    size = len(glob.glob(os.path.join(dirtarget_wcs, '*.fits')))
    err = np.empty(size)
    err[:] = np.nan
    date_obs = np.empty(size)
    date_obs[:] = np.nan
    altitudes = np.empty(size)
    altitudes[:] = np.nan
    total_sum = []
    saturated = []
    exposure_times = []
    raw_loc = []
    centroid_corrections = []
    image_stack = []
    im = None
    p1 = None
    p2 = None
    circ = None
    figure = None
    axis = None

    for ra, dec in zip(rightascension, declination):
        cent_ind = np.linspace(0, size-1, 9).astype(int)
        if centroid_plt:
            figure, axis = plt.subplots(nrows=3, ncols=3, figsize=(10,8))
            axis = axis.flatten()
        aper_sum = np.empty(size)
        aper_sum[:] = np.nan
        for i, item in enumerate(sorted(glob.glob(os.path.join(dirtarget_wcs,
                                                               '*.fits')))):
            o_file = os.path.join(dirtarget_wcs, item)
            hdulist = fits.open(o_file)
            if hdulist[0].header['WCSMATCH'] < 10:
                print('Less than 10 stars matched in WCS calculation.\n')
                continue
            # Determine the exposure time for item.
            exposure_times.append(hdulist[0].header['EXPTIME'])

            # Determine what physical position corresponds to the right
            # ascension and declination that are defined as SkyCoord objects
            # with units of hourangles and degrees.
            header = fits.getheader(o_file)
            w = WCS(header)
            coords = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            px, py = w.wcs_world2pix(coords.ra.deg, coords.dec.deg, 1)
            px_int = int(px)
            py_int = int(py)

            # Get the image array.
            image_array = fits.getdata(o_file)

            # Read in a 40 x 40 square centered at the star for which the
            # counts are being summed.
            star = image_array[(py_int - 19):(py_int + 21),
                               (px_int - 19):(px_int + 21)]

            # Ensure that the square is entirely on the image.
            if ((py_int - 19) < 0) or ((py_int + 21) > 2084):
                continue
            if ((px_int - 19) < 0) or ((px_int + 21) > 3072):
                continue

            # Flatten the array so that it is one-dimensional.
            star_flat = star.reshape(1, len(star) * len(star[0]))

            # Determine if the image is saturated at the star's position using
            # the expected saturation level. If saturated, the loop will move
            # on to the next image.
            print('Saturation value: {}'.format(hdulist[0].header['SATLEVEL']))
            print('Max aperture value: {}'.format(int(np.amax(star_flat))))
            if np.amax(star_flat) >= hdulist[0].header['SATLEVEL']:
                saturated.append(item)
                continue

            arr_x_centroid, arr_y_centroid = c.centroid_2dg(star)
            pix_x_centroid = int((px_int - 19) + arr_x_centroid)
            pix_y_centroid = int((py_int - 19) + arr_y_centroid)
            pix_centroid_coords = (pix_x_centroid, pix_y_centroid)
            print('\nPlate-solution (x, y) object location: ({}, {})'.format(px_int, py_int))
            print('\nCentroided (x, y) object location: ({}, {})'.format(pix_x_centroid, pix_y_centroid))
            centroid_coords = SkyCoord.from_pixel(pix_x_centroid,
                                                  pix_y_centroid, w)

            print('\nCentroided R.A. and declination: {}'.format(centroid_coords))

            i_plot = np.where(cent_ind == i)[0]

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
            aperture = SkyCircularAperture(centroid_coords, radius)
            annulus = SkyCircularAnnulus(centroid_coords, r_in=r_in,
                                         r_out=r_out)

            secpix1 = abs(hdulist[0].header['SECPIX1'])

            # Determine the area of the aperture and annulus using the
            # arcseconds per pixel in the horizontal dimension header keyword.
            aper_area = np.pi * (radius / secpix1) ** 2
            area_out = np.pi * (r_out / secpix1) ** 2
            area_in = np.pi * (r_in / secpix1) ** 2
            annulus_area = area_out - area_in

            if centroid_plt:
                if i in cent_ind:
                    p1 = axis[i_plot[0]].scatter([px_int+1], [py_int+1],
                              c="thistle", label="Original",
                              edgecolors='mistyrose')
                    p2 = axis[i_plot[0]].scatter([pix_x_centroid+1],
                              [pix_y_centroid+1], c="rebeccapurple",
                              label="Corrected", edgecolors='mistyrose')
                    im = axis[i_plot[0]].imshow(star, extent=(px_int - 19,
                                                              px_int + 21,
                                                              py_int - 19,
                                                              py_int + 21),
                                                cmap='magma', origin='lower')
                    axis[i_plot[0]].set_title('Image {}'.format(i), size=10)
                    circ = Circle((pix_x_centroid+1, pix_y_centroid+1),
                                  radius.value/secpix1, fill=False,
                                  label='Aperture', ls='-', color='mistyrose')
                    axis[i_plot[0]].add_patch(circ)
                    for label in (axis[i_plot[0]].get_xticklabels() +
                                  axis[i_plot[0]].get_yticklabels()):
                        label.set_fontsize(8)

            apers = (aperture, annulus)

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

            aper_sum[i] = phot_table['residual_aperture_sum'][0]
            err[i] = np.sqrt(source + bkg_sum)

            # Determine the time at which the image was taken in Julian Days
            # and the altitude and write them to an array.
            date = hdulist[0].header['DATE-OBS']
            t = Time(date)
            time = t.jd

            date_obs[i] = (time)
            altitudes[i] = (hdulist[0].header['OBJCTALT'])

            hdulist.close()

        total_sum.append(aper_sum)

        if centroid_plt:
            figure.subplots_adjust(right=0.8)
            cbar_ax = figure.add_axes([0.85, 0.15, 0.05, 0.7])
            cb = figure.colorbar(im, cax=cbar_ax)
            plt.setp(cb.ax.get_yticklabels(), fontsize=8)

            plt.figlegend([p1, p2, circ], ['Original', 'Corrected', 'Aperture'],
                          fontsize=14)
            figure.suptitle('Aperture Centroiding on {}, {}'.format(ra, dec),
                            fontsize=16)
            figure.text(0.5, 0.04, 'x [pixel]', ha='center')
            figure.text(0.04, 0.5, 'y [pixel]', va='center',
                        rotation='vertical')
            plt.savefig(os.path.join(dirtarget, 'ISR_Images', fil, 'WCS',
                                     'output/centroid_{}.pdf'.format(name)))
        else:
            continue

    return total_sum, err, date_obs, altitudes, saturated, exposure_times


def multi_filter_analysis(dirtarget, date, target, filters):
    """Produces color light curve for multi-filter observations.

    If the observation was done in more than one filter, the user is prompted
    for whether or not they would like to plot a color light curve. If they do,
    they are then prompted for which color they would like to plot by asking for
    the filters they would like to subtract. Input is comma-delimited, with
    order corresponding to order of subtraction, e.g. if the input is "B,V" the
    resultant plot will be B-V color over time. The outputted lightcurve is
    saved to dirtarget. Assumes user would only like to calculate one color
    light curve and that each color observation has the same amount of data
    points.

    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    date : str
        Date of observation.
    target : str
        Name of target.
    filters : list

    Returns
    -------
    None
    """
    # Sets answer to 'N' if only one filter was used for observation so that the
    # function does not plot a color light curve and returns None.
    answer = 'n'
    if len(filters) > 1:
        answer = input('\nWould you like to plot a color light curve? (Y/N): ')
    if answer.lower() == 'y':
        # Determine color user would like to plot in light curve.
        fil_color = input('\nInput magnitude filters to subtract (comma-delimited): ').split(',')
        dates = []
        mags = []
        err = []
        for n, fil in enumerate(filters):
            date_fil = []
            mag_fil = []
            err_fil = []
            os.chdir(os.path.join(dirtarget, 'ISR_Images', fil, 'WCS'))
            # Read in magnitudes of target in fil from output file generated by
            # SIA for analysis corresponding to that filter.
            with open('output_{}_{}.txt'.format(date, fil), "r") as f:
                for line in f:
                    if line.startswith(target):
                        line = line.split(',')
                        # Saves time and magnitude values for given filter.
                        date_fil.append(float(line[1]))
                        mag_fil.append(float(line[2]))
                        err_fil.append(float(line[3]))
            dates.append(date_fil)
            mags.append(mag_fil)
            err.append(err_fil)

        # Creates individual arrays for each filter of magnitudes and takes the
        # time of each color to be the average of the observation times of each
        # filter for a given point.
        mags1 = np.array(mags[0], dtype='float')
        mags2 = np.array(mags[1], dtype='float')
        err1 = np.array(err[0], dtype='float')
        err2 = np.array(err[1], dtype='float')
        dates = np.average(dates, axis=0)

        # Calculates color of target for each observation time.
        colors = mags1 - mags2
        color_err = np.sqrt(err1 ** 2 + err2 ** 2)

        # Plots and saves color light curve.
        fig = plt.figure(figsize=(10,8))
        plt.errorbar(dates, colors, yerr=color_err, fmt='o', c='cadetblue')
        plt.title('Color Light Curve of {}, {}'.format(target, date))
        plt.xlabel('Time [JD]')
        plt.ylabel('{}-{} Color'.format(filters[0], filters[1]))
        plt.gca().invert_yaxis()
        plt.savefig(os.path.join(dirtarget,
                                 'color_lightcurve_{}.pdf'.format(date)))

    else:
        pass
