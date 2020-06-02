import glob
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
from matplotlib import use
import numpy as np
import os
import sys
import warnings

sys.path.insert(0, 'Calibration')
from get_counts import get_counts

use('agg')
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


def perform_photometry(target, dirtarget, filters, date, coords, comp_ra,
                       comp_dec, comp_mags, clabel, cra, cdec, aper_rad,
                       ann_in_rad, ann_out_rad, set_rad=False):
    """Extracts light curves from image data using differential photometry.

    Calls photometry, which calls get_counts to get the aperture sums for the
    target, comparison stars, and check star. Then these aperture sums are
    converted in ounts_to_mag from net count values to magnitudes for the
    target and check star.

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
    aper_rad : float
        User-specified aperture radius in arcseconds.
    ann_in_rad : float
        User-specified annulus inner radius in arcseconds.
    ann_out_rad : float
        User-specified annulus outer radius in arcseconds.
    set_rad : Boolean
        Determine whether user would like to use default aperture/annulus radii
        or specify their own.

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

        aper_sum, comp_aper_sums, check_aper_sum, err, check_err, date_obs, \
            altitudes, final_comp_mags, saturated, exposure_times, \
            centroid_coords, init_coord_list, \
            image_num = photometry(dirtarget, fil, coords, comp_ra, comp_dec,
                                   cra, cdec, comp_mags, aper_rad, ann_in_rad,
                                   ann_out_rad, date, set_rad)

        write_net_counts(dirtarget, fil, date, comp_aper_sums, aper_sum,
                         check_aper_sum, err, check_err, date_obs, altitudes,
                         target, clabel)

        target_mags, target_err, check_mags, check_err = \
            counts_to_mag(aper_sum, comp_aper_sums, err, check_err,
                          final_comp_mags, check_aper_sum, fil, date_obs,
                          date)

        mag_plot(target_mags, target_err, date_obs, target, date, fil,
                 dirtarget, check_mags, check_err)

        write_file(target_mags, target_err, date_obs, target, dirtarget, fil,
                   altitudes, clabel, check_mags, check_err, date, image_num)

        # Write output file to summarize target saturation levels.
        path = os.path.join(dirtarget, 'ISR_Images', fil, 'WCS',
                            'output/saturated_{}_{}.txt'.format(date, fil))

        with open(path, 'w+') as f:
            f.write('#Files in which the target star met or exceeded the ' +
                    'saturation level:\n')
            for o_file in saturated:
                f.write(str(o_file) + '\n')
            f.close()

        # Write output file to summarize target pixel locations before and
        # after centroiding for non-saturated images.
        path = os.path.join(dirtarget, 'ISR_Images', fil, 'WCS',
                            'output/centroid_{}_{}.txt'.format(date, fil))
        with open(path, 'w+') as f:
            f.write('Summary of target pixel locations before and after ' +
                    'centroiding for non-saturated images.\n')
            f.write('#IMAGE-NUM,#RAW-XPIX,#RAW-YPIX,#CENT-XPIX,#CENT-YPIX,' +
                    '#DATE[JD]\n')
            for im, coord_i, date_i, init_coord_i in zip(image_num,
                                                         centroid_coords,
                                                         date_obs,
                                                         init_coord_list):
                input_list = [im, init_coord_i, coord_i, date_i]
                input_string = ",".join(map(str, input_list))
                f.write(input_string + '\n')
            f.close()

        print('\nSaturated images ({}): {}'.format(fil, saturated))
        print('\nExposure times ({}): {}'.format(fil, exposure_times))

    multi_filter_analysis(dirtarget, date, target, filters)


def photometry(dirtarget, fil, coords, comp_ra, comp_dec, cra, cdec, comp_mags,
               aper_rad, ann_in_rad, ann_out_rad, date, set_rad=False):
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
    aper_rad : float
        User-specified aperture radius in arcseconds.
    ann_in_rad : float
        User-specified annulus inner radius in arcseconds.
    ann_out_rad : float
        User-specified annulus outer radius in arcseconds.
    date : str
        Date of observation.
    set_rad : Boolean
        Determine whether user would like to use default aperture/annulus radii
        or specify their own.

    Returns
    -------
    aper_sum : numpy.ndarray
        Array of floats corresponding to aperture sums of target star.
    comp_apers_return : numpy.ndarray
        Array of arrays of aperture sum floats for each comparison star.
    check_apers : numpy.ndarray
        Array of floats corresponding to aperture sums of check star.
    err : numpy.ndarray
        Array of floats corresponding to uncertainty of each aperture sum.
    date_obs : numpy.ndarray
        Array of floats corresponding to Julian Date of each image.
    altitudes : numpy.ndarray
        Array of floats corresponding to the altitude of each image.
    comp_mags : numpy.ndarray
        Array of floats of comparison star magnitudes.
    saturated : list
        List of strings containing the file path of any image whose source
        meets or exceeds the expected saturation level.
    exposure_times : numpy.ndarray
        Array of floats corresponding to exposure time of each image.
    centroid_coords : numpy.ndarray
        Array of strings corresponding to the central pixel coordinates (x,y)
        of the centroided aperture. If the centroid routine failed, the string
        'init' is returned for that image.
    init_coords : numpy.ndarray
        Array of strings corresponding to the pixel coordinates (x,y) of the
        R.A. and dec. of the image's source, according to its WCS solution.
    image_num : numpy.ndarray
        Array of integers containing the number of each image in dirtarget.
    """
    # Get aperture sum, error of aperture sum, times of data collection,
    # and altitudes for target.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        aper_sum, err, date_obs, altitudes, saturated, exposure_times, \
            init_coord_list, centroid_coords, image_num, sat_qual, cent_qual \
            = get_counts(dirtarget, coords[0], coords[1], fil, aper_rad,
                         ann_in_rad, ann_out_rad, "target", date, set_rad,
                         True)
    aper_sum = aper_sum[0]
    err = err[0]
    date_obs = date_obs[0]
    altitudes = altitudes[0]
    saturated = saturated[0]
    exposure_times = exposure_times[0]
    centroid_coords = centroid_coords[0]
    init_coord_list = init_coord_list[0]
    image_num = image_num[0]
    sat_qual = sat_qual[0]
    cent_qual = cent_qual[0]
    aper_sum = np.array(aper_sum, dtype=float)

    # Get aperture sums for each somparison star.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        comp_apers, comp_err, comp_date_obs, comp_altitudes, comp_saturated, \
            comp_exposure_times, comp_init_coord_list, comp_centroid_coords, \
            comp_image_num, comp_sat_qual, comp_cent_qual = \
            get_counts(dirtarget, comp_ra, comp_dec, fil, aper_rad, ann_in_rad,
                       ann_out_rad, "comp", date, set_rad, True)
    comp_apers = np.array(comp_apers, dtype=float)

    # Get aperture sum of the check star.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        check_apers, check_err, check_date_obs, check_altitudes, \
            check_saturated, check_exposure_times, check_init_coord_list, \
            check_centroid_coords, check_image_num, check_sat_qual, \
            check_cent_qual = get_counts(dirtarget, cra, cdec, fil, aper_rad,
                                         ann_in_rad, ann_out_rad, "check",
                                         date, set_rad, True)
    check_err = check_err[0]
    check_date_obs = check_date_obs[0]
    check_altitudes = check_altitudes[0]
    check_saturated = check_saturated[0]
    check_exposure_times = check_exposure_times[0]
    check_centroid_coords = check_centroid_coords[0]
    check_init_coord_list = check_init_coord_list[0]
    check_image_num = check_image_num[0]
    check_sat_qual = check_sat_qual[0]
    check_cent_qual = check_cent_qual[0]
    check_apers = np.array(check_apers, dtype=float)[0]

    comp_sat_qual_tot = np.sum(comp_sat_qual, axis=0)
    comp_cent_qual_tot = np.sum(comp_cent_qual, axis=0)

    data_qual = np.sum([sat_qual, cent_qual, check_sat_qual, check_cent_qual,
                        comp_sat_qual_tot, comp_cent_qual_tot], axis=0)

    good_im = np.where(data_qual == 0)[0]

    bool_im = (data_qual == 0)

    path = os.path.join(dirtarget, 'ISR_Images', fil, 'WCS',
                        'output/data_quality_{}_{}.txt'.format(date, fil))

    with open(path, 'w+') as f:
        comp_n = len(comp_sat_qual)
        key_str_0 = '#0=passed data quality check.\n'
        key_str_1 = '#1=failed data quality check.\n'
        header_str = str('#IMAGE NUMBER,IMAGE USED,TARGET [SAT, CENT],CHECK ' +
                         '[SAT, CENT],C1,...,C{} [SAT, CENT]\n'.format(comp_n))
        f.write(key_str_0)
        f.write(key_str_1)
        f.write(header_str)
        comp_sat_qual_w = comp_sat_qual.astype(str)
        comp_sat_qual_w = list(zip(*comp_sat_qual_w))
        comp_cent_qual_w = comp_cent_qual.astype(str)
        comp_cent_qual_w = list(zip(*comp_cent_qual_w))
        for num, bool, t_sat, t_cent, c_sat, c_cent, comp_sat, comp_cent in \
            zip(image_num, bool_im, sat_qual, cent_qual, check_sat_qual,
                check_cent_qual, comp_sat_qual_w, comp_cent_qual_w):
            targ_qual = [t_sat, t_cent]
            check_qual = [c_sat, c_cent]
            comp_qual = [comp_sat, comp_cent]
            input_list = [num, bool, targ_qual, check_qual, comp_qual]
            input_string = ",".join(map(str, input_list))
            f.write(input_string + '\n')
        f.close()

    aper_sum = aper_sum[good_im]
    comp_apers_return = []
    for ap in comp_apers:
        comp_apers_return.append(ap[good_im])
    comp_apers_return = np.array(comp_apers_return, dtype=float)
    check_apers = check_apers[good_im]
    check_err = check_err[good_im]
    err = err[good_im]
    date_obs = date_obs[good_im]
    altitudes = altitudes[good_im]
    exposure_times = exposure_times[good_im]
    centroid_coords = centroid_coords[good_im]
    init_coord_list = init_coord_list[good_im]
    image_num = image_num[good_im]

    print('\nImages passed data quality check ({}): {}'.format(fil, image_num))

    """
    try:
        os.mkdir(os.path.join(dirtarget, 'ISR_Images', fil, 'WCS', 'output',
                 'centroid'))
    except FileExistsError:
        pass

    for i, im in enumerate(image_arr):
        figure, axis = plt.subplots(nrows=1, ncols=1, figsize=(10,8))
        axis.set_title('Image {}'.format(i+1), size=10)
        for label in (axis.get_xticklabels() + axis.get_yticklabels()):
            label.set_fontsize(8)

            p1 = axis.scatter(init_coord_list[i][0], init_coord_list[i][1],
                              c="thistle", label="Original",
                              edgecolors='mistyrose')
            p2 = axis.scatter(centroid_coords[i][0], centroid_coords[i][1],
                              c="rebeccapurple", label="Corrected",
                              edgecolors='mistyrose')
            print(im)
            im = axis.imshow(im, cmap='jet', origin='lower')
            c_targ = Circle((centroid_coords[i][0]+1, centroid_coords[i][1]+1),
                             5*pix_radius, fill=False, label='Aperture',
                             ls='-', color='mistyrose')
            axis.add_patch(c_targ)

            figure.subplots_adjust(right=0.8)
            cbar_ax = figure.add_axes([0.85, 0.15, 0.05, 0.7])
            cb = figure.colorbar(np.log10(im), cax=cbar_ax)
            plt.setp(cb.ax.get_yticklabels(), fontsize=8)

            plt.figlegend([p1, p2, c_targ], ['Original', 'Corrected',
                          'Aperture'], fontsize=14)
            figure.text(0.5, 0.04, 'x [pixel]', ha='center')
            figure.text(0.04, 0.5, 'y [pixel]', va='center',
                        rotation='vertical')
            plt.savefig(os.path.join(dirtarget, 'ISR_Images', fil, 'WCS',
                                     'output',
                                     'centroid/imcent_{}.pdf'.format(i)))
    """

    return aper_sum, comp_apers_return, check_apers, err, check_err, \
        date_obs, altitudes, comp_mags, saturated, exposure_times, \
        centroid_coords, init_coord_list, image_num


def write_net_counts(dirtarget, fil, date, comp_aper_sums, aper_sum,
                     check_aper_sum, t_err, check_err, date_obs, altitudes,
                     target, clabel):
    """Save output file with net count values.

    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    fil : str
        Name of filter used for images which are currently being processed.
    date : str
        Date of observation.
    comp_aper_sum : numpy.ndarray
        Filtered array of comparison star magnitudes.
    aper_sum : numpy.ndarray
        Array of aperture sums for target star in counts.
    check_aper_sum : numpy.ndarray
        Filtered array of aperture sum for check star in counts.
    t_err : numpy.ndarray
        Array of error values for each aperture sum of the target star in
        counts.
    check_err : numpy.ndarray
        Array of error values for each aperture sum of the check star in
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
                        'output/net_counts_{}_{}.txt'.format(date, fil))

    with open(path, 'w+') as f:
        f.write('#SOFTWARE=STEPUP Image Analysis\n#DELIM=,\n#DATE=JD\n' +
                '#OBSTYPE=CCD\n')
        comp_n = len(comp_aper_sums)
        header_str = str('#TARGET NAME,DATE,TARGET COUNTS,ERR,FILTER,' +
                         'C1,...,C{} COUNTS,CHECK LABEL,'.format(comp_n) +
                         'CHECK COUNTS,CHECK ERR,AIRMASS\n')
        f.write(header_str)
        comp_aper_sums = comp_aper_sums.astype(str)
        comp_sums = list(zip(*comp_aper_sums))
        zip_n = zip(date_obs, aper_sum, t_err, check_err, check_aper_sum,
                    altitudes)
        for n, (date_i, tsum, err, csum, cerr, alt) in enumerate(zip_n):
            csums = ",".join(comp_sums[n])
            zenith = np.deg2rad(90 - alt)
            airmass = 1 / np.cos(zenith)
            input_list = [target, date_i, tsum, err, fil, csums, clabel, csum,
                          airmass]
            input_string = ",".join(map(str, input_list))
            f.write(input_string + '\n')
        f.close()


def counts_to_mag(aper_sum, comp_aper_sums, err, check_err, comp_mags,
                  check_aper_sum, fil, date_obs, date):
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
    date : str
        Date of observation.

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
        figure = plt.figure(figsize=(10, 8))
        plt.plot(date_obs, 2.5 * np.log10(obj), "o", c='cadetblue',
                 label='{}'.format(fil))
        plt.ylabel("Instrumental Magnitude")
        plt.xlabel("Time [JD]")
        plt.gca().invert_yaxis()
        plt.title("Comp {}, {} = {} mag".format(i + 1, fil, mag))
        plt.savefig("output/comp{}_{}_{}.pdf".format(i + 1, date, fil))
        plt.legend()

    for i, (mag, obj) in enumerate(zip(comp_mags, comp_aper_sums)):
        # Using magnitude value of comparison star (mag) and aperture sum
        # of comparison star (obj), each image's target count value
        # (aper_sum) is determined.
        scaled_mags[i] = mag - 2.5 * np.log10(aper_sum / obj)

        # If the check star is in the image, the magnitudes of each star are
        # determined for each image.
        check_mags[i] = mag - 2.5 * np.log10(check_aper_sum / obj)

    target_err = (2.5 * np.log10((aper_sum + err) / aper_sum))
    c_err = (2.5 * np.log10((check_aper_sum + check_err) / check_aper_sum))

    # For each image, the scaled magnitude value for each comparison star is
    # averaged.
    target_mags = np.average(scaled_mags, axis=0)

    # If the check star is in the image, the calculated magnitudes are averaged
    # for each image.
    check_mags = np.array(check_mags)
    check_mags_f = np.average(check_mags, axis=0)

    print('\nTarget mags ({}): '.format(fil), target_mags)
    print('\nTarget error ({}): '.format(fil), target_err)
    print('\nCheck mags ({}): '.format(fil), check_mags_f)
    print('\nCheck error ({}): '.format(fil), c_err)

    return target_mags, target_err, check_mags_f, c_err


def mag_plot(target_mags, target_err, date_obs, target, date, fil, dirtarget,
             check_mags, check_err):
    """Plots magnitudes of target star with error and check star over time.

    Creates a subplot with the target magntidues and its error over time in
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
    axarr[1].errorbar(date_obs, check_mags, yerr=check_err, fmt='o',
                      c='cadetblue', label='{}'.format(fil))
    axarr[1].set_ylim(axarr[1].get_ylim()[::-1])
    axarr[1].invert_yaxis
    axarr[1].legend()
    axarr[1].set_title('Check Star')
    axarr[1].set_ylabel('Magnitude')
    axarr[1].set_xlabel('Time [JD]')
    f.savefig(os.path.join(dirtarget, 'ISR_Images', fil, 'WCS',
                           'output/lightcurve_{}_{}.pdf'.format(date, fil)))


def write_file(target_mags, target_err, date_obs, target, dirtarget, fil,
               altitudes, clabel, cmags, cerrs, date, image_num):
    """Writes file of target scaled magnitudes.

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
        f.write('TARGET,IMAGE NUMBER,DATE,TARGET MAG,ERROR,FILTER,CHECK ' +
                'LABEL,CHECK MAG,CHECK ERR,AIRMASS\n')
        zip_num = zip(image_num, date_obs, target_mags, target_err, cmags,
                      cerrs, altitudes)
        for num, date_i, mag, err, cmag, cerr, alt in zip_num:
            zenith = np.deg2rad(90 - alt)
            airmass = 1 / np.cos(zenith)
            input_list = [target, num, date_i, mag, err, fil, clabel, cmag,
                          cerr, airmass]
            input_string = ",".join(map(str, input_list))
            f.write(input_string + '\n')
        f.close()


def multi_filter_analysis(dirtarget, date, target, filters):
    """Produces color light curve for multi-filter observations.

    If the observation was done in more than one filter, the user is prompted
    for whether or not they would like to plot a color light curve. If they do,
    they are then prompted for which color they would like to plot by asking
    for the filters they would like to subtract. Input is comma-delimited,
    with order corresponding to order of subtraction, e.g. if the input is
    "B,V" the resultant plot will be B-V color over time. The outputted
    lightcurve is saved to dirtarget. Assumes user would only like to calculate
    one color light curve and that each color observation has the same amount
    of data points.

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
    # Sets answer to 'N' if only one filter was used for observation so that
    # the function does not plot a color light curve and returns None.
    answer = 'n'
    if len(filters) > 1:
        answer = input('\nWould you like to plot a color light curve? (Y/N): ')
    if answer.lower() == 'y':
        # Determine color user would like to plot in light curve.
        fil_color = input('\nInput magnitude filters to subtract ' +
                          '(comma-delimited): ').split(',')
        dates = []
        mags = []
        err = []
        for n, fil in enumerate(filters):
            date_fil = []
            mag_fil = []
            err_fil = []
            os.chdir(os.path.join(dirtarget, 'ISR_Images', fil, 'WCS',
                                  'output'))
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

        # Creates individual arrays for each filter of magnitudes and takes
        # the time of each color to be the average of the observation times of
        # each filter for a given point.
        mags1 = np.array(mags[0], dtype='float')
        mags2 = np.array(mags[1], dtype='float')
        err1 = np.array(err[0], dtype='float')
        err2 = np.array(err[1], dtype='float')
        dates = np.average(dates, axis=0)

        # Calculates color of target for each observation time.
        colors = mags1 - mags2
        color_err = np.sqrt(err1 ** 2 + err2 ** 2)

        # Plots and saves color light curve.
        fig = plt.figure(figsize=(10, 8))
        plt.errorbar(dates, colors, yerr=color_err, fmt='o', c='cadetblue')
        plt.title('Color Light Curve of {}, {}'.format(target, date))
        plt.xlabel('Time [JD]')
        plt.ylabel('{}-{} Color'.format(filters[0], filters[1]))
        plt.gca().invert_yaxis()
        plt.savefig(os.path.join(dirtarget,
                                 'color_lightcurve_{}.pdf'.format(date)))

    else:
        pass
