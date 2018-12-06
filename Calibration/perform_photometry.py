import os
import glob
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy import units as u
from photutils import SkyCircularAperture, SkyCircularAnnulus
import numpy as np
from photutils import aperture_photometry
from astropy.time import Time
import matplotlib.pyplot as plt
from matplotlib import gridspec


def perform_photometry(target, dirtarget, filters, date, coords, comp_ra,
                       comp_dec, comp_mags, vsp_code, cname, check_ra,
                       check_dec, rname, ref_ra, ref_dec, verbose=False):
    """Run photometry part of image analysis routine.

    Call photometry, which calls get_counts to get the aperture sums for the
    target, comparison stars, check star, and ref star as well as the error on
    the target aperture sum. get_counts also returns a list of all of the times
    and altitudes for each image. Then these aperture sums are converted in
    counts_to_mag from count values to magnitude values using the aperture sums
    of the comparison stars and their mangitudes. This process is repeated for
    the target, check, and reference aperture sums as well as the target error.
    Then, mag_plot is called, which plots the target magnitudes with error over
    time and check magnitudes over time. Then write_file is called to write an
    output file in the AAVSO extended file format to be submitted to AAVSO
    using WebObs.

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
        List of strings of comparison stars' right ascension.
    comp_mags : list
        List of floats representing the magnitudes of the comparison stars.
    vsp_code : str
        Code to indentify AAVSO photometry table used for analysis.
    cname : str
        Label of check star corresponding to vsp_code.
    check_ra : list
        List of string of right ascension of check star.
    check_dec : list
        List of string of delination of check star.
    rname : str
        Label of reference star corresponding to vsp_code.
    ref_ra : list
        List of string of right ascension of reference star.
    ref_dec : str
        List of string of declination of reference star.
    verbose : Boolean
        Specifies whether or not user would like to print more information
        about the status of the code.

    Returns
    -------
    None
    """
    aper_sum, comp_aper_sums, check_aper_sum, ref_aper_sum, err, date_obs, altitudes, final_comp_mags, saturated, exposure_times = photometry(dirtarget, filters, coords, comp_ra, comp_dec, check_ra, check_dec, ref_ra, ref_dec, comp_mags)

    target_mags, target_err, check_mags, ref_mags = counts_to_mag(aper_sum, comp_aper_sums, err, final_comp_mags, check_aper_sum, ref_aper_sum)

    mag_plot(target_mags, target_err, date_obs, target, date, filters,
             dirtarget, check_mags)

    write_file(target_mags, target_err, date_obs, target, vsp_code, dirtarget,
               filters, altitudes, cname, check_mags, rname, ref_mags, date)

    print('Saturated images: {}'.format(saturated))
    print('Exposure times: {}'.format(exposure_times))


def photometry(dirtarget, filters, coords, comp_ra, comp_dec, check_ra,
               check_dec, ref_ra, ref_dec, comp_mags):
    """Get aperture sums for target, comparison, check, and reference stars.

    Calls get_counts for the list of right ascension(s) and declination(s) for
    the target star, each comparison star, the check star, and the reference
    star. The err, date_obs, and altitudes are defined when get_counts is
    called for the target star. Then, if any of the comparison stars are not
    image (i.e. they have either nan or negative values in their aper_sum
    in the arrays), then they are removed from the comp_aper_sum array and the
    corresponding magnitude is removed from comp_mags. The same process is
    repreated for the check and reference star, except there is no
    corresponding magnitude to delete.
    
    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    filters : list
        List containing string of each filter keyword found in header of flat 
        field and light frame images.
    coords : list
        List of list of string of target right ascension and declination.
    comp_ra : list
        List of strings of comparison stars' right ascension.
    comp_dec : list
        List of strings of comparison stars' right ascension.
    check_ra : list
        List of string of right ascension of check star.
    check_dec : list
        List of string of delination of check star.
    ref_ra : list
        List of string of right ascension of reference star.
    ref_dec : str
        List of string of declination of reference star.
    comp_mags : list
        List of floats representing the magnitudes of the comparison stars.
    
    Returns
    -------
    aper_sum : numpy.ndarray
        Array of aperture sums for target star in counts.
    final_comp_apers : numpy.ndarray
        Filtered array of aperture sums for comparison star in counts.
    check_aper_sum : numpy.ndarray
        Filtered array of aperture sum for check star in counts.
    ref_aper_sum : numpy.ndarray
        Filtered array of aperture sum for reference star in counts.
    err : numpy.ndarray
        Array of error values for each aperture sum of the target star in
        counts.
    date_obs : np.ndarray
        Array of times each image was taken in Julian Days.
    altitudes : np.ndarray
        Array of object altitudes for each image.
    final_comp_mags : numpy.ndarray
        Filtered array of comparison star magnitudes.
    """
    for fil in filters:
        # Get aperture sum, error of aperture sum, times of data collection,
        # and altitudes for target.
        aper_sum, err, date_obs, altitudes, saturated_dudes, exposure_times = get_counts(dirtarget, coords[0],
                                                                                         coords[1], fil)
        aper_sum = np.array(aper_sum, dtype=float)

        # Get aperture sums for each somparison star.
        comp_apers = get_counts(dirtarget, comp_ra, comp_dec, fil)[0]
        comp_apers = np.array(comp_apers, dtype=float)

        # Get aperture sum of the check and reference stars.
        check_aper_sum = (get_counts(dirtarget, check_ra, check_dec, fil))[0]
        ref_aper_sum = (get_counts(dirtarget, ref_ra, ref_dec, fil))[0]
        check_aper_sum = np.array(check_aper_sum, dtype=float)
        ref_aper_sum = np.array(ref_aper_sum, dtype=float)

        # Determine if any comparison stars are not in the image by checking
        # for the presence of nan values.
        bad_index = []
        for i, aper in enumerate(comp_apers):
            for row in aper:
                if np.any(np.isnan(row)):
                    bad_index.append(i)

        # Remove aperture sum(s) of comparison star(s) that are not in the
        # image.
        new_comp_apers = np.delete(comp_apers, bad_index, 0)
        new_comp_mags = np.delete(comp_mags, bad_index)


        # Determine if any comparison stars are not in the image by checking
        # for the presence of negative values.
        bad_index = []
        for i, aper in enumerate(new_comp_apers):
            for row in aper:
                if np.any(row <= 0):
                    bad_index.append(i)

        # Remove aperture sum(s) and magntidues of comparison star(s) that are
        # not in the image.
        final_comp_apers = np.delete(new_comp_apers, bad_index, 0)
        final_comp_mags = np.delete(new_comp_mags, bad_index)

        # Determine if the check star is not in the image by checking for the
        # presence of nan or negative values.
        bad_index = []
        for i, row in enumerate(check_aper_sum):
            if np.any(np.isnan(row)):
                bad_index.append(i)

        for i, row in enumerate(check_aper_sum):
            if np.any(row <= 0):
                bad_index.append(i)

        # If check star is not in image, let check_aper_sum = None.
        if len(bad_index) != 0:
            print('Check star either contains nan or non-positive values.')
            check_aper_sum = None

        # Determine if the reference star is not in the image by checking for
        # the presence of nan or negative values.
        bad_index = []
        for i, row in enumerate(ref_aper_sum):
            if np.any(np.isnan(row)):
                bad_index.append(i)

        for i, row in enumerate(ref_aper_sum):
            if np.any(row <= 0):
                bad_index.append(i)

        # If check star is not in image, let ref_aper_sum = None.
        if len(bad_index) != 0:
            print('Reference star either contains nan or non-positive values.')
            ref_aper_sum = None

    return aper_sum, final_comp_apers, check_aper_sum, ref_aper_sum, err, date_obs, altitudes, final_comp_mags, saturated_dudes, exposure_times


def counts_to_mag(aper_sum, comp_aper_sums, err, comp_mags, check_aper_sum,
                  ref_aper_sum):
    """Scales the aperture sums from counts to magnitudes.

    Initializes arrays for target

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
    ref_aper_sum : numpy.ndarray
        Filtered array of aperture sum for reference star in counts.

    Returns
    -------
    target_mags : numpy.ndarray
        Array of magnitude values for the target.
    target_err : numpy.ndarray
        Array of error values for target magnitudes.
    check_mags_f : numpy.ndarray
        Array of magnitude values for the check star.
    ref_mags_f : numpy.ndarray
        Array of magnitude values for the reference star.
    """
    # Initialize arrays for target magnitudes, check magnitudes, reference
    # magnitudes, and errors in target magnitudes.
    scaled_mags = np.empty(comp_aper_sums.shape)
    scaled_mags[:] = np.nan
    check_mags = np.empty(comp_aper_sums.shape)
    check_mags[:] = np.nan
    ref_mags = np.empty(comp_aper_sums.shape)
    ref_mags[:] = np.nan

    for i, (mag, obj) in enumerate(zip(comp_mags, comp_aper_sums)):
        # Using magnitude value of comparison star (mag) and aperture sum 
        # of comparison star (obj), each image's target count value 
        # (aper_sum) is determined.
        scaled_mags[i] = mag - 2.5 * np.log10(aper_sum / obj)

        # If the check star and reference star is in the image, the magnitudes
        # of each star are determined for each image.
        if np.all(check_aper_sum != None):
            check_mags[i] = mag - 2.5 * np.log10(check_aper_sum / obj)
        if np.all(ref_aper_sum != None):
            ref_mags[i] = mag - 2.5 * np.log10(ref_aper_sum / obj)

    target_err = (2.5 * np.log10((aper_sum + err) / aper_sum))[0]

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

    # If the reference star is in the image, the calculated magnitudes are
    # averaged for each image.
    ref_mags_f = None
    if np.all(ref_aper_sum != None):
        ref_mags = np.array(ref_mags)
        ref_mags_f = np.average(ref_mags, axis=0)
    else:
        ref_mags_f = 0

    print('\nTarget mags: ', target_mags)
    print('\nTarget error: ', target_err)
    print('\nCheck mags: ', check_mags_f)
    print('\nRef mags: ', ref_mags_f)

    return target_mags, target_err, check_mags_f, ref_mags_f


def mag_plot(target_mags, target_err, date_obs, target, date, filters,
             dirtarget, check_mags):
    """Plot magnitudes of target star with error and check star over time.

    Creates a subplot with the target magntidues and their error over time in
    Julian Days on top and the check star magnitudes on the bottom with the
    same x-axis as the target star. A PDF of the lightcurve is written to
    /home/depot/STEPUP/raw/<name-of-target>/<date-of-observation>/ISR_Images/
    <filter-name>/WCS/accurate_WCS.

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
    date_obs : numpy.ndarray
        Array of times each image was taken in Julian Days.
    date : str
        Date of observation.
    filters : list
        List containing string of each filter keyword found in header of flat 
        field and light frame images.
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    check_mags : numpy.ndarray
        Array of magnitude values for the check star.

    Returns
    -------
    None
    """
    for fil in filters:
        f, axarr = plt.subplots(2, sharex=True,
                                gridspec_kw = {'height_ratios':[3, 1]})
        axarr[0].errorbar(date_obs, target_mags, yerr=target_err, fmt='o')
        axarr[0].set_title('Light Curve of {}, {}'.format(target, date))
        axarr[0].set_ylabel('Magnitude')
        axarr[0].invert_yaxis
        axarr[0].set_ylim(axarr[0].get_ylim()[::-1])
        axarr[1].scatter(date_obs, check_mags)
        axarr[1].set_ylim(axarr[1].get_ylim()[::-1])
        axarr[1].invert_yaxis
        axarr[1].set_title('Check Star')
        axarr[1].set_ylabel('Magnitude')
        axarr[1].set_xlabel('Time (JD)')
        f.savefig(os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS',
                               'lightcurve_{}.pdf'.format(date)))
        plt.show()


def write_file(target_mags, target_err, date_obs, target, vsp_code, dirtarget,
               filters, altitudes, cname, cmags, rname, rmags, date):
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
    comp_codes : tuple
        Identifier for each comparison star within the VSP.
    vsp_code : str
        Identifier for AAVSO Variable Star Plotter chart.
    target : str
        Name of target.
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    filters : list
        List containing string of each filter keyword found in header of flat
        field and light frame images.
    altitudes : numpy.ndarray
        Array of altitudes of each image.
    cname : str
        AAVSO Unique Identifier (AUID) for the comparison star.
    cmags : numpy.ndarray
        Array of count values of aperture sum of comparison star of closest
        location to the target star.
    kmags : numpy.ndarray
        Array of count values of apeture sum of check star.
    kname : str
        AAVSO Unique Identifier (AUID) for the check star.
    kmags : numpy.ndarray
        Array of count values of apeture sum of check star.
        
    Returns
    -------
    None
    """
    for fil in filters:
        path = os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS', 
                            'output_{}.txt'.format(date))
        
        with open(path, 'w+') as f:
            f.write('#TYPE=Extended\n#OBSCODE=PIT01\n#SOFTWARE=STEPUP ' +
                    'Image Analysis\n#DELIM=,\n#DATE=JD\n#OBSTYPE=CCD\n')
            for date_i, mag, err, cmag, rmag, alt in zip(date_obs, target_mags,
                                                       target_err, cmags, rmags,
                                                       altitudes):
                zenith = np.deg2rad(90 - alt)
                airmass = 1 / np.cos(zenith)
                input_list = [target, date_i, mag, err, fil, 'NO', 'STD', cname,
                              cmag, rname, rmag, airmass, 'na', vsp_code, 'na']
                input_string = ",".join(map(str, input_list))
                f.write(input_string + '\n')
        f = open(os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS',
                              'output_{}.txt'.format(date)), "r")
        lines = f.readlines()
        f.close()                         
        f = open(os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS',
                              'output_{}.txt'.format(date)), "w")
        f.write(lines[0])
        f.write(lines[1])
        f.write(lines[2])
        f.write(lines[3])
        f.write(lines[4])
        f.write(lines[5])
        for line in lines:
            if line.startswith(target):
                line_list = line.split(",")
                if line_list[1] != "nan":
                    f.write(line)
        f.close()

def get_counts(dirtarget, rightascension, declination, fil):
    """Determine count values for star(s).

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
    """
    dirtarget_wcs = os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS')
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
    
    for ra, dec in zip(rightascension, declination):
        aper_sum = np.empty(size)
        aper_sum[:] = np.nan
        for i, item in enumerate(sorted(glob.glob(os.path.join(dirtarget_wcs,
                                                               '*.fits')))):
            o_file = os.path.join(dirtarget_wcs, item)
            hdulist = fits.open(o_file)
            if hdulist[0].header['WCSMATCH'] < 20:
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

            # Read in a 30 x 30 square centered at the star for which the
            # counts are being summed.
            star = image_array[(py_int - 14):(py_int + 16),
                               (px_int - 14):(px_int + 16)]

            # Ensure that the square is entirely on the image.
            if ((py_int - 14) < 0) or ((py_int + 16) > 2084):
                continue
            if ((px_int - 14) < 0) or ((px_int + 16) > 3072):
                continue

            # Flatten the array so that it is one-dimensional.
            star_flat = star.reshape(1, 900)

            # Determine if the image is saturated at the star's position using
            # the expected saturation level. If saturated, the loop will move
            # on to the next image.
            if np.amax(star_flat) >= hdulist[0].header['SATLEVEL']:
                saturated.append(item)
                continue

            # Define aperture and annulus radii.
            radius = 5.1669 * u.arcsec
            r_in = 6.3151 * u.arcsec
            r_out = 8.6115 * u.arcsec

            # Create SkyCircularAperture and SkyCircularAnnulus objects
            # centered at the position of the star whose counts are being
            # summed.
            aperture = SkyCircularAperture(coords, radius)
            annulus = SkyCircularAnnulus(coords, r_in=r_in,
                                         r_out=r_out)

            secpix1 = abs(hdulist[0].header['SECPIX1'])

            # Determine the area of the aperture and annulus using the
            # arcseconds per pixel in the horizontal dimension header keyword.
            aper_area = np.pi * (radius / secpix1) **2
            area_out = np.pi * (r_out / secpix1) ** 2
            area_in = np.pi * (r_in /secpix1) ** 2
            annulus_area = area_out - area_in

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
            source_err = np.sqrt(phot_table['residual_aperture_sum'])

            bkg_err = np.sqrt(bkg_sum)

            aper_sum[i] = phot_table['residual_aperture_sum'][0]
            err[i] = np.sqrt(source_err ** 2  + bkg_err ** 2)

            # Determine the time at which the image was taken in Julian Days
            # and the altitude and write them to an array.
            date = hdulist[0].header['DATE-OBS']
            t = Time(date)
            time = t.jd

            date_obs[i] = (time)
            altitudes[i] = (hdulist[0].header['OBJCTALT'])

        total_sum.append(aper_sum)

    return total_sum, err, date_obs, altitudes, saturated, exposure_times
