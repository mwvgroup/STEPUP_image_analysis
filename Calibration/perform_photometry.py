import os
import glob
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from photutils import SkyCircularAperture, SkyCircularAnnulus
import numpy as np
from photutils import aperture_photometry
from astropy.time import Time
import matplotlib.pyplot as plt
from matplotlib import gridspec


def perform_photometry(target, dirtarget, filters, date, coords, comp_ra,
                       comp_dec, comp_mags, vsp_code, cname, check_ra, check_dec,
                       rname, ref_ra, ref_dec, verbose=False):
    """Run photometry part of image analysis routine.
    """
    aper_sum, comp_aper_sums, check_aper_sum, ref_aper_sum, err, date_obs, altitudes, final_comp_mags = photometry(dirtarget, filters, coords, comp_ra, comp_dec, check_ra, check_dec, ref_ra, ref_dec, comp_mags)
    target_mags, target_err, check_mags, ref_mags, date_obs = counts_to_mag(aper_sum, comp_aper_sums, err, final_comp_mags, check_aper_sum, ref_aper_sum, date_obs)
    mag_plot(target_mags, target_err, date_obs, target, date, filters,
             dirtarget, check_mags)
    write_file(target_mags, target_err, date_obs, target, vsp_code, dirtarget,
               filters, altitudes, cname, check_mags, rname, ref_mags)


def photometry(dirtarget, filters, coords, comp_ra, comp_dec, check_ra, check_dec,
               ref_ra, ref_dec, comp_mags):
    for fil in filters:
        aper_sum, err, date_obs, altitudes = get_counts(dirtarget, coords[0], coords[1], fil)
        aper_sum = np.array(aper_sum, dtype=float)

        comp_apers = get_counts(dirtarget, comp_ra, comp_dec, fil)[0]
        comp_apers = np.array(comp_apers, dtype=float)

        check_aper_sum = (get_counts(dirtarget, check_ra, check_dec, fil))[0]
        ref_aper_sum = (get_counts(dirtarget, ref_ra, ref_dec, fil))[0]
        
        check_aper_sum = np.array(check_aper_sum, dtype=float)
        ref_aper_sum = np.array(ref_aper_sum, dtype=float)

        bad_index = []
        for i, aper in enumerate(comp_apers):
            for row in aper:
                if np.any(np.isnan(row)):
                    bad_index.append(i)

        new_comp_apers = np.delete(comp_apers, bad_index, 0)
        new_comp_mags = np.delete(comp_mags, bad_index)

        bad_index = []
        for i, aper in enumerate(new_comp_apers):
            for row in aper:
                if np.any(row <= 0):
                    bad_index.append(i)

        final_comp_apers = np.delete(new_comp_apers, bad_index, 0)
        final_comp_mags = np.delete(new_comp_mags, bad_index)

        bad_index = []
        for i, row in enumerate(check_aper_sum):
            if np.any(np.isnan(row)):
                bad_index.append(i)

        for i, row in enumerate(check_aper_sum):
            if np.any(row <= 0):
                bad_index.append(i)

        if len(bad_index) != 0:
            print('Check star either contains nan or non-positive values.')
            check_aper_sum = None

        print(ref_aper_sum)
        bad_index = []
        for i, row in enumerate(ref_aper_sum):
            if np.any(np.isnan(row)):
                bad_index.append(i)

        for i, row in enumerate(ref_aper_sum):
            if np.any(row <= 0):
                bad_index.append(i)

        if len(bad_index) != 0:
            print('Reference star either contains nan or non-positive values.')
            ref_aper_sum = None

        print(check_aper_sum, ref_aper_sum)

    return aper_sum, final_comp_apers, check_aper_sum, ref_aper_sum, err, date_obs, altitudes, final_comp_mags


def counts_to_mag(aper_sum, comp_aper_sums, err, comp_mags, check_aper_sum, ref_aper_sum, date_obs):
    scaled_mags = np.empty(comp_aper_sums.shape)
    scaled_mags[:] = np.nan
    ref_mags = np.empty(comp_aper_sums.shape)
    ref_mags[:] = np.nan
    check_mags = np.empty(comp_aper_sums.shape)
    check_mags[:] = np.nan
    scaled_err = np.empty(comp_aper_sums.shape)
    scaled_err[:] = np.nan
    for i, (mag, obj) in enumerate(zip(comp_mags, comp_aper_sums)):
        # Using magnitude value of comparison star (mag) and aperture sum 
        # of comparison star (obj), each image's target count value 
        # (aper_sum) is determined.
        scaled_mags[i] = mag - 2.5 * np.log10(aper_sum / obj)

        # Using magnitude value of comparison star (mag) and aperture sum 
        # of comparison star (obj), each image's target error count value 
        # (err) is determined. 
        scaled_err[i] = mag * (err / obj)
        if np.any(check_aper_sum != None):
            check_mags[i] = mag - 2.5 * np.log10(check_aper_sum / obj)
        if np.any(ref_aper_sum != None):
            ref_mags[i] = mag - 2.5 * np.log10(ref_aper_sum / obj)

    # For each image, the scaled magnitude value for each comparison star is 
    # averaged.
    target_mags = np.average(scaled_mags, axis=0)
    target_err = np.average(scaled_err, axis=0)

    check_mags_f = None
    if np.any(check_aper_sum != None):
        check_mags = np.array(check_mags)
        check_mags_f = np.average(check_mags, axis=0)
    else:
        check_mags_f = 0

    ref_mags_f = None
    if np.any(ref_aper_sum != None):
        ref_mags = np.array(ref_mags)
        ref_mags_f = np.average(ref_mags, axis=0)
    else:
        ref_mags_f = 0

    print('\nTarget mags: ', target_mags)
    print('\nTarget error: ', target_err)
    print('\nCheck mags: ', check_mags_f)
    print('\nRef mags: ', ref_mags_f)

    return target_mags, target_err, check_mags_f, ref_mags_f, date_obs


def mag_plot(target_mags, target_err, date_obs, target, date, filters,
             dirtarget, check_mags):
    for fil in filters:
        f, axarr = plt.subplots(2, sharex=True, gridspec_kw = {'height_ratios':[3, 1]})
        axarr[0].errorbar(date_obs, target_mags, yerr=target_err, fmt='o')
        axarr[0].set_title('Light Curve of {}, {}'.format(target,date))
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
                               'lightcurve.pdf'))
        plt.show()


def write_file(target_mags, target_err, date_obs, target, vsp_code, dirtarget,
               filters, altitudes, cname, cmags, rname, rmags):
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
                            'output.txt')
        
        with open(path, 'w+') as f:
            f.write('#TYPE=Extended\n#OBSCODE=NTHC\n#SOFTWARE=STEPUP ' +
                    'Image Analysis\n#DELIM=,\n#DATE=JD\n#OBSTYPE=CCD\n')
            for date, mag, err, cmag, rmag, alt in zip(date_obs, target_mags,
                                                       target_err, cmags, rmags,
                                                       altitudes):
                zenith = np.deg2rad(90 - alt)
                airmass = 1 / np.cos(zenith)
                input_list = [target, date, mag, err, fil, 'NO', 'STD', cname,
                              cmag, rname, rmag, airmass, 'na', vsp_code, 'na']
                input_string = ",".join(map(str, input_list))
                f.write(input_string + '\n')


def get_counts(dirtarget, rightascension, declination, fil):
    """Determine count values for star(s).
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
    good_indices = []
    
    for ra, dec in zip(rightascension, declination):
        aper_sum = np.empty(size)
        aper_sum[:] = np.nan
        for i, item in enumerate(glob.glob(os.path.join(dirtarget_wcs, '*.fits'))):
            o_file = os.path.join(dirtarget_wcs, item)
            hdulist = fits.open(o_file)
            if hdulist[0].header['WCSMATCH'] >= 20:
                good_indices.append(i)
                coords = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
                radius = 9 * u.arcsec
                r_in = 11 * u.arcsec
                r_out = 15 * u.arcsec

                aperture = SkyCircularAperture(coords, radius)
                annulus = SkyCircularAnnulus(coords, r_in=r_in,
                                             r_out=r_out)

                secpix1 = abs(hdulist[0].header['SECPIX1'])

                aper_area = np.pi * (radius / secpix1) **2
                area_out = np.pi * (r_out / secpix1) ** 2
                area_in = np.pi * (r_in /secpix1) ** 2
                annulus_area = area_out - area_in

                apers = (aperture, annulus)

                phot_table = aperture_photometry(hdulist, apers)

                source_err = np.sqrt(phot_table['aperture_sum_0'])

                bkg_mean = phot_table['aperture_sum_1'] / annulus_area
                bkg_sum = bkg_mean * aper_area
                final_sum = phot_table['aperture_sum_0'] - bkg_sum
                phot_table['residual_aperture_sum'] = final_sum

                bkg_err = np.sqrt(bkg_sum)

                aper_sum[i] = (phot_table['residual_aperture_sum'][0])
                err[i] = (np.sqrt(source_err + bkg_err))

                date = hdulist[0].header['DATE-OBS']
                t = Time(date)
                time = t.jd
                
                date_obs[i] = (time)
                altitudes[i] = (hdulist[0].header['OBJCTALT'])

        total_sum.append(aper_sum)

        print(total_sum)

    return total_sum, err, date_obs, altitudes
