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
                       comp_dec, comp_mags, vsp_code, rname, ref_ra, ref_dec,
                       cname, check_ra, check_dec, verbose=False):
    """Run photometry part of image analysis routine.
    """
    aper_sum, comp_aper_sums, check_aper_sum, ref_aper_sum, err, date_obs, altitudes = photometry(dirtarget, filters, coords, comp_ra, comp_dec, ref_ra, ref_dec, check_ra, check_dec)
    print('Final length of date_obs: ', len(date_obs))
    target_mags, target_err, check_mags, ref_mags, date_obs = counts_to_mag(aper_sum, comp_aper_sums, err, comp_mags, check_aper_sum, ref_aper_sum, date_obs)
    print('Target magnitudes: ', target_mags, 'Target error: ', target_err, 'Observation dates: ', date_obs)
    mag_plot(target_mags, target_err, date_obs, target, date, filters,
             dirtarget, check_mags, cname)

    write_file(target_mags, target_err, date_obs, target, vsp_code, dirtarget, filters,
               altitudes, rname, ref_mags, cname, check_mags)


def photometry(dirtarget, filters, coords, comp_ra, comp_dec, ref_ra, ref_dec,
               check_ra, check_dec):
    for fil in filters:
        aper_sum, err, date_obs, altitudes = get_counts(dirtarget, coords[0], coords[1], fil)
        aper_sum = np.array(aper_sum, dtype=float)

        comp_apers = get_counts(dirtarget, comp_ra, comp_dec, fil)[0]
        comp_apers = np.array(comp_apers, dtype=float)

        check_aper_sum = (get_counts(dirtarget, check_ra, check_dec, fil))[0]
        ref_aper_sum = (get_counts(dirtarget, ref_ra, ref_dec, fil))[0]
        
        check_aper_sum = np.array(check_aper_sum, dtype=float)
        ref_aper_sum = np.array(ref_aper_sum, dtype=float)

    return aper_sum, comp_apers, check_aper_sum, ref_aper_sum, err, date_obs, altitudes


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
        if np.all(obj != np.nan) and np.all(obj > 0):
            print('This object contains only positive numbers.', obj)
            scaled_mags[i] = mag - 2.5 * np.log10(aper_sum / obj)

            # Using magnitude value of comparison star (mag) and aperture sum 
            # of comparison star (obj), each image's target error count value 
            # (err) is determined. 
            # Is this the right way to calculate this?
            scaled_err[i] = mag * (err / obj)
            if np.all(check_aper_sum != np.nan) and np.all(check_aper_sum > 0):
                print('Check star is in the image.')
                check_mags[i] = mag - 2.5 * np.log10(check_aper_sum / obj)
            else:
                print('Check star is not in the image.')
            if np.all(ref_aper_sum != np.nan) and np.all(ref_aper_sum > 0):
                print('Reference star is in the image.')
                ref_mags[i] = mag - 2.5 * np.log10(ref_aper_sum / obj)
            else:
                print('Reference star is not in the image.')
        else:
            print('This object contains either nan or negative values.', obj)
            continue

    # For each image, the scaled magnitude value for each comparison star is 
    # averaged.
    good_mags = []
    for index in scaled_mags:
        nantest = np.isnan(index)
        if np.any(nantest):
            continue
        else:
            good_mags.append(index)
    good_err = []
    for index in scaled_err:
        nantest = np.isnan(index)
        if np.any(nantest):
            continue
        else:
            good_err.append(index)
    good_mags = np.array(good_mags)
    good_err = np.array(good_err)
    target_mags = np.average(good_mags, axis=0)
    target_err = np.average(good_err, axis=0)
    check_mag = np.array(check_mags)
    ref_mags = np.array(ref_mags)
    check_mags_f = np.average(check_mags, axis=0)
    ref_mags_f = np.average(ref_mags, axis=0)

    print(target_err)

    return target_mags, target_err, check_mags_f, ref_mags_f, date_obs


def mag_plot(target_mags, target_err, date_obs, target, date, filters,
             dirtarget, scaled_refmags, kname):
        f, axarr = plt.subplots(2, sharex=True, gridspec_kw = {'height_ratios':[3, 1]})
        axarr[0].errorbar(date_obs, target_mags, yerr=target_err, fmt='o')
        axarr[0].set_title('Light Curve of {}, {}'.format(target,date))
        axarr[0].set_ylabel('Magnitude')
        axarr[0].invert_yaxis
        axarr[0].set_ylim(axarr[0].get_ylim()[::-1])
        axarr[1].scatter(date_obs, scaled_refmags)
        axarr[1].set_ylim(axarr[1].get_ylim()[::-1])
        axarr[1].invert_yaxis
        axarr[1].set_title('Check Star')
        axarr[1].set_ylabel('Magnitude')
        axarr[1].set_xlabel('Time (JD)')
        f.savefig(os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS',
                               'lightcurve.pdf'))
        plt.show()


def write_file(target_mags, target_err, date_obs, target, vsp_code, dirtarget,
               filters, altitudes, cname, cmags, kname, kmags):
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
            for date, mag, err, cmag, kmag, alt in zip(date_obs, target_mags,
                                                       target_err, cmags, kmags,
                                                       altitudes):
                zenith = np.deg2rad(90 - alt)
                airmass = 1 / np.cos(zenith)
                input_list = [target, date, mag, err, fil, 'NO', 'STD', cname,
                              cmag, kname, kmag, airmass, 'na', vsp_code, 'na']
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
                radius = 90 * u.arcsec
                r_in = 110 * u.arcsec
                r_out = 150 * u.arcsec

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

        good_date_obs = np.take(date_obs, good_indices)
        print(date_obs, len(date_obs))
        print(good_date_obs, len(good_date_obs))
        print(err, len(err))

    return total_sum, err, good_date_obs, altitudes
