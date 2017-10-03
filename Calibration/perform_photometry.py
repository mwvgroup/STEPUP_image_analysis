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


def perform_photometry(target, dirtarget, filters, date, coords, comp_coords,
                       comp_mag, vsp_code, comp_codes, verbose=False):
    """Perform photometry on dataset and plot lightcurve.

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
        Date of observation (YYYY-MM-DD).
    coords : tuple
        Tuple containing strings of RA and Dec of target star
        ('HH:MM:SS', '+/-DD:MM:SS).
    comp_coords : tuple
        Tuple containing tuples of string(s) of RA and Dec of comparison
        star(s) (('HH:MM:SS', '+/-DD:MM:SS')).
    comp_mag : list
        List of float(s) of magnitude(s) of comparison star(s).
    verbose : boolean, optional
        Print more information about status of program.

    Returns
    -------
    None
    """
    aper_sum, err, date_obs, comp_aper_sum, comp_err = photometry(dirtarget,
                                                                  filters,
                                                                  coords,
                                                                  comp_coords)

    target_mags, target_err = counts_to_mag(aper_sum, comp_aper_sum, err,
                                            comp_err, comp_mag)

    mag_plot(target_mags, target_err, date_obs, target, date, filters,
             dirtarget)

    write_file(target_mags, target_err, date_obs, target, comp_codes, vsp_code,
               dirtarget, filters)


def photometry(dirtarget, filters, coords, comp_coords):
    """Perform aperture photometry on dataset.

    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    filters : list
        List containing string of each filter keyword found in header of flat
        field and light frame images.
    coords : tuple
        Tuple containing strings of RA and Dec of target star ('HH:MM:SS',
        '+/-DD:MM:SS).
    comp_coords : tuple
        Tuple containing tuples of string(s) of RA and Dec of comparison
        star(s) (('HH:MM:SS', '+/-DD:MM:SS')).

    Returns
    -------
    aper_sum : list
        List of count values for each image's target aperture sum.
    comp_aper_sum : list
        List of lists containing count values for the aperture sum of each
        comparison star in each image.
    err : list
        List of count values of error of aperture sum.
    date_obs : list
        List of time of observation in Julian Days.
    """
    mid_ind = len(comp_coords)

    for fil in filters:
        last_in = len(os.listdir(os.path.join(dirtarget, fil, 'WCS',
                                              'accurate_WCS')))
        aper_sum = np.empty((last_in))
        aper_sum[:] = np.nan
        err = np.empty((last_in))
        err[:] = np.nan
        date_obs = np.empty((last_in))
        date_obs[:] = np.nan

        for i, path in enumerate(os.listdir(os.path.join(dirtarget, fil, 'WCS',
                                                         'accurate_WCS'))):
            if path.endswith('.fits'):
                o_file = os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS',
                                      path)
                hdulist = fits.open(o_file)

                coordinates = SkyCoord(coords[0], coords[1], unit=(u.hourangle,
                                                                   u.deg))
                radius = 6 * u.arcsec
                r_in = 8 * u.arcsec
                r_out = 12 * u.arcsec
                aperture = SkyCircularAperture(coordinates, radius)
                annulus = SkyCircularAnnulus(coordinates, r_in=r_in,
                                             r_out=r_out)

                secpix1 = abs(hdulist[0].header['SECPIX1'])
                # Calculate area of aperture.
                aperture_area = np.pi * ((secpix1 / 6) ** (-1)) ** 2
                # Calculate area of annulus.
                outer_area = np.pi * ((secpix1 / 12) ** (-1)) ** 2
                inner_area = np.pi * ((secpix1 / 8) ** (-1)) ** 2
                annulus_area = outer_area - inner_area
                apers = (aperture, annulus)

                # Create photometry table including aperture sum for target
                # aperture and annulus sum for target annulus.
                phot_table = aperture_photometry(hdulist, apers)
                source_err = np.sqrt(phot_table['aperture_sum_0'][0])
                bkg_mean = phot_table['aperture_sum_1'] / annulus_area
                bkg_sum = bkg_mean * aperture_area
                final_sum = phot_table['aperture_sum_0'] - bkg_sum
                phot_table['residual_aperture_sum'] = final_sum
                bkg_err = np.sqrt(bkg_sum)

                date = hdulist[0].header['DATE-OBS']
                t = Time(date)
                time = t.jd
                date_obs[i] = time

                aper_sum[i] = phot_table['residual_aperture_sum'][0]

                err[i] = np.sqrt((source_err)**2 + (bkg_err)**2)
    for fil in filters:
        last_in = len(os.listdir(os.path.join(dirtarget, fil, 'WCS',
                                              'accurate_WCS')))
        mid_ind = len(comp_coords)
        comp_aper_sum = np.empty((mid_ind, last_in))
        comp_aper_sum[:] = np.nan
        comp_err = np.empty((mid_ind, last_in))
        comp_err[:] = np.nan
        for k, coord in enumerate(comp_coords):
            for j, path in enumerate(os.listdir(os.path.join(dirtarget, fil,
                                                             'WCS',
                                                             'accurate_WCS'))):
                if path.endswith('.fits'):
                    o_file = os.path.join(dirtarget, fil, 'WCS',
                                          'accurate_WCS', path)
                    hdulist = fits.open(o_file)
                    comp_coordinates = SkyCoord(coord[0], coord[1],
                                                unit=(u.hourangle, u.deg))
                    radius = 6 * u.arcsec
                    r_in = 8 * u.arcsec
                    r_out = 12 * u.arcsec
                    comp_aperture = SkyCircularAperture(comp_coordinates,
                                                        radius)
                    comp_annulus = SkyCircularAnnulus(comp_coordinates,
                                                      r_in=r_in, r_out=r_out)
                    comp_apers = (comp_aperture, comp_annulus)
                    comp_phot_table = aperture_photometry(hdulist, comp_apers)
                    comp_src_err = np.sqrt(comp_phot_table['aperture_sum_0'][0])
                    comp_bkg_mean = comp_phot_table['aperture_sum_1'] / annulus_area
                    comp_bkg_sum = comp_bkg_mean * aperture_area
                    bkg_err = np.sqrt(comp_bkg_sum)
                    # Remove background value from aperture sum and add residual
                    # aperture sum to photometry table.
                    comp_final_sum = comp_phot_table['aperture_sum_0'] - comp_bkg_sum
                    comp_phot_table['residual_aperture_sum'] = comp_final_sum

                    comp_aper_sum[k, j] = comp_phot_table['residual_aperture_sum'][0]
                    comp_err[k, j] = np.sqrt((bkg_err)**2 + (comp_src_err)**2)

    return aper_sum, err, date_obs, comp_aper_sum, comp_err


def counts_to_mag(aper_sum, comp_aper_sum, err, comp_err, comp_mag):
    """Converts instrumental measurements of brightness to magnitude values.
    """
    scaled_mags = np.empty(comp_aper_sum.shape)
    scaled_mags[:] = np.nan
    scaled_err = np.empty(comp_aper_sum.shape)
    scaled_err[:] = np.nan
    for i, mag in enumerate(comp_mag):
        for j, obj in enumerate(comp_aper_sum):
            scaled_mags[i] = mag - 2.5 * np.log10(aper_sum / obj)
            scaled_err[i] = mag*(err/obj)

    target_mags = np.average(scaled_mags, axis=0)
    target_err = np.average(scaled_err, axis=0)

    return target_mags, target_err


def mag_plot(target_mags, target_err, date_obs, target, date, filters,
             dirtarget):
    """Plots and saves figure of magnitude and error as a function of time.
    """
    for fil in filters:
        fig1 = plt.gcf()
        x = date_obs
        y = target_mags
        err = target_err
        plt.errorbar(x, y, yerr=err, fmt='o')
        plt.title('Light Curve of {}, {}'.format(target, date))
        plt.ylabel('Magnitude')
        plt.xlabel('JD')
        plt.gca().invert_yaxis()
        plt.show()
        fig1.savefig(os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS',
                                  'lightcurve.pdf'))


def write_file(target_mags, target_err, date_obs, target, comp_codes, vsp_code,
               dirtarget, filters):
    for fil in filters:
        path = os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS',
                            'output.txt')
        with open(path, 'w+') as f:
            f.write('#TYPE=Visual\n#OBSCODE=NTHC\n#SOFTWARE=STEPUP Image ' +
                    'Analysis\n#DELIM=,\n#DATE=JD\n#OBSTYPE=Visual\n')
            for mag, date in zip(target_mags, date_obs):
                    input_list = ['#' + target, mag, date, 'na', comp_codes[0],
                                  comp_codes[1], comp_codes[2], comp_codes[3],
                                  comp_codes[4], 'na', vsp_code, 'na']
                    input_string = ",".join(map(str, input_list))
                    f.write(input_string + '\n')

        f.close()


perform_photometry('AGDra', '/Users/helenarichie/tests2/ISR_Images', ['R'],
                   '2017-08-09', ('16:01:41.00', '66:48:10.0'),
                   [('16:02:54.40', '66:41:33.9'), ('16:00:56.46', '66:42:57.5'),
                    ('16:00:24.08', '66:49:29.6'), ('16:00:08.77', '66:49:20.0'),
                    ('16:01:08.41', '66:55:21.4')], [10.708, 11.644, 11.980,
                                                     12.555, 12.900],
                   'X21126DBA', (111, 120, 123, 129, 132), verbose=False)

