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
                       comp_mag, verbose=False):
    """Peforms photometry on dataset.

    Use "verbose=True" to print information about program status. Then calls
    photometry, which returns lists of aperture sums for the target and
    comparison stars as well as the error of the target aperture sums and the
    time that each aperture sum corresponds with. It then calls counts_to_mag,
    which returns a list of pixel values scaled to magnitude using the aperture
    sums of the comparison stars and their known magnitudes. It does the same
    for the error values passed into the function. It then calls mag_plot,
    which plots the each point (magnitude, time) with an error bar for each
    magnitude value.

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
        Tuple containing strings of RA and Dec of target star ((HH:MM:SS),
        (+/-DD:MM:SS))."
    comp_coords : tuple
        Tuple containing strings of RA and Dec of comparison stars
        ((HH:MM:SS), (+/-DD:MM:SS))."
    comp_mag : float
        Magnitude of comparison star.
    verbose : boolean, optional
        Print information about status of program.

    Returns
    -------
    None
    """
    aper_sum, aper_sum_c, aper_error, date_obs = photometry(filters, dirtarget,
                                                            coords, comp_coords)

    target_mag = counts_to_mag(comp_mag, aper_sum, aper_sum_c, aper_error)

    mag_plot(target_mag, date_obs, target, date, aper_error)


def photometry(filters, dirtarget, coords, comp_coords):
    """Perform aperture photometry on dataset.

    Creates aperture and annulus to measure the counts from the star and
    background level, respectively, for both the target star and for the
    comparison star(s) The residual aperture sum for the target and for the
    comparison stars are returned.

    Inputs
    ------
    filters : list
        List containing string of each filter keyword found in header of flat
        field and light frame images.
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    coords : tuple
        Tuple containing strings of RA and Dec of target star ((HH:MM:SS),
        (+/-DD:MM:SS))."
    comp_coords : tuple
        Tuple containing strings of RA and Dec of comparison stars
        ((HH:MM:SS), (+/-DD:MM:SS))."

    Returns
    -------
    aper_sum : list
        List of pixel values of target star in each image with background value
        removed.
    aper_sum_c : list
        List of pixel values of comparison stars in each images with background
        value removed.
    aper_error : list
        Error values returned on photometry table for aperture.
    date_obs : list
        List of values stored in FITS headers for keyword "DATE-OBS".
    """
    aper_sum = []
    aper_error = []
    date_obs = []

    for fil in filters:
        for path in os.listdir(os.path.join(dirtarget, fil, 'WCS',
                                            'accurate_WCS')):
            if path.endswith('.fits'):
                o_file = os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS',
                                      path)

                hdulist = fits.open(o_file)

                ra = coords[0]
                dec = coords[1]
                # Create SkyCoord object for target.
                coordinates = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

                # Define radius of aperture and inner/outer radius of annulus.
                radius = 6 * u.arcsec
                r_inner = 8 * u.arcsec
                r_outer = 12 * u.arcsec

                # Create SkyCircularAperture object for target.
                aperture = SkyCircularAperture(coordinates, radius)

                # Create SkyCircularAnnulus  object for target.
                annulus = SkyCircularAnnulus(coordinates, r_in=r_inner,
                                             r_out=r_outer)

                # Retrieve arcseconds per pixel in right ascensions value.
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
                error = np.sqrt(phot_table['aperture_sum_0'][0])


                # Calculate background value.
                bkg_mean = phot_table['aperture_sum_1'] / annulus_area
                bkg_sum = bkg_mean * aperture_area
                # Remove background value from aperture sum and add residual
                # aperture sum to photometry table.
                final_sum = phot_table['aperture_sum_0'] - bkg_sum
                phot_table['residual_aperture_sum'] = final_sum
                bkg_error = np.sqrt(final_sum)

                # Add residual aperture sum, aperture error, and observation
                # date to lists to be returned by function.
                aper_sum.append(phot_table['residual_aperture_sum'][0])
                aper_error.append((error - bkg_error)) 
                date_obs.append(hdulist[0].header['DATE-OBS'])

        aper_sum_c = []

        for coords in comp_coords:
            aper_sum_1 = []
            for path in os.listdir(os.path.join(dirtarget, fil, 'WCS',
                                                'accurate_WCS')):
                if path.endswith('.fits'):
                    o_file = os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS',
                                          path)

                    hdulist_c = fits.open(o_file)

                    ra_c = coords[0]
                    dec_c = coords[1]
                    radius_c = 6 * u.arcsec
                    r_inner_c = 8 * u.arcsec
                    r_outer_c = 12 * u.arcsec

                    # Retrieve arcseconds per pixel in right ascension value.
                    secpix1 = abs(hdulist[0].header['SECPIX1'])

                    # Calculate area of aperture.
                    aperture_area = np.pi * ((secpix1 / 6) ** (-1)) ** 2
                    # Calculate area of annulus.
                    outer_area = np.pi * ((secpix1 / 12) ** (-1)) ** 2
                    inner_area = np.pi * ((secpix1 / 8) ** (-1)) ** 2
                    annulus_area = outer_area - inner_area

                    # Create SkyCoord object for comparison star.
                    coordinates_c = SkyCoord(ra_c, dec_c, unit=(u.hourangle,
                                                                u.deg))

                    # Create SkyCircularAperture object for comparison star.
                    aperture_c = SkyCircularAperture(coordinates_c, radius_c)

                    # Create SkyCircularAnnulus object for comparison star.
                    annulus_c = SkyCircularAnnulus(coordinates_c, r_in=r_inner_c,
                                                   r_out=r_outer_c)

                    apers = (aperture_c, annulus_c)

                    # Create photometry table including aperture sum for target
                    # aperture and annulus sum for target annulus.
                    phot_table_c = aperture_photometry(hdulist_c, apers)

                    # Calculate background value.
                    bkg_mean = phot_table_c['aperture_sum_1'] / annulus_area
                    bkg_sum = bkg_mean * aperture_area
                    # Remove background value from aperture sum and add residual
                    # aperture sum to photometry table.
                    final_sum = phot_table_c['aperture_sum_0'] - bkg_sum
                    phot_table_c['residual_aperture_sum_c'] = final_sum

                    # Add residual aperture sum to list to be returned by function.
                    aper_sum_1.append(phot_table_c['residual_aperture_sum_c'][0])
            aper_sum_c.append(aper_sum_1)

    return aper_sum, aper_sum_c, aper_error, date_obs


def counts_to_mag(comp_mag, aper_sum, aper_sum_c, aper_error):
    target_mag = []
    error_mag = []
    final_mag = []
    for mag in comp_mag:
        for comp_sum in aper_sum_c:
            scaled_mag = []
            for i in range(0, len(aper_sum)):
                comp_scale = mag - 2.5 * np.log10(aper_sum[i] / comp_sum[i])
                scaled_mag.append(comp_scale)
            target_mag.append(scaled_mag)
    for j in range(0, len(target_mag[0])):
        final_mag.append(np.mean([target_mag[0][j], target_mag[1][j]]))
                
    return final_mag


def mag_plot(target_mag, date_obs, target, date, error_mag):
    t = Time(date_obs)
    times = t.jd
    x = times
    y = target_mag
    plt.title('Light Curve of {}, {}'.format(target, date))
    plt.ylabel('Magnitude')
    plt.xlabel('JD')
    # plt.errorbar(x, y, yerr=error_mag)
    plt.plot(x, y, 'o')
    plt.gcf().autofmt_xdate()
    plt.gca().invert_yaxis()
    plt.show()
