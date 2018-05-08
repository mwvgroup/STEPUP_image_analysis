import os
import glob
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from photutils import SkyCircularAperture, SkyCircularAnnulus
import numpy as np
from photutils import aperture_photometry
from astropy.time import Time


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

        good_date_obs = np.take(date_obs, good_indices)
        print(date_obs, len(date_obs))
        print(good_date_obs, len(good_date_obs))

    return total_sum, err, good_date_obs, altitudes
