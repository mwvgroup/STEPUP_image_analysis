import os
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from photutils import SkyCircularAperture, SkyCircularAnnulus, aperture_photometry
from astropy.io import fits


def perform_photometry(filters, dirtarget, coords, comparison_coords):
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
    comparison_coords : tuple 
        Tuple containing strings of RA and Dec of comparison stars
        ((HH:MM:SS), (+/-DD:MM:SS))."
        
    Returns
    -------
    residual_aperture_sum : list
        List of pixel values of target star in each image with background value
        removed.
    residual_aperture_sum_c : list
        List of pixel values of comparison stars in each images with background
        value removed.
    """
    residual_aperture_sum = []
    for fil in filters:
        for path in os.listdir(os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS')):
            if path.endswith('.fits'):
                o_file = os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS', path)
                hdulist = fits.open(o_file)
                ra = coords[0]
                dec = coords[1]
                # Create SkyCoord object for target.
                coordinates = SkyCoord(ra, dec, unit = (u.hourangle, u.deg))
                # Define radius of aperture and inner/outer radius of annulus.
                radius = 6 * u.arcsec
                r_inner = 8 * u.arcsec
                r_outer = 12 * u.arcsec
                # Create SkyCircularAperture object for target.
                aperture = SkyCircularAperture(coordinates, radius)
                # Create SkyCircularAnnulus  object for target.
                annulus = SkyCircularAnnulus(coordinates, r_in=r_inner, r_out=r_outer)
                # Retrieve arcseconds per pixel in right ascensions value.
                secpix1 = abs(hdulist[0].header['SECPIX1'])
                # Calculate area of aperture.
                aperture_area = np.pi * ((secpix1/6)**(-1))**2
                # Calculate area of annulus.
                annulus_area = np.pi * ((secpix1/12)**(-1))**2 - np.pi * ((secpix1/8)**(-1))**2
                apers = (aperture, annulus)
                # Create photometry table including aperture sum for target
                # aperture and annulus sum for target annulus.
                phot_table = aperture_photometry(hdulist, apers)
                # Calculate background value.
                bkg_mean = phot_table['aperture_sum_1'] / annulus_area
                bkg_sum = bkg_mean * aperture_area
                # Remove background value from aperture sum and add residual
                # aperture sum to photometry table.
                final_sum = phot_table['aperture_sum_0'] - bkg_sum
                phot_table['residual_aperture_sum'] = final_sum
                # Add residual aperture sum to list to be returned by function.
                residual_aperture_sum.append(phot_table['residual_aperture_sum'][0])
                
        residual_aperture_sum_c = []
        for path in os.listdir(os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS')):
            if path.endswith('.fits'):
                o_file = os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS', path)
                hdulist_c = fits.open(o_file)
                ra_c = comparison_coords[0]
                dec_c = comparison_coords[1]
                radius_c = 6 * u.arcsec
                r_inner_c = 8 * u.arcsec
                r_outer_c = 12 * u.arcsec
                # Retrieve arcseconds per pixel in right ascension value.
                secpix1 = abs(hdulist[0].header['SECPIX1'])
                # Calculate area of aperture.
                aperture_area = np.pi * ((secpix1/6)**(-1))**2
                # Calculate area of annulus.
                annulus_area = np.pi * ((secpix1/12)**(-1))**2 - np.pi * ((secpix1/8)**(-1))**2
                # Create SkyCoord object for comparison star.
                coordinates_c = SkyCoord(ra_c, dec_c, unit = (u.hourangle, u.deg))
                # Create SkyCircularAperture object for comparison star.
                aperture_c = SkyCircularAperture(coordinates_c, radius_c)
                # Create SkyCircularAnnulus object for comparison star.
                annulus_c = SkyCircularAnnulus(coordinates_c, r_in=r_inner_c, r_out=r_outer_c)
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
                residual_aperture_sum_c.append(phot_table_c['residual_aperture_sum_c'][0])
                
    return residual_aperture_sum, residual_aperture_sum_c
