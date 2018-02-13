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
                       comp_mags, vsp_code, cname, c_coords, kname,
                       k_coords, verbose=False):
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
    comp_coords : list
        List containing tuples of string(s) of RA and Dec of comparison
        star(s) [('HH:MM:SS', '+/-DD:MM:SS')].
    comp_mags : list
        List of float(s) of magnitude(s) of comparison star(s).
    vsp_code : str
        Identifier for AAVSO Variable Star Plotter chart.
    cname : str
        AAVSO Unique Identifier (AUID) for the comparison star.
    c_coords : tuple
        Tuple containing strings of RA and Dec of comparison star of closest
        location to target star ('HH:MM:SS', '+/-DD:MM:SS).
    kname : str
        AAVSO Unique Identifier (AUID) for the check star.
    k_coords : tuple
        Tuple containing strings of RA and Dec of check star
        ('HH:MM:SS', '+/-DD:MM:SS).
    verbose : boolean, optional
        Print more information about status of program.
    Returns
    -------
    None
    """
    aper_sum, err, date_obs, comp_aper_sum, altitudes, cmags, kmags = photometry(dirtarget, filters, coords, comp_coords, cname, c_coords, kname, k_coords)

    target_mags, target_err, scaled_cmags, scaled_kmags = counts_to_mag(aper_sum, comp_aper_sum, err, comp_mags, kmags, cmags)

    mag_plot(target_mags, target_err, date_obs, target, date, filters,
             dirtarget)

    write_file(target_mags, target_err, date_obs, target, vsp_code, dirtarget, filters,
               altitudes, cname, cmags, kname, kmags)


def photometry(dirtarget, filters, coords, comp_coords, cname, c_coords, kname,
               k_coords):
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
    comp_coords : list
        List containing tuples of string(s) of RA and Dec of comparison
        star(s) [('HH:MM:SS', '+/-DD:MM:SS')].
    cname : str
        AAVSO Unique Identifier (AUID) for the comparison star.
    c_coords : tuple
        Tuple containing strings of RA and Dec of comparison star of closest
        location to target star ('HH:MM:SS', '+/-DD:MM:SS).
    kname : str
        AAVSO Unique Identifier (AUID) for the check star.
    k_coords : tuple
        Tuple containing strings of RA and Dec of check star
        ('HH:MM:SS', '+/-DD:MM:SS). 

    Returns
    -------
    aper_sum : numpy.ndarray
        Array of count values for each image's target aperture sum.
    comp_aper_sum : numpy.ndarray
        Array containing count values for the aperture sum of each comparison
        star in each image.
    err : numpy.ndarray
        Array of count values of error of aperture sum.
    date_obs : numpy.ndarray
        Array of time of observation in Julian Days.
    altitudes : numpy.ndarray
        Array of altitudes of each image.
    cmags : numpy.ndarray
        Array of count values of aperture sum of comparison star of closest
        location to the target star.
    kmags : numpy.ndarray
        Array of count values of apeture sum of check star.
    """
    for fil in filters:
        path = os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS')
        # Get number of FITS files in path.
        last_in = len(glob.glob(os.path.join(path, '*.fits')))
        # Initialize numpy array of nan vaules with size that is equal to the
        # number of images in the dataset for the aperture sums.
        aper_sum = np.empty((last_in))
        aper_sum[:] = np.nan
        # Initialize numpy array of nan vaules with size that is equal to the
        # number of images in the dataset for the err values for aperture sums.
        err = np.empty((last_in))
        err[:] = np.nan
        
        # Initialize numpy array of nan vaules with size that is equal to the
        # number of images in the dataset for the time that each image was
        # taken.
        date_obs = np.empty((last_in))
        date_obs[:] = np.nan

        # Initialize numpy array of nan values with size equal to the number
        # of images in dataset for aperture sum of comparison star of closest
        # proximity to target star.
        cmags = np.empty((last_in))
        cmags[:] = np.nan

        # Initialize numpy array of nan values with size equal to the number
        # of images in dataset for aperture sum of check star.
        kmags = np.empty((last_in))
        kmags[:] = np.nan

        # Initialize numpy array of nan values with size that is equal to the
        # number of images in the dataset for the altitude of each image.
        altitudes = np.empty((last_in))
        altitudes[:] = np.nan

        # Open all files in path that are FITS files.
        for i, item in enumerate(glob.glob(os.path.join(path, '*.fits'))):
            if item.endswith('.fits'):
                o_file = os.path.join(path, item)
                hdulist = fits.open(o_file)
                if hdulist[0].header['WCSMATCH'] >= 20:

                    # Create SkyCoord object at target, comparison star, and check
                    # star position.
                    coordinates = SkyCoord(coords[0], coords[1], unit=(u.hourangle, 
                                                                           u.deg))
                    c_coord = SkyCoord(c_coords[0], c_coords[1],
                                       unit=(u.hourangle, u.deg))
                    k_coord = SkyCoord(k_coords[0], k_coords[1],
                                       unit=(u.hourangle, u.deg))
                        
                    # Set radius of aperture.
                    radius = 6 * u.arcsec
                    # Set inner/outer radius of annulus.
                    r_in = 8 * u.arcsec
                    r_out = 12 * u.arcsec
                        
                    # Create SkyCircularAperture object at target, comp star, and
                    # check star position.
                    aperture = SkyCircularAperture(coordinates, radius)
                    c_aperture = SkyCircularAperture(c_coord, radius)
                    k_aperture = SkyCircularAperture(k_coord, radius)
                        
                    # Create SkyCircularAnnulus object at target comp star, and
                    # check star position.
                    annulus = SkyCircularAnnulus(coordinates, r_in=r_in,
                                                 r_out=r_out)
                    c_annulus = SkyCircularAnnulus(c_coord, r_in=r_in,
                                                   r_out=r_out)
                    k_annulus = SkyCircularAnnulus(k_coord, r_in=r_in,
                                                   r_out=r_out)

                    # Get value of SECPIX1 from FITS headers.
                    secpix1 = abs(hdulist[0].header['SECPIX1'])
                    # Calculate area of aperture.
                    aperture_area = np.pi * ((secpix1 / 6) ** (-1)) ** 2
                    # Calculate area of annulus.
                    outer_area = np.pi * ((secpix1 / 12) ** (-1)) ** 2
                    inner_area = np.pi * ((secpix1 / 8) ** (-1)) ** 2
                    annulus_area = outer_area - inner_area
                        
                    apers = (aperture, annulus)
                    c_apers = (c_aperture, c_annulus)
                    k_apers = (k_aperture, k_annulus)

                    # Create photometry table for target, comp star, and check
                    # star aperture and annulus sum.
                    phot_table = aperture_photometry(hdulist, apers)
                    c_phot_table = aperture_photometry(hdulist, c_apers)
                    k_phot_table = aperture_photometry(hdulist, k_apers)
                        
                    # Calculate error value for aperture sum.
                    source_err = np.sqrt(phot_table['aperture_sum_0'][0])
                        
                    # Calculate background value for target.
                    bkg_mean = phot_table['aperture_sum_1'] / annulus_area
                    bkg_sum = bkg_mean * aperture_area
                    final_sum = phot_table['aperture_sum_0'] - bkg_sum
                    phot_table['residual_aperture_sum'] = final_sum

                    # Calculate background value for comp star.
                    c_bkg_mean = c_phot_table['aperture_sum_1'] / annulus_area
                    c_bkg_sum = c_bkg_mean * aperture_area
                    c_final_sum = c_phot_table['aperture_sum_0'] - c_bkg_sum
                    c_phot_table['residual_aperture_sum'] = c_final_sum

                    # Calculate background value for check star.
                    k_bkg_mean = k_phot_table['aperture_sum_1'] / annulus_area
                    k_bkg_sum = k_bkg_mean * aperture_area
                    k_final_sum = k_phot_table['aperture_sum_0'] - k_bkg_sum
                    k_phot_table['residual_aperture_sum'] = k_final_sum
                    
                    # Calculate error value for target background level.
                    bkg_err = np.sqrt(bkg_sum)

                    # Add aperture sum for target, comp star, and check star to
                    # their respective arrays.
                    aper_sum[i] = phot_table['residual_aperture_sum'][0]
                    cmags[i] = c_phot_table['residual_aperture_sum'][0]
                    kmags[i] = k_phot_table['residual_aperture_sum'][0]

                    # Add error to err.
                    err[i] = np.sqrt(source_err + bkg_err)
                    
                    # Retrieve time that image was taken and convert value to 
                    # Julian Days.
                    date = hdulist[0].header['DATE-OBS']
                    t = Time(date)
                    time = t.jd
                    # Add time to date_obs.ß
                    date_obs[i] = time

                    # Add altitude to altitudes.
                    altitudes[i] = hdulist[0].header['OBJCTALT']
                
    for fil in filters:
        path = os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS')
        # Get last index number to initialize comp_aper_sum and comp_err numpy 
        # arrays.
        last_in = len(glob.glob(os.path.join(path, '*.fits')))
        # Get middle index number to initialize comp_aper_sum and comp_err numpy 
        # arrays.
        mid_ind = len(comp_coords)
        
        # Initialize numpy array of nan vaules with size that is equal to the
        # numpber of comparison stars by the number of images in the dataset 
        # for the comparison aperture sums.
        comp_aper_sum = np.empty((mid_ind, last_in))
        comp_aper_sum[:] = np.nan
        
        # Open all files in path that are FITS files.
        for k, coord in enumerate(comp_coords):
            for j, item in enumerate(glob.glob(os.path.join(path, '*.fits'))):
                if item.endswith('.fits'):
                    o_file = os.path.join(path, item)
                    hdulist = fits.open(o_file)
                    if hdulist[0].header['WCSMATCH'] >= 20:
                
                        # Create SkyCoords object at position of k'th comparison
                        # star.
                        comp_coordinates = SkyCoord(coord[0], coord[1],
                                                    unit=(u.hourangle, u.deg))
                        
                        # Set radius of aperture.
                        radius = 6 * u.arcsec
                        # Set inner/outer radius of annulus.
                        r_in = 8 * u.arcsec
                        r_out = 12 * u.arcsec
                        
                        # Create SkyCircularAperture object for comparison star.
                        comp_aperture = SkyCircularAperture(comp_coordinates,
                                                            radius)
                        # Create SkyCircularAnnulus object for comparison star.
                        comp_annulus = SkyCircularAnnulus(comp_coordinates,
                                                          r_in=r_in, r_out=r_out)
                        
                        comp_apers = (comp_aperture, comp_annulus)
                        
                        # Create photometry table for comparison aperture and 
                        # annulus sum.
                        comp_phot_table = aperture_photometry(hdulist, comp_apers)
                        
                        # How do I shorten this line?
                        # Calculate background value.
                        # How do I shorten this line?
                        comp_bkg_mean = comp_phot_table['aperture_sum_1'] / annulus_area
                        comp_bkg_sum = comp_bkg_mean * aperture_area
                        # How do I shorten this line?
                        comp_final_sum = comp_phot_table['aperture_sum_0'] - comp_bkg_sum
                        comp_phot_table['residual_aperture_sum'] = comp_final_sum
                    

                        # Add aperture sum for k'th comparison star and i'th image
                        # to comp_aper_sum.
                        # How do I shorten this line.
                        comp_aper_sum[k, j] = comp_phot_table['residual_aperture_sum'][0]

    print(aper_sum, comp_aper_sums, cmags, kmags)

    return aper_sum, err, date_obs, comp_aper_sum, altitudes, cmags, kmags


def counts_to_mag(aper_sum, comp_aper_sum, err, comp_mags, kmags, cmags):
    """Converts instrumental measurements of brightness to magnitude values.
    
    Using the aperture sums of the target star and the comparison stars
    determined by photometry and the known magnitudes of the comparison stars, 
    counts_to_mag converts the count value of the target aperture sum to a
    magnitude value. That value is then averaged for each image and returned 
    Also returns error values that have been converted from count to magnitude
    values.
    
    Parameters
    ----------
    aper_sum : numpy.ndarray
        Array of count values for each image's target aperture sum.
    comp_aper_sum : numpy.ndarray
        Array containing count values for the aperture sum of each comparison
        star in each image.  
    err : numpy.ndarray
        Array of count values of error of aperture sum.
    comp_mag : list
        List of float(s) of magnitude(s) of comparison star(s).
    cmags : numpy.ndarray
        Array of count values of aperture sum of comparison star of closest
        location to the target star.
    kmags : numpy.ndarray
        Array of count values of apeture sum of check star.
    
    Returns
    -------
    target_mags : numpy.ndarray
        Array of magnitude values of target star for each image.
    target_err : numpy.ndarray 
        Array of error values of magntidue of target star for each image.
    cmags : numpy.ndarray
        Array of magnitude values of comparison star of closest location to
        the target star.
    kmags : numpy.ndarray
        Array of magnitude values of check star.
    """
    # Initialize numpy array for scaled magntidue values of nan values with 
    # size of comp_aper_sum, which is the number of comparison stars by the 
    # number of images.
    scaled_mags = np.empty(comp_aper_sum.shape)
    scaled_mags[:] = np.nan
    scaled_cmags = np.empty(comp_aper_sum.shape)
    scaled_cmags[:] = np.nan
    scaled_kmags = np.empty(comp_aper_sum.shape)
    scaled_kmags[:] = np.nan
    # Initialize numpy array for scaled error values of nan values with size
    # of comp_aper_sum, which is the number of comparison stars by the number
    # of images.
    scaled_err = np.empty(comp_aper_sum.shape)
    scaled_err[:] = np.nan
    for i, (mag, obj) in enumerate(zip(comp_mags, comp_aper_sum)):
        # Using magnitude value of comparison star (mag) and aperture sum 
        # of comparison star (obj), each image's target count value 
        # (aper_sum) is determined.
        if np.all(obj != np.nan) and np.all(obj > 0):
            scaled_mags[i] = mag - 2.5 * np.log10(aper_sum / obj)

            # Using magnitude value of comparison star (mag) and aperture sum 
            # of comparison star (obj), each image's target error count value 
            # (err) is determined. 
            # Is this the right way to calculate this?
            scaled_err[i] = mag * (err / obj)
        else:
            continue

        print(cmags)
        print(kmags)

        scaled_cmags[i] = - 2.5 * np.log10(cmags)
        scaled_kmags[i] = - 2.5 * np.log10(kmags)

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

    return target_mags, target_err, scaled_cmags, scaled_kmags


def mag_plot(target_mags, target_err, date_obs, target, date, filters,
             dirtarget):
    """Plots and saves figure of magnitude and error as a function of time.
    
    Plots on the x-axis the time in Julian Days that each image was taken and 
    on the y-axis the magnitude value of the target star at that time. It also
    plots an error bar for each y-value, which is the error on magnitude 
    values. It then saves the figure to the directory containing FITS files 
    with accurate WCS information for each filter.
    
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
    date : str
        Date of observation (YYYY-MM-DD).
    filters : list
        List containing string of each filter keyword found in header of flat
        field and light frame images.
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
        
    Returns
    -------
    None
    """
    for fil in filters:
        fig1 = plt.gcf()
        plt.errorbar(date_obs, target_mags, yerr=target_err, fmt='o')
        plt.title('Light Curve of {}, {}'.format(target, date))
        plt.ylabel('Magnitude')
        plt.xlabel('JD')
        plt.gca().invert_yaxis()
        fig1.savefig(os.path.join(dirtarget, fil, 'WCS', 'accurate_WCS',
                                  'lightcurve.pdf'))


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
