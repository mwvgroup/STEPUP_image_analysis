from astropy.io import fits
import numpy as np
import os
import glob

date = input('Enter date of observation (MM/DD/YYYY): ')
target = input('Enter name of target: ')
dirdark = '/home/depot/STEPUP/raw/calibration/Dark/default'
dirstar = '/home/depot/STEPUP/raw/' + date


def get_images(dirstar, dirdark):
    """Retrieves all flat, bias, dark, and target images.

    Calls get_biases, get_flats, get_darks, and get_star. Each function (except
    get_darks, which looks in dirdark) looks at all FITS file headers at the
    "IMAGETYP" keyword and adds each type to its respective array. It then
    returns an array of bias, dark, flat, and target images in a tuple of
    that order.

    Parameters
    ----------
    dirstar : str
        Directory in which bias, flat, and target images are stored.
    dirdark : str
        Directory in which dark images are stored.

    Returns
    -------
    biases : numpy array
        3D numpy array containing all bias images found in dirstar.
    darks : numpy array
        3D numpy array containing all dark images found in dirdark.
    flats : numpy array
        3D numpy array containing all flat images found in dirstar.
    star : numpy array
        3D numpy array containing all target images found in dirstar.
    """
    
    def get_biases(dirstar):
        """Retrieves all bias images from dirstar.
    
        Searches in primary HDU headers of all files in dirstar for files with
        keyword "IMAGETYP" that point to "Bias Frame". It then puts all biases
        found into an array which is returned to the caller.
    
        Parameters
        ----------
        dirstar : string
            module variable which is the directory that
            image, flat, and bias files are saved in
            /home/depot/STEPUP/raw/<date>
    
        Returns
        -------
        biases : numpy array
            3D array of all bias images stored in dirstar
        """
        
        files = glob.glob(os.path.join(dirstar, '*.fit'))
        biases = []
        for file in files:
            hdulist = fits.open(file)
            if hdulist[0].header['IMAGETYP'] == 'Bias Frame':
                image = fits.getdata(file)
                biases.append(image)
            hdulist.close()
        biases = np.array(biases, dtype = float)
        return biases
    
    def get_darks(dirdark):
        """Retrieves all dark images from dirdark.
    
        Searches in primary HDU headers of all files in dirdark for files with
        keyword "IMAGETYP" that point to "Dark Frame". It then puts all darks
        found into an array which is returned to the caller.
    
        Parameters
        ----------
        dirdark : str
            Module variable that is the directory in which all darks are
            stored. /home/depot/STEPUP/raw/caliration/Dark/default
    
        Returns
        -------
        darks : numpy array
            3D numpy array of all dark images located within dirdark
        """
        
        files = glob.glob(os.path.join(dirdark, '*.fit'))
        darks = []
        for file in files:
            hdulist = fits.open(file)
            if hdulist[0].header['IMAGETYP'] == 'Dark Frame':
                image = fits.getdata(file)
                darks.append(image)
            hdulist.close()
        darks = np.array(darks, dtype = float)
        return darks
    
    
    def get_flats(dirstar):
        """Retrieves all flat files from dirstar.
    
        Searches in primary HDU headers of all files in dirstar for files with
        keyword "IMAGETYP" that point to "Flat Frame". It then puts all flats
        found into an array which is returned to the caller.
    
        Parameters
        ----------
        dirstar : str
            Module variable that is the directory in which all flat, bias, and
            target images are located. /home/depot/STEPUP/raw/<date>
    
        Returns
        -------
        flats : numpy array
            3D array of all flat images stored in dirstar
    
        """
        
        files = glob.glob(os.path.join(dirstar, '*.fit'))
        flats = []
        for file in files:
            hdulist = fits.open(file)
            if hdulist[0].header['IMAGETYP'] == 'Flat Field':
                image = fits.getdata(file)
                flats.append(image)
            hdulist.close()
        flats = np.array(flats, dtype = float)
        return flats
    
    def get_star(dirstar):
        """Retrieves all raw target images from dirstar.
    
        Searches in primary HDU headers of all files in dirstar for files with
        keyword "IMAGETYP" that point to "Light Frame". It then puts all target
        images found into an array which is returned to the caller.
    
        Parameters
        ----------
        dirstar : str
            Module variable that is the directory in which all flat, bias, and
            target images are located. /home/depot/STEPUP/raw/<date>
    
        Returns
        -------
        star : numpy array
            3D array of all raw target images stored in dirstar
        """
        
        files = glob.glob(os.path.join(dirstar, '*.fit'))
        star = []
        for file in files:
            hdulist = fits.open(file)
            if hdulist[0].header['IMAGETYP'] == 'Light Frame':
                image = fits.getdata(file)
                star.append(image)
            hdulist.close()
        star = np.array(star, dtype = float)
        return star
    
    biases = get_biases(dirstar)
    darks = get_darks(dirdark)
    flats = get_flats(dirstar)
    star = get_star(dirstar)
    
    return (biases, darks, flats, star)

(biases, darks, flats, star) = get_images(dirstar, dirdark)
