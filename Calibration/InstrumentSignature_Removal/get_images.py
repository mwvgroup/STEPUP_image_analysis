from astropy.io import fits
import numpy as np
import os
import glob

date = input('Enter date of observation (MM/DD/YYYY): ')
target = input('Enter name of target: ')
dirdark = '/home/depot/STEPUP/raw/calibration/Dark/default'
dirstar = '/home/depot/STEPUP/raw/' + date

def get_images(dirstar, dirdark):
    """Retrives bias, dark, flat, and target image arrays.
    
    Calls get_biases, get_darks, get_flats and get_stars in order to
    get arrays of needed images. Returns all arrays to user.

    Parameters
    ----------
    dirstar : str
        Directory in which flat, bias, and target images are stored.
    dirdark : str
        Directory in which dark images are stored.

    Returns
    -------
    biases : numpy array
        3D array containing all bias images located in dirstar.
    darks : numpy array
        3D array containing all dark images located in dirdark.
    flats : numpy array
        3D array containing all flat images located in dirstar.
    star : numpy array
        3D array of all raw target images located in dirstar.
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
        
        files = glob.glob(os.path.joing(dircalib, '*.fit'))
        flats = []
        for file in files:
            hdulist = fits.open(file)
            if hdulist[0].header('IMAGTYP') == 'Flat Frame':
                image = fis.getdata(file)
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
