# STEPUP New Code:
  
  # Module variables:
    dirtarget = str, /home/depot/STEPUP/raw/<date>, retrieved from user input
    dirdark = str, /home/depot/STEPUP/raw/calibration/Dark/default
    target = str, <starname>, retrieved from user input
    date = str, <MM/DD/YYYY>, retrieved user input
    
# ISR_main:

#   ISR.py

#       ISR_main:
Creates ISR FITS files by executing a preliminary calibration sequence. Imports and calls get_unfiltered_calibimages, get_filtered_calibimages, and instrument_signature_removal. Saves ISR science images to dirtarget/ISR_Images/<filter-name>.

Paramaters: dirtarget, dirdark, target

Returns: None

#       get_unfiltered_calibimages:
Creates and saves master bias and dark by searching in dirtarget for all bias frames and creates master bias by taking the median of the array of those images. It then saves the master bias todirtarget/mcalib. It then searches in dirdark for all dark frames and subtracts the mbias and time-corrects each image. It then creates the master dark by taking the median of the array of those images. It then saves the master dark to dirtarget/mcalib.

Paramaters: dirtarget, dirdark

Returns: exptime

#       get_filtered_calibimages:
Creates and saves master flat for each filter. Creates a list of all unique filter names on flat and light images. It then loops through that list looking for all flatswith the same filter name for each iteration of the loop. It then subtracts the bias and normalizes each individual flat and averages that array. The master flat is then saved as "<filter-name>_mflat.fits" to dirtarget/mcalib.

Parameters: dirtarget

Returns: image_filters

#       instrument_signature_removal:
Removes instrument signatures from raw science images. Retrieves all FITS files with IMAGETYP keyword, "Light Frame" and all master calibration images and puts them into arrays. It then calculates the expected saturation, which is later added to the header. It then removes the master bias, dark, and flat by filer from each image in the light frame array and saves it as a FITS file with the light frame header as dirtarget/ ISR_Images/<filter-name>/<target-name>_<filter-keyword>_*.fits

Parameters: dirtarget, target, exptime, image_filters

Returns: None
  
# Calibration_main:
  
#       Calibration_main:
  
   # Calibration:
  
#       mmm:
Estimates sky background in stellar contaminated field.
 
 Parameters: science_images, minsky, highbad
 
 Returns: skymod, sigma, skew
