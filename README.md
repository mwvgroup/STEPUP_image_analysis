# STEPUP New Code:
  
  # Module variables:
    dirtarget = str, /home/depot/STEPUP/raw/<date>, retrieved from user input
    dirdark = str, /home/depot/STEPUP/raw/calibration/Dark/default, retrieved from user input
    target = str, <starname>, retrieved from user input
    date = str, <MM/DD/YYYY>, retrieved from program
    dark_exptime = float, exposure time for darks, retrieved by calling get_darks and looking in FITS header
    mbias = numpy array, 2D array containing master flat image, retrieved by calling create_mbias
    mdark = numpy array, 2D array containing master dark image, retrieved by calling create_mdark
    mflat = numpy array, 2D array containing master flat image, retrieved by calling create_mflat

# main:
  
      main:
Overview - Imports and uses modules get_calibimages, get_science_images, create_mcalib, and target_reduction.

Parameters: dirtarget, dirdark

Returns: will alwasy return reduced_target, but can return others if needed.

# Calibration:

  # Instrument-Signature Removal:
Overview - includes: get_star, get_flat, get_bias, get_dark functions, create_mbias, create_mflat, create_mdark, and target_reduction.

      get_star, get_flat, get_bias get_dark:
Each function loops through all files in dirtarget (except get_dark, which uses files in dirdark) by opening up FITS header and looking at the keyword 'IMAGETYP'. 

Paramaters: dirtarget or dirdark 

Returns: 3D numpy array of flat, bias, dark, or target images. 

      create_mbias, create_mdark, create_mflat:
Each function calls the get_<calibration_image> of their respective type and assigns the 3D numpy array to a variable flats, biases, or darks. Performs necessary calculations (this includes taking the median along the 3rd axis) to create the master images. Each function returns a 2D numpy array of this data. Each function will be assigned to module variable mbias, mdark, and mflat that will be called by target_reduction and calibration.

Parameters: biases, darks, flats, mbias, mdark

Returns: 2D arrays of master calibration images.

      target_reduction:
Performs necessary calculations on star array using mbias, mdark, and mflat.

Paramaters: target, mbias, mdark, mflat, and dark_exptime. 

Returns: 3D numpy array of instrument-signature-removed images.

  # Saturation/Astrometry:
Overview - includes: calibration
    
      calibration:
Subtracts the medians of all master calibration images from the initial saturation level. Then, for each image in calibstar, opens up primary HDU and adds new 'SATLEVEL' to header. 

Parameters: mbias, mdark, mflat, and dark_exptime. 

Returns: 3D numpy array of reduced images whose headers include saturation level.
