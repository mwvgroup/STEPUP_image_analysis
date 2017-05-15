# STEPUP New Code:
  
  # Module variables retrieved from user input:
    dirstar = /home/depot/STEPUP/raw/<date>
    dirdark = /home/depot/STEPUP/raw/calibration/Dark/default
    target = <starname>
    date = <MM/DD/YYYY>

# Calibration:

  # Instrument-Signature Removal:
Overview - includes: get_star, get_flat, get_bias, get_dark functions, create_mbias, create_mflat, create_mdark, and reduction.

   # get_star, get_flat, get_bias, get_dark:
Each function loops through all files in dirstar (except get_dark, which uses files in dirdark) by opening up FITS header and looking at the keyword 'IMAGETYP'. Each function returns a 3D numpy array of flat, bias, dark, or target images. 

   # create_mbias, create_mdark, create_mflat:
Each function calls the get_<calibration_image> of their respective type of calibration image and assigns the 3D numpy array to a variable flats, biases, or darks. Performs necessary calculations (this includes taking the median along the 3rd axis) to create the master copies. Each function returns a 2D numpy array of this data.

  # reduction:
Calls create_mbias, create_mdark, create_mflat and assigns each of their outputs to the appropriate variable name. Calls get_star and assigns to variable name star. Opens HDU of dark images to get exposure time, assigns time to exptime in order to time correct master dark. Performs necessary calculations on star array. Returns 3D numpy array of instrument-signature-removed images. 

  # Saturation/Astrometry:
