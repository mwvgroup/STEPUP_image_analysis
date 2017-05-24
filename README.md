# STEPUP New Code:
  
  # Module variables:
    dirtarget = str, /home/depot/STEPUP/raw/<date>, retrieved from user input
    dirdark = str, /home/depot/STEPUP/raw/calibration/Dark/default
    target = str, <starname>, retrieved from user input
    date = str, <MM/DD/YYYY>, retrieved user input
    dark_exptime = float, exposure time for darks, retrieved by calling in FITS header of dark images.
    mbias = numpy array, 2D array containing master flat image, retrieved by calling create_mbias
    mdark = numpy array, 2D array containing master dark image, retrieved by calling create_mdark
    mflat = numpy array, 2D array containing master flat image, retrieved by calling create_mflat

# ISR_main:
  
      main:
Overview - Imports and uses modules get_calibimages, create_mcalib, and instrument_signature_removal. Saves ISR reduced images to new directory.

Parameters: dirtarget, dirdark

Returns: None
    
# ISR (Instrument Signature Removal):
     
      get_mcalib:
Loops through all FITS files in dirtarget. If header index 'IMAGETYP' == 'Flat Field', file is appended to flats and 'Bias Frame', to biases. Then loops through dirdark and does the same thing for 'Dark Frame' image types. 

Paramaters: dirtarget, dirdark 

Returns: 3D numpy array of flat, bias, or dark images. 

      create_mcalib :: create_mbias, create_mdark, create_mflat:
Each function has their respective type of calibration image passed in as a parameter. Performs necessary calculations (this includes taking the median along the first axis) to create the master images. Each function returns a 2D numpy array of this data.

Parameters: biases, darks, flats, mbias, mdark

Returns: 2D arrays of master calibration images.

      instrument_signature_removal:
Loops through FITS file in dirtarget for 'Light Frame' image types and reduces them using mbias, time-corrected mdark, and mflat. Then adds expected saturation to header. Saves instrument signature removed FITS files to new directory.

Paramaters: target, mbias, mdark, mflat, dark_exptime, and target. 

Returns: None
