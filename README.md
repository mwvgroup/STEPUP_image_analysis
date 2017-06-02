# STEPUP New Code:
  
  # Module variables:
    dirtarget = str, /home/depot/STEPUP/raw/<date>, retrieved from user input
    dirdark = str, /home/depot/STEPUP/raw/calibration/Dark/default
    target = str, <starname>, retrieved from user input
    date = str, <MM/DD/YYYY>, retrieved user input

# ISR_main:
  
#       ISR_main:
Overview - Imports and uses modules get_calibimages, create_mcalib, and instrument_signature_removal. Saves ISR reduced images to new directory.

Parameters: dirtarget, dirdark

Returns: science_images
    
  # ISR:

#       get_mcalib:
Loops through all FITS files in dirtarget. If header index 'IMAGETYP' == 'Flat Field', file is appended to flats and 'Bias Frame', to biases. Then loops through dirdark and does the same thing for 'Dark Frame' image types. 

Paramaters: dirtarget, dirdark 

Returns: mbias, mdark, mfalt, 3D numpy array of bias, dark, or flat images. 

#       create_mcalib
 Module containing functions to create master bias, dark, and flat images.
        
        create_mbias:
Takes median of darks along first axis of bias images.
 
Parameters: biases
  
Returns: mbias
  
        create_mdark:
Takes median of darks along first axis of dark images. Then subtracts mbias.

Parameters: darks, mbias

Returns: mdark
  
          create_mflat:
Takes median of flats along first axis. Then subtracts mbias and mdark. Then divides by the median, taken along the first axis, of flats.

Parameters: flats, mbias, mdark

Returns: mflat

#       instrument_signature_removal:
Loops through FITS file in dirtarget for 'Light Frame' image types and reduces them using mbias, time-corrected mdark, and mflat. Then adds expected saturation to header. Saves instrument signature removed FITS files to new directory.

Paramaters: target, mbias, mdark, mflat, dark_exptime, and target. 

Returns: isr_scimages
  
# Calibration_main:
  
#       Calibration_main:
  
   # Calibration:
  
#       mmm:
Estimates sky background in stellar contaminated field.
 
 Parameters: science_images, minsky, highbad
 
 Returns: skymod, sigma, skew
