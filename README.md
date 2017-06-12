# STEPUP New Code:
  
  # Module variables:
    dirtarget = str, /home/depot/STEPUP/raw/<date>, retrieved from user input
    dirdark = str, /home/depot/STEPUP/raw/calibration/Dark/default
    target = str, <starname>, retrieved from user input
    date = str, <MM/DD/YYYY>, retrieved user input

# ISR_main:
  
#       ISR_main:
Overview - Imports and uses modules get_calibimages, create_mcalib, and instrument_signature_removal. Saves ISR reduced images to new directory.

Parameters: dirtarget, dirdark, filters

Returns: science_images
    
  # ISR:

#       get_mcalib:
Looks at headers of all FITS files in dirtarget and dirdark. Returns arrays of each type of calibration image. Flat images are separated by filters into their own arrays.

Paramaters: dirtarget, dirdark, target, filters

Returns: biases, darks, r_flats, b_flats, v_flats, dark_exptime, exptime, bias_prihdr, dark_prihdr, r_flat_prihdr, b_flat_prihdr, v_flat_prihdr, r_filter, b_filter, v_filter 

      # Currently working on restructuring the filter separation process by creating a filter object.

#       create_mcalib
 Module containing functions to create master bias, dark, and flat images.
        
        create_mbias:
Takes median of darks along first axis of bias images. Then saves mbias.fits to dirtarget in "mcalib" folder.
 
Parameters: biases, bias_prihdr, dirtarget
  
Returns: mbias
  
        create_mdark:
Subtracts mbias from and time-corrects each dark image. Then takes median along first axis of all dark images. Then saves mdark.fits to dirtarget in "mcalib" folder.

Parameters: darks, mbias, dark_prihdr, dirtarget, dark_exptime, exptime

Returns: mdark
  
        create_mflat:
Subtracts mbias from and normalizes each flat image. Then takes average along first axis of all flat images. Then saves <filter>_mflat.fits to dirtarget in "mcalib" folder. (This algorithm is executed for each filter of flat images.)

Parameters: r_flats, b_flats, v_flats, mbias, r_flat_prihdr, b_flat_prihdr, v_flat_prihdr, r_filter, b_filter, v_filter, dirtarget

Returns: r_mflat, b_mflat, v_mflat

#       instrument_signature_removal:
Loops through FITS file in dirtarget for 'Light Frame' image types and reduces them using mbias, mdark, and <filter>_mflat for each different filter. Then adds expected saturation to header. Saves instrument signature removed FITS files to new directory.

Paramaters: irtarget, target, mbias, mdark, r_mflat, b_mflat, v_mflat, dark_exptime, exptime, r_filter, b_filter, v_filter

Returns: r_isr_scimages, b_isr_scimages, v_isr_scimages
  
# Calibration_main:
  
#       Calibration_main:
  
   # Calibration:
  
#       mmm:
Estimates sky background in stellar contaminated field.
 
 Parameters: science_images, minsky, highbad
 
 Returns: skymod, sigma, skew
