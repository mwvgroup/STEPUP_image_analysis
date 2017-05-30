import numpy as np

def mmm(science_images, highbad = False, minsky = 20, maxiter = 50,
        readnoise = False, nsky = False, integer = "discrete", debug = False):
    """Estimate sky background in a stellar contaminated field.
    
    MMM assumes that contaminated sky pixel values overwhelminlgy
    display positive departures from the true value.
    
    Parameters
    ----------
    science_images : numpy array
        3D array containing ISR science images. 
    highbad
        (Optional) Scalar value of the (lowest) "bad" pizel level (e.g. 
        cosmic rays or saturated pixels). If not supplied, then there is
        assumed to be no bad pixels. 
    minsky : int
        (Optional) Minimum number of sky values to be used. MMM will 
        return an error if fewer sky elements are supplied. 
        Default = 20.
    maxiter : int
        (Optional) Maximum number of iterations allowed, default = 50.
    readnoise 
        (Optional) The read noise (or minimum noise for any pixel). 
        Normally, mmm determines the (robust) median by averaging the 
        central 20% of the sky values. In some cases where the noise is 
        low, and pizel values are quantized, a larger fraction may be 
        needed. By supplying the optional read noise parameter, mmm 
        is better able to adjust the fraction of pixels used to determine
        the median.
    integer : str
        (Optional) Set this keyword if the input science_images array only 
        contains discrete integer values. This keyword is only needed if 
        the array is of type float, but contains only discrete integer 
        values. 
    debug 
        If this keyword is set and non-zero, then additional information is
        displayed at the terminal.
        
    Returns
    -------
    skymod 
        Scalar giving estimated mode of the sky values.
    sigma 
        Scalar giving standard deviation of the peak in 
        the sky histogram. If for some reason it is impossible
        to derive skymore, then sigma = -1.0
    skew : Scalar giving skewness of the peak in the sky histogram.
    """
 
    if np.nan: science_images = science_images[np.where(science_images == 
                                                     science_images)]
    nsky = len(science_images)
    
    if nsky < minsky:
        sigma -= 1.0
        skew = 0.0
        skymod = np.nan
        print('ERROR - Input array must contain at least ' + str(minksy) + 
              ' elements.')
        return(skymod, sigma, skew)
    
    nlast = nsky - 1
    
    if debug:
        print('Processing ' + str(nsky) + ' element array.')
        
    sz_sky = np.shape(science_images)
    
    sky = np.sort(science_images)  # Sort in ascending values
    
    skymid = 0.5 * sky[int((nsky - 1)/2)] + 0.5 * sky[int(nsky/2)] # Median value (per pixel?)
    
    cut1 = np.min([skymid - sky[0], sky[nsky - 1] - skymid])
    
    if highbad:
        cut1[np.where(cut1 > highbad - skymid)[0]] = highbad - skymid
    cut2 = skymid + cut1  # this creates a range of acceptable sky values
    cut1 = skymid - cut1
        
    good = np.where((sky <= cut2) & (sky >= cut1))[0]
    ngood = len(good)
    
    if ngood == 0:
        sigma -= 1.0
        skew = 0.0
        skymod = 0.0
        print('ERROR - No sky values fall within ' + str(cut1) + \
              ' and ' + str(cut2))
        return(skymod, sigma, skew)
    
    delta = sky[good] - skymid # deviation of sky from median
    sumd = np.sum(delta.astype('float64'))        
    sumdsq = np.sum(delta.astype('float64')**2)

    maximm = np.max(good) # Highest value accepted at upper end of vector
    minimm = np.min(good) 
    minimm -= 1 # Highest value rejected at lower end of vector
    
    skymed = 0.5 * sky[(minimm + maximm + 1)/2] + 0.5 * sky[(minimm + maximm)/2 +1] # Median
    skymn = sumd / (maximm - minimm) - skumn**2 # Mean
    sigma = np.sqrt(sumdsq / (maximm - minimm) - skymn**2) # Sigma
    skymn = skymn + skymid # Add median which was subtracted off earlier

    # If mean is less than mode, then the contamination is lsight, and the
    # mean value is what we really want.
    # skymod = (skymed < skymn) ? 3.*skymed - 2.*skymn : skymn
    
    if skymed < skymn:
        skymod = 3. * skymed - 2. * skymn
    
    else: 
        skymod = skymn
