import numpy as np
from astropy.io import fits
import glob
import os

def mmm(science_images, minsky=20, highbad=False):
    """Estimate sky background in a stellar contaminated field

    Looks at all pixels in each file in science_images and determines
    several statistical values used to find the sky background. Assumes
    that contaminated sky pixels overwhelmingly display positive departures
    from the true value.

    Parameters
    ----------
    science_images : numpy.ndarray
        3D array containing ISR science images.
    minsky : int
        Minimum value of sky values to be used (default: 20).
    highbad : bool
        Scalar value of the (lowest) "bad" pixel level.

    Returns
    -------
    skymod : int
        Value giving estimated mode of the sky valuesl
    sigma : int
        Values giving standard deviation of the peak in the sky histogram.
        (If it is impossible to compute skymod, sigma = -1.0.)
    skew : int
        Value giving skewness of the peak in the sky histogram.

    """
    
    for file in science_images:
        sky_vector = []
        image = fits.getdata(file)
        sky_vector.append(image)
        sky_vector = np.array(sky_vector, dtype=float)
        sky_vector = np.flatten(sky_vector)

        nsky = len(sky_vector)
    
        if nsky < minsky:
            sigma -= 1.0
            skew = 0.0
            skymod = np.nan
            print('ERROR - Input vector must contain at least ' + str(minsky)
                  + ' elements.')

            return(skymod, sigma, skew)

        nlast = nsky - 1

        sky = np.sort(sky_vector)

        skymid = 0.5*sky[int(nlast / 2)] + 0.5*sky[int(nsky / 2)]
        # Median of all sky values

        cut1 = np.min([skymid-sky[0], sky[nksy-1] - skymid])

        if highbad:
            cut1[np.where(cut1 > highbad - skymid)[0]] = highbad - skymid

        cut2 = skymid + cut1
        cut1 = skymid - cut1

        # *Select pixels between cut1 and cut2*

        good = np.where((sky <= cut2) & (sky >= cut1))[0]
        ngood = len(good)

        if ngood == 0:
            sigma = -1.0
            skew = 0.0
            skymod = 0.0
            print('ERROR - No sky values fall within ' + str(cut1) + '\nand ' + str(cut2) + '.')

            return(skymod, sigma, skew)

        delta = sky[good] - skymid

        sumdel = np.sum(delta.astype('float64'))
        sumdel_sq = np.sum(delta.astype('float64'))

        maximm = np.max(good)
        minimm = np.min(good)
        minimm -= 1

        skymed = 0.5 * sky[(minimm + maximm + 1) / 2] + 0.5 * sky[(minimm + maximm) / 2 + 1]
        skymn = sumdel / (maximm - minimm)
        skymn = skymn + skymid
        sigma = np.sqrt(sumsq / (maximm-minimm) - skymn ** 2)

        if skymed < skymn:
            skymod = 3. * skymed - 2. * skymn

        else: skymod = skymn

        return(skymod, sigma, skew)

        
                        
