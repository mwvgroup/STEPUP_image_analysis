import sys
sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/ISR_main')
import ISR_main
sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/Calibration')
import mmm

dirtarget = dirdark = '/Users/helenarichie/yo'
target = 'HAT-P-3'

def main(dirtaget, dirdark, target):
    science_images = ISR_main.ISR_main(dirtarget, dirdark, target)

    (skymod, sigma, skew) = mmm.mmm(science_images, minsky=20, highbad=False)

    return(skymod, sigma, skew)

(skymod, sigma, skew) = main(dirtarget, dirdark, target)
