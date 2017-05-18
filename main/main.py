#from astropy.io import fits
import sys
sys.path.insert(0, 'directory where code is stored')
import get_calibimages 
import get_science_images
import create_mcalib
import calibtarget

date = input('Enter the date of observation (MM/DD/YYYY): ')
target = input('Enter the name of the target: ')
dirtarget = '/home/depot/STEPUP/raw/' + date
dirdark = '/home/depot/STEPUP/raw/calibration/Dark/default'

def main(dirtarget, dirdark):
    biases = get_calibimages.get_biases(dirtarget)
    flats = get_calibimages.get_flats(dirtarget)
    darks = get_calibimages.get_darks(dirdark)
    target = get_science_images.get_science_images(dirtarget)

    
    #hdulist_dark = fits.open(dark)
    #dark_exptime = hdulist_dark[0].header['EXPTIME']
    #this part is not currently working
    
    mbias = create_mcalib.create_mbias(biases)
    mdark = create_mcalib.create_mdark(darks, mbias)
    mflat = create_mcalib.create_mflat(flats, mbias, mdark)

    calibrated_target = calibtarget.calibtarget(target, mbias, mdark, mflat, dark_exptime)
    
    pass

    
    
    
