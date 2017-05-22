import sys
sys.path.insert(0, 'directory where code is stored')
import get_calibimages 
import create_mcalib
import isr_removal

date = input('Enter the date of observation (MM/DD/YYYY): ')
target = input('Enter the name of the target: ')
dirtarget = '/home/depot/STEPUP/raw/' + date
dirdark = '/home/depot/STEPUP/raw/calibration/Dark/default'

def main(dirtarget, dirdark):
    
    (biases, darks, flats, dark_exptime) = get_calibimages.get_calibimages(dirtarget, dirdark)
    
    mbias = create_mcalib.create_mbias(biases)
    mdark = create_mcalib.create_mdark(darks, mbias)
    mflat = create_mcalib.create_mflat(flats, mbias, mdark)

    isr_removal.isr_removal(dirtarget, mbias, mdark, mflat, dark_exptime)
    
main(dirtarget, dirdark)
    
    
    
