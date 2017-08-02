import sys
sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/ISR_main')
import ISR
sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/Calibration')
import perform_astrometry

date = input('Input date of observation (YYYY-MM-DD): ')
target = input('Input target name: ')
dirtarget = '/Users/helenarichie/tests'
dirdark = dirtarget

filters = ISR.ISR_main(dirtarget, dirdark, target)

print('\nInstrument signature removal completed.\nAstrometry in progress....\n')

add_WCS_info.add_WCS_info(dirtarget, filters)
