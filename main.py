import sys
sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/ISR')
import ISR
sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/Calibration')
import add_WCS_info

date = input('Input date of observation (YYYY-MM-DD): ')
target = input('Input target name: ')
dirtarget = '/Users/helenarichie/tests2'
dirdark = dirtarget
coords = input('Input RA and Dec of target :')

filters = ISR.ISR_main(dirtarget, dirdark, target)

dirtarget += '/ISR_Images'

print('\nInstrument signature removal completed.\nPhotometry in progress....\n')

add_WCS_info.add_WCS_info(dirtarget, filters)
