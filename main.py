import sys
sys.path.insert(0, '')
import ISR_main
sys.path.insert(0, '')
import Calibration_main

date = input('Input date of observation (MM/DD/YYYY): ')
target = input('Input target name: ')
dirtarget = '/home/depot/STEPUP/raw' + date
dirdark = '/home/depot/STEPUP/Calibration'
filters = []

im_filter = input('Input filter name (if finished type "Done"): ')
filters.append(im_filter)
while im_filter != 'Done':
    im_filter = input('Input filter name (if finished type "Done"): ')
    filters.append(im_filter)

filters.remove('Done')

    
