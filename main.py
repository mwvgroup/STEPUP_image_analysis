import sys
sys.path.insert(0, '')
import ISR
sys.path.insert(0, '')
import Calibration_main

date = input('Input date of observation (MM/DD/YYYY): ')
target = input('Input target name: ')
dirtarget = '/home/depot/STEPUP/raw' + date
dirdark = '/home/depot/STEPUP/Calibration'

ISR.ISR_main(dirtarget, dirdark, target)
