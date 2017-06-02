import sys
sys.path.insert(0, '')
import ISR_main
sys.path.insert(0, '')
import Calibration_main

date = input('Input date of observation: ')
target = input('Input target name: ')
dirtarget = '/home/depot/STEPUP/raw' + date
dirdark = '/home/depot/STEPUP/Calibration'

def main(dirtaget, dirdark, target):
    science_images = ISR_main.ISR_main(dirtarget, dirdark, target)

    
