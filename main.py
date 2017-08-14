import sys
sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/ISR')
import ISR
sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/Calibration')
import add_WCS_info

date = input('Input date of observation (YYYY-MM-DD): ')
target = input('Input target name: ')
dirtarget = '/Users/helenarichie/tests2'
dirdark = dirtarget
ra = input('Input RA of target (HH:MM:SS): ')
dec = input('Input Dec of target (+/-DD:MM:SS): ')
coords = (ra, dec)
ra_c = input('Input RA of comparison star (HH:MM:SS): ')
dec_c = input('Inpur Dec of comparison star (+/-DD:MM:SS): ')
comparison_coords = (ra_c, dec_c)

filters = ISR.ISR_main(dirtarget, dirdark, target)

answer = print('\nInstrument signature removal completed. Continue? (Y/N): ')
if answer == 'Y':
    print('\nPhotometry in progress...')
else:
    return None

dirtarget += '/ISR_Images'

add_WCS_info.add_WCS_info(target, dirtarget, filters, verbose=True)

perform_photometry.perform_photometry(target, dirtarget, filters, date, coords,
                                      comparison_coords, comparison_mag,
                                      verbose=False)
