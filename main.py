import sys
sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/ISR')
import ISR
sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/Calibration')
import add_WCS_info, perform_photometry

date = input('Input date of observation (YYYY-MM-DD): ')
target = input('Input target name: ')
dirtarget = '/Users/helenarichie/tests3'
dirdark = dirtarget
ra = input('Input RA of target (HH:MM:SS): ')
dec = input('Input Dec of target (+/-DD:MM:SS): ')
coords = (ra, dec)
ra_c = input('Input RA of comparison star (HH:MM:SS): ')
dec_c = input('Inpur Dec of comparison star (+/-DD:MM:SS): ')
comparison_coords = (ra_c, dec_c)
filters = ISR.ISR_main(dirtarget, dirdark, target)

dirtarget += '/ISR_Images'

print('\nInstrument signature removal completed.\nPhotometry in progress....\n')

add_WCS_info.add_WCS_info(dirtarget, filters, verbose=True)

(residual_aperture_sum,
 residual_aperture_sum_c) = perform_photometry.perform_photometry(filters,
                                                                  dirtarget,
                                                                  coords,
                                                                  comparions_coords)
