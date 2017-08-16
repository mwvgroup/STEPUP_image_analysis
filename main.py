import sys
sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/ISR')
import ISR
sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/Calibration')
import perform_astrometry
import perform_photometry


def main(verbose=False):
    date = input('Input date of observation (YYYY-MM-DD): ')
    target = input('Input target name: ')
    dirtarget = '/Users/helenarichie/tests2'
    dirdark = dirtarget
    ra = input('Input RA of target (HH:MM:SS): ')
    dec = input('Input Dec of target (+/-DD:MM:SS): ')
    coords = (ra, dec)
    ra_c = input('Input RA of comparison star (HH:MM:SS): ')
    dec_c = input('Inpur Dec of comparison star (+/-DD:MM:SS): ')
    comp_coords = (ra_c, dec_c)

    filters = ISR.ISR_main(dirtarget, dirdark, target)

    answer = print('\nInstrument signature removal completed.',
                   'Continue to astrometry (Y/N): ')
    answer = answer.upper()
    if answer == 'Y':
        print('\nAstrometry in progress...')
    else:
        return None

    dirtarget += '/ISR_Images'

    perform_astrometry.perform_astrometry(target, dirtarget, filters, verbose=False)

    answer = print('\nAstrometry completed.',
                   'Continue to photometry? (Y/N): ')
    answer = answer.upper()
    if answer == 'Y':
        print('\nAstrometry in progress...')
    else:
        return None

    perform_photometry.perform_photometry(target, dirtarget, filters, date, coords,
                                          comp_coords, comp_mag, verbose=False)
