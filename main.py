import os
import sys

computer = input("\nWork computer or Helena's computer? (W/H): ")
if computer == 'W':
    sys.path.insert(0, '/home/depot/STEPUP/STEPUP_image_analysis/ISR')
    import ISR
    sys.path.insert(0, '/home/depot/STEPUP/STEPUP_image_analysis/Calibration')
    import perform_astrometry
    import perform_photometry
if computer == 'H':
    sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/ISR')
    import ISR
    sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/Calibration')
    import perform_photometry
    import perform_astrometry

def main(verbose=False):
    target = input('\nInput target name: ')
    date = input('\nInput date of observation: ')
    dirtarget = None
    dirdark = None
    if computer == 'W':
        dirtarget = os.path.join('/home/depot/STEPUP/raw', target, date)
        dirdark = '/home/depot/STEPUP/raw/Calibration/Dark/Default'
    if computer == 'H':
        dirtarget = '/Users/helenarichie/tests2'
        dirdark = dirtarget
    analysis = input('\nWould you like to run the whole image analyis routine? (Y/N): ')
    if analysis == 'Y':
        print('\nInstrument signature removal in progress...')
        filters = ISR.ISR_main(dirtarget, dirdark, target)
        print('I\nnstrument signature removal completed.')
        answer = input('\nHave you saved a new-image.fits file to the appropriate directory? (Y/N): ')
        if answer == 'Y':
            print('\nAstrometry in progress...')
            dirtarget += '/ISR_Images'
            perform_astrometry.perform_astrometry(target, dirtarget, filters,
                                                  verbose=False)
        else:
            return None
        print('\nAstrometry completed.\nPhotometry in progress...')

        if computer == 'W':
            os.chdir(os.path.join('/home/depot/STEPUP/raw', target, date))
        if computer == 'H':
            os.chdir('/Users/helenarichie/tests2')
        ra = []
        dec = []
        comp_ra = []
        comp_dec = []
        comp_mag = []
        vsp_code = None
        rname = None
        ref_ra = []
        ref_dec = []
        cname = None
        check_ra = []
        check_dec = []

        with open('input-file.txt') as f:
            for line in f:
                if line.startswith('#RA='):
                    ra.append(line[4:].strip('\n'))
                if line.startswith('#DEC='):
                    dec.append(line[5:].strip('\n'))
                if line.startswith('#VSPCODE='):
                    vsp_code = line[9:].strip('\n')
                if line.startswith('#COMPMAGS='):
                    comp_mag = line[10:].strip('\n').split(',')
                if line.startswith('#COMPRA='):
                    comp_ra = line[8:].strip('\n').split(',')
                if line.startswith('#COMPDEC='):
                    comp_dec = line[9:].strip('\n').split(',')
                if line.startswith('#CLABEL='):
                    cname = line[8:].strip('\n')
                if line.startswith('#CRA='):
                    ref_ra.append(line[5:].strip('\n'))
                if line.startswith('#CDEC='):
                    ref_dec.append(line[6:].strip('\n'))
                if line.startswith('#KLABEL='):
                    kname = line[8:].strip('\n')
                if line.startswith('#KRA'):
                    check_ra.append(line[5:].strip('\n')) 
                if line.startswith('#KDEC='):
                    check_dec.append(line[6:].strip('\n'))
        comp_mags = []
        for mag in comp_mag:
            comp_mags.append(float(mag))

        coords = []
        coords.append(ra)
        coords.append(dec)

        for fil in filters:

            if computer == 'W':
                os.chdir(os.path.join('/home/depot/STEPUP/raw', target, date,
                                      'ISR_Images', fil, 'WCS/accurate_WCS'))
            if computer == 'H':
                os.chdir('/Users/helenarichie/tests2/ISR_Images/R/WCS/accurate_WCS')

            perform_photometry.perform_photometry(target, dirtarget, filters,
                                                  date, coords, comp_ra,
                                                  comp_dec, comp_mags, vsp_code,
                                                  rname, ref_ra, ref_dec, cname,
                                                  check_ra, check_dec,
                                                  verbose=False)

    else:
        cont_analysis = input('\nWould you still like to perform a function? (Y/N): ')
        while cont_analysis == 'Y':
            answer = input('\nWhich would you like to run? (ISR, ASTROM, PHOT): ')
            which_analysis(answer, dirtarget, dirdark, target, date, computer)
            cont_analysis = input('\nWould you still like to perform a function? (Y/N): ')
            
        
def which_analysis(answer, dirtarget, dirdark, target, date, computer):
    if answer == 'ISR':
        print('\nInstrument signature removal in progress...')
        filters = ISR.ISR_main(dirtarget, dirdark, target)
        print('\nInstrument signature removal completed.')

    if answer == 'ASTROM':
        im = input('\nHave you saved a new-image.fits file to the appropriate directory? (Y/N): ')
        if im == 'Y':
            filters = list(input('\nEnter filter of observation: '))
            print('\nAstrometry in progress...')
            dirtarget += '/ISR_Images'
            perform_astrometry.perform_astrometry(target, dirtarget,
                                                  filters, verbose=False)
            print('\nAstrometry completed.')
        else:
            return None
        
    if answer == 'PHOT':
        filters = list(input('\nEnter filter of observation: '))
        print('\nPhotometry in progress...')
        dirtarget += '/ISR_Images'
        path = None
        if computer == 'W':
            path = os.path.join('/home/depot/STEPUP/raw', target, date)
        if computer == 'H':
            path = '/Users/helenarichie/tests2'

        os.chdir(path)

        ra = []
        dec = []
        comp_ra = []
        comp_dec = []
        comp_mag = []
        vsp_code = None
        rname = None
        ref_ra = []
        ref_dec = []
        cname = None
        check_ra = []
        check_dec = []

        with open('input-file.txt') as f:
            for line in f:
                if line.startswith('#RA='):
                    ra.append(line[4:].strip('\n'))
                if line.startswith('#DEC='):
                    dec.append(line[5:].strip('\n'))
                if line.startswith('#VSPCODE='):
                    vsp_code = line[9:].strip('\n')
                if line.startswith('#COMPMAGS='):
                    comp_mag = line[10:].strip('\n').split(',')
                if line.startswith('#COMPRA='):
                    comp_ra = line[8:].strip('\n').split(',')
                if line.startswith('#COMPDEC='):
                    comp_dec = line[9:].strip('\n').split(',')
                if line.startswith('#CLABEL='):
                    cname = line[8:].strip('\n')
                if line.startswith('#CRA='):
                    ref_ra.append(line[5:].strip('\n'))
                if line.startswith('#CDEC='):
                    ref_dec.append(line[6:].strip('\n'))
                if line.startswith('#KLABEL='):
                    kname = line[8:].strip('\n')
                if line.startswith('#KRA'):
                    check_ra.append(line[5:].strip('\n')) 
                if line.startswith('#KDEC='):
                    check_dec.append(line[6:].strip('\n'))

        comp_mags = []
        for mag in comp_mag:
            comp_mags.append(float(mag))

        coords = []
        coords.append(ra)
        coords.append(dec)

        for fil in filters:
            path = None
            if computer == 'W':
                path = os.path.join('/home/depot/STEPUP/raw', target, date,
                                    'ISR_Images', fil, 'WCS/accurate_WCS')
            if computer == 'H':
                path = os.path.join('/Users/helenarichie/tests2/ISR_Images',
                                    fil, 'WCS/accurate_WCS')
            print(path)

            os.chdir(path)
            print(dirtarget)

            perform_photometry.perform_photometry(target, dirtarget,
                                                  filters, date, coords,
                                                  comp_ra, comp_dec,
                                                  comp_mags, vsp_code,
                                                  rname, ref_ra, ref_dec,
                                                  cname, check_ra,
                                                  check_dec, verbose=False)

main(verbose=False)
