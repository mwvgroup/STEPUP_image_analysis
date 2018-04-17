import os
import sys
sys.path.insert(0, '/home/depot/STEPUP/STEPUP_image_analysis/ISR')
# sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/ISR')
import ISR
sys.path.insert(0, '/home/depot/STEPUP/STEPUP_image_analysis/Calibration')
# sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/Calibration')
import perform_astrometry
import perform_photometry


def main(verbose=False):
    target = input('Input target name: ')
    date = input('Input date of observation: ')
    dirtarget = os.path.join('/home/depot/STEPUP/raw', target, date)
    dirdark = '/home/depot/STEPUP/raw/Calibration/Dark/Default'
    # dirtarget = '/Users/helenarichie/tests2'
    # dirdark = dirtarget

    # filters = ISR.ISR_main(dirtarget, dirdark, target)
    # filters = ["R"]

    answer = input('\nInstrument signature removal completed.\nContinue to astrometry (Y/N): ')
    if answer == 'Y':
        im = input('\nHave you saved a new-image.fits file to the appropriate directory? (Y/N): ')
        if im == 'Y':
            print('\nAstrometry in progress...')
        
    else:
        return None

    dirtarget += '/ISR_Images'

    perform_astrometry.perform_astrometry(target, dirtarget, filters, verbose=False)

    answer = input('\nAstrometry completed.\nContinue to photometry? (Y/N): ')
    
    if answer == 'Y':
        if im == 'Y':
            print('\nPhotometry in progress...')
    else:
        return None

    os.chdir(os.path.join('/home/depot/STEPUP/raw', target, date))
    # os.chdir(os.path.join('/Users/helenarichie/tests2'))

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

    os.chdir(os.path.join('/home/depot/STEPUP/raw/WDra/2018-01-30/ISR_Images/R/WCS/accurate_WCS'))

    # os.chdir(os.path.join("/Users/helenarichie/tests2/ISR_Images/r'/WCS/accurate_WCS"))
    
    perform_photometry.perform_photometry(target, dirtarget, filters, date, coords,
                                          comp_ra, comp_dec, comp_mags, vsp_code,
                                          rname, ref_ra, ref_dec, cname, check_ra,
                                          check_dec, verbose=False)

main(verbose=True)
