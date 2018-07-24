import os
import sys


# Determine if code is being ran from Helena's personal computer or a computer
# in the Astro Lab and import functions ISR, perform_astrometry, and perform_
# photometry to the corresponding directories.
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
    """Runs image analysis routine or specified steps of routine on dataset.

    Determines what steps of STEPUP_image_analysis the user would like to run.
    If the user specifies that they would like to run the entire routine, the
    code runs without pausing, except for to check that a "new-image.fits" file
    has been saved to ISR_Images. Otherwise, the user can run individual parts
    of the routine as many times as they would like (as long as the individual
    functions allow to be repeatdely ran).

    Parameters
    ----------
    verbose : Boolean
        Specifies whether or not user would like to print more information
        about the status of the code.

    Returns
    -------
    None
    """
    # Retrieves target name and date of observation from user input and
    # initializes dirtarget and dirdark variables.
    target = input('\nInput target name: ')
    date = input('\nInput date of observation: ')
    dirtarget = None
    dirdark = None

    # Assigns variables according to which computer the code is being ran from.
    if computer == 'W':
        dirtarget = os.path.join('/home/depot/STEPUP/raw', target, date)
        dirdark = '/home/depot/STEPUP/raw/Calibration/Dark/Default'
    if computer == 'H':
        dirtarget = '/Users/helenarichie/tests2'
        dirdark = dirtarget

    # Determines if user would like to run entire image analysis routine.
    analysis = input('\nWould you like to run the whole image analyis routine? (Y/N): ')
    if analysis == 'Y':
        print('\nInstrument signature removal in progress...')
        # Removes instrument signatures from dataset.
        filters = ISR.ISR_main(dirtarget, dirdark, target)
        print('\nInstrument signature removal completed.')
        # Determines if user has saved new-image.fits WCS calibration file to
        # ISR_Images directory that was created in ISR function.
        answer = input('\nHave you saved a new-image.fits file to the appropriate directory? (Y/N): ')
        if answer == 'Y':
            print('\nAstrometry in progress...')
            dirtarget += '/ISR_Images'
            # Calculates WCS information for dataset.
            perform_astrometry.perform_astrometry(target, dirtarget, filters,
                                                  verbose=False)
        else:
            return None
        print('\nAstrometry completed.\nPhotometry in progress...')
        # Changes current working directory according to which computer the
        # code is being ran from.
        if computer == 'W':
            os.chdir(os.path.join('/home/depot/STEPUP/raw', target, date))
        if computer == 'H':
            os.chdir('/Users/helenarichie/tests2')

        # Initialize variables to obtain information from input-file.txt file.
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

        # Read in information from input-file.txt and append/assign it to
        # corresponding list/variable.
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
                    check_ra.append(line[5:].strip('\n'))
                if line.startswith('#CDEC='):
                    check_dec.append(line[6:].strip('\n'))
                if line.startswith('#RLABEL='):
                    rname = line[8:].strip('\n')
                if line.startswith('#RRA'):
                    ref_ra.append(line[5:].strip('\n'))
                if line.startswith('#RDEC='):
                    ref_dec.append(line[6:].strip('\n'))

        # Ensure that all magntidue values are floats.
        comp_mags = []
        for mag in comp_mag:
            comp_mags.append(float(mag))

        # Create list of two different lists of either right ascensions or
        # declinations.
        coords = []
        coords.append(ra)
        coords.append(dec)

        for fil in filters:
            # Change current working directory according to which computer the
            # code is being ran from.
            if computer == 'W':
                os.chdir(os.path.join('/home/depot/STEPUP/raw', target, date,
                                      'ISR_Images', fil, 'WCS/accurate_WCS'))
            if computer == 'H':
                os.chdir(os.path.join('/Users/helenarichie/tests2/ISR_Images/',
                                      fil, '/WCS/accurate_WCS/'))

            # Perform absolute relative photometry on, create light curve from,
            # and generate output file of dataset.
            perform_photometry.perform_photometry(target, dirtarget, filters,
                                                  date, coords, comp_ra,
                                                  comp_dec, comp_mags, vsp_code,
                                                  cname, check_ra, check_dec,
                                                  rname, ref_ra, ref_dec,
                                                  verbose=False)

    else:
        # Determine if user would like to run a certain function from routine.
        cont_analysis = input('\nWould you still like to perform a function? (Y/N): ')
        while cont_analysis == 'Y':
            answer = input('\nWhich would you like to run? (ISR, ASTROM, PHOT): ')
            # Run function specified by "answer" variable.
            which_analysis(answer, dirtarget, dirdark, target, date, computer)
            # Determine if user has finished running STEPUP_image_analysis.
            cont_analysis = input('\nWould you still like to perform a function? (Y/N): ')


def which_analysis(answer, dirtarget, dirdark, target, date, computer):
    """Run one of three functions in image analysis routine.

    Based on user input, "answer" variable specifies which function from
    STEPUP_image_analysis that the user would like to run. Once the function
    finishes running, the user may run another function, but note that ISR and
    perform_astrometry create new directories and thus will cause an error if
    the user attempts to run them twice in a row without removing the created
    directories. However, perform_photometry may be ran as many times in a row
    as the user wishes.

    Parameters
    ----------
    answer : str
        Specifies whether user would like to perform instrument signature
        removal, astrometry, or photometry.
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    dirdark : str
        Directory containing all dark images.
    target : str
        Name of target.
    date : str
        Date of observation.
    computer : str
        Either "W" or "H" depending on whether code is being ran from Helena's
        computer or Astro Lab computer.

    Returns
    -------
    None
    """
    if answer == 'ISR':
        print('\nInstrument signature removal in progress...')
        # Removes instrument signatures from dataset.
        filters = ISR.ISR_main(dirtarget, dirdark, target)
        print('\nInstrument signature removal completed.')

    if answer == 'ASTROM':
        im = input('\nHave you saved a new-image.fits file to the appropriate directory? (Y/N): ')
        # Determines if user has saved new-image.fits WCS calibration file to
        # ISR_Images directory that was created in ISR function.
        if im == 'Y':
            # Determines filters using in observation.
            filters = input('\nEnter filters of observation (comma-delimited): ').split(",")
            print('\nAstrometry in progress...')
            dirtarget += '/ISR_Images'
            # Calculates WCS information for dataset.
            perform_astrometry.perform_astrometry(target, dirtarget,
                                                  filters, verbose=False)
            print('\nAstrometry completed.')
        else:
            return None

    if answer == 'PHOT':
        # Determines filters using in observation.
        filters = input('\nEnter filters of observation (comma-delimited): ').split(",")
        print('\nPhotometry in progress...')
        dirtarget += '/ISR_Images'
        # Initialize variable for and determine path based on which computer
        # code is being ran from.
        path = None
        if computer == 'W':
            path = os.path.join('/home/depot/STEPUP/raw', target, date)
        if computer == 'H':
            path = '/Users/helenarichie/tests2'
        # Change current working directory to specified path.
        os.chdir(path)

        # Initialize variables to obtain information from input-file.txt file.
        ra = []
        dec = []
        comp_ra = []
        comp_dec = []
        comp_mag = []
        vsp_code = None
        cname = None
        check_ra = []
        check_dec = []
        rname = None
        ref_ra = []
        ref_dec = []

        # Read in information from input-file.txt and append/assign it to
        # corresponding list/variable.
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
                    check_ra.append(line[5:].strip('\n'))
                if line.startswith('#CDEC='):
                    check_dec.append(line[6:].strip('\n'))
                if line.startswith('#RLABEL='):
                    rname = line[8:].strip('\n')
                if line.startswith('#RRA'):
                    ref_ra.append(line[5:].strip('\n'))
                if line.startswith('#RDEC='):
                    ref_dec.append(line[6:].strip('\n'))

        # Ensure that all magntidue values are floats.
        comp_mags = []
        for mag in comp_mag:
            comp_mags.append(float(mag))

        # Create list of two different lists of either right ascensions or
        # declinations.
        coords = []
        coords.append(ra)
        coords.append(dec)

        for fil in filters:
            # Initialize variable for and determine path based on which computer
            # code is being ran from.
            path = None
            if computer == 'W':
                path = os.path.join('/home/depot/STEPUP/raw', target, date,
                                    'ISR_Images', fil, 'WCS/accurate_WCS')
            if computer == 'H':
                path = os.path.join('/Users/helenarichie/tests2/ISR_Images',
                                    fil, 'WCS/accurate_WCS')

            # Change current working directory to specified path.
            os.chdir(path)

            # Perform absolute relative photometry on, create light curve from,
            # and generate output file of dataset.
            perform_photometry.perform_photometry(target, dirtarget,
                                                  filters, date, coords,
                                                  comp_ra, comp_dec,
                                                  comp_mags, vsp_code,
                                                  cname, check_ra, check_dec,
                                                  rname, ref_ra,
                                                  ref_dec, verbose=False)

main(verbose=False)
