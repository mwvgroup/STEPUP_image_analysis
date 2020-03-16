import os
import sys
sys.path.insert(0, 'ISR')
import ISR
sys.path.insert(0, 'Calibration')
import perform_astrometry
import perform_photometry
from shutil import copyfile


def main():
    """Execute specified functions of STEPUP Image Analysis.

    Based on whether the use would like to run SIA interactively, the user is
    either prompted for which functions they would like to perform. If the user
    is not running SIA interactively, the functions to be ran are executed in
    the order that they are specified in the input-file.txt file, without
    pausing. If the user is running SIA interactively, they will enter one of
    the three function names, which then runs until completion. Once the
    funcion has completed, the user will be prompted to run another function if
    they wish, and if not SIA will finish running.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    # Initialize variables for assignment based on the contents of the
    # input-file.txt file.
    set_rad = False
    functions = None

    # Determine directory containing input-file.txt, target data, and some or
    # all calibration data as well as if user would like specify functions at
    # command line.
    dirtarget = input('\nInput target directory: ')
    while not os.path.exists((os.path.join(dirtarget, 'input-file.txt'))):
        dirtarget = input('\nThis directory does not contain an input file. Check to ensure that the file exists and is saved in your data directory.\n\nInput target directory or enter "Q" to quit: ')
        if dirtarget.lower() == 'q':
            return None
    interactive = input('\nWould you like to run SIA interatively? (Y/N): ').lower().strip(' ')

    # Specify keywords SIA should search for in input-file.txt.
    str_keywords = ['#TARGET=', '#DATE=', '#DIRDARK=', '#CLABEL=', '#APERRAD=',
                    '#ANNINRAD=', '#ANNOUTRAD=']
    list_keywords = ['#FILTERS=', '#RA=', '#DEC=', '#COMPMAGS=', '#COMPRA=',
                     '#COMPDEC=', '#CRA=', '#CDEC=']

    # Read in information from input-file.txt and assign it to corresponding
    # variable.
    str_input_values = []
    list_input_values = []
    with open(os.path.join(dirtarget, 'input-file.txt')) as f:
        for line in f:
            # Look for string variables.
            for keyword in str_keywords:
                if line.startswith(keyword):
                    str_input_values.append(line[len(keyword):].strip('\n'))

            for keyword in list_keywords:
                # Look for list variables.
                if line.startswith(keyword):
                    list_input_values.append(line[len(keyword):].strip('\n').split(','))

            if line.startswith('#APERRAD='):
                if len(line.strip(' ').strip('\n')) == 9:
                    pass
                else:
                    set_rad = True
            if line.startswith('#FUNCTIONS='):
                if interactive == 'n':
                    functions = line[11:].strip('\n').split(',')

    target = str_input_values[0]
    date = str_input_values[1]
    dirdark = str_input_values[2]
    clabel = str_input_values[3]
    aper_rad = str_input_values[4]
    ann_in_rad = str_input_values[5]
    ann_out_rad = str_input_values[6]
    filters = list_input_values[0]
    ra = list_input_values[1]
    dec = list_input_values[2]
    comp_mag = list_input_values[3]
    comp_ra = list_input_values[4]
    comp_dec = list_input_values[5]
    cra = list_input_values[6]
    cdec = list_input_values[7]

    # If dark files are stored in same directory as input-file.txt, target
    # data, and other calibration files, dirdark is assigned to same string as
    # dirtarget.
    if dirdark == '':
        dirdark = dirtarget

    # Ensure that all magntidue values are floats.
    comp_mags = []
    for mag in comp_mag:
        comp_mags.append(float(mag))

    # Create list of two different lists of either right ascensions or
    # declinations.
    coords = []
    coords.append(ra)
    coords.append(dec)

    if interactive == 'y':
        # Allows user to specify functions to be ran at command line.
        cont_analysis = interactive
        while cont_analysis == 'y':
            answer = input('\nWhich function would you like to run? (ISR, ASTROM, PHOT): ').lower().strip(' ')
            # Run function specified by "answer" variable.
            which_analysis(interactive, answer, target, date, filters, coords,
                           dirtarget, dirdark, comp_mags, comp_ra, comp_dec,
                           clabel, cra, cdec, set_rad,
                           aper_rad, ann_in_rad, ann_out_rad)
            # Determine if user has finished running STEPUP Image Analysis.
            cont_analysis = input('\nWould you still like to perform a function? (Y/N): ').lower().strip(' ')
        print('\nGoodbye.')

    else:
        # Runs which_analysis for function in order of specification in the
        # input-file.txt file if user is not running SIA interactively.
        for function in functions:
            answer = function.lower()
            which_analysis(interactive, answer, target, date, filters, coords,
                           dirtarget, dirdark, comp_mags, comp_ra, comp_dec,
                           clabel, cra, cdec, set_rad, aper_rad, ann_in_rad,
                           ann_out_rad)


def which_analysis(interactive, answer, target, date, filters, coords,
                   dirtarget, dirdark, comp_mags, comp_ra, comp_dec, clabel,
                   cra, cdec, set_rad, aper_rad, ann_in_rad, ann_out_rad):
    """Run one of three functions in image analysis routine.

    Based on user input, "answer" variable specifies which function from
    STEPUP Image Analysis that the user would like to run. Once the function
    finishes running, the user may run another function, but note that ISR and
    perform_astrometry create new directories and thus will cause an error if
    the user attempts to run them twice in a row without removing the created
    directories. However, perform_photometry may be ran as many times in a row
    as the user wishes.

    Parameters
    ----------
    interactive : str
        User input specifying whether they would like to run SIA interactively.
    answer : str
        Specifies whether user would like to perform instrument signature
        removal, astrometry, or photometry.
    target : str
        Name of target.
    date : str
        Date of observation.
    filters : list
        List containing string of each filter keyword found in header of flat
        field and light frame images.
    coords : list
        List of list of string of target right ascension and declination.
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    dirdark : str
        Directory containing all dark images.
    comp_mags : list
        List of floats representing the magnitudes of the comparison stars.
    comp_ra : list
        List of strings of comparison stars' right ascension.
    comp_dec : list
        List of strings of comparison stars' declination.
    clabel : str
        Name of check star.
    cra : list
        List of string of right ascension of check star.
    cdec : list
        List of string of delination of check star.
    set_rad : Boolean
        Determine whether user would like to use default aperture/annulus radii
        or specify their own.
    aper_rad : float
        User-specified aperture radius in arcseconds.
    ann_in_rad : float
        User-specified annulus inner radius in arcseconds.
    ann_out_rad : float
        User-specified annulus outer radius in arcseconds.

    Returns
    -------
    None
    """
    if answer == 'isr':
        print('\nInstrument signature removal in progress...')
        # Removes instrument signatures from dataset.
        filters = ISR.ISR_main(dirtarget, dirdark, target)
        print('\nInstrument signature removal completed.')

    if answer == 'astrom':
        im = None
        try:
            copyfile(os.path.join(dirtarget, 'new-image.fits'),
                         os.path.join(dirtarget, 'ISR_Images/new-image.fits'))
            im = 'y'
        except FileNotFoundError:
            print('\nnew-image.fits not found in raw data directory. Try again.')

        # Determines if user has saved new-image.fits WCS calibration file to
        # ISR_Images directory that was created in ISR function.
        if im == 'y':
            print('\nAstrometry in progress...')
            # Calculates WCS information for dataset.
            perform_astrometry.perform_astrometry(target, dirtarget, filters,
                                                  verbose=False)
            print('\nAstrometry completed.')
        else:
            return None

    if answer == 'phot':
        print('\nPhotometry in progress...')
        # Perform absolute differential photometry on, create light curve(s)
        # from, and generate output file of dataset.
        perform_photometry.perform_photometry(target, dirtarget, filters, date,
                                              coords, comp_ra, comp_dec,
                                              comp_mags, clabel, cra, cdec, set_rad,
                                              aper_rad, ann_in_rad, ann_out_rad)
        print('\nPhotometry completed.')


main()
