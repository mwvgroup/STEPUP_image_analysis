import os
import subprocess
import glob
from shutil import move
from shutil import copyfile
from astropy.io import fits


def perform_astrometry(target, dirtarget, filters, verbose=False):
    """Adds accurate WCS information to all headers in dataset.

    Uses imstar command from WCSTools to generate a table of positions and
    coordinates of previous image that has accurate WCS information. It then
    adds the WCS information from previous image to the header of the current
    image in the dataset to get a preliminary estimate of the true WCS
    information. It writes out this file and then uses imwcs to correct the WCS
    information in the image beign processed to reflect the positions of its
    sources. The result is a complete set of instrument signature removed,
    astrometrically calibrated images saved in the WCS folder.

    Parameters
    ----------
    target : str
        Name of target.
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    filters : list
        List containing string of each filter keyword found in header of flat
        field and light frame images.
    verbose : boolean, optional
        Print information about status of program.

    Returns
    -------
    None
    """
    dirtarget = os.path.join(dirtarget, 'ISR_Images')
    os.chdir(dirtarget)

    for fil in filters:
        try:
            os.mkdir(os.path.join(dirtarget, fil, 'WCS'))
        except FileExistsError:
            pass

    # Generate list of strings of three-digit numbers from 0 to 999
    # used to name files that are written.
    numbers1 = list(range(0, 999, 1))
    numbers2 = []
    for i in numbers1:
        numbers2.append(str(i))
    numbers = []
    for i in numbers2:
        if len(i) == 1:
            numbers.append('00{}'.format(i))
        if len(i) == 2:
            numbers.append('0{}'.format(i))
        if len(i) == 3:
            numbers.append(i)

    for fil in filters:
        os.chdir(dirtarget)

        images = sorted(glob.glob(os.path.join(dirtarget, fil, '*.fits')))

        copyfile('new-image.fits', os.path.join(fil, 'new-image.fits'))

        im_name = 'new-image.fits'
        for n, image in enumerate(images):
            print('Start: {}'.format(im_name))
            # Creates table with positions of stars and their coordinates.
            os.chdir(os.path.join(dirtarget, fil))
            args = '-vhi' if verbose else '-hi'
            subprocess.call(['imstar', args, '100', '-tw', im_name])

            solved_hdu = fits.open(os.path.join(dirtarget, fil, im_name))
            solved_header = solved_hdu[0].header

            hdu = fits.open(os.path.join(dirtarget, fil, image))
            header = hdu[0].header
            data = hdu[0].data

            diff = fits.HeaderDiff(solved_header, header).diff_keywords[0]

            for keyword in diff:
                # Skips unneeded header keywords and adds uncommon keywords
                # and their value to images.
                if keyword not in ('COMMENT', 'HISTORY'):
                    header.set(keyword, solved_header[keyword])

            tab_name = im_name.rstrip('.fits') + '.tab'

            im_name = target + '-' + fil + '-{}c.fits'.format(numbers[n+1])

            out_path = os.path.join(dirtarget, fil, im_name)
            hdu.writeto(out_path, overwrite=True)

            print('End: {}'.format(im_name))
            print(tab_name)

            subprocess.call(['imwcs', '-w', '-i', '100', '-c', tab_name,
                             im_name])

            im_name = target + '-' + fil + '-{}cw.fits'.format(numbers[n+1])

        plate_solved = glob.glob(os.path.join(dirtarget, fil, '*cw.fits'))
        out_path = os.path.join(dirtarget, fil, 'WCS')
        for image in plate_solved:
            move(image, out_path)

        scratch = glob.glob(os.path.join(dirtarget, fil, '*c.fits'))
        for image in scratch:
            os.remove(image)

        scratch = glob.glob(os.path.join(dirtarget, fil, '*cw.tab'))
        for image in scratch:
            os.remove(image)
