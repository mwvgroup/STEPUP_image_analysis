import sys
sys.path.insert(0, 'directory where code is stored')
import get_calibimages
import create_mcalib
import sys
sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/ISR_main/ISR')
import get_calibimages
import create_mcalib
import instrument_signature_removal

def ISR_main(dirtarget, dirdark, target):
    """Creates ISR FITS files by executing a preliminary calibration sequence.

    Imports all ISR modules (get_calibimage, instrument_signature_removal, and calibration)
    which get all calibration images, create master calibration images, and
    subtract them from the raw science images. It then adds the expected
    saturation to the header and saves the ISR science images to the specified
    directory.

    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    dirdark : str
        Directory containing all dark images.

    Returns
    -------
    None
    """

    (biases, darks, flats, dark_exptime) = get_calibimages.get_calibimages(dirtarget, dirdark)

    mbias = create_mcalib.create_mbias(biases)
    mdark = create_mcalib.create_mdark(darks, mbias)
    mflat = create_mcalib.create_mflat(flats, mbias, mdark)

    science_images = instrument_signature_removal.instrument_signature_removal(dirtarget, mbias, mflat, mdark, dark_exptime, target)
    
    return science_images    

    
