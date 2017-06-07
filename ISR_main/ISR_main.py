import sys
sys.path.insert(0, '')
import get_calibimages
import create_mcalib
import instrument_signature_removal

def ISR_main(dirtarget, dirdark, target):
    """Creates ISR FITS files by executing a preliminary calibration sequence.

    Imports all ISR modules (get_calibimage, instrument_signature_removal, and
    calibration) which get all calibration images, create master calibration
    images, and subtracts them from the raw science images. It then adds the
    expected saturation to the header and saves the ISR science images to the
    specified directory.

    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    dirdark : str
        Directory containing all dark images.

    Returns
    -------
    science_images : numpy.ndarray
        3D array containing ISR science images.
    """

    (biases, darks, flats, dark_exptime, bias_prihdr, flat_prihdr, dark_prihdr, exptime) = get_calibimages.get_calibimages(dirtarget, dirdark)

    mbias = create_mcalib.create_mbias(biases, bias_prihdr, dirtarget)
    mdark = create_mcalib.create_mdark(darks, mbias, dark_prihdr, dirtarget, dark_exptime, exptime)
    mflat = create_mcalib.create_mflat(flats, mbias, flat_prihdr, dirtarget)

    science_images = instrument_signature_removal.instrument_signature_removal(dirtarget, mbias, mdark, mflat, dark_exptime, target)
    
    return science_images    

    
