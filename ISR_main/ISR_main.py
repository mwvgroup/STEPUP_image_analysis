import sys
sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/ISR_main/ISR')
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
    (biases, darks, r_flats, b_flats, v_flats, dark_exptime,
    exptime, bias_prihdr, dark_prihdr, r_flat_prihdr, b_flat_prihdr,
    v_flat_prihdr, r_filter, b_filter, v_filter) = get_calibimages.get_calibimages(dirtarget, dirdark)

    mbias = create_mcalib.create_mbias(biases, bias_prihdr, dirtarget)
    mdark = create_mcalib.create_mdark(darks, mbias, dark_prihdr, dirtarget, dark_exptime, exptime)
    (r_mflat, b_mflat, v_mflat) = create_mcalib.create_mflat(r_flats, b_flats, v_flats, mbias, r_flat_prihdr,
                                                             b_flat_prihdr, v_flat_prihdr, r_filter, b_filter,
                                                             v_filter, dirtarget)

    (r_science_images,
     b_scimages,
     v_scimages) = instrument_signature_removal.instrument_signature_removal(dirtarget, target, mbias, mdark,
                                                                               r_mflat, b_mflat, v_mflat, dark_exptime,
                                                                               exptime)
    
    return r_science_images, b_scimages, v_scimages

    
