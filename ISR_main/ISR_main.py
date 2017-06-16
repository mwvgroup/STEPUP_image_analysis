import sys
sys.path.insert(0, '/Users/helenarichie/GitHub/STEPUP_image_analysis/ISR_main/ISR')
import get_unfiltered_calibimages
import get_filtered_calibimages
import instrument_signature_removal

def ISR_main(dirtarget, dirdark, target):
    """Creates ISR FITS files by executing a preliminary calibration sequence.

    Imports and calls get_unfiltered_calibimages, get_filtered_calibimages, and
    instrument_signature_removal. Saves ISR science images to
    dirtarget/ISR_Images/<filter-name>.
    
    Parameters
    ----------
    dirtarget : str
        Directory containing all bias, flat, and raw science images.
    dirdark : str
        Directory containing all dark images.
    traget : str
        Name of target.

    Returns
    -------
    n/a
    """


    exptime = get_unfiltered_calibimages.get_unfiltered_calibimages(dirtarget, dirdark)
    image_filters = get_filtered_calibimages.get_filtered_calibimages(dirtarget)

    instrument_signature_removal.instrument_signature_removal(dirtarget, target, exptime, image_filters)

    
