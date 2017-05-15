import numpy as np

def create_mdark():
    mbias = create_mbias()
    darks = get_darks(dirdark)
    mdark = np.meadian(darks, 2) - mbias
    return mdark
    
