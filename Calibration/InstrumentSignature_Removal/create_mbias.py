import numpy as np

def create_mbias():
    biases = get_biases(dircalib)
    mbias = np.median(biases, 2)
    return mbias
