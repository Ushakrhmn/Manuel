# diag.py

import numpy as np
from scipy.linalg import eigh

def diagonalize(H_flat):
    H = np.array(H_flat, dtype=np.complex128).reshape(6, 6)
    evals, evecs = eigh(H)
    return evals.tolist(), evecs.flatten(order='F').tolist()  # column-major for C
