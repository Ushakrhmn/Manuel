import ctypes
import numpy as np

lib = ctypes.CDLL("./libchi2.so")



lib.compute_chi2.argtypes = [ctypes.c_double] * 7
lib.compute_chi2.restype = ctypes.c_double

lib.compute_chi2_ih.argtypes = [ctypes.c_double] * 7
lib.compute_chi2_ih.restype = ctypes.c_double

lib.free_globes.argtypes = []

# Test initialization and chi2 computation


def test_compute_chi2():
    params_nh = [np.arcsin(np.sqrt(0.02237)),   # theta13
                 np.arcsin(np.sqrt(0.52)),      # theta23
                 0.0,                           # delta_cp
                 2.528e-3,                      # dm31
                 5.0,                           # mu
                 10.0,                          # N
                 0.2]                           # m0
    
    params_nh_ctypes = [ctypes.c_double(param) for param in params_nh]
    
    try:
        chi2_nh = lib.compute_chi2(*params_nh_ctypes)
        print("Chi2 (NH):", chi2_nh)
    except Exception as e:
        print("Error in compute_chi2 for NH:", e)

# Run tests

test_compute_chi2()
