import os
import sys
from ctypes import cdll, c_int, c_double, POINTER

import numpy as np

# LIB = cdll.LoadLibrary('/Users/danielodell/6Li/cc/libk3cd.so')
LIB = np.ctypeslib.load_library('/Users/danielodell/6Li/cc/libk3cd.so', '')

k3cotdelta = LIB.k3cotdelta
k3cotdelta.argtypes = (
    c_double,
    POINTER(c_double),
    POINTER(c_double),
    POINTER(c_double),
    c_int,
    c_double,
    c_double
)
k3cotdelta.restype = c_double

def k3cotdelta_py(k, V, q, wq, qmax, mass):
    Vp = V.copy()
    qp = q.copy()
    wqp = wq.copy()
    return k3cotdelta(k,
            Vp.ctypes.data_as(POINTER(c_double)), 
            qp.ctypes.data_as(POINTER(c_double)),
            wqp.ctypes.data_as(POINTER(c_double)),
            qp.size, qmax, mass)