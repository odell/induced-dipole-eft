'''
Functions for computing binding energies.
'''
from pyexpat.errors import XML_ERROR_ATTRIBUTE_EXTERNAL_ENTITY_REF
import numpy as np
from scipy import optimize
import scatter
from utility import *
import constants as const

def kernel_matrix_gen(en, v_matrix, q, wq):
    nq = np.size(q)
    kernel = np.zeros((nq, nq))
    for i in range(nq):
        qi = q[i]
        G0 = 1/(2*const.MU*en-qi**2)
        for j in range(nq):
            wj = wq[j]
            qj = q[j]
            mwq2 = 2*const.MU*wj*qj**2
            kernel[i, j] = float(i == j) - mwq2*v_matrix[i, j]*G0
    return kernel


def det_kernel(en, v_matrix, x_matrix_bare, gi, q, wq):
    v = v_matrix + gi*x_matrix_bare
    k = kernel_matrix_gen(en, v, q, wq)
    return np.linalg.det(k)


def diagonalize(v_bare_matrix, x_bare_matrix, g, p, wp):
    nq = np.size(p)
    v_matrix = v_bare_matrix+g*x_bare_matrix
    ham = np.zeros((nq, nq))
    for (i, pi) in enumerate(p):
        for (j, pj) in enumerate(p):
            ham[i, j] = (i == j)*pi*pj/(2*const.MU) + wp[j]*pj**2*v_matrix[i, j]
    
    return np.linalg.eig(ham)


def bound_states(v_bare_matrix, x_bare_matrix, g, p, wp):
    evals, evecs = diagonalize(v_bare_matrix, x_bare_matrix, g, p, wp)
    ii = np.sort(np.intersect1d(
        np.where(np.real(evals) < 0)[0],
        np.where(np.imag(evals) == 0)[0]
    ))
    return evals[ii], evecs[:, ii]


def spectrum(v_bare_matrix, x_bare_matrix, g, p, wp):
    evals, evecs = diagonalize(v_bare_matrix, x_bare_matrix, g, p, wp)
    bound_states = np.real(np.array(sorted(list(filter(lambda j: np.imag(j) == 0,
        filter(lambda i: i < 0, evals)))))[::-1])
    return bound_states
