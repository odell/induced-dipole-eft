'''
Defines common constants and functions.
'''

import numpy as np
from scipy.special import roots_legendre, spherical_jn

def kd(i, j):
    '''
    Kronecker delta
    '''
    return 1.0 if i == j else 0.0


def log_mesh(a, b, n):
    '''
    Returns a logarithmicaly spaced mesh: nodes (x) and weights (w), from a to
    b (noninclusive).
    '''
    x, w = roots_legendre(n)
    x *= (np.log(b+1)-np.log(a+1))/2
    x += (np.log(b+1)+np.log(a+1))/2
    w *= (np.log(b+1)-np.log(a+1))/2
    w = np.exp(x)*w
    x = np.exp(x)-1
    return x, w


def ft_matrix_gen(pot, l, lp, p, r, wr):
    '''
    Given a coordinate-space potential (pot), returns the partial-wave
    projected, momentum-space matrix elements.
    '''
    nq = np.size(p)
    nr = np.size(r)
    A = np.zeros((nq, nr))
    B = np.zeros((nr, nq))
    for (i, pi) in enumerate(p):
        for ((j, rj), wj) in zip(enumerate(r), wr):
            A[i, j] = 2/np.pi * wj * rj**2 * spherical_jn(l, pi*rj)
            B[j, i] = pot(rj) * spherical_jn(lp, pi*rj)
    return np.matmul(A, B)


# @jit
# def ft_matrix_gen_jit(pot, l, lp, p, r, wr, F=1):
#     '''
#     Given a coordinate-space potential (pot), returns the partial-wave
#     projected, momentum-space matrix elements.
#     '''
#     nq = p.size
#     nr = r.size
#     A = np.zeros((nq, nr))
#     B = np.zeros((nr, nq))
#     for (i, pi) in enumerate(p):
#         for ((j, rj), wj) in zip(enumerate(r), wr):
#             A[i, j] = F * wj * rj**2 * spherical_jn(l, pi*rj)
#             B[j, i] = pot(rj) * spherical_jn(lp, pi*rj)
#     return np.matmul(A, B)
# 

def v_matrices_gen(pot, l, lp, p, r, wr):
    '''
    Generates coupled-channel potential matrices.
    '''
    v00_matrix = ft_matrix_gen(lambda ri: pot(ri, l, l), l, l, p, r, wr)
    v02_matrix = ft_matrix_gen(lambda ri: pot(ri, l, lp), l, lp, p, r, wr)
    v20_matrix = ft_matrix_gen(lambda ri: pot(ri, lp, l), lp, l, p, r, wr)
    v22_matrix = ft_matrix_gen(lambda ri: pot(ri, lp, lp), lp, lp, p, r, wr)
    return (v00_matrix, v02_matrix, v20_matrix, v22_matrix)


def progress_bar(count, tot):
    '''
    Prints a progress bar to output.
    count / tot * 100 = %
    '''
    numblocks = round(count/tot*50)
    numspaces = 50-numblocks
    s = '|' + u'\u2588'*numblocks + ' '*numspaces + '|' + ' {0:.0f}%'.format(count/tot*100)
    print(s)
    # clear_output(wait=True)

