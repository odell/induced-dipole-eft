import numpy as np
from constants import BETA4
from utility import ft_matrix_gen


def local_LO_counterterm(r, R):
    return np.exp(-(r/R)**4)


def local_NLO_counterterm(r, R):
    derivative_factor = (16*r**6 - 12*r**2*R**4) / R**8
    return derivative_factor * local_LO_counterterm(r, R)


class LocalCounterterm:
    def __init__(self, r, wr, q, R, ell):
        self.q = np.copy(q)
        self.Lambda = 2/R
        self.x_tilde_LO = ft_matrix_gen(lambda r: local_LO_counterterm(r, R),
                ell, ell, self.q, r, wr)
        self.x_tilde_NLO = ft_matrix_gen(lambda r: local_NLO_counterterm(r, R),
                ell, ell, self.q, r, wr)
    

    def gen(self, gLO, gNLO):
        return gLO*self.x_tilde_LO + gNLO*self.x_tilde_NLO


def h(q, L, n, l):
    return (q/L)**l*np.exp(-(q/L)**n)


def nlo_factor(p, k, L):
    return ((p/L)**2 + (k/L)**2)/2


class Counterterm1:
    n = 4
    def __init__(self, q_nodes, R, ell):
        self.q = np.copy(q_nodes)
        self.Lambda = 2/R
        self.l = ell

        self.lo = np.array(
            [[h(p, self.Lambda, self.n, self.l)*h(k, self.Lambda, self.n,
                self.l) for p in self.q] for k in self.q]
        )
        self.nlo = np.array(
            [[self.lo[i, j] * nlo_factor(p, k, self.Lambda) for (i, p) in
                enumerate(self.q)] for (j, k) in enumerate(self.q)]
        )

    
    def gen(self, gLO, gNLO):
        return gLO * self.lo + gNLO * self.nlo