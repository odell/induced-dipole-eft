import os
import sys
import numpy as np
import constants as const
import utility
import scatter
from bind import spectrum
from counterterm import Counterterm1, LocalCounterterm
import mu2

from scipy import optimize

def regulator(r, R):
    return (1 - np.exp(-(r/R)**2))**4


def potential(r, R):
    return regulator(r, R) * -const.C4/r**4


class LocalSystem:
    '''
    A system is characterized by a cutoff in coordinate space.
    '''
    def __init__(self,
                 r_c,
                 rmesh=None,
                 qmesh=None,
                 nq=200,
                 kmin=0.01/const.BETA4,
                 kmax=0.3/const.BETA4,
                 ell=0
    ):
        '''
        Upon instantiation, several things are generated:
            * momentum mesh (q, wq)
            * momentum-space matrix elements
                * v_vdW (v_tilde)
                * LO counterterm (xlo_tilde)
                * NLO counterterm (xnlo_tilde)
        '''
        self.r_c = r_c
        if rmesh == None:
            r, wr = utility.log_mesh(0, 10*const.BETA4, 2000)
        else:
            r, wr = rmesh
        self.r_nodes = np.copy(r)
        self.r_weights = np.copy(wr)
        if qmesh is None:
            self.q, self.wq = utility.log_mesh(0, 10*2/self.r_c, nq)
        else:
            self.q, self.wq = qmesh
        self.ks = np.linspace(kmin, kmax, 30)
        self.ell = ell
        self.counterterm = LocalCounterterm(r, wr, self.q, self.r_c, ell)
        self.v_tilde = utility.ft_matrix_gen(lambda r: potential(r, self.r_c),
                ell, ell, self.q, r, wr)
    
    
    def phase_shifts(self, glo, gnlo):
        '''
        Calculates phase shifts (in radians).
        '''
        xterm = self.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        return np.array(
            [scatter.phase_shift(ki, v, self.q, self.wq, 10*2/self.r_c,
                2*const.MU, degrees=False) for ki in self.ks]
        )

        
    def cross_sections(self, glo, gnlo, l):
        '''
        Calculates the lth partial-wave term in the cross section.
        '''
        deltas = self.phase_shifts(glo, gnlo)
        return np.array(
            [4*np.pi/k**2 * (2*l+1) * np.sin(delta)**2 for (k, delta) in zip(self.ks, deltas)]
        )


    def effective_range_expansion(self, glo, gnlo, l, use_c=False):
        '''
        Calculates k^(2l+1) cot(delta).
        '''
        xterm = self.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        if use_c:
            return self.k3cotd_gen_fast(self.ks, glo, gnlo)
        else:
            return np.array(
                [ki**(2*l)*scatter.kcotdelta(ki, v, self.q, self.wq, 10*2/self.r_c, 2*const.MU) for ki in self.ks]
            )
    
    
    def k3cotd(self, p, glo, gnlo):    
        '''
        Calculates k^3 cot(delta) as at a specificied momentum, p.
        '''
        xterm = self.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        pp = p*p
        print(p)
        scatter_temp = scatter.kcotdelta(p, v, self.q, self.wq, 10*2/self.r_c, 2*const.MU)
        return pp*scatter_temp
    
    
    def k3cotd_gen(self, ks, glo, gnlo):    
        '''
        Calculates k^3 cot(delta) as function of ks (array) so that we can curve_fit the phaseshifts.
        '''
        xterm = self.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        return np.array(
            [ki**2*scatter.kcotdelta(ki, v, self.q, self.wq, 10*2/self.r_c, 2*const.MU) for ki in ks]
        )
    
    
    def k3cotd_gen_fast(self, ks, glo, gnlo):    
        '''
        Calculates k^3 cot(delta) as function of ks (array) so that we can curve_fit the phaseshifts.
        '''
        xterm = self.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        return np.array(
            [cscatter.k3cotdelta_py(ki, v, self.q, self.wq, 10*2/self.r_c, 2*const.MU) for ki in ks]
        )
    
    
    def effective_range_parameters(self, glo, gnlo, p0=[-1, -1, 1, 1], use_c=False):
        '''
        Returns the $P$-wave effective range parameters for given LO (glo) and NLO (gnlo) strengths.
        Fits up to and including c_3 (cubic term).
        '''
        k3cds = self.effective_range_expansion(glo, gnlo, use_c=use_c)
        pars, cov = optimize.curve_fit(lambda x, c0, c1, c2, c3: c0 + c1*x + c2*x**2 + c3*x**3,
                                       self.ks, k3cds,
                                       p0=p0, maxfev=20000)
        sig = np.sqrt(np.diag(cov))
        return pars, sig
    
    
    def a1_and_r1(self, glo, gnlo, p0=None, use_c=False):
        '''
        Returns a_1 and r_1 after fitting the effective range parameters.
        '''
        pars, sig = self.effective_range_parameters(glo, gnlo, p0=p0, use_c=use_c)
        return -1/pars[0], 2*pars[2]
    
    
    def diff_a1(self, glo):
        a1, _ = self.a1_and_r1(glo, 0)
        return a1_tune - a1
    
    
    def bound_state_spectrum(self, glo, gnlo):
        xterm = self.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        return spectrum(v, 0, 0, self.q, self.wq)


class NonlocalSystem:
    '''
    A system is characterized by a cutoff in coordinate space.
    '''
    def __init__(self,
                 r_c,
                 rmesh=None,
                 nq=200,
                 kmin=0.01/const.BETA4,
                 kmax=0.3/const.BETA4,
                 ell=0
    ):
        '''
        Upon instantiation, several things are generated:
            * momentum mesh (q, wq)
            * momentum-space matrix elements
                * v_vdW (v_tilde)
                * LO counterterm (xlo_tilde)
                * NLO counterterm (xnlo_tilde)
        '''
        self.r_c = r_c
        self.l = ell
        if rmesh == None:
            r, wr = utility.log_mesh(0, 10*const.BETA4, 2000)
        else:
            r, wr = rmesh
        self.r_nodes = np.copy(r)
        self.r_weights = np.copy(wr)
        self.q, self.wq = utility.log_mesh(0, 10*2/self.r_c, nq)
        self.ks = np.linspace(kmin, kmax, 30)
        self.v_tilde = utility.ft_matrix_gen(lambda r: potential(r, self.r_c/10),
                self.l, self.l, self.q, r, wr)
        xreg = np.array([[np.exp(-(p*self.r_c/2)**4)*np.exp(-(k*self.r_c/2)**4)
            for p in self.q] for k in self.q])
        self.v_tilde *= xreg
        self.counterterm = Counterterm1(self.q, self.r_c, self.l)
    
    
    def phase_shifts(self, glo, gnlo):
        '''
        Calculates phase shifts (in radians).
        '''
        xterm = self.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        return np.array(
            [scatter.phase_shift(ki, v, self.q, self.wq, 10*2/self.r_c,
                2*const.MU, degrees=False) for ki in self.ks]
        )

        
    def cross_sections(self, glo, gnlo):
        '''
        Calculates the P-wave term in the cross section.
        '''
        deltas = self.phase_shifts(glo, gnlo)
        return np.array(
            [4*np.pi/k**2 * 3 * np.sin(delta)**2 for (k, delta) in zip(self.ks, deltas)]
        )


    def effective_range_expansion(self, glo, gnlo, use_c=False):
        '''
        Calculates k^3 cot(delta).
        '''
        xterm = self.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        if use_c:
            return self.k3cotd_gen_fast(self.ks, glo, gnlo)
        else:
            return np.array(
                [ki**(2*self.l)*scatter.kcotdelta(ki, v, self.q, self.wq, 10*2/self.r_c, 2*const.MU) for ki in self.ks]
            )
    
    
    def k3cotd(self, p, glo, gnlo):    
        '''
        Calculates k^3 cot(delta) as at a specificied momentum, p.
        '''
        xterm = self.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        scatter_temp = scatter.kcotdelta(p, v, self.q, self.wq, 10*2/self.r_c, 2*const.MU)
        return p**(2*self.l)*scatter_temp
    
    
    def k3cotd_gen(self, ks, glo, gnlo):    
        '''
        Calculates k^3 cot(delta) as function of ks (array) so that we can curve_fit the phaseshifts.
        '''
        xterm = self.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        return np.array(
            [ki**(2*self.l)*scatter.kcotdelta(ki, v, self.q, self.wq, 10*2/self.r_c, 2*const.MU) for ki in ks]
        )
    
    
    def k3cotd_gen_fast(self, ks, glo, gnlo):    
        '''
        Calculates k^3 cot(delta) as function of ks (array) so that we can curve_fit the phaseshifts.
        '''
        xterm = self.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        return np.array(
            [cscatter.k3cotdelta_py(ki, v, self.q, self.wq, 10*2/self.r_c, 2*const.MU) for ki in ks]
        )
    
    
    def effective_range_parameters(self, glo, gnlo, p0=[-1, -1, 1, 1], use_c=False):
        '''
        Returns the $P$-wave effective range parameters for given LO (glo) and NLO (gnlo) strengths.
        Fits up to and including c_3 (cubic term).
        '''
        k3cds = self.effective_range_expansion(glo, gnlo, use_c=use_c)
        pars, cov = optimize.curve_fit(lambda x, c0, c1, c2, c3: c0 + c1*x + c2*x**2 + c3*x**3,
                                       self.ks, k3cds,
                                       p0=p0, maxfev=20000)
        sig = np.sqrt(np.diag(cov))
        return pars, sig
    
    
    def a1_and_r1(self, glo, gnlo, p0=None, use_c=False):
        '''
        Returns a_1 and r_1 after fitting the effective range parameters.
        '''
        pars, sig = self.effective_range_parameters(glo, gnlo, p0=p0, use_c=use_c)
        return -1/pars[0], 2*pars[2]
    
    
    def diff_a1(self, glo):
        a1, _ = self.a1_and_r1(glo, 0)
        return a1_tune - a1
    
    
    def bound_state_spectrum(self, glo, gnlo):
        xterm = self.counterterm.gen(glo, gnlo)
        v = self.v_tilde + xterm
        return spectrum(v, 0, 0, self.q, self.wq)
