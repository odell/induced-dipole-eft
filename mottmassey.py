'''
Defines the Mott-Massey (MM) potential.
Provides an easy and consistent way of accessing the MM matrix elements.
'''
import numpy as np
from interaction import LocalSystem
import bind
import utility

def alpha(r):
    '''
    Numerator of the MM potential.
    '''
    return 9/2 - 2/3*np.exp(-2*r) * (r**5 + 9/2*r**4 + 9*r**3 + 27/2*r**2 + 27/2*r + 27/4)


def mm_potential(r):
    '''
    Full MM potential.
    '''
    return -1/2 * alpha(r) / r**4