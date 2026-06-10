'''Lorentzian model functions for power-spectrum fitting.'''

import numpy as np


def lore(r, D, nu0, nu):
    '''
    Calculates a Lorentzian.

    There is a + 1.e-20 in the denominator to avoid division by zero.
    The full-line integral of this form equals r**2, so r is the rms
    amplitude of the component.

    Takes:

    r  : amplitude
    D  : half width
    nu0: centroid
    nu : np array where to evaluate

    Returns:

    The values
    '''
    return (r**2 * D / np.pi) / (D**2 + (nu - nu0)**2 + 1.e-20)


def nu_char(D, nu0):
    '''Characteristic frequency sqrt(nu0^2 + D^2).'''
    return np.sqrt(np.asarray(D)**2 + np.asarray(nu0)**2)


def q_factor(D, nu0):
    '''Quality factor Q = nu0 / (2 D).'''
    return np.asarray(nu0) / np.asarray(D) / 2.
