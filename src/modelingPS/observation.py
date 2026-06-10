'''Observation container for power spectra.

Loads the 5-column .dat files used by this project
(frequency, meaningless zeros, power, power error, effective noise
level) and evaluates sums of Lorentzians against them. The effective
noise level is the uncertainty floor left by the Poisson-noise
subtraction: signal below it is unreliable.
'''

import numpy as np

from .models import lore

NOISE_MODES = ('none', 'quadrature', 'mask')


class Observation:
    '''
    Stores an observation PSD and evaluates Lorentzian models on it.

    Attributes:

    infile     : path of the data file
    freq       : frequencies of the power spectrum
    power      : powers (Poisson noise already subtracted)
    err_power  : errors on the powers
    noise_level: effective noise level per bin
    noise_mode : how the noise level enters the chi^2 weights:
                 'none'      : weigh = 1 / err^2          (legacy behaviour)
                 'quadrature': weigh = 1 / (err^2 + noise_level^2)
                 'mask'      : points with power < noise_level get weight 0
    weigh      : the chi^2 weights
    '''

    def __init__(self, infile, noise_mode='none'):
        if noise_mode not in NOISE_MODES:
            raise ValueError(
                f'noise_mode must be one of {NOISE_MODES}, got {noise_mode!r}')

        self.infile = infile
        self.noise_mode = noise_mode

        dat = np.genfromtxt(infile, skip_header=1).T

        self.freq = dat[0]
        self.power = dat[2]
        self.err_power = dat[3]

        if dat.shape[0] > 4:
            self.noise_level = dat[4]
        else:
            self.noise_level = np.zeros_like(self.freq)

        if noise_mode == 'quadrature':
            self.weigh = 1. / (self.err_power**2 + self.noise_level**2)
        elif noise_mode == 'mask':
            self.weigh = 1. / self.err_power**2
            self.weigh[self.power < self.noise_level] = 0.
        else:
            self.weigh = 1. / self.err_power**2

    def n_points(self):
        '''Number of data points contributing to the chi^2.'''
        return int(np.count_nonzero(self.weigh))

    def model(self, par):
        '''Sum of Lorentzians; their number is len(par) // 3.'''
        n = len(par) // 3
        tot = np.zeros(len(self.power))
        for i in range(n):
            r, D, nu0 = par[3 * i:3 * (i + 1)]
            tot += lore(r, D, nu0, self.freq)
        return tot

    def fun(self, par, full=False):
        '''Chi^2 (or residuals if full=True) of n Lorentzians.'''
        chi = self.power - self.model(par)
        if full:
            return chi
        return np.sum(chi**2 * self.weigh)

    def funQPO(self, par, res, fr, full=False):
        '''
        Chi^2 with the solution res frozen plus one extra Lorentzian.

        par: (r, D) of the extra (QPO) Lorentzian, whose frequency is
        frozen at fr.
        '''
        chi = self.power - self.model(res) - lore(par[0], par[1], fr, self.freq)
        if full:
            return chi
        return np.sum(chi**2 * self.weigh)

    def funQPO2(self, par, fr, full=False):
        '''
        Chi^2 of n free Lorentzians plus a QPO at frozen frequency fr.

        par: [r, D, nu0] * n for the background, then (r, D) of the QPO.
        '''
        chi = (self.power - self.model(par[:-2])
               - lore(par[-2], par[-1], fr, self.freq))
        if full:
            return chi
        return np.sum(chi**2 * self.weigh)
