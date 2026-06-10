'''Fitting machinery: incremental Lorentzian fits with an F-test stopping
criterion, plus the two QPO search strategies.

Ported from src/crawler_fix.py keeping the algorithm identical (same
optimizer, constraints, initial guesses and F-test). The only
behavioural change is the corrected Q computation in qpo_try_search:
the QPO width is the last fitted parameter (x[-1]), not the width of
the first background Lorentzian (x[1]).
'''

import numpy as np
from scipy.optimize import LinearConstraint as LC
from scipy.optimize import minimize as mini
from scipy.stats import f as ftest

DEFAULT_GUESS = (0.3, 0.5, 0.0)
QPO_GUESS = (0.085, 1.e-2)
MAXITER = 20000


def make_args(method, const, bounds):
    '''Route constraints or bounds to the optimizer depending on method.'''
    kws = {}
    if method in ('COBYLA', 'SLSQP', 'trust-constr'):
        kws['constraints'] = const
    if method in ('Nelder-Mead', 'Powell', 'L-BFGS-B', 'TNC'):
        kws['bounds'] = bounds
    return kws


def lorentzian_constraints(num):
    '''
    Constraints and bounds for num Lorentzians.

    0 <= r <= 1, 0 <= D <= 100, 0 <= nu0 <= 100 for each component.

    Returns (LinearConstraint, bounds array of shape (3 num, 2)).
    '''
    A = np.identity(3 * num)
    lb = np.zeros(3 * num)
    ub = np.full(3 * num, 100.)
    for i in range(num):
        ub[3 * i] = 1.
    return LC(A, lb=lb, ub=ub), np.array([lb, ub]).T


def fit_components(obs, x0_given=None, ftestlim=0.05, method='COBYLA'):
    '''
    Adds Lorentzians one at a time until an F-test flags no improvement.

    Takes:

    obs     : Observation to fit
    x0_given: optional initial guess, flat [r, D, nu0] * n
    ftestlim: probability used as the stopping condition
    method  : scipy.optimize.minimize method

    Returns:

    (params, num, chi2) of the last accepted step.
    '''
    leno = 0 if x0_given is None else len(x0_given) // 3
    npts = obs.n_points()

    num = 0
    cow = 0
    x0_old = []
    fit_use = None
    conv = False
    chis = []
    while not conv:

        num += 1
        cow += 1

        const, bounds = lorentzian_constraints(num)

        # previous solution plus the next given triplet, or the default
        if x0_given is not None and num <= leno:
            ini = 3 * (num - 1)
            x0 = np.append(x0_old, x0_given[ini:ini + 3])
        else:
            x0 = np.append(x0_old, DEFAULT_GUESS)

        min_fit = mini(obs.fun, x0=x0, **make_args(method, const, bounds),
                       method=method, options={'maxiter': MAXITER})

        chis.append(min_fit.fun / (npts - 3 * num))

        if cow == 1:
            conv = False
        else:
            prob = ftest.cdf(
                (chis[-2] - chis[-1]) * (npts - 3 * num) / chis[-1] / 3.,
                npts - 3 * num + 3,
                npts - 3 * num)
            conv = prob < ftestlim

        if not conv:
            # convergence means: this step was useless
            fit_use = min_fit

        x0_old = min_fit.x

    num -= 1
    return fit_use.x, num, fit_use.fun


def qpo_try_search(obs, solution, method='COBYLA', freqs=None):
    '''
    Fits a QPO Lorentzian over a grid of frozen trial frequencies.

    The previous solution seeds the background, which is refit together
    with the QPO (r, D) at each trial frequency; the solution with the
    minimum chi^2 wins.

    Takes:

    obs     : Observation to fit
    solution: parameters of the solution so far
    method  : scipy.optimize.minimize method
    freqs   : trial frequencies (default: 100 points in 1e-3 .. 10 Hz)

    Returns:

    sola : parameters including the QPO, [background..., r, D, nu0]
    freco: the trial frequencies
    chi2 : chi^2 at each trial frequency
    Qs   : Q of the fitted QPO at each trial frequency
    numQ : number of Lorentzians including the QPO
    loc  : index of the winning trial
    '''
    numQ = len(solution) // 3 + 1

    # background: 0 <= r <= 1, 0 <= D, nu0 <= 100
    # QPO: 0.01 <= r <= 0.1, 0 <= D <= 1
    numo = 3 * numQ - 1
    A = np.identity(numo)
    lb = np.zeros(numo)
    ub = np.full(numo, 100.)
    lb[-2] = 0.01
    for i in range(numQ - 1):
        ub[3 * i] = 1.
    ub[-2] = 0.1
    ub[-1] = 1.0
    const = LC(A, lb=lb, ub=ub)
    bounds = np.array([lb, ub]).T

    if freqs is None:
        freco = np.linspace(1.e-3, 10., 100, endpoint=True)
    else:
        freco = np.asarray(freqs, dtype=float)

    sols = []
    chi2 = []
    Qs = []
    for fr in freco:

        x0 = np.append(solution, QPO_GUESS)
        sols.append(mini(obs.funQPO2, x0=x0, args=(fr,),
                         **make_args(method, const, bounds),
                         method=method, options={'maxiter': MAXITER}))

        # Q of the trial QPO; its width is the last fitted parameter
        Qs.append(fr / sols[-1].x[-1] / 2.)
        chi2.append(obs.funQPO2(sols[-1].x, fr))

    loc = int(np.argmin(chi2))
    sola = np.append(sols[loc].x, [freco[loc]])

    return sola, freco, chi2, Qs, numQ, loc


def qpo_try_max(obs, solution, method='COBYLA'):
    '''
    Fits a QPO at the maximum of nu f(nu), then refits everything freely.

    First step: the given solution and the QPO frequency (the max of
    nu f(nu)) are frozen; only the QPO (r, D) are free, with
    0 <= r <= 0.1 and 0 <= D <= 100.

    Second step: everything is refit without constraints, seeded by the
    given solution, the first-step result and the implied centroid
    nu0 = sqrt(max(fr^2 - D^2, 1e-5)).

    Returns the same tuple shape as qpo_try_search (1-element
    diagnostics, loc = 0).
    '''
    imax = np.argmax(obs.freq * obs.power)
    fr = obs.freq[imax]

    A = np.identity(2)
    lb = np.zeros(2)
    ub = np.array([0.1, 100.])
    const = LC(A, lb=lb, ub=ub)
    bounds = np.array([lb, ub]).T

    parz = mini(obs.funQPO, x0=np.array(QPO_GUESS), args=(solution, fr),
                **make_args(method, const, bounds),
                method=method, options={'maxiter': MAXITER})

    nu0 = np.sqrt(max(fr**2 - parz.x[1]**2, 1.e-5))

    x0 = np.concatenate((solution, parz.x, [nu0]))

    refit = mini(obs.fun, x0=x0, method=method, options={'maxiter': MAXITER})
    sola = refit.x

    freco = [fr]
    chi2 = [obs.fun(sola)]
    Qs = [sola[-1] / sola[-2] / 2.]
    numQ = len(sola) // 3

    return sola, freco, chi2, Qs, numQ, 0


def fit(obs, x0_given=None, searchQPO=2, ftestlim=0.05, method='COBYLA'):
    '''
    Full fit of an observation: component search plus optional QPO step.

    Takes:

    obs      : Observation to fit
    x0_given : optional initial guess
    searchQPO: 0 none, 1 qpo_try_search, 2 qpo_try_max
    ftestlim : F-test stopping probability
    method   : scipy.optimize.minimize method

    Returns a dict:

    params       : best-fit parameters (QPO last, if searched)
    n_lorentzians: number of Lorentzians including the QPO
    chi2, dof, red_chi2
    search       : None, or dict(freqs, chi2, Qs, loc) from the QPO step
    '''
    if searchQPO not in (0, 1, 2):
        raise ValueError(f'searchQPO must be 0, 1 or 2, got {searchQPO!r}')

    base, num, base_chi2 = fit_components(obs, x0_given=x0_given,
                                          ftestlim=ftestlim, method=method)

    if searchQPO == 0:
        sola = np.asarray(base)
        numQ = num
        chi2 = float(base_chi2)
        search = None
    else:
        searcher = qpo_try_search if searchQPO == 1 else qpo_try_max
        sola, freco, chi2s, Qs, numQ, loc = searcher(obs, base, method=method)
        sola = np.asarray(sola)
        chi2 = float(chi2s[loc])
        search = {'freqs': np.asarray(freco), 'chi2': np.asarray(chi2s),
                  'Qs': np.asarray(Qs), 'loc': loc}

    dof = obs.n_points() - 3 * numQ
    return {
        'params': sola,
        'n_lorentzians': numQ,
        'chi2': chi2,
        'dof': dof,
        'red_chi2': chi2 / dof,
        'search': search,
    }
