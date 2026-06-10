'''Parameter uncertainties from the chi^2 curvature.

For chi^2 = sum((data - model)^2 / err^2) the parameter covariance at
the minimum is cov = 2 H^-1, with H the Hessian of the chi^2. Derived
quantities use the delta method: var(g) = grad(g) . cov . grad(g).

COBYLA stops slightly off the exact minimum, so the Hessian can come
out non-positive along some directions; those parameters get NaN
errors rather than misleading numbers.
'''

import numpy as np

DEFAULT_STEP = 1.e-4


def hessian(func, x, step=DEFAULT_STEP):
    '''
    Central-difference Hessian of func at x.

    The step for parameter i is step * |x_i| with an absolute floor of
    step for parameters near zero.
    '''
    x = np.asarray(x, dtype=float)
    n = x.size
    h = np.maximum(np.abs(x) * step, step)

    H = np.empty((n, n))
    f0 = func(x)

    for i in range(n):
        xp = x.copy()
        xm = x.copy()
        xp[i] += h[i]
        xm[i] -= h[i]
        H[i, i] = (func(xp) - 2. * f0 + func(xm)) / h[i]**2

        for j in range(i + 1, n):
            xpp = x.copy()
            xpm = x.copy()
            xmp = x.copy()
            xmm = x.copy()
            xpp[i] += h[i]
            xpp[j] += h[j]
            xpm[i] += h[i]
            xpm[j] -= h[j]
            xmp[i] -= h[i]
            xmp[j] += h[j]
            xmm[i] -= h[i]
            xmm[j] -= h[j]
            H[i, j] = ((func(xpp) - func(xpm) - func(xmp) + func(xmm))
                       / (4. * h[i] * h[j]))
            H[j, i] = H[i, j]

    return H


def covariance_from_chi2(chi2_func, params, step=DEFAULT_STEP):
    '''Covariance matrix at a chi^2 minimum: 2 * pinv(H).'''
    H = hessian(chi2_func, params, step=step)
    if not np.all(np.isfinite(H)):
        n = len(params)
        return np.full((n, n), np.nan)
    return 2. * np.linalg.pinv(H)


def parameter_errors(chi2_func, params, step=DEFAULT_STEP):
    '''
    1-sigma parameter errors from the chi^2 curvature.

    Returns (errors, cov); errors are NaN where the curvature is flat,
    negative or non-finite.
    '''
    cov = covariance_from_chi2(chi2_func, params, step=step)
    diag = np.array(np.diag(cov), dtype=float)
    diag[~(diag > 0.)] = np.nan
    return np.sqrt(diag), cov


def derived_error(grad, cov):
    '''Delta-method 1-sigma error of a scalar derived quantity.'''
    grad = np.asarray(grad, dtype=float)
    var = float(grad @ cov @ grad)
    return np.sqrt(var) if var > 0. else np.nan


def nu_char_error(D, nu0, cov):
    '''Error on nu_char = sqrt(D^2 + nu0^2); cov is the (D, nu0) block.'''
    nc = np.hypot(D, nu0)
    if nc == 0.:
        return np.nan
    return derived_error([D / nc, nu0 / nc], cov)


def q_factor_error(D, nu0, cov):
    '''Error on Q = nu0 / (2 D); cov is the (D, nu0) block.'''
    if D == 0.:
        return np.nan
    return derived_error([-nu0 / (2. * D * D), 1. / (2. * D)], cov)
