'''Tests for psfit.uncertainties.'''

import numpy as np
import pytest

from psfit.uncertainties import (covariance_from_chi2, derived_error,
                                 hessian, nu_char_error, parameter_errors,
                                 q_factor_error)


def test_hessian_of_quadratic():
    # f(x, y) = 2 x^2 + 3 y^2 + x y  ->  H = [[4, 1], [1, 6]] everywhere
    def func(p):
        return 2. * p[0]**2 + 3. * p[1]**2 + p[0] * p[1]

    H = hessian(func, np.array([0.7, -0.3]))
    assert np.allclose(H, [[4., 1.], [1., 6.]], atol=1.e-4)


def test_linear_model_error_is_exact():
    # chi2(a) = sum((y - a x)^2 / err^2)  ->  sigma_a = 1/sqrt(sum(x^2/err^2))
    rng = np.random.default_rng(0)
    x = np.linspace(1., 10., 50)
    err = np.full_like(x, 0.5)
    y = 2.0 * x + rng.normal(0., 0.5, x.size)

    def chi2(p):
        return np.sum((y - p[0] * x)**2 / err**2)

    a_best = np.sum(x * y / err**2) / np.sum(x**2 / err**2)
    errors, cov = parameter_errors(chi2, np.array([a_best]))
    assert errors[0] == pytest.approx(1. / np.sqrt(np.sum(x**2 / err**2)),
                                      rel=1.e-3)


def test_parameter_errors_flag_flat_directions():
    # a parameter the chi2 does not depend on must give NaN, not a number
    def chi2(p):
        return (p[0] - 1.)**2

    errors, cov = parameter_errors(chi2, np.array([1., 5.]))
    assert np.isfinite(errors[0])
    assert np.isnan(errors[1])


def test_covariance_nan_on_nonfinite_hessian():
    def chi2(p):
        return np.inf

    cov = covariance_from_chi2(chi2, np.array([1.0]))
    assert np.isnan(cov).all()


def test_derived_error_simple_sum():
    # z = x + y with independent unit variances -> sigma_z = sqrt(2)
    assert derived_error([1., 1.], np.identity(2)) == pytest.approx(np.sqrt(2.))


def test_nu_char_and_q_errors():
    cov = np.identity(2)  # uncorrelated unit errors on (D, nu0)
    # nu_char = 5 for (D, nu0) = (3, 4); grad = (3/5, 4/5) -> sigma = 1
    assert nu_char_error(3., 4., cov) == pytest.approx(1.0)
    # Q = nu0/(2D); grad = (-nu0/(2 D^2), 1/(2D)) = (-2/9, 1/6)
    assert q_factor_error(3., 4., cov) == pytest.approx(
        np.sqrt((2. / 9.)**2 + (1. / 6.)**2))
