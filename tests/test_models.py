'''Tests for modelingPS.models.'''

import numpy as np

from modelingPS.models import lore, nu_char, q_factor


def test_lore_peak_at_centroid():
    nu = np.linspace(0.5, 1.5, 1001)
    vals = lore(0.1, 0.02, 1.0, nu)
    assert np.isclose(nu[np.argmax(vals)], 1.0, atol=2.e-3)


def test_lore_peak_value():
    # at nu = nu0 the Lorentzian is r^2 / (pi * D)
    val = lore(0.1, 0.02, 1.0, np.array([1.0]))[0]
    assert np.isclose(val, 0.1**2 / (np.pi * 0.02))


def test_lore_integral_is_r_squared():
    # the full-line integral equals r^2; integrate wide around the peak
    nu = np.linspace(-200., 200., 2_000_001)
    vals = lore(0.3, 0.05, 1.0, nu)
    assert np.isclose(np.trapezoid(vals, nu), 0.3**2, rtol=1.e-3)


def test_lore_no_division_by_zero():
    # D = 0 evaluated exactly at nu = nu0 must stay finite
    val = lore(0.1, 0.0, 1.0, np.array([1.0]))
    assert np.isfinite(val).all()


def test_nu_char():
    assert np.isclose(nu_char(3.0, 4.0), 5.0)


def test_q_factor():
    assert np.isclose(q_factor(0.05, 1.0), 10.0)
