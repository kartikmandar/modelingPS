'''Tests for modelingPS.fitting.'''

import numpy as np
import pytest

from modelingPS.fitting import (fit, fit_components, lorentzian_constraints,
                           make_args, qpo_try_max, qpo_try_search)
from modelingPS.observation import Observation


@pytest.fixture
def synth_obs(synth_file):
    path, params = synth_file
    return Observation(str(path)), params


def test_make_args_routes_constraints_and_bounds():
    const, bounds = lorentzian_constraints(2)
    assert make_args('COBYLA', const, bounds) == {'constraints': const}
    assert make_args('SLSQP', const, bounds) == {'constraints': const}
    assert list(make_args('L-BFGS-B', const, bounds)) == ['bounds']
    assert list(make_args('Nelder-Mead', const, bounds)) == ['bounds']


def test_lorentzian_constraints_bounds():
    const, bounds = lorentzian_constraints(2)
    assert bounds.shape == (6, 2)
    assert np.allclose(bounds[:, 0], 0.)
    # r upper bounds are 1; D and nu0 upper bounds are 100
    assert np.allclose(bounds[[0, 3], 1], 1.)
    assert np.allclose(bounds[[1, 2, 4, 5], 1], 100.)


def test_fit_components_structure(synth_obs):
    # NOTE: the incremental fitter models the broadband continuum; narrow
    # features are intentionally left to the QPO step, so no absolute
    # goodness-of-fit is asserted here (see the legacy-equivalence test).
    obs, _ = synth_obs
    params, num, chi2 = fit_components(obs)
    assert num >= 1
    assert len(params) == 3 * num
    # the returned chi2 belongs to the returned parameters
    assert chi2 == pytest.approx(obs.fun(params))


def test_fit_components_matches_legacy_crawler(synth_file):
    # behaviour-preservation regression: the ported incremental fitter
    # must reproduce legacy/crawler_fix.py exactly on identical input
    import crawler_fix

    path, _ = synth_file
    obs = Observation(str(path))
    params, num, chi2 = fit_components(obs)

    legacy = crawler_fix.obs(str(path))
    legacy_params = legacy.fit(None, plot=False, searchQPO=0)

    assert num == len(legacy_params) // 3
    assert np.array_equal(params, legacy_params)
    assert chi2 == pytest.approx(legacy.fun(legacy_params))


def test_qpo_try_max_q_uses_qpo_width(synth_obs):
    obs, _ = synth_obs
    base = np.array([0.3, 0.1, 0.0])
    sola, freco, chi2, Qs, numQ, loc = qpo_try_max(obs, base)
    assert numQ == len(sola) // 3
    assert loc == 0
    assert len(freco) == len(chi2) == len(Qs) == 1
    # Q must be built from the QPO parameters at the end of sola
    assert Qs[0] == pytest.approx(sola[-1] / sola[-2] / 2.)


def test_qpo_try_search_q_uses_qpo_width(synth_obs):
    # regression test for the legacy bug: Q used x[1] (width of the FIRST
    # background Lorentzian) instead of x[-1] (the QPO width)
    obs, _ = synth_obs
    base = np.array([0.3, 0.1, 0.0])
    sola, freco, chi2, Qs, numQ, loc = qpo_try_search(obs, base,
                                                      freqs=[0.9, 1.0])
    assert numQ == 2
    assert len(freco) == len(chi2) == len(Qs) == 2
    # solution layout: [background..., r_qpo, D_qpo, fr_qpo]
    assert Qs[loc] == pytest.approx(freco[loc] / sola[-2] / 2.)


def test_fit_full_run(synth_obs):
    obs, _ = synth_obs
    out = fit(obs, searchQPO=2)
    params = out['params']
    assert out['n_lorentzians'] == len(params) // 3
    assert out['dof'] == obs.n_points() - len(params)
    assert out['red_chi2'] == pytest.approx(out['chi2'] / out['dof'])
    assert out['search'] is not None
    # some component should sit near the synthetic QPO at 1 Hz
    centroids = params[2::3]
    assert np.any(np.abs(centroids - 1.0) < 0.2)


def test_fit_no_qpo_search(synth_obs):
    obs, _ = synth_obs
    out = fit(obs, searchQPO=0)
    assert out['search'] is None
    assert len(out['params']) == 3 * out['n_lorentzians']


def test_fit_invalid_searchqpo(synth_obs):
    obs, _ = synth_obs
    with pytest.raises(ValueError):
        fit(obs, searchQPO=3)
