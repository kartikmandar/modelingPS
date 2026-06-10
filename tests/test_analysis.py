'''Tests for modelingPS.analysis.'''

import numpy as np
import pytest

from modelingPS.analysis import (fit_quality_figure, identify_qpos,
                            qpo_evolution_figure)
from modelingPS.results import FitResult


def qpo_dict(Q=5.0, r=0.1, D=0.05, nu0=0.5):
    return {'r': r, 'D': D, 'nu0': nu0,
            'r_err': 0.01, 'D_err': 0.01, 'nu0_err': 0.01,
            'nu_char': float(np.hypot(D, nu0)), 'nu_char_err': 0.01,
            'Q': Q, 'Q_err': 0.5}


def make_result(obs_num, red_chi2=1.0, n_lor=3, qpo=None):
    params = [0.3, 0.5, 0.1] * n_lor
    return FitResult(
        infile=f'1200120{obs_num}pow.dat', obs_num=obs_num,
        params=params, errors=[0.01] * len(params),
        n_lorentzians=n_lor, chi2=red_chi2 * 100., dof=100,
        red_chi2=red_chi2, noise_mode='none', search_qpo=2, qpo=qpo)


def test_identify_qpos_filters():
    results = [
        make_result(101, qpo=qpo_dict(Q=5.0)),             # keep
        make_result(102, qpo=qpo_dict(Q=0.5)),             # too broad
        make_result(103, qpo=qpo_dict(Q=float('nan'))),    # not finite
        make_result(104, qpo=qpo_dict(Q=8.0, nu0=-0.1)),   # negative centroid
        make_result(105, qpo=qpo_dict(Q=8.0, r=-0.01)),    # negative amplitude
        make_result(106, qpo=None),                        # no QPO searched
    ]
    kept = identify_qpos(results)
    assert [r.obs_num for r in kept] == [101]


def test_identify_qpos_resolvability_and_significance():
    # sub-resolution spike: D far below the frequency bin width
    spike = qpo_dict(Q=30.0, D=3.e-4, nu0=0.02)
    # resolved but amplitude only 2 sigma (r=0.02, r_err=0.01)
    weak = qpo_dict(Q=5.0, r=0.02)
    # good detection
    good = qpo_dict(Q=5.0)
    # resolved, no usable error -> the width cut alone must not kill it
    no_err = qpo_dict(Q=5.0, D=0.03, nu0=0.3)
    no_err['r_err'] = float('nan')

    results = [
        make_result(101, qpo=spike),
        make_result(102, qpo=weak),
        make_result(103, qpo=good),
        make_result(104, qpo=no_err),
    ]
    kept = identify_qpos(results, min_width=1.e-3, min_significance=3.0)
    assert [r.obs_num for r in kept] == [103, 104]
    # with the cuts disabled (defaults) everything Q > 2 survives
    assert len(identify_qpos(results)) == 4


def test_qpo_evolution_figure(tmp_path):
    results = [make_result(100 + i, qpo=qpo_dict(Q=3.0 + i))
               for i in range(5)]
    out = tmp_path / 'qpo_evolution.pdf'
    kept = qpo_evolution_figure(results, str(out))
    assert out.exists() and out.stat().st_size > 0
    assert len(kept) == 5


def test_qpo_evolution_figure_no_detections(tmp_path):
    results = [make_result(101, qpo=qpo_dict(Q=0.5))]
    with pytest.raises(ValueError):
        qpo_evolution_figure(results, str(tmp_path / 'qpo.pdf'))


def test_fit_quality_figure(tmp_path):
    results = [make_result(100 + i, red_chi2=0.8 + 0.3 * i,
                           n_lor=3 + (i % 3)) for i in range(8)]
    out = tmp_path / 'fit_quality_analysis.pdf'
    fit_quality_figure(results, str(out))
    assert out.exists() and out.stat().st_size > 0
