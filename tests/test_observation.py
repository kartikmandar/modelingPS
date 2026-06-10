'''Tests for psfit.observation.'''

import numpy as np
import pytest

from conftest import DATA_DIR, write_dat
from psfit.models import lore
from psfit.observation import Observation

REAL_FILE = str(DATA_DIR / '1200120101pow.dat')


def test_load_real_file():
    # the real files have 227 lines: 1 header ('READ SERR 1 2') + 226 rows
    obs = Observation(REAL_FILE)
    assert len(obs.freq) == 226
    assert len(obs.power) == 226
    assert len(obs.err_power) == 226
    assert len(obs.noise_level) == 226
    assert obs.freq[0] == pytest.approx(0.0025)
    assert np.allclose(obs.weigh, 1. / obs.err_power**2)
    assert obs.n_points() == 226


def test_noise_mode_quadrature():
    obs = Observation(REAL_FILE, noise_mode='quadrature')
    assert np.allclose(obs.weigh, 1. / (obs.err_power**2 + obs.noise_level**2))
    plain = Observation(REAL_FILE)
    # inflated errors can only lower the weights
    assert np.all(obs.weigh <= plain.weigh)


def test_noise_mode_mask():
    obs = Observation(REAL_FILE, noise_mode='mask')
    below = obs.power < obs.noise_level
    assert below.any()  # the high-frequency tail dips below the floor
    assert np.all(obs.weigh[below] == 0.)
    assert obs.n_points() == int((~below).sum())


def test_noise_mode_invalid():
    with pytest.raises(ValueError):
        Observation(REAL_FILE, noise_mode='bogus')


def test_fun_zero_on_noiseless_model(tmp_path, freq_grid):
    params = np.array([0.3, 0.1, 0.0])
    model = lore(*params, freq_grid)
    path = tmp_path / 'noiseless.dat'
    write_dat(path, freq_grid, model,
              np.full_like(freq_grid, 0.01), np.zeros_like(freq_grid))
    obs = Observation(str(path))
    assert obs.fun(params) == pytest.approx(0.0, abs=1.e-12)
    assert np.allclose(obs.fun(params, full=True), 0., atol=1.e-12)


def test_funqpo2_equals_fun_with_frozen_frequency(synth_file):
    path, _ = synth_file
    obs = Observation(str(path))
    par_qpo = np.array([0.3, 0.1, 0.0, 0.12, 0.05])
    fr = 1.0
    assert obs.funQPO2(par_qpo, fr) == pytest.approx(obs.fun(np.append(par_qpo, fr)))


def test_funqpo_equals_fun_with_frozen_solution(synth_file):
    path, _ = synth_file
    obs = Observation(str(path))
    base = np.array([0.3, 0.1, 0.0])
    qpo = np.array([0.12, 0.05])
    fr = 1.0
    full = np.concatenate([base, qpo, [fr]])
    assert obs.funQPO(qpo, base, fr) == pytest.approx(obs.fun(full))
