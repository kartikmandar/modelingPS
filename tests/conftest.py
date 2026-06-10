'''Shared test fixtures: synthetic power spectra with known Lorentzian content.'''

import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC = REPO_ROOT / 'src'
sys.path.insert(0, str(SRC))

from psfit.models import lore  # noqa: E402

DATA_DIR = REPO_ROOT / 'maxij1820_powspecs'


def make_spectrum(params, freq, noise_frac=0.05, seed=42):
    '''
    Synthetic PSD: sum of Lorentzians plus Gaussian noise.

    params: flat [r, D, nu0] * n
    freq  : frequency grid

    Returns (power, err).
    '''
    rng = np.random.default_rng(seed)

    model = np.zeros_like(freq)
    for i in range(len(params) // 3):
        r, D, nu0 = params[3 * i:3 * i + 3]
        model += lore(r, D, nu0, freq)

    err = np.maximum(noise_frac * model, noise_frac * model.max() * 1.e-2)
    power = model + rng.normal(0., err)

    return power, err


def write_dat(path, freq, power, err, noise_level):
    '''Write a 5-column .dat file in the real data format.'''
    rows = np.column_stack([freq, np.zeros_like(freq), power, err, noise_level])
    np.savetxt(path, rows, header='READ SERR 1 2', comments='')


@pytest.fixture
def freq_grid():
    return np.geomspace(1.e-3, 100., 200)


@pytest.fixture
def synth_file(tmp_path, freq_grid):
    '''A synthetic spectrum file: broad noise Lorentzian + QPO at 1 Hz.'''
    params = np.array([0.3, 0.1, 0.0, 0.12, 0.05, 1.0])
    power, err = make_spectrum(params, freq_grid)
    path = tmp_path / '9999999101pow.dat'
    write_dat(path, freq_grid, power, err, np.full_like(freq_grid, 1.e-4))
    return path, params
