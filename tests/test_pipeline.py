'''End-to-end pipeline tests on real data files (slow: real COBYLA fits).'''

import numpy as np
import pytest

from conftest import DATA_DIR
from modelingPS.pipeline import fit_directory, fit_observation
from modelingPS.results import load_results

pytestmark = pytest.mark.slow


def test_fit_observation_on_real_file():
    infile = str(DATA_DIR / '1200120105pow.dat')
    obs, result, search = fit_observation(infile, noise_mode='quadrature')

    assert result.obs_num == 105
    assert result.noise_mode == 'quadrature'
    assert result.n_lorentzians >= 1
    assert len(result.params) == 3 * result.n_lorentzians
    assert len(result.errors) == len(result.params)
    assert result.dof == obs.n_points() - len(result.params)
    assert result.red_chi2 == pytest.approx(result.chi2 / result.dof)
    # searchQPO=2 -> one trial frequency, QPO summary present
    assert search is not None and len(search['freqs']) == 1
    assert result.qpo is not None
    assert result.qpo['nu_char'] == pytest.approx(
        np.hypot(result.params[-2], result.params[-1]))


def test_fit_directory_writes_products(tmp_path):
    out = tmp_path / 'out'
    results = fit_directory(str(DATA_DIR), str(out),
                            pattern='120012010[56]pow.dat',
                            noise_mode='quadrature', verbose=False)

    assert len(results) == 2
    assert (out / 'fit_results.json').exists()
    assert (out / 'summary.pdf').exists()
    assert (out / 'plots' / '1200120105pow.dat.pdf').exists()
    assert (out / 'plots' / '1200120106pow.dat.pdf').exists()

    loaded, meta = load_results(str(out / 'fit_results.json'))
    assert len(loaded) == 2
    assert meta['noise_mode'] == 'quadrature'
    assert meta['search_qpo'] == 2


def test_fit_directory_no_match(tmp_path):
    with pytest.raises(FileNotFoundError):
        fit_directory(str(DATA_DIR), str(tmp_path / 'out'),
                      pattern='nope*.dat', verbose=False)
