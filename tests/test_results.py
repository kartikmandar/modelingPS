'''Tests for psfit.results.'''

import math

from psfit.results import (FitResult, load_results, obs_num_from_path,
                           save_results, write_legacy_txt)


def make_result(obs_num=101, params=(0.42, 1.4, -0.00025)):
    params = list(params)
    return FitResult(
        infile=f'./maxij1820_powspecs/1200120{obs_num}pow.dat',
        obs_num=obs_num,
        params=params,
        errors=[0.01, 0.1, float('nan')],
        n_lorentzians=len(params) // 3,
        chi2=250.0,
        dof=224,
        red_chi2=250.0 / 224,
        noise_mode='quadrature',
        search_qpo=2,
        qpo={'r': 0.42, 'D': 1.4, 'nu0': -0.00025,
             'r_err': 0.01, 'D_err': 0.1, 'nu0_err': float('nan'),
             'nu_char': 1.4, 'nu_char_err': 0.1,
             'Q': -8.9e-5, 'Q_err': 0.01},
    )


def test_obs_num_from_path():
    assert obs_num_from_path('./maxij1820_powspecs/1200120101pow.dat') == 101
    assert obs_num_from_path('1200120197pow.dat') == 197


def test_round_trip(tmp_path):
    results = [make_result(101), make_result(105)]
    path = tmp_path / 'fit_results.json'
    save_results(results, str(path), meta={'noise_mode': 'quadrature'})

    loaded, meta = load_results(str(path))
    assert meta == {'noise_mode': 'quadrature'}
    assert len(loaded) == 2
    # NaN != NaN breaks plain equality; compare field by field
    assert loaded[0].infile == results[0].infile
    assert loaded[0].obs_num == 101
    assert loaded[0].params == results[0].params
    assert loaded[0].errors[:2] == results[0].errors[:2]
    assert math.isnan(loaded[0].errors[2])
    assert loaded[0].n_lorentzians == results[0].n_lorentzians
    assert loaded[0].chi2 == results[0].chi2
    assert loaded[0].dof == results[0].dof
    assert loaded[0].red_chi2 == results[0].red_chi2
    assert loaded[0].noise_mode == 'quadrature'
    assert loaded[0].search_qpo == 2
    assert loaded[0].qpo['nu_char'] == results[0].qpo['nu_char']
    assert math.isnan(loaded[0].qpo['nu0_err'])
    assert loaded[1].obs_num == 105


def test_write_legacy_txt(tmp_path):
    results = [make_result(101)]
    path = tmp_path / 'resu.txt'
    write_legacy_txt(results, str(path))
    line = path.read_text().rstrip('\n')
    # exact legacy format: path then '  {0: 8.2}' per parameter
    expected = results[0].infile
    for v in results[0].params:
        expected += '  {0: 8.2}'.format(v)
    assert line == expected
