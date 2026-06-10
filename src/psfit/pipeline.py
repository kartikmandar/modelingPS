'''Batch fitting pipeline: fit every power spectrum in a directory,
write per-observation plots, a JSON results file, an optional
legacy-format txt and a summary figure.'''

import glob
import os

import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt

from .fitting import fit
from .models import lore, nu_char, q_factor
from .observation import Observation
from .results import (FitResult, obs_num_from_path, save_results,
                      write_legacy_txt)
from .uncertainties import nu_char_error, parameter_errors, q_factor_error


def fit_observation(infile, noise_mode='none', search_qpo=2,
                    method='COBYLA', x0_given=None):
    '''
    Fit one file and package everything into a FitResult.

    Returns (obs, result, search) where search holds the QPO-search
    diagnostics needed by plot_fit (None when search_qpo == 0).
    '''
    obs = Observation(infile, noise_mode=noise_mode)
    out = fit(obs, x0_given=x0_given, searchQPO=search_qpo, method=method)

    params = out['params']
    errors, cov = parameter_errors(obs.fun, params)

    qpo = None
    if search_qpo > 0:
        r, D, nu0 = params[-3], params[-2], params[-1]
        block = cov[-2:, -2:]
        qpo = {
            'r': float(r), 'D': float(D), 'nu0': float(nu0),
            'r_err': float(errors[-3]),
            'D_err': float(errors[-2]),
            'nu0_err': float(errors[-1]),
            'nu_char': float(nu_char(D, nu0)),
            'nu_char_err': float(nu_char_error(D, nu0, block)),
            'Q': float(q_factor(D, nu0)),
            'Q_err': float(q_factor_error(D, nu0, block)),
        }

    result = FitResult(
        infile=infile,
        obs_num=obs_num_from_path(infile),
        params=[float(v) for v in params],
        errors=[float(e) for e in errors],
        n_lorentzians=out['n_lorentzians'],
        chi2=out['chi2'],
        dof=out['dof'],
        red_chi2=out['red_chi2'],
        noise_mode=noise_mode,
        search_qpo=search_qpo,
        qpo=qpo,
    )
    return obs, result, out['search']


def plot_fit(obs, result, search, outpath):
    '''Per-observation figure (ported from crawler_fix.py): data with
    each Lorentzian, weighted residuals and QPO-search diagnostics.'''
    params = np.asarray(result.params)
    numQ = result.n_lorentzians

    fs = [lore(*params[3 * i:3 * i + 3], obs.freq) for i in range(numQ)]

    if search is not None:
        fig, axi = plt.subplots(3, 1, sharex=True)
    else:
        fig, axi = plt.subplots(1, 1, sharex=True)
        axi = [axi]

    axi[0].errorbar(obs.freq, obs.freq * obs.power,
                    yerr=obs.freq * obs.err_power, marker='o')

    ft = np.sum(fs, axis=0)
    axi[0].plot(obs.freq, obs.freq * ft, '-')
    for f in fs:
        axi[0].plot(obs.freq, obs.freq * f, ':')

    if obs.noise_mode != 'none':
        axi[0].plot(obs.freq, obs.freq * obs.noise_level, '--',
                    color='gray', lw=1, label='effective noise level')
        axi[0].legend(loc='lower left', fontsize=8)

    axi[0].set_xscale('log')
    axi[0].set_yscale('log')
    axi[0].set_ylabel(r'$\nu f(\nu)$')
    axi[-1].set_xlabel(r'$\nu$')
    axi[0].set_ylim([1.e-3, axi[0].get_ylim()[1]])

    if search is not None:
        axi[1].plot(obs.freq, obs.fun(params, full=True) * np.sqrt(obs.weigh))
        axi[1].axhline(0, ls='-', color='C1', zorder=-5)
        axi[1].set_xscale('log')
        axi[1].set_ylabel('weighted residuals')

        axi[2].plot(search['freqs'], search['chi2'], 'r')
        axi[2].set_ylabel(r'$\chi^2$', color='r')
        ax2 = axi[2].twinx()
        ax2.plot(search['freqs'], search['Qs'], 'C0')
        ax2.set_ylabel('Q', color='C0')

        for ax in axi:
            ax.axvline(search['freqs'][search['loc']], ls='--', color='k')

    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def summary_figure(results, outpath):
    '''Characteristic frequency of every component vs ObsID suffix,
    coloured by |Q| clipped to [0, 2] (ported summary plot).'''
    fig, axi = plt.subplots(1, 1)
    axi.set_yscale('log')

    all_q = []
    for res in results:
        par = np.asarray(res.params)
        all_q.append(np.abs(q_factor(par[1::3], par[2::3])))
    mino = max(0., min(q.min() for q in all_q))
    maxo = min(2., max(q.max() for q in all_q))

    sc = None
    for res in results:
        par = np.asarray(res.params)
        cen = nu_char(par[1::3], par[2::3])
        Qs = np.abs(q_factor(par[1::3], par[2::3]))
        sc = axi.scatter([res.obs_num] * len(cen), cen, marker='o', s=20,
                         c=Qs, cmap='autumn', vmin=mino, vmax=maxo)

    cb = plt.colorbar(sc, ax=axi)
    axi.set_xlabel('ObsID suffix')
    axi.set_ylabel(r'$\nu_{\rm char}$')
    cb.set_label('Q')

    fig.savefig(outpath)
    plt.close(fig)


def fit_directory(data_dir, out_dir, pattern='*pow.dat', noise_mode='none',
                  search_qpo=2, method='COBYLA', plot=True,
                  legacy_txt=None, verbose=True):
    '''
    Fit every file matching pattern in data_dir.

    Writes out_dir/fit_results.json, out_dir/summary.pdf, per-obs plots
    under out_dir/plots/ and (optionally) out_dir/<legacy_txt>.

    Returns the list of FitResults.
    '''
    filelist = sorted(glob.glob(os.path.join(data_dir, pattern)))
    if not filelist:
        raise FileNotFoundError(
            f'no files matching {pattern!r} in {data_dir!r}')

    os.makedirs(out_dir, exist_ok=True)
    plot_dir = os.path.join(out_dir, 'plots')
    if plot:
        os.makedirs(plot_dir, exist_ok=True)

    results = []
    for infile in filelist:
        if verbose:
            print('\n' + infile)
        obs, result, search = fit_observation(
            infile, noise_mode=noise_mode, search_qpo=search_qpo,
            method=method)
        results.append(result)
        if verbose:
            print(f'  n_lor={result.n_lorentzians}'
                  f'  red_chi2={result.red_chi2:.2f}')
        if plot:
            outpath = os.path.join(plot_dir,
                                   os.path.basename(infile) + '.pdf')
            plot_fit(obs, result, search, outpath)

    meta = {'data_dir': data_dir, 'pattern': pattern,
            'noise_mode': noise_mode, 'search_qpo': search_qpo,
            'method': method}
    save_results(results, os.path.join(out_dir, 'fit_results.json'),
                 meta=meta)

    if legacy_txt:
        write_legacy_txt(results, os.path.join(out_dir, legacy_txt))

    summary_figure(results, os.path.join(out_dir, 'summary.pdf'))

    return results
