'''Post-fit analysis: QPO identification/evolution and fit-quality
figures. Rebuilds qpo_evolution.pdf and fit_quality_analysis.pdf from a
fit_results.json produced by the pipeline (the original scripts for
those figures were lost).'''

import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt

QPO_MIN_Q = 2.0
RED_CHI2_GOOD = 1.5
RED_CHI2_MARGINAL = 2.5


def identify_qpos(results, min_q=QPO_MIN_Q, min_width=0., min_significance=0.):
    '''
    Keep fits whose QPO component looks like a real QPO.

    Criteria: a QPO was searched, its Q is finite and above min_q, and
    its amplitude, width and centroid are all positive. Optionally:

    min_width       : reject widths below this (in Hz). Widths far below
                      the frequency resolution are unresolved single-bin
                      noise spikes, not physical QPOs.
    min_significance: reject amplitudes below this many sigma, judged
                      only when the amplitude error is usable (finite
                      and positive) - detections whose Hessian failed
                      are not penalized for the missing error.
    '''
    out = []
    for res in results:
        qpo = res.qpo
        if qpo is None:
            continue
        if not np.isfinite(qpo['Q']):
            continue
        if qpo['r'] <= 0. or qpo['D'] <= 0. or qpo['nu0'] <= 0.:
            continue
        if qpo['Q'] <= min_q:
            continue
        if qpo['D'] < min_width:
            continue
        if min_significance > 0.:
            r_err = qpo['r_err']
            if (np.isfinite(r_err) and r_err > 0.
                    and qpo['r'] / r_err < min_significance):
                continue
        out.append(res)
    return out


def _yerr(values):
    '''Error bars for matplotlib: NaNs (unavailable errors) become 0.'''
    arr = np.asarray(values, dtype=float)
    return np.where(np.isfinite(arr), arr, 0.)


def qpo_evolution_figure(results, outpath, min_q=QPO_MIN_Q,
                         min_width=0., min_significance=0.):
    '''
    Three panels vs ObsID suffix: QPO characteristic frequency (log),
    quality factor and amplitude, all with 1-sigma error bars.

    Returns the list of QPO detections used.
    '''
    qpos = identify_qpos(results, min_q=min_q, min_width=min_width,
                         min_significance=min_significance)
    if not qpos:
        raise ValueError('no QPO detections to plot')

    nums = np.array([r.obs_num for r in qpos])
    nu = np.array([r.qpo['nu_char'] for r in qpos])
    nu_err = _yerr([r.qpo['nu_char_err'] for r in qpos])
    Q = np.array([r.qpo['Q'] for r in qpos])
    Q_err = _yerr([r.qpo['Q_err'] for r in qpos])
    amp = np.array([r.qpo['r'] for r in qpos])
    amp_err = _yerr([r.qpo['r_err'] for r in qpos])

    fig, axi = plt.subplots(3, 1, sharex=True, figsize=(10, 8))

    axi[0].errorbar(nums, nu, yerr=nu_err, fmt='o', color='C0')
    axi[0].set_yscale('log')
    axi[0].set_ylabel(r'$\nu_{\rm char}$ (Hz)')
    axi[0].set_title('QPO Frequency Evolution in MAXI J1820')

    axi[1].errorbar(nums, Q, yerr=Q_err, fmt='o', color='C1')
    axi[1].axhline(min_q, ls='--', color='gray', lw=1)
    axi[1].set_ylabel('Q factor')

    axi[2].errorbar(nums, amp, yerr=amp_err, fmt='o', color='C2')
    axi[2].set_ylabel('Amplitude (r)')
    axi[2].set_xlabel('ObsID suffix')

    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    return qpos


def fit_quality_figure(results, outpath):
    '''2x2 rebuild of fit_quality_analysis.pdf: reduced chi^2 vs obs,
    model complexity vs obs, chi^2 histogram, chi^2 vs complexity.'''
    nums = np.array([r.obs_num for r in results])
    red = np.array([r.red_chi2 for r in results])
    nlor = np.array([r.n_lorentzians for r in results])

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    ax = axes[0, 0]
    colors = np.where(red <= RED_CHI2_GOOD, 'green',
                      np.where(red <= RED_CHI2_MARGINAL, 'orange', 'red'))
    ax.scatter(nums, red, c=colors, s=25)
    ax.axhline(1.0, ls='--', color='blue', label=r'Perfect fit ($\chi^2$=1)')
    ax.axhline(RED_CHI2_GOOD, ls=':', color='green', label='Good threshold')
    ax.axhline(RED_CHI2_MARGINAL, ls=':', color='orange',
               label='Marginal threshold')
    ax.set_xlabel('ObsID suffix')
    ax.set_ylabel(r'Reduced $\chi^2$')
    ax.set_title('Fit Quality vs Observation')
    ax.legend()

    ax = axes[0, 1]
    ax.scatter(nums, nlor, s=25)
    ax.set_xlabel('ObsID suffix')
    ax.set_ylabel('Number of Lorentzians')
    ax.set_title('Model Complexity vs Observation')

    ax = axes[1, 0]
    ax.hist(red, bins=20, edgecolor='black', alpha=0.7)
    ax.axvline(1.0, ls='--', color='blue', label='Perfect fit')
    ax.axvline(RED_CHI2_GOOD, ls=':', color='green', label='Good threshold')
    med = float(np.median(red))
    ax.axvline(med, color='red', label=f'Median={med:.2f}')
    ax.set_xlabel(r'Reduced $\chi^2$')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Fit Quality')
    ax.legend()

    ax = axes[1, 1]
    sc = ax.scatter(nlor, red, c=nums, cmap='viridis', s=30)
    fig.colorbar(sc, ax=ax, label='ObsID suffix')
    ax.set_xlabel('Number of Lorentzians')
    ax.set_ylabel(r'Reduced $\chi^2$')
    ax.set_title('Fit Quality vs Model Complexity')

    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
