'''Command-line entry points for modelingPS.'''

import argparse
import os
import sys

from .analysis import fit_quality_figure, identify_qpos, qpo_evolution_figure
from .pipeline import fit_directory
from .results import load_results


DEFAULT_DATA_DIR = 'data/raw/maxij1820_powspecs'
DEFAULT_OUT_DIR = 'outputs/current'


def build_fit_parser():
    '''Build the parser for fitting a directory of power spectra.'''
    parser = argparse.ArgumentParser(
        description='Run the Lorentzian fitting pipeline over power spectra.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--data-dir', default=DEFAULT_DATA_DIR,
                        help='directory with the *pow.dat power spectra')
    parser.add_argument('--pattern', default='*pow.dat',
                        help='glob pattern for the data files')
    parser.add_argument('--out', default=DEFAULT_OUT_DIR,
                        help='output directory')
    parser.add_argument('--noise-mode', default='quadrature',
                        choices=['none', 'quadrature', 'mask'],
                        help='how the effective noise level enters the weights')
    parser.add_argument('--search-qpo', type=int, default=2,
                        choices=[0, 1, 2],
                        help='0: no QPO, 1: frequency scan, 2: max of nu f(nu)')
    parser.add_argument('--method', default='COBYLA',
                        help='scipy.optimize.minimize method')
    parser.add_argument('--no-plots', action='store_true',
                        help='skip the per-observation plots')
    parser.add_argument('--legacy-txt', default='resu_quadrature.txt',
                        help='legacy-format results filename ("" to skip)')
    return parser


def fit_main(argv=None):
    '''CLI entry point for fitting all spectra in a directory.'''
    parser = build_fit_parser()
    args = parser.parse_args(argv)

    results = fit_directory(args.data_dir, args.out, pattern=args.pattern,
                            noise_mode=args.noise_mode,
                            search_qpo=args.search_qpo, method=args.method,
                            plot=not args.no_plots,
                            legacy_txt=args.legacy_txt or None)
    print(f'\n{len(results)} observations fitted -> {args.out}/')


def build_analyze_parser():
    '''Build the parser for rebuilding post-fit analysis figures.'''
    parser = argparse.ArgumentParser(
        description='Rebuild QPO-evolution and fit-quality figures.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--results',
                        default=os.path.join(DEFAULT_OUT_DIR, 'fit_results.json'),
                        help='fit_results.json from modelingps-fit')
    parser.add_argument('--out-dir', default=DEFAULT_OUT_DIR,
                        help='where to write the figures')
    parser.add_argument('--min-q', type=float, default=2.0,
                        help='minimum Q for a component to count as a QPO')
    parser.add_argument('--min-width', type=float, default=1.e-3,
                        help='minimum QPO width in Hz; sub-resolution '
                             'widths (data bins start at 2.5e-3 Hz) are '
                             'unresolved noise spikes (0 to disable)')
    parser.add_argument('--min-significance', type=float, default=3.0,
                        help='minimum amplitude significance r/r_err when '
                             'the error is available (0 to disable)')
    return parser


def analyze_main(argv=None):
    '''CLI entry point for rebuilding post-fit analysis figures.'''
    parser = build_analyze_parser()
    args = parser.parse_args(argv)

    results, meta = load_results(args.results)
    print(f'{len(results)} fits loaded (noise_mode={meta.get("noise_mode")})')

    qpos = qpo_evolution_figure(
        results, os.path.join(args.out_dir, 'qpo_evolution.pdf'),
        min_q=args.min_q, min_width=args.min_width,
        min_significance=args.min_significance)
    loose = identify_qpos(results, min_q=args.min_q)
    dropped = [r.obs_num for r in loose if r not in qpos]
    print(f'{len(qpos)} QPO detections (Q > {args.min_q}, '
          f'width >= {args.min_width} Hz, '
          f'>= {args.min_significance} sigma where measurable)')
    if dropped:
        print(f'dropped as unresolved/insignificant: obs {dropped}')

    fit_quality_figure(
        results, os.path.join(args.out_dir, 'fit_quality_analysis.pdf'))
    print('figures written')


def main(argv=None):
    '''Dispatch subcommands for python -m modelingPS.cli.'''
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in {'-h', '--help'}:
        print('usage: python -m modelingPS.cli {fit,analyze} ...')
        print()
        print('subcommands:')
        print('  fit      run the Lorentzian fitting pipeline')
        print('  analyze  rebuild post-fit analysis figures')
        return

    command, rest = argv[0], argv[1:]
    if command == 'fit':
        fit_main(rest)
    elif command == 'analyze':
        analyze_main(rest)
    else:
        raise SystemExit(f'unknown command {command!r}; expected fit or analyze')


if __name__ == '__main__':
    main()
