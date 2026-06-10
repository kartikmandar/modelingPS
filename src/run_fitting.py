#!/usr/bin/env python3
'''Run the Lorentzian fitting pipeline over a directory of power spectra.

Example:
    python3 src/run_fitting.py --data-dir maxij1820_powspecs \
        --out results --noise-mode quadrature
'''

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from psfit.pipeline import fit_directory  # noqa: E402


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--data-dir', default='maxij1820_powspecs',
                        help='directory with the *pow.dat power spectra')
    parser.add_argument('--pattern', default='*pow.dat',
                        help='glob pattern for the data files')
    parser.add_argument('--out', default='results',
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
    parser.add_argument('--legacy-txt', default='resu.txt',
                        help='legacy-format results filename ("" to skip)')
    args = parser.parse_args()

    results = fit_directory(args.data_dir, args.out, pattern=args.pattern,
                            noise_mode=args.noise_mode,
                            search_qpo=args.search_qpo, method=args.method,
                            plot=not args.no_plots,
                            legacy_txt=args.legacy_txt or None)
    print(f'\n{len(results)} observations fitted -> {args.out}/')


if __name__ == '__main__':
    main()
