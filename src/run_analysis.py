#!/usr/bin/env python3
'''Rebuild the QPO-evolution and fit-quality figures from fit_results.json.

Example:
    python3 src/run_analysis.py --results results/fit_results.json \
        --out-dir results
'''

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from psfit.analysis import (fit_quality_figure, identify_qpos,  # noqa: E402
                            qpo_evolution_figure)
from psfit.results import load_results  # noqa: E402


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--results', default='results/fit_results.json',
                        help='fit_results.json from run_fitting.py')
    parser.add_argument('--out-dir', default='results',
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
    args = parser.parse_args()

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


if __name__ == '__main__':
    main()
