'''Reading and writing fit results.

The JSON files use Python's json module, which emits NaN literals for
missing uncertainties; they are read back by the same module, so the
files are self-consistent within this package.
'''

import json
import os
from dataclasses import asdict, dataclass


@dataclass
class FitResult:
    '''Everything the pipeline knows about one fitted observation.

    params/errors are flat [r, D, nu0] * n_lorentzians; when a QPO was
    searched it is the LAST triplet, summarized in the qpo dict
    together with the delta-method nu_char and Q errors.
    '''

    infile: str
    obs_num: int
    params: list
    errors: list
    n_lorentzians: int
    chi2: float
    dof: int
    red_chi2: float
    noise_mode: str
    search_qpo: int
    qpo: dict | None = None


def obs_num_from_path(path):
    '''1200120101pow.dat -> 101 (the last three digits of the ObsID).'''
    digits = ''.join(ch for ch in os.path.basename(path) if ch.isdigit())
    return int(digits[-3:])


def save_results(results, path, meta=None):
    '''Write a list of FitResults (and optional meta dict) as JSON.'''
    payload = {'meta': meta or {}, 'results': [asdict(r) for r in results]}
    with open(path, 'w') as fh:
        json.dump(payload, fh, indent=2)


def load_results(path):
    '''Read results written by save_results. Returns (results, meta).'''
    with open(path) as fh:
        payload = json.load(fh)
    results = [FitResult(**entry) for entry in payload['results']]
    return results, payload.get('meta', {})


def write_legacy_txt(results, path):
    '''resu_reb.txt-compatible output: one line per file, parameters at
    2 significant figures (kept identical to the original crawler).'''
    with open(path, 'w') as fh:
        for res in results:
            txt = res.infile
            for v in res.params:
                txt += '  {0: 8.2}'.format(v)
            fh.write(txt + '\n')
