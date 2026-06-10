# modelingPS — Lorentzian decomposition of MAXI J1820+070 power spectra

Fits sums of Lorentzians to NICER power density spectra of the black-hole
X-ray binary MAXI J1820+070 (2018 outburst, hard state -> HIMS) and tracks
the type-C QPO through the outburst. This package is a tested refactor of
Yuri Cavecchi's Barcelona crawler, with the legacy scripts retained for
regression comparison.

## Data

`data/raw/maxij1820_powspecs/1200120NNNpow.dat` contains 88 spectra
(NNN = 101-197 with gaps at 128, 129, 151, 153-155, 181, 184, 188).
Most files have 226 bins after the `READ SERR 1 2` header; 10 have 210
bins. Columns are frequency, zeros (meaningless), power, power error and
effective noise level. Normalization is (rms/mean)^2, 2-11 keV, with
Poisson noise already subtracted. Signal below the effective noise level
is unreliable because of uncertainty in the subtracted Poisson level.

## Usage

Use Pixi for the reproducible project environment:

    # tests (fast)
    pixi run test-fast
    # full suite, incl. real fits (~minutes)
    pixi run test

    # fit all spectra (noise-aware weights) and write outputs/current products
    pixi run fit

    # rebuild the QPO-evolution and fit-quality figures
    pixi run analyze

The package also exposes console scripts when installed in an equivalent
environment:

    modelingps-fit --data-dir data/raw/maxij1820_powspecs \
        --out outputs/current --noise-mode quadrature \
        --legacy-txt resu_quadrature.txt

    modelingps-analyze --results outputs/current/fit_results.json \
        --out-dir outputs/current

    # tests (fast)
    python -m pytest -m "not slow" -q
    # full suite, incl. real fits (~minutes)
    python -m pytest -q

`--noise-mode none` reproduces the legacy weighting (1/err²);
`quadrature` adds the effective noise level in quadrature to the errors;
`mask` zero-weights bins whose power is below the noise level.

## Method

1. Lorentzians are added one at a time (COBYLA, bounds 0 ≤ r ≤ 1,
   0 ≤ D, ν₀ ≤ 100); an F-test on the reduced-χ² ratio stops the loop.
2. A QPO Lorentzian is then added: at the max of νf(ν) with everything
   else frozen, followed by a free unconstrained refit (`--search-qpo 2`,
   default), or via a frequency scan (`--search-qpo 1`).
3. 1σ parameter errors come from the finite-difference Hessian of χ²
   (cov = 2H⁻¹); ν_char = √(ν₀²+D²) and Q = ν₀/2D errors via the delta
   method. Directions where the curvature is unusable get NaN.
4. A component counts as a QPO when Q > 2 with positive r, D, ν₀. The
   production figures additionally require the width to be resolvable
   (≥ 10⁻³ Hz; the data bins start at 2.5×10⁻³ Hz — narrower "QPOs" are
   single-bin noise spikes) and ≥ 3σ amplitude significance where the
   error is measurable (`--min-width` / `--min-significance`, 0 disables).

The incremental fitter is a verbatim port of the original crawler; a
regression test asserts bit-identical results against `legacy/crawler_fix.py`,
so `--noise-mode none --search-qpo 0` reproduces the legacy behaviour
exactly.

## Layout

- `src/modelingPS/` — the package: `models` (Lorentzian, ν_char, Q),
  `observation` (data loading + noise modes + χ²), `fitting`
  (incremental fit, F-test, QPO searches), `uncertainties` (Hessian +
  delta method), `results` (JSON + legacy txt), `pipeline` (batch run +
  plots), `analysis` (QPO evolution + fit quality), `cli` (console entry
  points)
- `pyproject.toml` — package metadata and `modelingps-fit` /
  `modelingps-analyze` entry points
- `legacy/crawler_fix.py` — legacy monolith (kept working; only the
  `qpo_try_search` Q bug was fixed: the QPO width is `x[-1]`, not `x[1]`)
- `legacy/crawler_fix_final.py` — untouched collaborators' reference copy
- `data/raw/maxij1820_powspecs/` — tracked input spectra (`*.dat`)
- `data/external/` — original compressed data bundle
- `tests/` — pytest suite (`slow` marker = real COBYLA fits)
- `pixi.toml`, `pixi.lock` — reproducible Pixi environment and task lockfile
- `outputs/current/` — current products: `fit_results.json` (full precision +
  uncertainties), `resu_quadrature.txt` (legacy format), `summary.pdf`,
  `plots/*.pdf`, `qpo_evolution.pdf`, `fit_quality_analysis.pdf`
- `outputs/legacy_weights/` — same pipeline with `--noise-mode none`
  for comparison
- `outputs/reference_legacy/` — frozen February 2026 reference run from the
  pre-refactor crawler, including the old per-observation PDFs
- `docs/` — local ignored notes/research, intentionally not part of git
