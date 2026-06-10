# modelingPS — Lorentzian decomposition of MAXI J1820+070 power spectra

Fits sums of Lorentzians to NICER power density spectra of the black-hole
X-ray binary MAXI J1820+070 (2018 outburst, hard state → HIMS) and tracks
the type-C QPO through the outburst. Continuation of Yuri Cavecchi's
Barcelona code; collaboration context in `docs/emailConvo.txt`.

## Data

`maxij1820_powspecs/1200120NNNpow.dat` — 88 spectra (NNN = 101–197 with
gaps at 128, 129, 151, 153–155, 181, 184, 188), 226 bins each after the
`READ SERR 1 2` header. Columns: frequency, zeros (meaningless), power,
power error, effective noise level. Normalization (rms/mean)², 2–11 keV,
Poisson noise subtracted. Signal below the effective noise level is
unreliable (uncertainty of the subtracted Poisson level).

## Usage

Use a Python with numpy/scipy/matplotlib (e.g. the `pspec` conda env):

    # fit all spectra (noise-aware weights) and write results/ products
    python src/run_fitting.py --data-dir maxij1820_powspecs \
        --out results --noise-mode quadrature

    # rebuild the QPO-evolution and fit-quality figures
    python src/run_analysis.py --results results/fit_results.json \
        --out-dir results

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

The incremental fitter is a verbatim port of the original crawler — a
regression test asserts bit-identical results against `crawler_fix.py` —
so `--noise-mode none --search-qpo 0` reproduces the legacy behaviour
exactly.

## Layout

- `src/psfit/` — the package: `models` (Lorentzian, ν_char, Q),
  `observation` (data loading + noise modes + χ²), `fitting`
  (incremental fit, F-test, QPO searches), `uncertainties` (Hessian +
  delta method), `results` (JSON + legacy txt), `pipeline` (batch run +
  plots), `analysis` (QPO evolution + fit quality)
- `src/run_fitting.py`, `src/run_analysis.py` — CLI runners
- `src/crawler_fix.py` — legacy monolith (kept working; only the
  `qpo_try_search` Q bug was fixed: the QPO width is `x[-1]`, not `x[1]`)
- `src/crawler_fix _final.py` — untouched collaborators' reference copy
- `tests/` — pytest suite (`slow` marker = real COBYLA fits)
- `results/` — current products: `fit_results.json` (full precision +
  uncertainties), `resu_quadrature.txt` (legacy format), `summary.pdf`,
  `plots/*.pdf`, `qpo_evolution.pdf`, `fit_quality_analysis.pdf`
- `results_legacy_weights/` — same pipeline with `--noise-mode none`
  for comparison
- `resu_reb.txt`, `summa_reb.pdf`, `maxij1820_powspecs/*.pdf` and the
  root-level `qpo_evolution.pdf` / `fit_quality_analysis.pdf` — frozen
  February 2026 reference run from the pre-refactor crawler
