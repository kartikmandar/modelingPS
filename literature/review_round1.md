# Literature review: Lorentzian timing of NICER data for MAXI J1820+070 (2018 rise, ObsIDs 1200120101–197)

*For the author of `/home/user/modelingPS`. The review is built from the upstream web-search and verification passes plus a critic pass; I ran no new searches while writing it. Each verification status below is the one those passes reported. A source marked **[UNVERIFIED]** must be checked on ADS before you cite it. Figures noted "from memory" in the upstream notes are flagged where they appear.*

---

## 0. Executive summary

- **Coverage.** Your sample runs from about MJD 58189 to 58305.8. That spans the rising hard state (LHS), the long plateau, the hard-intermediate state (HIMS), and the type-C→type-B switch with the jet ejection on 2018-07-06. **ObsID 1200120197 contains the switch.** Unless you split it in time, it mixes type-C, the extra non-harmonic QPO and type-B intervals.
- **Validated choices.** Three choices match published practice: the Belloni–Psaltis–van der Klis (BPK02) parametrisation with ν_char = ν_max = √(ν0²+Δ²) and Q = ν0/2Δ, the (rms/mean)² normalisation, and a single 2–11 keV band. The QPO centroid frequency does not depend on energy from 0.2 to 200 keV, but Q and rms do. Most hard-state NICER spectra need about 4 Lorentzians, and 5–6 when a harmonic or subharmonic is present.
- **Main weaknesses** of the current statistics, from most to least serious:
  - The F-test for adding a component is not calibrated, because the null hypothesis puts the extra normalisation on a boundary (Protassov et al. 2002).
  - χ² with error bars taken from the data biases the QPO amplitude and the Poisson level low (Barret & Vaughan 2012).
  - COBYLA gives no Hessian. Hessian errors fail for weak components or parameters near a bound.
  - A fixed 2/⟨rate⟩ Poisson level is wrong by a few per cent at count rates of ~10⁴ c/s because of dead time.
- **Biggest publishable gap.** No paper applies time-dependent Comptonisation (vKompth) rms/lag modelling to the whole LHS→HIMS type-C sequence of this source. The only such work covers the C→B transition (Ma R. et al. 2023). A clean, statistically calibrated Lorentzian catalogue across all ObsIDs, compared quantitatively with the frequency correlations and with Li et al. (2024), would also be publishable on its own.

---

## 1. MAXI J1820+070 timing results to compare against

### 1.1 Outburst timeline (verified sources)

| Phase | MJD range | Source |
|---|---|---|
| Optical onset (ASASSN-18ey), 7.20 ± 0.97 d before the X-ray onset; first MAXI detection 32 ± 9 mCrab | ~58180–58189 | Tucker et al. 2018 |
| Rising hard state (LHS), including the plateau | ~58189 – ~58290 | Homan et al. 2020; Li et al. 2024 |
| HIMS | ~58290 – 58305.6 | Li et al. 2024 |
| Ejection launch (radio) | 58305.60 ± 0.04 (Li); ~2 h before the C→B switch (Wood 2021) | Li et al. 2024; Wood et al. 2021 |
| Last type-C detection | 58305.61 | Li et al. 2024 |
| Type-C→type-B switch | **58305.6774** (2018-07-06 16:15 UTC) | Homan et al. 2020 |
| Radio flare onset | 58305.773 ± 0.006 | Homan et al. 2020 |
| Type-B no longer seen | by ~58305.81 | Homan et al. 2020 |
| SIMS | ~58305.6 – 58310 | Li et al. 2024 |
| Soft state | 58310 – ~58383/58393 | Li et al. 2024 |
| Decaying hard state (outside your sample) | ~58393/58400 – 58443 | Homan et al. 2020 |

The hard→soft window of 58303.5–58310.7 attributed to Shidatsu et al. (2019) comes from a secondary citation. **[UNVERIFIED numbers]**

**ObsID mapping.**

| ObsID | Date / MJD | Used by |
|---|---|---|
| 103 | 2018-03-13 (MJD 58190) | Dziełak et al. 2021 |
| 104 | 2018-03-15 (MJD 58192) | Dziełak et al. 2021 |
| 106 | 03-21 (MJD 58198) | Kara et al. 2019 |
| 120 | 04-04 (MJD 58212) | Kara et al. 2019 |
| 130 | 04-16 (MJD 58224) | Kara et al. 2019 |
| 134/135 | 04-21/23 (MJD 58229–58231) | Kara et al. 2019 |
| 142 | 05-02 (MJD 58240) | Kara et al. 2019 |
| 148 | 05-08 (MJD 58246) | Kara et al. 2019 |
| 192–197 | MJD 58301–58305 (the IMS) | Li et al. 2024 |

Li et al. cut ObsIDs 196 and 197 into 12 segments. Ma R. et al. (2023) cut ObsID 197 into 15 orbit-length segments: 10 type-C and 3 type-B.

### 1.2 Published type-C QPO frequencies vs date

| MJD (approx.) | Date (2018) | ν_QPO (Hz) | Instrument | Source / notes |
|---|---|---|---|---|
| ~58190–58197 | March (INTEGRAL revolutions 1932–1934) | 0.030 → 0.035 → 0.042 (JEM-X 0.04) | INTEGRAL | Mereminskiy et al., ATel 11488 **[UNVERIFIED in the final pass]** |
| 58193.42 | Mar 16 | first type-C detection; full range 0.01–8.05 over the outburst | HXMT+NICER | Li et al. 2024 |
| early LHS | March (exact date not in the abstract) | 0.0477 | AstroSat/LAXPC | Mudambi et al. 2020 |
| 58190–58192 | Mar 13–15 | highest Lorentzian ν_max ≈ 2.5 (broadband, not the QPO) | NICER | Dziełak et al. 2021 |
| ~58198 (day 10) | Mar 21 | 0.03 and 0.08 (two related QPOs) | NICER/XRT | Stiele & Kong 2020 (T0 = MJD 58188) |
| ~58200–58290 | plateau | ~0.04–0.3 | NICER/XRT | Stiele & Kong 2020 |
| any t | — | log₁₀ f = 0.01869 t − 2.9158 (t in days from MJD 58119), e-folding ≈ 23 d; gives ≈0.04 Hz at MJD 58199 and ≈0.26 Hz at MJD 58289 | NuSTAR | Buisson et al. 2019 / ATel **[fit from a search snippet; UNVERIFIED]** |
| ~58304 (day 116) | Jul 5 | 1.57 and 3.17 | NICER | Stiele & Kong 2020 |
| ~58304.6 → 58305.6 | last 24.8 h | 1.57 → 3.94 | NICER | Stiele & Kong 2020 |
| 58305 | Jul 6 | 4.5 (+ 9) rising to ~8 (+ ~16) before 58305.677 | NICER | Homan et al. 2020 |
| 58305 | Jul 6 | type-C 4.4 → 7.7 | NICER | Ma R. et al. 2023 |
| ~58305.1 (≈0.5 d before the switch) | Jul 6 | extra **non-harmonic** QPO at 3.5–5.9 that turns into the type-B | NICER, 2–12 keV | Jin et al. 2026 |
| 58305.61 | Jul 6 | last type-C, ≈8.05 (top of the range) | HXMT+NICER | Li et al. 2024 |
| 58305.677–58305.81 | Jul 6 | **type-B** 4.5 → ~3 (Homan); 4.2 → 3.1 (Ma R.) | NICER | Homan et al. 2020; Ma R. et al. 2023 |

**Piecewise e-folding times of the type-C frequency** (Li et al. 2024; a direct quantitative benchmark):

| Phase | e-folding time |
|---|---|
| LHS rise | 7.63 ± 0.78 d |
| Plateau | 23.99 ± 0.65 d |
| Bright decline | −50.27 ± 7.86 d |
| IMS, first part | 8.28 ± 0.76 d |
| IMS, second part | 1.02 ± 0.05 d |

**Characteristic frequencies of the soft-reverberation lag** (not the QPO). Kara et al. (2019) found 2.8, 5.6, 7.9, 11.1, 15.7 and 22.0 Hz across ObsIDs 106, 120, 130, 134/135, 142 and 148, while the Fe line profile stayed constant. De Marco et al. (2021) report the same rise, then a sudden drop in lag frequency with a jump in amplitude near the transition.

**Outside your range, for context only.** Xu et al. (2020) followed the decay QPO from 11.6 ± 0.2 to 2.8 ± 0.1 mHz. Bellavita et al. (2025) found an "imaginary QPO" falling from about 5 to 1 Hz after MJD ~58393.

### 1.3 Published Lorentzian component counts

| Reference | Data | Components |
|---|---|---|
| Dziełak et al. 2021 | NICER ObsIDs 103/104, 2–10 keV | 4 Lorentzians (the 4th has zero normalisation at 0.3–2 keV); Poisson level fitted above 300 Hz over 0.005–1250 Hz |
| Stiele & Kong 2020 | NICER/XRT | ≤2 zero-centred BLN + ≤2 QPOs + a 3–8 Hz peaked noise |
| Li et al. 2024 | NICER | 1–3 BLN + 1–2 QPO |
| Méndez et al. 2024 | MAXI J1820+070 | PDS alone needs 4; the joint PDS + real/imaginary cross-spectrum fit needs **7** |
| Bhargava et al. 2021 | NICER, 0.01–12 keV | Belloni-style ν_max Lorentzians, used for a Relativistic Precession Model (RPM) spin of a = 0.799 (+0.016/−0.015); kept only observations with total rms ≥ ~10% |

The fourth line matters directly for your stopping rule: an F-test on the PDS alone undercounts components.

### 1.4 Specific traps for the tracking code

1. **Harmonic confusion.**
   - Stiele & Kong and Homan both see pairs of QPOs at about a 1:2 ratio.
   - The second harmonic is sometimes stronger than the fundamental (BPK02).
   - Fiori et al. (2025) find optical QPOs at about half the main X-ray QPO frequency at two hard-state epochs. On that reading the "X-ray type-C" could itself be a harmonic of the precession frequency.
   - Shui et al. (2023) show that the phase difference between harmonics changes with time.
   - Track the fundamental consistently and record your labelling rule.
2. **The extra non-harmonic QPO on the last day** (Jin et al. 2026, 3.5–5.9 Hz). An F-test can add this component, and a harmonic-only interpretation would label it wrongly.
3. **Mixed QPO types inside ObsIDs 196/197.** See Section 3.
4. **Short-timescale variability.**
   - The QPO is intermittent, with coherent stretches of about 5 cycles on average (Zhang P. et al. 2023).
   - QPO lags change over tens of seconds below ~0.3 Hz (Bollemeijer & Uttley 2025).
   - Averaging a long ObsID therefore smears ν0 and lowers Q.

---

## 2. Methodology and statistics: best practice vs the current pipeline

### 2.1 Parametrisation and conventions (keep these)

**Profile.** BPK02 write each component as P(ν) = (r²Δ/π)/[Δ² + (ν−ν0)²], with r the integrated fractional rms and Δ the half width at half maximum (HWHM). They define:
- ν_max = √(ν0²+Δ²), the peak of νP_ν. This equals your ν_char.
- Q = ν0/2Δ.

Two consequences:
- Using Lorentzians instead of broken power laws shifts the break frequency by ≲15%.
- No 1/f component is needed.

**Nowak (2000) convention.** Nowak uses Q = f0/FWHM in a different functional form. Convert before comparing numbers with Nowak or Pottschmidt et al. (2003).

**Q > 2 threshold.** It is a convention: van der Klis (2006) calls Q > 2 features QPOs and Q < 2 features peaked noise, and BPK02 call the threshold "arbitrary". Cite van der Klis 2006 or BPK02 for it, not Remillard & McClintock 2006, where the wording was not found.

**Labelling (Klein-Wolt & van der Klis 2008).**
- Find L_LF, the strongest narrow low-frequency QPO.
- Call the broad component below it L_b and the one above it L_h.
- Call a single high-frequency broad component L_l. Use L_u only when two high-frequency components are seen.
- Label weaker QPOs at about ×2 and ×½ of ν_LF as L_LF+ and L_LF−.
- When a broad component's ν0 runs negative, fix ν0 = 0.
- Keep a component if its integrated power is significant at >3σ or an F-test gives >3σ improvement.

**Correlation diagnostics.**
- Wijnands & van der Klis (1999): the hump/LF QPO sits about 2–8 times above ν_b, with the ratio falling as ν_b rises (van der Klis 2006).
- Psaltis, Belloni & van der Klis (1999): ν_LF (or ν_h) correlates with ν_l over about 3 decades.
- Plot your ν_max values over the BPK02 correlation tracks to catch misidentified components.

**Harmonic ratios.** van Doesburgh & van der Klis (2020) measure QPO2/QPO1 = 2.011 (GX 339-4), 1.975 (XTE J1550-564) and 1.988 (H1743-322), and QPO3/QPO1 ≈ 3.01. The offsets from integers are formally significant, largely because a Lorentzian does not match the true QPO profile. Do not force an exact 2:1 tie; fit the ratio, or tie it only as a test. Subharmonic ratios can be about 1.5.

### 2.2 How the current approach compares

| Element | Current | Issue | Best practice / fix |
|---|---|---|---|
| Objective | χ² with errors P/√(MW) | Data-derived weights bias the QPO rms and the Poisson level low. The bias is still ~2.5% in R at M = 64, against <0.3% for maximum likelihood at any M (Barret & Vaughan 2012; Humphrey et al. 2009, the general Poisson case). | Use the Whittle/Gamma likelihood for averaged spectra: S = 2 Σ_j M_j W_j [P_j/S_j + ln S_j]. If you keep χ², use model variances S_j/√(M_j W_j) and require M·W ≳ 50 in every bin, including the lowest log bins (van der Klis 1989; Papadakis & Lawrence 1993). |
| Optimiser | COBYLA (derivative-free, local) | No gradient and no Hessian. Multi-Lorentzian likelihood surfaces have several modes (component swaps, harmonic aliases). | Fit in log(r), log(Δ), log(ν0) or ν0 ≥ 0 with L-BFGS-B or another bounded optimiser. Use several starts, seeded by the previous ObsID's solution. For final runs, consider nested sampling (UltraNest, Buchner 2021). |
| Stopping rule | Analytic F-test | Not calibrated: under the null r = 0 sits on the boundary and ν0, Δ are undefined. Protassov et al. (2002) found real false-positive rates of 2.6% and 1.5% for nominal 5% tests, and the error can go either way. The PDS alone also undercounts components (4 vs 7, Méndez et al. 2024). | (a) Calibrate ΔS by parametric bootstrap from the null fit, with Gamma(MW) draws around the model (Vaughan 2010; Hurkett et al. 2008; Huppenkothen et al. 2013). Use 500–1000 simulations per decision, on a subset or as tabulated calibration curves. (b) As a cheap screen, AICc/BIC from the Whittle likelihood (Liddle 2007). (c) Bayesian evidence via nested sampling (Buchner et al. 2014). In every case add the Klein-Wolt criterion r/σ_r > 3, with σ_r from a profile or posterior lower bound, and Q > 2 for calling a QPO. |
| Goodness of fit | reduced χ² | Approximate for nonlinear models, because degrees of freedom are ill-defined (Andrae et al. 2010). | Keep it as a diagnostic only. |
| Errors | Hessian | Agrees with ΔS = 1 to within a few % for well-detected QPOs (Barret & Vaughan 2012). Fails for weak components, at bounds (r→0, Δ→0) and for degenerate pairs. | Finite-difference Hessian in log space. Compare with ΔS = 1 profiles or emcee (Hogg & Foreman-Mackey 2018) for the QPO's ν_max, Q and r, and flag disagreement above ~10%. |
| Poisson level | Subtracted | At ~10⁴ c/s dead time pushes it below 2/⟨rate⟩, and the size of the effect depends on flux. | Fit a free constant per ObsID (Bhargava et al. 2021), or take the mean power above 100–300 Hz (Ma R. et al. 2023; Dziełak et al. 2021). Propagate its uncertainty. |
| Stationarity | Averaging over an ObsID | A drifting QPO broadens the peak and inflates significance (Hübner et al. 2022, ApJS 259, 32). | Fit per GTI or orbit where S/N allows, and always on the last-day ObsIDs. |
| QPO rms | Additive sum of Lorentzians | On HXMT data a multiplicative (convolution) QPO–noise coupling is preferred in the time domain. ν0 and width barely change between the models, but the QPO rms does (Zhou et al. 2022). | Quote rms with this model caveat. ν0 and Q are robust. |

### 2.3 Log-space caveat

If you ever fit the logarithm of the power, raw-periodogram log powers sit 0.25068 dex below the true log power (Vaughan 2005). Fitting averaged powers linearly with the Whittle likelihood avoids this.

---

## 3. NICER reduction recipe for Phase 1 (nicerl2 with CALDB xti20240206)

**Software version.** Use HEASoft ≥ 6.32.1; 6.34 or later is better. 6.32.0 has a per-MPU normalisation bug in `nifpmsel`. Check for the `NIMPUFIX` keyword in the products (HEASARC 6.32 release notes).

**Current nicerl2 defaults** (filtcolumns = NICERV6):

| Parameter | Default |
|---|---|
| underonly_range | 0–500 |
| overonly_range | 0–30 |
| cor_range | \*-\* (no COR_SAX cut) |
| elv | 15 |
| br_earth | 30 |
| min_fpm | 7 |
| detlist | launch |

Geomagnetic data are read from `$GEOMAG_PATH`. `niautoscreen` and `nifpmsel` run by default since 6.32.

### 3.1 Suggested pipeline

Check the exact parameter names with `plist nicerl2` and `fhelp nicerl2` in your installed version.

```bash
export CALDB=...            # with xti20240206 installed
export GEOMAG_PATH=/path/to/geomag   # populated with nigeodown
for o in $(seq 1200120101 1200120197); do
  nicerl2 indir=$o clobber=YES \
     geomag_path=$GEOMAG_PATH \
     detlist=launch,-14,-34 \        # 50-FPM convention used in every MAXI J1820 paper found
     filtcolumns=NICERV6 \
     # leave under/overonly at defaults; run a second pass with cor_range=4-* (Kara et al. 2019 used COR_SAX>4)
     > logs/nicerl2_$o.log 2>&1
done
```

Then, for each ObsID:

1. **Energy band.** Select 2–11 keV, which is PI 200–1099 (NICER PI channel = 10 eV). This band avoids optical loading and noise problems below 0.5 keV and maximises the type-C rms, which rises with energy from ~1% to ~10% (Ma R. et al. 2023). **Also keep a 0.5–12 keV and a 0.5–2 keV product** for the follow-ups in Section 5. Barycentring with `barycorr` is optional for power spectra but needed for any cross-mission comparison.
2. **Segments within continuous GTIs only.**
   - Never FFT across a gap. The HEASARC Tips page uses MAXI J1820+070 ObsID 1200120107 as its example: 5912 GTIs from telemetry saturation, and including the gaps gives a spurious "QPO" near 55 Hz.
   - Use an integer number of segments per GTI (Bollemeijer et al. 2024 used 64 s).
   - Time resolution: dt ≤ 1/2048 s (Homan used 16 s × 1/2048 s; Bhargava 26.2144 s × 0.4 ms; Ma R. 16 s × 0.25 ms). This gives a Nyquist frequency of ≥ 1024 Hz for fitting the Poisson level.
   - Segment length: use **≥ 128–256 s for ObsIDs ~101–110**, where the QPO is at 0.01–0.05 Hz and the broad noise sits lower still. 16–64 s is enough later.
   - Keep the same segment length and dt as the 88 existing spectra in one run, so you can compare old and new spectra one by one.
3. **Normalisation and Poisson level.**
   - Use (rms/mean)² Hz⁻¹ normalisation (Belloni & Hasinger 1990; Miyamoto et al. 1991).
   - Fit the Poisson level per ObsID rather than assuming 2/⟨rate⟩. Typical per-event dead time is ~22 µs per FPM, rising to ~40 µs for very large events, and it freezes only the triggered FPM (Webbe & Young 2024; Remillard et al. 2022). Dead-time-correcting the light curves before the FFT brings no significant gain (Webbe & Young 2024). Stevens et al. (2018) fitted a level only 1% below 2/r at >16,000 c/s. A rough estimate (mine): at ~500 c/s per FPM, 2·r_det·τ ≈ 2%.
   - Note that dead time comes from all events in an FPM, including soft photons and over/undershoots. That is another reason to fit the level in-band.
   - As a cross-check on a few bright ObsIDs, compute the cospectrum or Fourier Amplitude Difference (FAD) between two disjoint MPU groups (Bachetti et al. 2015; Bachetti & Huppenkothen 2018).
4. **Background.** It is negligible: median 0.87 c/s (0.4–12 keV, 50 FPMs), >10 c/s in 5% of GTIs (Remillard et al. 2022), against 10³–2.5×10⁴ c/s from the source.
   - Still compute 3C50 (`nibackgen3c50`) or SCORPEON rates per ObsID and apply rms_src = rms_obs·(S+B)/S, with Poisson level 2(S+B)/S² (Uttley et al. 2014, Eq. 6). That lets you show the change is ≪1%.
   - `nicerl3-lc` offers only SCORPEON or Space Weather background, and none by default. The Space Weather model is still "Gendreau et al. in prep." and uses COR_SAX and Kp from the geomag files.
5. **Segment-level filters.**
   - Background flares: reject segments where the 13–15 keV rate exceeds 1 c/s per 52 FPMs (Homan et al. 2020). Dziełak used 2 c/s.
   - FPM outliers: reject a segment when any FPM's rate lies more than 4 × (IQR/2) from the median FPM rate (Stevens et al. 2018).
   - Drop GTIs shorter than ~10 s.
6. **ObsID-specific handling.**
   - **130 and 131** (2018-04-16/17): the ISS robotic arm crossed NICER's field of view, causing flux drops that appear as a shadow moving across the FPMs. Mask those intervals (Bollemeijer et al. 2024).
   - **Brightest phases:** the number of active FPMs was cut, in some cases to 17 (Homan et al. 2020), and some detectors were off from **ObsID 196 onward** (Wang et al. 2021). Record the active-FPM count per segment. Make sure no segment straddles a change in the FPM set. Normalise per FPM for any comparison of count rates.
   - **197** (and check 196): split by orbit or GTI. Label type-C before MJD 58305.677 and type-B for 58305.677–58305.81, and look for the extra QPO in the ~0.5 d before the switch (Homan et al. 2020; Ma R. et al. 2023; Jin et al. 2026).
7. **Robustness runs.** Repeat the component counting and QPO tracking with and without FPMs 14 and 34, and with and without COR_SAX > 4. HEASARC now says FPMs 14 and 34 "still produce good data most of the time" and recommends autoscreening, so the choice should be tested and stated.

---

## 4. Tools

| Tool | Use here | Reference |
|---|---|---|
| **Stingray 2** | Reads NICER event files. Builds averaged PDS and cross spectra with GTI handling, lags, coherence and rms/covariance spectra (~1 s per bright NICER observation). Its `modeling` module has a PSD log-likelihood that takes the number of averaged segments m, emcee sampling and posterior-predictive LRT calibration. Use it to validate or replace the home-made FFT and objective. | Huppenkothen et al. 2019; Bachetti et al. 2024 (JOSS) |
| HENDRICS | Command-line quick look (HENreadevents, HENfspec, HENlags); a check on your own spectra | Bachetti 2018, ascl:1805.019 |
| Cospectrum / FAD | Poisson-free checks of the noise level | Bachetti et al. 2015; Bachetti & Huppenkothen 2018 |
| UltraNest | Evidence-based choice of component number; full posteriors on ν_max and Q | Buchner 2021 |
| emcee recipes | MCMC errors and convergence checks | Hogg & Foreman-Mackey 2018 |
| nDspec | Joint spectral-timing modelling with a PSD-Lorentzian tutorial; A&A (accepted), alpha release | Lucchini et al. 2025, arXiv:2512.10615 |
| vKompth (XSPEC) | Fits QPO rms and lag spectra with time-dependent Comptonisation | Bellavita et al. 2022, ascl:2207.020 |
| Gaussian-process QPO search | Time-domain cross-check for short or gappy GTIs; mainly built for transients | Hübner et al. 2022 (ApJ 936, 17) |
| QPOML | Machine-learning prediction of QPO properties from spectra. Labels come from Lorentzian fits, so use it at most as an initial-guess provider | Kiker et al. 2023 |
| Cross-spectral errors | Coherence bias correction; energy-dependent cross-spectrum and rms errors | Vaughan & Nowak 1997; Ingram 2019; Uttley et al. 2014 |

A minimal Stingray sketch (check the exact signatures against your installed version):

```python
from stingray import EventList, AveragedPowerspectrum
ev = EventList.read("ni1200120150_0mpu7_cl.evt", fmt="hea")
ev = ev.filter_energy_range([2.0, 11.0])
ps = AveragedPowerspectrum.from_events(ev, segment_size=64, dt=1/4096, norm="frac")
```

---

## 5. Physical interpretation and publishable follow-ups

### 5.1 Framework

- **Geometric origin.** Type-C QPOs are read as Lense–Thirring precession of a hot inner flow inside a truncated disc (Ingram, Done & Fragile 2009; review by Ingram & Motta 2019). In the joint model of Ingram & Done (2011), the low-frequency broadband break is the viscous frequency at the truncation radius and the QPO is precession of the same flow, so **ν_b and ν_QPO should move together**. Your Lorentzian catalogue tests this directly. Ingram & van der Klis (2013) give the exact analytic power and lag spectra; Rapisarda et al. (2016) fit the model jointly to PDS and cross spectra.
- **Radiative layer.** vKompth (Karpouzas et al. 2020; García et al. 2021; Bellavita et al. 2022) explains the energy dependence of rms and lag through corona size, feedback onto the disc and the oscillation of kT_e. It complements precession: it describes how the oscillation is radiated, not what sets its frequency.
- **Inclination.** The source is at high inclination: jet axis 63 ± 3° (Atri et al. 2020) or 64 ± 5° (Wood et al. 2021). That predicts strong type-C rms (Motta et al. 2015) and QPO lags that change from hard to soft as ν rises (van den Eijnden et al. 2017). QPO amplitude correlates with jet inclination, and the amplitudes are reproduced only if the spin–orbit misalignment is ≳10–15° (Vincentelli et al. 2026). The measured misalignment is >40° (Poutanen et al. 2022), the natural setup for Lense–Thirring precession.
- **Constraints specific to MAXI J1820+070 that the discussion must address:**
  1. A contracting corona with a fixed inner disc, from reverberation and reflection (Kara et al. 2019; Buisson et al. 2019), against truncation-radius readings of the rising frequencies. Axelsson & Veledina (2021) and Kawamura et al. (2022) instead favour a truncated disc at a few tens of Rg.
  2. A type-C QPO detected above 200 keV with soft lags reaching ~0.9 s, read as a precessing small-scale jet (Ma X. et al. 2021).
  3. ν_QPO correlates with photon index Γ and anticorrelates with Compton optical depth (Li et al. 2024; compare the saturation in Titarchuk & Shaposhnikov 2009).
  4. Variable QPO lags and reduced disc–corona coherence below ~0.3 Hz (Bollemeijer & Uttley 2025).
  5. The corona grows in height about 5 d before the radio flare (Wang et al. 2021). A dual corona grows to ~10⁵ km horizontally, then collapses at the C→B switch, with feedback falling from ~0.2–0.3 to 0 (Ma R. et al. 2023).
  6. Timing-based spin estimates disagree: a ≈ 0.80 from RPM (Bhargava et al. 2021) against a ≲ 0.15 from Lense–Thirring with the optical subharmonic (Fiori et al. 2025).
- **System parameters for converting frequency to radius.**
  - Distance D = 2.96 ± 0.33 kpc; outburst peak ~15 ± 3% L_Edd (Atri et al. 2020).
  - M_BH = 5.73–8.34 M☉ at 95% for an orbital inclination of 66–81°, or 8.48 (+0.79/−0.72) M☉ if i = 63° (Torres et al. 2020).
  - P_orb = 0.68549 d and f(M) = 5.18 ± 0.15 M☉ (Torres et al. 2019).

### 5.2 Publishable follow-ups (most to least feasible)

1. **A statistically calibrated Lorentzian catalogue of the full LHS→HIMS rise.** Report ν0, Δ, r, Q and ν_max with MCMC or profile errors, as a machine-readable table, with:
   - piecewise e-folding fits of log ν_QPO(t) compared with Li et al. (2024);
   - correlation plots ν_LF vs ν_b (Wijnands–van der Klis) and ν_LF vs ν_l (Psaltis–Belloni–van der Klis);
   - QPO rms vs ν on the Motta et al. (2015) and Vincentelli et al. (2026) tracks.
   Stiele & Kong (2020) had only about six QPOs above 3σ (as stated by Li et al. 2024), so this is a real improvement.
2. **Segment-resolved analysis of the final ~1.5 days** (ObsIDs 192–197). Test the non-harmonic QPO of Jin et al. (2026) with a calibrated significance, and map the 1–2 Hz peaked noise that appears before the type-B (Homan et al. 2020).
3. **Joint fits of the PDS with the real and imaginary cross spectra** (Méndez et al. 2024). These recover components the PDS misses (7 vs 4 for this source) and test whether your F-test count is complete.
4. **Energy-resolved QPO rms and lag spectra for every type-C ObsID** (Ingram 2019 errors), extended below 2 keV, then vKompth dual-corona fits jointly with the time-averaged spectra. This fills the gap noted above, following the MAXI J1535-571 template of Zhang Y. et al. (2022) and the method of Ma R. et al. (2023).
5. **Spectral correlations.** Fit the time-averaged NICER spectra (Γ, kT_in, R_in) and test ν_QPO vs Γ and vs optical depth.
6. **Cheap extras.** rms–flux per ObsID (Uttley, McHardy & Vaughan 2005; Wang & Zhang 2024), harmonic phase differences (de Ruiter et al. 2019) and harmonic lags.

---

## 6. Ranked recommendations

1. **Split ObsID 1200120197, and check 196, by orbit or GTI.** Stop type-C tracking at MJD 58305.677 and label later segments type-B or post-transition. Mask the robotic-arm intervals in ObsIDs 130/131.
2. **Never FFT across GTI gaps.** Use integer segments within GTIs and record the active-FPM count. Watch for the spurious 55 Hz feature from ObsID 107's fragmented GTIs.
3. **Fit the Poisson level per ObsID** (a free constant, with the PDS extended to ≥1 kHz) instead of 2/⟨rate⟩, and propagate its uncertainty.
4. **Replace χ² with data-derived errors by the Whittle/Gamma(M·W) likelihood**, or at least use model variances and enforce M·W ≳ 50 per bin.
5. **Replace the analytic F-test** with a bootstrap-calibrated ΔS from the null model, or BIC/AICc from the Whittle likelihood, plus r/σ_r > 3 from a profile or posterior and Q > 2. Measure the false-positive rate on QPO-free simulations.
6. **Errors:** finite-difference Hessian in log parameters, cross-checked against ΔS = 1 or emcee for the QPO parameters. Flag disagreement above ~10% and parameters at bounds.
7. **Optimiser:** use multiple or warm starts from the neighbouring ObsID, a bounded gradient-based method, and log parametrisation. Fix ν0 = 0 for broad components that run negative.
8. **Labelling:** adopt the Klein-Wolt & van der Klis (2008) L_b/L_LF/L_h/L_l/L_u scheme. Flag pairs of Q > 2 components with a frequency ratio of about 2, 0.5 or 1.5. Do not force exact integer ratios. Account for Jin et al. (2026) on the last day.
9. **Validate the reprocessed spectra** against the published frequencies in Section 1.2 and against Stingray, and compare the new spectra with the current 88 one by one.
10. **Robustness:** with and without FPMs 14 and 34, and with and without COR_SAX > 4. Apply and report the background correction (S+B)/S on rms.
11. **Longer segments (≥128–256 s)** for the earliest ObsIDs, where the QPO is at ~0.01–0.05 Hz.
12. **Keep 0.5–2 keV and 0.5–12 keV products** from the reprocessing now, for the rms/lag and vKompth work later.

---

## Annotated bibliography

### MAXI J1820+070: timing, states and geometry
- **Stiele, H. & Kong, A. K. H. (2020)**, "A timing study of MAXI J1820+070 based on Swift/XRT and NICER monitoring in 2018/19", ApJ (889, 142 from memory). arXiv:1912.07625; doi:10.3847/1538-4357/ab64ef. The closest prior NICER QPO catalogue: two related QPOs rising from 0.03 to 1.57 Hz and from 0.08 to 3.17 Hz, ~116 d in the hard state, HIMS→soft in under 1 d.
- **Homan, J. et al. (2020)**, "A rapid change in X-ray variability and a jet ejection in the black hole transient MAXI J1820+070", ApJL 891, L29. arXiv:2003.01012; doi:10.3847/2041-8213/ab7932. Type-C at 4.5/9 Hz rising to ~8/16 Hz; type-B (4.5→3 Hz) from MJD 58305.677; radio flare at 58305.773; FPMs cut to as few as 17.
- **Li, Y., Yan, Z., Gao, C., Yu, W. (2024)**, "X-ray spectral and timing evolution during the 2018 outburst of MAXI J1820+070", MNRAS (accepted). arXiv:2407.08421. Type-C from MJD 58193.42 to 58305.61 over 0.01–8.05 Hz; piecewise e-folding times; state boundaries.
- **Dziełak, M. A., De Marco, B., Zdziarski, A. A. (2021)**, "A spectrally stratified hot accretion flow in the hard state of MAXI J1820+070", MNRAS (506, 2020 from memory). arXiv:2102.11635; doi:10.1093/mnras/stab1700. Four Lorentzians in NICER ObsIDs 103/104, added one at a time; 50 FPMs.
- **De Marco, B. et al. (2021)**, "The inner flow geometry in MAXI J1820+070 during hard and hard-intermediate states", A&A 654, A14. arXiv:2102.07811; doi:10.1051/0004-6361/202140567. Reverberation-lag frequency rises through the outburst, then drops near the transition; R_in > ~10 Rg until close to it.
- **Kara, E. et al. (2019)**, "The corona contracts in a black-hole transient", Nature 565, 198. arXiv:1901.03877; doi:10.1038/s41586-018-0803-x. Lag turnover frequency rises from 2.8 to 22 Hz across ObsIDs 106–148 with a constant Fe line; screening with COR_SAX > 4.
- **Buisson, D. J. K. et al. (2019)**, "MAXI J1820+070 with NuSTAR I…", MNRAS 490, 1350. arXiv:1909.04688; doi:10.1093/mnras/stz2681. Break and QPO frequencies rise while the reflection spectrum stays stable. The log-frequency fit is **[UNVERIFIED]**.
- **Wang, J. et al. (2021)**, "Disk, Corona, Jet Connection in the Intermediate State of MAXI J1820+070…", ApJL 910, L3. arXiv:2103.05616; doi:10.3847/2041-8213/abec79. Coronal height grows about 5 d before the radio flare; some detectors off from ObsID 196.
- **Wang, J. et al. (2022)**, "The NICER Reverberation Machine…", ApJ 930, 18. arXiv:2205.00928; doi:10.3847/1538-4357/ac6262. 88 automated groups (a coincidence with your 88 spectra), soft lags in 84 of them.
- **Ma, R., Méndez, M., García, F. et al. (2023)**, "A variable corona during the transition from type-C to type-B QPOs in MAXI J1820+070", MNRAS. arXiv:2307.12728; doi:10.1093/mnras/stad2284. Type-C 4.4→7.7 Hz, type-B 4.2→3.1 Hz; ObsID 197 split into 15 segments; dual-corona vKompth.
- **Jin, P., Méndez, M., García, F., Altamirano, D., Ma, R. (2026)**, "Discovery of a QPO non-harmonically related to the Type-C QPO in the HIMS of MAXI J1820+070", A&A (accepted). arXiv:2606.18004. Extra QPO at 3.5–5.9 Hz ~0.5 d before the HIMS→SIMS switch, evolving into the type-B.
- **Bellavita, C. et al. (2025)**, "The nature of an imaginary QPO in the soft-to-hard transition of MAXI J1820+070", A&A 696, A128. arXiv:2502.12283; doi:10.1051/0004-6361/202453092. Seen only in the decay, outside your range.
- **Méndez, M. et al. (2024)**, "Unveiling hidden variability components… using both the Fourier power and cross spectra", MNRAS (527, 9405 from memory). arXiv:2312.03476. For this source, 4 Lorentzians in the PDS alone vs 7 in the joint fit.
- **Kawamura, T. et al. (2022)**, "A full spectral-timing model… the low/hard state of MAXI J1820+070", MNRAS 511, 536. arXiv:2107.12517; doi:10.1093/mnras/stac045. A propagating-fluctuation model with a disc truncated at tens of Rg reproduces the double-humped PSD.
- **Ma, X. et al. (2021)**, "Discovery of oscillations above 200 keV in a black hole X-ray binary with Insight-HXMT", Nat. Astron. 5, 94. arXiv:2009.10607; doi:10.1038/s41550-020-1192-2. QPO soft lag reaches ~0.9 s at 150–200 keV; precessing-jet interpretation.
- **Ma, X. et al. (2023)**, "A detailed view of LFQPO in the broadband 0.2–200 keV with Insight-HXMT and NICER", ApJ 948, 116. arXiv:2303.00481; doi:10.3847/1538-4357/acc4c3. Centroid independent of energy; FWHM and rms depend on energy; lag-energy turns over above 1 Hz.
- **Zhou, D.-K. et al. (2022)**, "Determination of QPO properties in the presence of strong broad-band noise… MAXI J1820+070", MNRAS. arXiv:2206.12905; doi:10.1093/mnras/stac1789. Multiplicative coupling preferred; QPO rms depends on the model.
- **Yang, Z.-X. et al. (2022)**, "The accretion flow geometry of MAXI J1820+070 through broadband noise research with Insight-HXMT", ApJ. arXiv:2204.00739; doi:10.3847/1538-4357/ac63af. Multi-Lorentzian noise up to ~150 keV.
- **Bhargava, Y., Belloni, T. et al. (2021)**, "A timing-based estimate of the spin of the black hole in MAXI J1820+070", MNRAS 508, 3104. arXiv:2109.14371; doi:10.1093/mnras/stab2848. RPM spin a = 0.799; Poisson level fitted per observation.
- **Mudambi, S. P. et al. (2020)**, "Unveiling the temporal properties of MAXI J1820+070 through AstroSat observations", ApJL. arXiv:2001.00755; doi:10.3847/2041-8213/ab66bc. QPO at 47.7 mHz.
- **Gao, C., Yan, Z., Yu, W. (2023)**, "LFQPO in MAXI J1820+070: Revealing distinct Compton and Reflection Contributions", MNRAS. arXiv:2302.00915; doi:10.1093/mnras/stad434. The QPO in a 2–11 keV-like band mixes reflection and Comptonised modulation.
- **Zhang, P. et al. (2023)**, "Intermittent QPO properties of MAXI J1820+070 revealed by Insight-HXMT", A&A 677, A178. arXiv:2307.07709; doi:10.1051/0004-6361/202346309. Coherent stretches of ~5 cycles.
- **Shui, Q. C. et al. (2023)**, "A Phase-resolved View of the LFQPOs from MAXI J1820+070", ApJ. arXiv:2311.03635; doi:10.3847/1538-4357/acfc42. Phase difference between harmonics varies with time.
- **Fiori, M. et al. (2025)**, "Optical and X-ray timing analysis of the 2018–2020 outburst… MAXI J1820+070", A&A 697, A222. arXiv:2504.06001; doi:10.1051/0004-6361/202553899. Optical QPO at ~½ ν_X; Lense–Thirring spin a ≲ 0.15.
- **Bollemeijer, N. et al. (2024)**, "Evidence for a dynamic corona in the short-term time lags of MAXI J1820+070", MNRAS 528, 558. arXiv:2312.09835; doi:10.1093/mnras/stad3912. Robotic-arm occultations in ObsIDs 130/131.
- **Bollemeijer, N. & Uttley, P. (2025)**, "Variable QPO lags and reduced coherence between the disc and corona in MAXI J1820+070", MNRAS. arXiv:2507.18437; doi:10.1093/mnras/staf1238.
- **Bollemeijer, N., Uttley, P., You, B. (2026)**, "Broad-band spectral-timing: simultaneous NICER and HXMT… anticorrelation…", MNRAS. arXiv:2609.27608; doi:10.1093/mnras/stag1829. Fluxes below 1 keV and above 100 keV anticorrelated below ~0.02 Hz.
- **Axelsson, M. & Veledina, A. (2021)**, "Accretion geometry of MAXI J1820+070 probed by frequency-resolved spectroscopy", MNRAS 507, 2744. arXiv:2103.08795; doi:10.1093/mnras/stab2191.
- **Prabhakar, G. et al. (2022)**, "Accretion Scenario of MAXI J1820+070 during 2018 Outbursts with Multi-mission Observations", MNRAS. arXiv:2204.13363; doi:10.1093/mnras/stac1176. Some QPOs need a broad plus a narrow Lorentzian.
- **Paice, J. A. et al. (2019)**, "A black hole X-ray binary at ~100 Hz: multiwavelength timing of MAXI J1820+070 with HiPERCAM and NICER", MNRAS Letters (490, L62 from memory). arXiv:1910.04174; doi:10.1093/mnrasl/slz148. Optical lags X-rays by ~+165 ms.
- **Thomas, J. K. et al. (2022)**, "Synchronous X-ray/optical QPOs from MAXI J1820+070", MNRAS Letters. arXiv:2111.13642; doi:10.1093/mnrasl/slab132.
- **You, B. et al. (2021)**, "Insight-HXMT observations of jet-like corona in MAXI J1820+070", Nat. Commun. 12, 1025. arXiv:2102.07602; doi:10.1038/s41467-021-21169-5.
- **Wang, Y. & Zhang, S.-N. (2024)**, "RMS-flux slope in MAXI J1820+070…", ApJ (accepted). arXiv:2401.01515.
- **Xu, Y. et al. (2020)**, "Evidence for disk truncation at low accretion states of MAXI J1820+070…", ApJ. arXiv:2003.01778; doi:10.3847/1538-4357/ab7cdb.
- **Tucker, M. A. et al. (2018)**, "ASASSN-18ey: The Rise of a New Black-Hole X-ray Binary", ApJL 867, L9. arXiv:1808.07875; doi:10.3847/2041-8213/aae88a.
- **Shidatsu, M. et al. (2019)**, "X-ray and optical monitoring of state transitions in MAXI J1820+070", ApJ 874, 183. arXiv:1903.01686; doi:10.3847/1538-4357/ab09ff. The MJD window is **[UNVERIFIED]**.
- **Bright, J. S. et al. (2020)**, "An extremely powerful long-lived superluminal ejection from MAXI J1820+070", Nat. Astron. 4, 697. arXiv:2003.01083; doi:10.1038/s41550-020-1023-5.
- **Wood, C. M. et al. (2021)**, "The varying kinematics of multiple ejecta from MAXI J1820+070", MNRAS. arXiv:2105.09529; doi:10.1093/mnras/stab1479.
- **Atri, P. et al. (2020)**, "A radio parallax to the black hole X-ray binary MAXI J1820+070", MNRAS Letters 493, L81. arXiv:1912.04525; doi:10.1093/mnrasl/slaa010.
- **Torres, M. A. P. et al. (2019)**, "Dynamical confirmation of a black hole in MAXI J1820+070", ApJL 882, L21. arXiv:1907.00938; doi:10.3847/2041-8213/ab39df.
- **Torres, M. A. P. et al. (2020)**, "The binary mass ratio in the black hole transient MAXI J1820+070", ApJL. arXiv:2003.02360; doi:10.3847/2041-8213/ab863a.
- **Poutanen, J. et al. (2022)**, "Black hole spin-orbit misalignment in the X-ray binary MAXI J1820+070", Science 375, 874. arXiv:2109.07511; doi:10.1126/science.abl4679.
- **Mereminskiy, I. et al. (2018)**, ATel 11488. **[UNVERIFIED in the final pass]**. No refereed Mereminskiy paper on this source was found.

### Lorentzian methodology and QPO phenomenology
- **Belloni, Psaltis & van der Klis (2002)**, ApJ 572, 392. arXiv:astro-ph/0202213; doi:10.1086/340290. The parametrisation you use; 4–5 components per spectrum.
- **Nowak (2000)**, MNRAS 318, 361. arXiv:astro-ph/0005232; doi:10.1046/j.1365-8711.2000.03668.x.
- **Psaltis, Belloni & van der Klis (1999)**, ApJ 520, 262. arXiv:astro-ph/9902130.
- **Wijnands & van der Klis (1999)**, ApJ 514, 939 (volume and page from memory). arXiv:astro-ph/9810342; doi:10.1086/306993.
- **van der Klis (2006)**, "Rapid X-ray variability", in Compact Stellar X-ray Sources (CUP). arXiv:astro-ph/0410551.
- **Klein-Wolt & van der Klis (2008)**, ApJ 675, 1407. arXiv:0711.1274. The labelling recipe and the 3σ retention criterion.
- **van Straaten et al. (2002)**, ApJ 568, 912. arXiv:astro-ph/0107562.
- **Pottschmidt et al. (2003)**, A&A 407, 1039. arXiv:astro-ph/0202258; doi:10.1051/0004-6361:20030906. Four Lorentzians in the hard state of Cyg X-1.
- **van Doesburgh & van der Klis (2020)**, MNRAS 496, 5262. doi:10.1093/mnras/staa1867. Harmonic ratios near, but measurably off, 2 and 3.
- **Casella, Belloni & Stella (2005)**, ApJ 629, 403. arXiv:astro-ph/0504318. The A/B/C QPO classes.
- **de Ruiter et al. (2019)**, MNRAS 485, 3834. arXiv:1903.03135; doi:10.1093/mnras/stz665.
- **Ingram & Motta (2019)**, New Astron. Rev. 85, 101524. arXiv:2001.08758; doi:10.1016/j.newar.2020.101524.
- **Belloni & Motta (2016)**, "Transient Black Hole Binaries", ASSL 440, 61. arXiv:1603.07872; doi:10.1007/978-3-662-52859-4_2.
- **Motta et al. (2015)**, MNRAS 447, 2059. arXiv:1404.7293; doi:10.1093/mnras/stu2579.
- **van den Eijnden et al. (2017)**, MNRAS 464, 2643. arXiv:1610.03469.
- **Vincentelli et al. (2026)**, "The strength of Type-C QPOs… correlates with the jet inclination", MNRAS (accepted). arXiv:2605.27510.

### Statistics
- **Protassov et al. (2002)**, ApJ 571, 545. arXiv:astro-ph/0201547; doi:10.1086/339856. The F-test and LRT are invalid for components whose null value lies on a boundary.
- **Barret & Vaughan (2012)**, ApJ 746, 131. arXiv:1112.0535; doi:10.1088/0004-637X/746/2/131. Whittle likelihood vs χ² bias; Hessian vs ΔS errors.
- **Vaughan (2010)**, MNRAS 402, 307. arXiv:0910.2706; doi:10.1111/j.1365-2966.2009.15868.x. Posterior-predictive LRT.
- **Vaughan (2005)**, A&A 431, 391. arXiv:astro-ph/0412697; doi:10.1051/0004-6361:20041453. The 0.25068 dex bias of log powers.
- **Papadakis & Lawrence (1993)**, MNRAS 261, 612. doi:10.1093/mnras/261.3.612.
- **van der Klis (1989)**, "Fourier techniques in X-ray timing", NATO ASI C 262, 27. https://dare.uva.nl/record/218477
- **Humphrey, Liu & Buote (2009)**, ApJ 693, 822. arXiv:0811.2796; doi:10.1088/0004-637X/693/1/822.
- **Cash (1979)**, ApJ 228, 939. doi:10.1086/156922.
- **Hurkett et al. (2008)**, ApJ 679, 587. arXiv:0802.0441; doi:10.1086/586881.
- **Huppenkothen et al. (2013)**, ApJ 768, 87. doi:10.1088/0004-637X/768/1/87. The arXiv id is inconsistent between the upstream passes (1212.1011 vs 1303.6648); **[check the arXiv id]**.
- **Hübner et al. (2022a)**, "Pitfalls of Periodograms…", ApJS 259, 32. arXiv:2108.07418.
- **Hübner et al. (2022b)**, "Searching for QPOs… using Gaussian Processes", ApJ 936, 17. arXiv:2205.12716.
- **Liddle (2007)**, MNRAS 377, L74. arXiv:astro-ph/0701113; doi:10.1111/j.1745-3933.2007.00306.x.
- **Buchner et al. (2014)**, A&A 564, A125. arXiv:1402.0004; doi:10.1051/0004-6361/201322971.
- **Andrae, Schulze-Hartung & Melchior (2010)**, arXiv:1012.3754.
- **Hogg & Foreman-Mackey (2018)**, ApJS 236, 11. arXiv:1710.06068; doi:10.3847/1538-4365/aab76e.

### NICER reduction, dead time and normalisation
- **Remillard et al. (2022)**, "An Empirical Background Model for the NICER X-ray Timing Instrument" (3C50), AJ 163, 130. arXiv:2105.09901; doi:10.3847/1538-3881/ac4ae6.
- **Webbe & Young (2024)**, "The Effects of Instrumental Deadtime on NICER Timing Products", RASTI 3, 593. arXiv:2409.12574; doi:10.1093/rasti/rzae040.
- **Stevens et al. (2018)**, ApJL 865, L15. arXiv:1809.07556; doi:10.3847/2041-8213/aae1a4. FPM 14/34 exclusion and the IQR outlier filter.
- **Zhang W. et al. (1995)**, ApJ 449, 930. doi:10.1086/176111. Dead-time modification of FFT power spectra.
- **Bachetti et al. (2015)**, ApJ 800, 109. arXiv:1409.3248. The cospectrum.
- **Bachetti & Huppenkothen (2018)**, ApJL 853, L21. arXiv:1709.09700; doi:10.3847/2041-8213/aaa83b. The FAD method.
- **Uttley et al. (2014)**, A&ARv 22, 72. arXiv:1405.6575; doi:10.1007/s00159-014-0072-0. Use v2 or later.
- **Gendreau et al. (2016)**, Proc. SPIE 9905, 99051H. doi:10.1117/12.2231304.
- **HEASARC documentation**: nicerl2 help, NICER Analysis Tips & Caveats, HEASoft 6.32 and 6.31.1 release threads, FPM selection thread, background thread.
  - https://heasarc.gsfc.nasa.gov/lheasoft/ftools/headas/nicerl2.html
  - https://heasarc.gsfc.nasa.gov/docs/nicer/data_analysis/nicer_analysis_tips.html
  - https://heasarc.gsfc.nasa.gov/docs/nicer/analysis_threads/heasoft632
  - https://heasarc.gsfc.nasa.gov/docs/nicer/analysis_threads/fpmsel-using
  - https://heasarc.gsfc.nasa.gov/docs/nicer/analysis_threads/background
- **Belloni & Hasinger (1990)**, A&A 227, L33. Normalisation; verified only through citing papers.

### Tools and physical models
- **Huppenkothen et al. (2019)**, "Stingray", ApJ 881, 39. arXiv:1901.07681; doi:10.3847/1538-4357/ab258d.
- **Bachetti et al. (2024)**, "Stingray 2", JOSS 9, 7389. doi:10.21105/joss.07389.
- **Bachetti (2018)**, HENDRICS, ascl:1805.019.
- **Buchner (2021)**, UltraNest, JOSS 6, 3001. doi:10.21105/joss.03001; arXiv:2101.09604.
- **Lucchini et al. (2025)**, nDspec, A&A (accepted). arXiv:2512.10615.
- **Kiker et al. (2023)**, QPOML, MNRAS 524, 4801. arXiv:2306.04055; doi:10.1093/mnras/stad1643.
- **Vaughan & Nowak (1997)**, ApJL 474, L43. arXiv:astro-ph/9610257; doi:10.1086/310430.
- **Nowak et al. (1999)**, ApJ 510, 874. arXiv:astro-ph/9807278.
- **Ingram (2019)**, "Error formulae for the energy-dependent cross-spectrum", MNRAS 489, 3927. arXiv:1909.01385.
- **Vaughan et al. (2003)**, MNRAS 345, 1271. arXiv:astro-ph/0307420; doi:10.1046/j.1365-2966.2003.07042.x.
- **Stevens & Uttley (2016)**, MNRAS 460, 2796. arXiv:1605.01753; doi:10.1093/mnras/stw1093.
- **Ingram, Done & Fragile (2009)**, MNRAS 397, L101. arXiv:0901.1238.
- **Ingram & Done (2011)**, MNRAS 415, 2323. arXiv:1101.2336.
- **Ingram & van der Klis (2013)**, MNRAS 434, 1476. arXiv:1306.3823; doi:10.1093/mnras/stt1107.
- **Ingram et al. (2016)**, MNRAS 461, 1967. arXiv:1607.02866; doi:10.1093/mnras/stw1245. Fe-line centroid modulated over the QPO cycle.
- **Rapisarda et al. (2016)**, MNRAS 462, 4078. arXiv:1607.08178.
- **Uttley, McHardy & Vaughan (2005)**, MNRAS 359, 345. arXiv:astro-ph/0502112; doi:10.1111/j.1365-2966.2005.08886.x.
- **Karpouzas et al. (2020)**, MNRAS 492, 1399. arXiv:1912.05380; doi:10.1093/mnras/stz3502.
- **García et al. (2021)**, MNRAS 501, 3173. arXiv:2012.10354.
- **Bellavita et al. (2022)**, vKompth, MNRAS 515, 2099. arXiv:2206.13609; doi:10.1093/mnras/stac1922; ascl:2207.020.
- **Méndez et al. (2022)**, Nat. Astron. 6, 577. doi:10.1038/s41550-022-01617-y. This is about GRS 1915+105, not MAXI J1820+070.
- **Zhang, Y. et al. (2022)**, MAXI J1535-571 corona evolution with type-C QPOs, MNRAS 512, 2686. arXiv:2203.05308.
- **Titarchuk & Shaposhnikov (2009)**, ApJ 699, 453. arXiv:0902.2852. Author order is from memory.

### Flagged / unverified: do not cite without checking
- Remillard & McClintock (2006), ARA&A 44, 49, arXiv:astro-ph/0606352. The Q > 2 wording was not found in it.
- Zdziarski et al. (2021), ApJL 909, L9, arXiv:2101.04482. Authors and volume from memory.
- Belloni (2010), Lecture Notes in Physics 794, 53, arXiv:0909.2474. Id from memory.
- Belloni & Hasinger (1990b), A&A 230, 103.
- Miyamoto et al. (1991), ApJ 383, 784. Bibliographic record confirmed; authors and normalisation formula not checked.
- Prigozhin et al. (2016), SPIE 9905, 99051I. The DOI–paper pairing is unconfirmed.
- Mereminskiy ATel 11488.
- The Buisson ATel frequency fit.
- The Shidatsu et al. (2019) MJD window.
- No "Cavecchi et al." paper on MAXI J1820+070 was found. Cavecchi's relevant work is as a coauthor on Vincentelli et al. 2026.
- Two internal inconsistencies resolved here:
  - The Topic 3 notes say the C→B switch falls after ObsID 197. Homan (2020) and Ma R. (2023) place it **inside** ObsID 197.
  - The task's "Méndez et al. 2022 for MAXI J1820" is the GRS 1915+105 paper. The MAXI J1820+070 vKompth paper is Ma R. et al. (2023).