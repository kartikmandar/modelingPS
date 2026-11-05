# Literature review, round 3: gap-fill and novelty check

*For the author of `/home/user/modelingPS`. This round fills the open queries from round 2, Section 12.3, with six targeted search passes plus a synthesis and spot-check pass. It does not repeat sources that rounds 1 and 2 already cover in full; where a known source is mentioned again, it is because there is new detail or a correction. Verification status: "verified" means the paper or page was seen on a primary page (arXiv abstract or PDF, VizieR, HEASARC, a journal page) or on a search-result page by one of the passes. Section 8.2 lists everything that is partly verified or unverified. The synthesis pass re-checked four things directly:*
- *the VizieR J/ApJ/889/142 table list;*
- *the Bhargava et al. 2021 author list on arXiv;*
- *the NICER Science Nugget of 2018-03-22 (tentative 66 Hz);*
- *the ATel #11960 retraction title in the ATel search index (the ATel pages themselves are blocked).*

*It also checked the code facts quoted below (`pixi.toml` platforms; COBYLA, `scipy.stats.f` and the finite-difference Hessian in `src/modelingPS/`) against the repository.*

---

## 0. Summary: what changed compared with rounds 1 and 2

### 0.1 Novelty verdict

**A per-epoch NICER multi-Lorentzian catalogue of the 2018 rise already exists, and it is machine-readable.** It is Stiele & Kong 2020 (ApJ 889, 142), with all six tables on VizieR as J/ApJ/889/142.

The NICER part covers:
- Table 3: 196 noise-component rows.
- Table 4: 212 QPO rows.
- Table 5: 39 rows for an extra Lorentzian.
- Table 6: 6 rows for extra QPOs or peaked noise.

The rows run from day 1.58 to day 116.9 after MJD 58188, which is about 90 epochs and essentially our 88 ObsIDs. Rounds 1 and 2 cited this paper only as a QPO-frequency track. Round 2 concluded that "a full Lorentzian catalogue of the 2018 rise has not been published as an open … release". **That conclusion no longer holds and must not appear in any paper or proposal.**

Further partial per-ObsID tables now confirmed:
- Bollemeijer et al. 2024: Table D1 plus `qpofitsc.txt` in their Zenodo package, 63 ObsIDs from 126 to 189, no errors.
- Bhargava et al. 2021: Table 1, 21 ObsIDs from 130 to 191, with errors.
- Zhang Y. et al. 2026: Table A1, 27 ObsIDs from 138 to 197, QPO frequency and hue only.
- Bollemeijer & Uttley 2025: Table 1, 16 ObsIDs, ν_QPO and hue.
- Ma X. 2023: 6 epochs.
- Li et al. 2024: Table C1, ObsIDs 196/197 only.

**Not novel:** coverage of the rise; "the first per-ObsID Lorentzian/QPO table"; the ν_QPO(t) track; the ν_QPO–hue track.

**Still defensible (Section 1.3):**
1. The first ObsID- and MJD-keyed catalogue from a current-CALDB nicerl2 reprocessing, with an audit of the screening.
2. A restricted 2–11 keV band. No published catalogue uses it.
3. Full BPK02 parameters (ν0, Δ, r, Q, ν_max) with posterior or profile errors for every component, including the broadband noise. Stiele & Kong give only ν_max and rms for their noise components.
4. A calibrated, pre-registered rule for selecting components, published with its decision table.
5. A per-ObsID resolution of which QPO is the fundamental and which the harmonic. Published tables disagree by a factor of about 2, and the RPM spin depends on the choice.
6. A quantitative row-by-row validation against the three existing tables.
7. The hierarchical model across ObsIDs.
8. A documented FAIR release.
9. Two first results: the first HFQPO upper-limit table for this source, and the first placement of MAXI J1820 on the Körding 2007 variability plane. Neither has been published.

### 0.2 Other changes

- **The 55 Hz artefact now has a primary record.** ATel #11951 (Zhang L. et al., 14 Aug 2018) claimed a 55.12 ± 0.06 Hz QPO in ObsIDs 1200120107–120. ATel #11960, "The 55 Hz signal we detected in MAXI J1820+070 is not a QPO", retracted it as caused by small gaps in the GTI file. There is also a quantitative mechanism: 100-event packets at about 55 packets/s per MPU (HEASoft 6.32 thread).
- **Pile-up has a CALDB threshold.** The xti20240206 release notes say "Above ~2 Crab pile-up becomes significant". The brightest rise ObsIDs (around 107–112) are roughly 2–5 Crab in NICER counts. Webbe & Young 2024 validated dead-time behaviour only up to about 7,000 c/s.
- **New source papers.** Two new NICER timing papers in the 2024–2026 window: Kimura et al. 2024/25 (PASJ), a shot analysis with six hard-state phases, and Uttley & Malzac 2025 (MNRAS 536, 3284), whose Fig. 1 uses our ObsID 104. Also new: NuSTAR flaring and soft lags (Buisson 2021), HXMT evidence that the noise frequencies depend on energy (Gao, Yu & Yan 2024/25), and a missing pre-window NICER lag paper (Omama 2023).
- **Software.** No public package does incremental Lorentzian search plus QPO tracking across ObsIDs. Stingray's merged `gpmodeling` module (tinygp + jaxns) replaces the planned celerite2 GP cross-check. ISIS+SITAR is the preferred independent third code.
- **Comparisons.** Altamirano 2008 (VizieR J/ApJ/685/436) gives a ready neutron-star Lorentzian comparison set. The Körding 2007 recipe allows a first placement on the variability plane. Monitoring pipelines (MAXI Nova-Alert, Swift/BAT, WATCHDOG, OHMAN) do no automated Lorentzian or QPO tracking.

### 0.3 Corrections to rounds 1 and 2

| # | Where | What was said | Correction |
|---|---|---|---|
| C1 | Round 1, l. 283 (§5.2) | "Stiele & Kong (2020) had only about six QPOs above 3σ (as stated by Li et al. 2024)" | A misreading. Li et al. 2024 (Sect. 4.3) say the six significant QPOs fall in the **decaying IMS/HS**. In the rise, S&K Table 4 has about 155 QPO rows by day 117, roughly 76 of them at σ ≥ 3 (rough count from the PDF text). |
| C2 | Round 2, §10 verdict (l. ~421) | No full Lorentzian catalogue of the rise published as an open release | Superseded: S&K 2020 and VizieR J/ApJ/889/142 exist (Section 1.1). |
| C3 | Round 2, l. 97 and l. 407 | Fogantini 2025 "says the same features appear in MAXI J1820" / "report them in MAXI J1820" | The abstract analyses only NICER data of Cyg X-1. It says only that "similar features in MAXI J1348-630 and MAXI J1820+070 support this interpretation", probably referring to Méndez 2024 and Bellavita 2025. Reword. |
| C4 | Round 2, l. 414 and l. 598 | Mudambi et al. 2020 an "unconfirmed recollection; do not cite" | Verified in full: ApJL 889, L17, arXiv:2001.00755 (Section 3.1). Round 1, l. 331 had it right. |
| C5 | Round 1, l. 62 | Mudambi QPO date "not in the abstract" | Observed 2018-03-30/31 (about MJD 58207–58208), about 93 ks, LAXPC20. |
| C6 | Round 1, l. 440 | Prigozhin 2016 DOI pairing unconfirmed | Confirmed: SPIE 9905, 99051I, doi:10.1117/12.2231718. The full text was still not retrieved. |
| C7 | Round 2, l. 596 | Link between the 55 Hz feature and MPU packet limits is "an inference" | Now quantitatively consistent: 55 packets/s × 100 events ≈ 5.5×10³ ev/s per MPU, ≈ 3.8×10⁴ c/s for 7 MPUs, matching Enoto's 2021 throughput figure. The ATel #11960 authors attribute it to GTI gaps. No HEASARC page states the packet link outright, so keep the word "consistent with". |
| C8 | Round 1, l. 214 | "Dead-time-correcting the light curves before the FFT brings no significant gain (Webbe & Young 2024)" | Webbe & Young tested only 31–7024 c/s (about 1–135 c/s per FPM) on GX 339-4. They excluded telemetric dead time (≳ 20,000 c/s) and warn that their fit "will breakdown at very high count rates". The statement is unvalidated at MAXI J1820's 400–1000 c/s per FPM. |
| C9 | Round 1, §1.3 table and l. 97 | Bhargava 2021 implicitly treated as tracking the same QPO as the other papers | Bhargava model **only the higher-frequency member** of the harmonic pair. Their ν_LT = 0.276 Hz at ObsID 130 compares with 0.12 Hz (Bollemeijer 2024) and 0.120 Hz (Ma X. 2023) for the "fundamental" at the same ObsID. See Section 2.3. |
| C10 | Round 1, l. 66 / §1.2 | Early 66 Hz hint not recorded | The NICER Science Nugget of 2018-03-22 says NICER "marginally detected a 66 Hz Quasi-Periodic Oscillation (QPO) in the early data", "a tentative detection". It has no refereed follow-up and no retraction found. Treat it as instrumental until shown otherwise (Section 5.3). |

---

## 1. Novelty check

### 1.1 The prior catalogue: Stiele & Kong 2020 and VizieR J/ApJ/889/142

**Verification.** The upstream pass read the PDF appendix tables. The synthesis pass confirmed the six-table list on the CfA VizieR mirror.

**Tables:**
- Table 1: Swift/XRT broadband noise, 88 rows.
- Table 2: Swift/XRT QPOs, 101 rows.
- Table 3: NICER noise components, 196 rows. Columns are ν_max and rms for BLN1/BLN2, and ν_max, rms, Q and σ for the peaked noise PN1.
- Table 4: NICER QPOs, 212 rows. Columns are ν0 with asymmetric errors, Δ (HWHM) with errors or upper limits, rms, Q, σ, and a component number Num = 1 or 2.
- Table 5: an extra Lorentzian needed for the higher-frequency QPO, 39 rows, days 23–116.
- Table 6: extra QPOs or peaked noise, 6 rows.

**Coverage.** The NICER rows run from day 1.58 to day 116.9 (T0 = MJD 58188, 2018-03-11), about 90 distinct epochs. Day 116 is split into sub-intervals, for example:
- 116.27 and 116.28: 4.4–4.6 Hz.
- 116.91: 7.3–8.1 Hz.

Day 116.9 is MJD ≈ 58304.9. Li et al. 2024 put the last type-C at 8.05 Hz on MJD 58305.61 (day ≈ 117.6), so **check during the cross-match whether S&K cover the final hours of ObsID 197**. That is the type-C→B segment round 1 flagged.

**Method:**
- Pre-processed archive event files from 2018–19 (HEAsoft 6.23 era). Observations under 1 ks excluded.
- Full NICER band (0–15 keV), 6.1×10⁻³ to 50 Hz.
- Poisson level fitted as a constant to the high-frequency tail. No dead-time model.
- ISIS/GHATS following Belloni et al. 2006: at most 2 zero-centred noise components, at most 2 QPOs, and peaked noise.
- Significance is the normalisation divided by its error.

**Gaps that leave room for us:**
- No ObsID or MJD column; rows are indexed only by "days since 2018-03-11".
- The noise components are given as ν_max and rms only, not ν0 and Δ.
- No calibrated model-selection rule and no decision record.
- Some errors are pegged or non-physical, for example "+6.0685" in Table 5.
- Different band and different reduction.

**Community status.** Bollemeijer & Uttley 2025 (Table 1 caption) send readers to S&K 2020 "for a more detailed account of power-spectral fits". It is the de facto reference catalogue.

### 1.2 Other per-ObsID tables (partial precedents)

| Source | ObsIDs | Content | Errors | Band / processing |
|---|---|---|---|---|
| Bollemeijer et al. 2024, Table D1 | 126–189 (63) | Mean QPO fundamental 0.077 → 0.51 (ObsID 155) → 0.22 Hz (182/185); exposure; 0.5–1 and 3–10 keV rates | None | HEAsoft 6.28 |
| Bollemeijer 2024 Zenodo (10.5281/zenodo.10391399), `datafiles.zip/lightcurves_J1820_6.28/qpofitsc.txt` | 126–189 (63 rows × 16 numeric columns) | ObsID suffix, ν_QPO, more Lorentzian parameters | None; no header or documentation | 64 s × 1/128 s light curves |
| Bollemeijer & Uttley 2025, Table 1 | 16 rise ObsIDs: 104, 108, 112, 120, 130, 131, 135–137, 141–143, 145, 147, 148, 156, 187 | ν_QPO 0.034–0.50 Hz; hue −12 to 138°; HR(4–10/2–4); 52-FPM rate | None | Stated Zenodo DOI 10.5281/zenodo.14284350 returned 404 |
| Zhang Y. et al. 2026, Table A1 | 27 ObsIDs: 138, 140–147, 149–153, 156, 157, 161–163, 165, 171, 172, 193–197 | ν_QPO, hue, vKompth corona size, feedback fraction, χ²_red | QPO only | 0.3–12 keV; QPO required at > 3σ in every narrow band, which excludes ObsIDs before 138 |
| Bhargava et al. 2021, Table 1 | 21 ObsIDs, 130–191 | ν_LT, ν_L, ν_U, e.g. 130: 0.276 ± 0.021, 2.79, 29.2 Hz | 1σ | 0.01–12 keV, PDS 0.03–1250 Hz; higher QPO member only |
| Ma X. et al. 2023, Table 1 | 6 NICER epochs: 107 (0.038), 130 (0.120; printed "1200210130"), 189 (0.419), 194 (0.775), 196a (2.073), 196b (3.010 Hz); 2 HXMT matches | ν_QPO only; other parameters in figures | — | Energy-resolved |
| Li et al. 2024, Table C1 | 12 segments of 196/197 | QPO ν, rms, Q, σ | Yes (MCMC) | nicerl2 + SCORPEON; full HXMT track in figures only |

**No outburst-long Lorentzian table** in:
- Basak 2025 (ObsIDs 120, 123, 133, 134; two hue bins).
- Mudambi 2020 (one AstroSat epoch).
- Prabhakar 2022 (spectral tables only).
- Yang Z.-X. 2022 (HXMT, selected observations).
- Gao C. 2024/25 (HXMT, near the peak).
- Dziełak 2021 (ObsIDs 103/104).
- Kawamura 2023 (ObsIDs 106/108).
- Fiori 2025: reprocessed ObsIDs 101–278 with nicerl2, HEAsoft 6.29, CALDB 20210707, 1–12 keV, but the NICER QPO track is shown only in Fig. 13.
- Jin 2026 (6 intervals).
- Ma R. 2023 (ObsID 197).
- Méndez 2024 (ObsID 120).

A Zenodo API search for "MAXI J1820+070" (25 hits) found no other timing-fit release. **A keyword search of VizieR was not possible** because the CDS main site reset the connection, so per-exposure HXMT tables on VizieR cannot be excluded.

### 1.3 What remains novel

1. **The reduction is ObsID-keyed, current, and audited.**
   - All 88 ObsIDs 101–197 reprocessed with nicerl2 under HEASoft ≥ 6.35.2 and CALDB xti20240206.
   - A per-ObsID screening audit: exposure, NUM_FPM and surviving MPUs under default, relaxed and NONE settings.
   - Robotic-arm masking for ObsIDs 130/131.
   - A Poisson level that accounts for dead time and pile-up.
   - Nearest precedent: Fiori 2025, which used a 2021 CALDB and published no table.
2. **The 2–11 keV band.** No catalogue uses it: S&K 0–15, Zhang 0.3–12, Fiori 1–12, Bhargava 0.01–12 keV.
3. **Full BPK02 parameters for every component**, with posterior or profile errors.
4. **A calibrated, pre-registered stopping rule** (SDDR or calibrated likelihood-ratio test in place of the F-test), with a public decision table.
5. **Resolution of fundamental against harmonic** across the rise, and how much the RPM spin depends on it (Section 2.3).
6. **A validation table** comparing our values with S&K (VizieR), Bollemeijer's `qpofitsc.txt` and Zhang's Table A1, including the systematic offsets from band, reduction and Poisson treatment.
7. **The hierarchical model across ObsIDs.** Round 2 found no precedent; none was found this round either.
8. **A FAIR release:** ECSV on Zenodo with ObsID, MJD start/mid, exposure, band and component labels, documented columns and units, asymmetric errors and a test statistic for every component. No existing release has this.
9. **Two products with no published counterpart:**
   - An HFQPO upper-limit table (Section 3.3).
   - Placement on the Körding 2007 variability plane, including the θ(f) transition in the HIMS (Section 6.2).

**Framing for the paper.** Present items 1–8 as a "re-analysis and extension of the S&K 2020 catalogue". Results that depend only on ν_QPO(t) are **confirmations**.

---

## 2. Recent MAXI J1820 NICER papers (2024–2026; new to this review)

**How the search was done.** arXiv listings were read on 2026-11-05: title "MAXI J1820", 23 hits for 2024–26; all fields "MAXI J1820+070" NICER, 12 hits. About 15 searches were run. Sixteen arXiv IDs already in rounds 1 and 2 were skipped.

### 2.1 New and in scope

**Kimura, Negoro, Yamada, Iwakiri, Sako & Ohsawa 2024/25** (PASJ, doi:10.1093/pasj/psae099, arXiv:2411.03602).

Data and method:
- NICER light curves at 0.015 s, 0.3–7 keV, split into 0.3–2 and 2–7 keV.
- Simultaneous Tomo-e Gozen optical at about 67 fps.
- MJD ≈ 58189.5–58303.5.
- Shot superposition: each day is detrended with LOWESS, shots are cut over ±5 s, and kept if R_near > 0.3. Each averaged shot is fitted with exponential rise and decay.

Results:
- Shot timescale about 0.2 s.
- Mean duration 0.37 s in the soft band and 0.31 s in the hard band, so hard shots are about 20% narrower.
- Mean amplitude 0.68 in X-rays against 0.33 in the optical.
- X-ray shot amplitude peaks at onset and fades towards the intermediate state.

Six hard-state phases (MJD converted from BJD):

| Phase | MJD |
|---|---|
| H1 | 58189.5–58194.5 |
| H2 | 58197.5–58215.5 |
| H3 | 58217.5–58236.0 |
| H4 | 58237.5–58252.5 |
| H5 | 58258.5–58289.5 |
| H6 | 58289.5–58303.5 |

The paper does no PDS fitting. It notes "wavy" QPO modulation in the 10 s shot profiles.

**Uttley & Malzac 2025** (MNRAS 536, 3284–3307, doi:10.1093/mnras/stae2514, arXiv:2312.08302).

- Model: fluctuations propagate through the disc into a compact corona with r_cor = 10 R_g and Γ ≈ 1.47.
- The lag between disc and power law changes from hard to soft above about 1 Hz.
- The hard lag between two power-law bands scales as τ ∝ ν^−0.7.
- Fig. 1 shows lags and PDS from **NICER ObsID 1200120104** (about 6 ks, 2018-03-15; processed as in Wang 2022) in 0.5–0.9, 2–3 and 5–10 keV. The model is not fitted to these data.
- The second half of the full text was not read.

### 2.2 Partial matches

- **Davidson et al. 2025** (arXiv:2505.07962; submitted to ApJL, v2 Nov 2025). NICER, NuSTAR and Swift spectra of MAXI J1348-630, MAXI J1535-571 and MAXI J1820+070. The corona moves outward near jet quenching and at transient-jet launch. Spectral only; context for the late-HIMS ObsIDs.
- **Yang P. et al. 2026** (A&A accepted, arXiv:2606.17976). An imaginary QPO in NICER data of AT2019wey falls from about 5 to 1 Hz. Across sources, the QPO lag rises with frequency. A search snippet suggests MAXI J1820 (i = 66–81°) is a comparison source in the full text; **not confirmed**.
- **Omama et al. 2023** (ApJ 945, 92, arXiv:2302.11263). NICER lags from differential cross-correlation. It predates the window but is missing from both earlier reviews.

### 2.3 The fundamental-versus-harmonic problem across published tables

At ObsID 130 (2018-04-16, S&K day ≈ 36.07):

| Source | Reported value(s) | Labelled as |
|---|---|---|
| Bollemeijer 2024 | 0.12 Hz | fundamental |
| Ma X. 2023 | 0.120 Hz | fundamental |
| Bhargava 2021 | 0.276 ± 0.021 Hz | ν_LT |
| S&K 2020 | 0.130 and 0.229 Hz QPOs, plus a broad 0.31 Hz component | two QPOs |

At ObsID 152, Bhargava give 1.023 Hz against Bollemeijer's 0.46 Hz.

Bhargava state that they model only the higher member of the pair. Note that the ratios are not exactly 2 (0.276/0.12 ≈ 2.3; 0.229/0.130 ≈ 1.76), so the "pair" may be partly blended with a broad component.

**Consequence.** The RPM spin a = 0.799 (+0.016/−0.015) assumes the higher peak is the Lense–Thirring frequency. Fiori 2025 obtain a ≲ 0.15 using the optical subharmonic. Round 1 recorded this disagreement. **This round finds that the published tables themselves disagree on the identification, and nobody has resolved it.**

---

## 3. Other-mission MAXI J1820 timing

### 3.1 AstroSat

**Mudambi et al. 2020** (ApJL 889, L17), now verified in full.
- Data: about 93 ks over 15 orbits from 2018-03-30. LAXPC20 only, because LAXPC10 had a gain change and LAXPC30 was off. PDS over 0.004–30 Hz.
- Five Lorentzians (χ²/dof = 123/114):
  - QPO at 47.7 (+1.6/−2.0) mHz, width 11.9 (+5.0/−3.6) mHz, so Q ≈ 4 (derived).
  - Weak QPO at 109.4 (+2.8/−1.2) mHz, width < 1.7 mHz. The ratio to the main QPO is 2.29, so it is not a clean harmonic.
  - Zero-centred humps of width 0.24 ± 0.02 Hz and 10.7 Hz.
  - A hump at 1.04 Hz, width 1.6 Hz.
- Hard lags (reference band 4.15–5.37 keV): about 100 ms below about 2 Hz and about 10 ms above.
- Only Lorentzian norms are given, not fractional rms.
- The spectrum has Γ = 1.61 and kT_in = 0.22 keV.
- The upstream pass maps the date to about NICER ObsIDs 113–115. **Check this against the ObsID MJDs.**

**Prabhakar, Mandal, Athulya & Nandi 2022** (MNRAS, arXiv:2204.13363), already cited; new detail.
- LAXPC ObsIDs 9000001994 (MJD 58208), 9000002216 (MJD 58309) and 9000002324 (MJD 58355).
- A QPO appears only at MJD 58208: 0.05 Hz, Q = 4.1 ± 0.7, 9.1σ, rms 9.1 ± 0.3%. This is an independent confirmation of Mudambi's 47.7 mHz one day later.
- Their NICER PDS (0.01–55.55 Hz) include, for example, 0.42 Hz with Q = 3.19 at only 2.3σ on day 56.43.
- Total rms is above 40% early, about 20% on the plateau and about 30% in Phase III.
- No Swift/XRT timing and no HFQPO search.

**Spectral only, new:**
- Banerjee S. et al. 2024 (ApJ 964, 189; arXiv:2402.08237). Two hard-state coronae with R_in ≳ 51 R_g; soft-state a = 0.77 ± 0.21 at i = 64 ± 5°.
- Chakraborty et al. 2020 (arXiv:2009.02465). AstroSat plus NuSTAR, hard state.

**Not found:** any AstroSat timing paper on MAXI J1820 from the Chand, Jana, Mondal, Bhuvana, Dutta, Sreehari or Debnath/Chakrabarti (TCAF) groups. Treat any recalled paper of this kind as nonexistent until it is found on ADS.

### 3.2 NuSTAR, XMM-Newton, Swift/XRT, HXMT

**Buisson et al. 2021** (MNRAS 500, 3976; arXiv:2011.02472), new.
- NuSTAR observation of the hard-to-soft transition.
- About 10 ks of variability between two flux plateaus: a dip phase, then a flare phase, with F_var ≈ 10%, strongest at high energy.
- The disc lags the power law by 20 (+1.6/−1.2) s, which they read as the disc refilling with α ≈ 10⁻³.
- Implication for us: ks-scale red noise and dips near the HIMS end (ObsIDs about 189–197) can bias the lowest Lorentzian and the total rms.

**Stiele & Kong 2020 (XRT part).**
- Two noise components up to day 22.4.
- A 6.19 Hz QPO (Q = 18.4, 2.97σ) on day 374, outside our range.

**Gao, Yu & Yan 2024/25** (ApJ, to appear; arXiv:2409.06414), HXMT, new.
- Noise components below 0.1 Hz shift to lower frequency as photon energy rises; those above 10 Hz shift higher.
- In their words, "the higher the photon energy is, the wider the plateau and the smaller the fractional variability."
- Table 2 gives ν_max and rms against energy in 9 bands.
- **Our 2–11 keV characteristic frequencies therefore depend on the band**, and comparisons across missions (LAXPC 3–80 keV, HXMT, Bhargava 0.01–12 keV) must be treated as band-dependent.

**Other HXMT papers, no NICER data, listed for completeness:**
- Fan et al. 2024 (ApJ 969, 61; arXiv:2404.12161).
- Peng et al. 2025 (arXiv:2502.14508).

### 3.3 High-frequency QPOs: no genuine detection, no published upper limit

- **ATel #11951** (Zhang L. et al., 14 Aug 2018): a 55.12 ± 0.06 Hz feature with FWHM 2.00 ± 0.20 Hz and about 1% rms, clear in ObsIDs 1200120107–120.
- **ATel #11960** retracted it: "The 55 Hz signal we detected in MAXI J1820+070 is not a QPO". The cause was small gaps in the GTI file, as the NICER team pointed out.
  - The ATel pages are blocked to our fetcher. The content comes from search-result text and the ATel index title.
- **NICER Science Nugget, 2018-03-22:** a tentative, marginal 66 Hz QPO in early data. There is no follow-up. It plausibly comes from the same GTI or telemetry mechanism, but that is an inference.
- A snippet links this hint to ATel #11423 ("…rapidly-brightening black hole X-ray binary in the hard state"). The ATel page did not load, so the 66 Hz content of that ATel is unverified.
- **Bhargava 2021** fit up to 1250 Hz in 21 ObsIDs:
  - ν_LT 0.276–1.143 Hz.
  - ν_L 2.79–7.22 Hz.
  - ν_U 22.0–46.3 Hz, highest at ObsID 190.

  They report no HFQPO search, no limit and no mention of 55 Hz. RPM predictions appear only as a dashed line in their Fig. 3.
- **Result:** a per-ObsID HFQPO upper-limit table for 101–197 would be new in the literature.

---

## 4. Software and pipelines

**Bottom line.** No public, maintained package does what modelingPS does: incremental multi-Lorentzian fits to averaged PDS, a stopping rule, QPO identification, and tracking across ObsIDs. Every outburst-sample study found used private scripts. Keep the modelingPS driver and release it; that release fills a gap in itself.

The current code was checked:
- `src/modelingPS/fitting.py` uses scipy COBYLA and `scipy.stats.f`, with `ftestlim=0.05`.
- `src/modelingPS/uncertainties.py` uses a finite-difference Hessian, cov = 2H⁻¹.
- `pixi.toml` targets `osx-arm64` only.

| Tool | Language / licence | Use for modelingPS | Verdict |
|---|---|---|---|
| Stingray `modeling.gpmodeling` (PR #739, G. Joshi, merged 2023-09-29) | Python, MIT; tinygp + jaxns | Nested-sampling GP QPO test against red noise on **light curves** (a JAX version of the Hübner+2022 method) | **Adopt** for the GP cross-check (round-2 item 14) on a subset of GTIs. Needs Phase-1 event files. |
| tinygp (JOSS doi:10.21105/joss.04455) | JAX, MIT | Backend for the above; no inference routines | Indirect |
| jaxspec (Dupourqué et al. 2024, A&A 690, A317) | JAX/NumPyro, MIT | Energy-spectrum fitter with Poisson/Gaussian likelihoods, not Whittle | Template only for a ~50-line JAX Whittle/Gamma(MW) Lorentzian likelihood |
| SITAR v1.2.8 (Nowak; ISIS/S-Lang; licence not stated) plus Nowak et al. 2008 (arXiv:0806.0560) | S-Lang | Averaged PSD, rebinning, cross-spectra; ISIS accepts **user-defined fit statistics** | **Preferred third code** for the three-code validation |
| XSPEC + flx2xsp (Ingram & Done 2012 route) | — | Common practice (GX 339-4 arXiv:2111.13403 uses an incremental `ftest`; RX J0440.9+4431 arXiv:2402.14226; GRS 1915+105 arXiv:2403.09341; Swift J1727 arXiv:2410.07333); no pipelines released | Fallback cross-check only: χ² only, heavy rebinning needed |
| BXA (Buchner) | Python, licence not seen | UltraNest evidence through XSPEC; non-Whittle statistics | Only if the team standardises on XSPEC |
| DIAMONDS (Corsaro & De Ridder 2014) / FAMED (Corsaro et al. 2020) | C++ / mixed, licences not confirmed | Automated evidence-based Lorentzian peak-bagging across tens of thousands of stars; assumes χ²₂ statistics and comb-like spacing | Design template for replacing the F-test loop |
| Gruberbauer et al. 2009 (A&A 506, 1043) | Methods only | "Ignorance" priors for Bayesian Lorentzian detection | Basis for SDDR priors (compare Malacaria's prior-volume sensitivity, log B from −30.3 to −23.4) |
| QPOEstimation (Hübner) | Python, no licence | Original Hübner+2022 code; about 2 years with no push | Do not depend on it |
| mind_the_gaps (Gúrpide & Middleton 2025) | Python, MIT | GP likelihood-ratio test for unevenly sampled data | Little benefit for nearly evenly sampled GTIs; context only |
| Pioran.jl (Lefkir et al. 2025) | Julia, MIT | Bending power laws only; no QPO component | Ruled out |
| EzTao (ascl:2201.001) | Python | CARMA; superseded by Stingray gpmodeling here | Low priority |

**Automated trackers.**
- No public outburst QPO tracker exists.
- Published recipes:
  - Belloni+2012: a grid stepped in 1 Hz, with a norm/error > 3 cut.
  - Marcel et al. 2020 (arXiv:2005.10359): a zero-centred Lorentzian plus extra Lorentzians at the largest residuals, keeping Q > 2 and > 5σ.
  - de Ruiter+2019: an RXTE pipeline with no released code.
- No HMM, Kalman or state-space QPO tracker was found for X-ray binaries.

**Not found; do not cite:** "BayesianPSD", "pyPSD", public code for Malacaria+2024. The `logpsplinepsd` package page did not load.

**Hardware.** JAX GPU support on `osx-arm64` (Metal) is experimental. 88 spectra × ~226 bins is cheap on CPU, so plan for CPU-only JAX.

---

## 5. NICER bright-source systematics and HEASoft changes

This section adds to round 2, Section 7, without repeating it.

### 5.1 Pile-up

- The CALDB release notes for xti20200722 (Markwardt & Steiner) and for **xti20240206, the CALDB this project uses** (Markwardt & Hare), both say: "Above ~2 Crab pile-up becomes significant and will distort the spectrum in a complicated and non-linear fashion."
- Both notes also give an approximately 8% dead-time correction in the Crab calibration spectra.
- Rate conversion:
  - Crab ≈ 10,531 ± 60 c/s (0.3–12 keV, 50 FPM, not dead-time corrected). Remillard et al. 2022; **search snippet only**.
  - MAXI J1820 ObsID 107 is at about 20,000 to more than 50,000 c/s (HEASARC Tips).
  - The brightest rise ObsIDs are therefore about 2–5 Crab in NICER counts. This ratio depends on spectral shape.
- No NICER timing paper that models pile-up was found.
- Enoto's 2021 workshop slide says NICER is "free from pileups, dead time, and data transfer loss". It conflicts with the CALDB notes and refers to Crab-pulsar rates.
- Inference (not stated by any source):
  - Pile-up lowers the Poisson level in the same direction as dead time.
  - It moves summed soft-photon events into 2–11 keV by an amount that depends on flux.
  - This would bias the measured rate and fractional rms around ObsIDs 107–112.

### 5.2 Dead time

- The NICERDAS Calibration Guide (2018-08-14) says:
  - Dead time is recorded per event in the ufa files and as MPU_DEADTIME in the .mkf file.
  - "Each MPU is parallel".
  - All events cause dead time: resets, overshoots, and forced triggers at 5 /s/FPM.
- HEASARC Tips: dead time matters from about 20,000 c/s with 52 FPMs.
- Webbe & Young 2024 (RASTI): about 22–23 µs per event, rising 2–3% with rate. Their limits are listed in C8.
- No per-MPU model of dead time or telemetry exists in the literature.

### 5.3 Telemetry saturation and the 55 Hz features

**Mechanism (2018 data).**
- 100-event packets at about 55 packets/s per MPU, so "tens of gaps per second" (HEASoft 6.32 thread).
- In ObsID 107 this gives 5912 GTIs, with GTI spacings peaking at 1/100 to 1/10 s.
- HEASARC's fix is to use only continuous segments of 16 s or more.
- Primary record: ATel 11951/11960 (Section 3.3).

**A second, unrelated mechanism (2023 onward).**
- A 55 Hz feature plus a broad 18–30 Hz feature in low-Sun-angle data, on MPUs 1 and 3 only, independent of energy (Tips page).
- It does not affect 2018 data but should be cited when the artefact is discussed.
- O'Connor 2023 (GCN 34549) reports that the 55 Hz feature disappears in orbit-night data.

**Bollemeijer, Uttley & You 2026, Appendix C** (arXiv:2609.27608), for ObsIDs 107–112:
- They fill inter-GTI gaps under 1 s in each MPU light curve with Poisson draws. This removes spurious variability around 20 Hz and leaves low frequencies visibly unchanged.
- They keep a fixed FPM set of 14–28 of 52 per ObsID, normalised to 52. Without this, the coherence drops spuriously.
- They fit the Poisson level above 200 Hz (256 s segments × 1/512 s).

### 5.4 Autoscreening risk for shredded 2018 data (current niautoscreen help)

- `mpugtiscr = "0.0:0.003:0:0.5"`: an MPU is cut if it has more than 0.5 GTI transitions per second.
- `roundrobbinscr = "0:0.0015:0:0.5"`.
- `lowmemscr = "0:0:0:200"`, while the help text says 500. This is an internal inconsistency; record which value actually applies.
- Both the nicerl2 and niautoscreen help pages say the team recommends `lowmemfrscr` over `lowmemscr`.
- With 5912 GTIs in ObsID 107, whole MPUs (or all of them) could be removed.

### 5.5 HEASoft 6.33–6.37.1

| Version (date) | NICER-relevant change |
|---|---|
| 6.33.0 (2024-02-23), NICERDAS v012 | `lowmemscr` per MPU (default at most 200 lost events per pointing GTI); `max_lowmem` disabled; `threshfilter` |
| 6.35.2 (2025-06-25), NICERDAS v14 | niautoscreen v1.6 fixes lowmem screening that "would only screen out MPU0". Also adds `lowmemfrscr`, `autoscreen_tseg` (day default 100 s), `timerange` (drops TIME < 100 s after an MPU reboot), NICERV6 columns (restoring MPU_*_COUNT and adding SUN_ELV), nifpmsel v1.10, and the new task `nicerl2-screen` |
| 6.36.0 (2025-09-26) | Minor: a nicerl2 stray-file fix, a nitargdet CALDB fix, barycorr support for FPM_SEL and filter files, a caldbinfo fix |
| 6.37.0 (2026-08-13) / 6.37.1 (2026-09-02) | No NICERDAS changes. **xselect v3.0 (C++) calls the new `ftextractor` in place of `extractor` for NICER.** Round 2 found that extractor ≤ 5.37 drops TIMEZERO, so re-check TIMEZERO handling (untested) |

Current nicerl2 defaults: `mingti=5.0`, `erodedilate=5.0` (applied only to telemetry-selected criteria, not to user GTIs), `underonly_range 0-500`, `overonly_range 0-30`, `min_fpm=7`, `filtcolumns=NICERV6`, `tasks=ALL`.

### 5.6 Absolute timing (matters only for cross-mission or optical lags)

- Better than 100 ns after L2 calibration. Instrumental delays of about 440 ns (fast chain) and 950 ns (slow-only) are corrected by nicertimecal. A +1 s offset is corrected through TIMEZERO.
- Archived data carry mixed calibration levels; "user must apply most recent nicerl2 and CALDB" (Calibration Guide).
- Timing chain 86 ns rms (Hamaguchi, IACHEC 2019).
- Tuo et al. 2022 (ApJS) use NICER Crab data as the reference, because the Jodrell Bank ephemeris has about 40 µs systematics.
  - HXMT offsets relative to NICER: HE 24.7 µs, ME 10.1 µs, LE 864.7 µs.
  - The NICER systematic of 5.5 µs is from a snippet only.
- Enoto et al. 2021 (Science 372, 187): the X-ray main pulse leads the radio by phase 0.9915 ± 0.00004 (about 304 µs), which is source-intrinsic.
- No dedicated NICER Crab absolute-offset paper was found. Bogdanov et al. 2019 (ApJL 887, L25) report "no timing or Fourier anomalies" in the MSP data, with no number given.
- None of this matters for 0.01–10 Hz Lorentzian fits.

---

## 6. NS, AGN and ULX comparisons, and monitoring

### 6.1 Neutron-star Lorentzian catalogues

**Altamirano et al. 2008** (ApJ 685, 436; arXiv:0806.0962; VizieR J/ApJ/685/436).
- 149 RXTE observations of 4U 1636-53 across the island and banana states.
- Follows the universal atoll correlations.
- The hectohertz QPO sets its frequency by a different mechanism from the other components.
- Harmonic mode switching does not explain the narrow low-frequency QPOs.
- **Its public Lorentzian tables make it a ready comparison set for ν_LF–ν_b and ν_LF–ν_l.**

**van Straaten, van der Klis & Wijnands 2005** (ApJ 619, 455).
- Accreting millisecond pulsars follow the atoll correlations.
- In SAX J1808.4-3658 and XTE J0929-314 the correlations are shifted by a factor of about 1.5 in the kHz QPO frequencies.
- Lesson: **fit a multiplicative shift factor** between tracks rather than claim they are "consistent".

**van Straaten, van der Klis & Méndez 2003** (ApJ 596, 1155; exact title not confirmed).
- 4U 1608-52 fitted with Lorentzians plus a power law for the very-low-frequency noise. Timing is nearly identical to 4U 0614+09 and 4U 1728-34.
- A companion paper (astro-ph/0307041) ties the components to position on the colour-colour track.

**kHz QPO coherence.**
- Boutelier, Barret & Miller 2009 (arXiv:0907.3223) and Barret et al. 2011 (arXiv:1011.4601, XTE J1701-462): the Q and rms of the lower kHz QPO drop above about 700–850 Hz, a proposed ISCO signature.
- In 4U 1636-53, Q rises to about 200 near 850 Hz and then falls.
- Sanna et al. 2010 dispute the drop; seen only through Barret 2011's description.

### 6.2 Mass scaling and the variability plane

**McHardy et al. 2006** (Nature 444, 730), coefficients verified from the full text this round.
- log T_B = 2.10(±0.15) log M₆ − 0.98(±0.15) log L₄₄ − 2.32(±0.2), with T_B in days.
- Fitted to soft or radio-quiet black-hole states only.
- Within Cyg X-1, T_B ∝ L^−1.3±0.2.
- In hard states, the Lorentzian nearest the break plays the role of T_B, but hard states were not fitted.

**Körding et al. 2007** (MNRAS 380, 301; arXiv:0706.2959).
- log ν = log Ṁ[g/s] − 2 log M[M☉] − 14.7 − 0.9 θ(f), where θ = 0 in the soft state (ν = ν_hb) and θ = 1 in the hard state (ν = ν_l).
- So hard-state ν_l sits **0.9 dex below** the soft-state break at fixed Ṁ and M.
- Hard-state Ṁ estimators:
  - Radio: Ṁ = 3×10¹⁷ (L_R/10³⁰)^(12/17) g/s.
  - X-ray: Ṁ ≈ 3.4×10¹⁷ (L_2–10/10³⁶)^0.5 (M/M_GX339)^0.43 g/s.
- Calibrated on GX 339-4, XTE J1550-564 and Cyg X-1.

**Extensions of the plane:**
- Scaringi et al. 2015 (Sci. Adv. 1, e1500686) extend it to white dwarfs and YSOs.
- González-Martín & Vaughan 2012 (A&A 544, A80): a bending power law is preferred in 15 of 104 AGN, ⟨ν_b⟩ ≈ 2×10⁻⁴ Hz, and only RE J1034+396 shows a QPO.

**MAXI J1820 on the plane.** No published paper places it there. The upstream pass's order-of-magnitude estimate, **not from any paper**, uses M = 8.48 M☉ and d = 2.96 kpc:
- Hard-state ν_l of about 2.3, 3.3 and 7.4 Hz for F(2–10) = 5×10⁻⁹, 10⁻⁸ and 5×10⁻⁸ erg cm⁻² s⁻¹.
- This is plausible against the observed hump of a few Hz; Mudambi's 1.04 Hz hump of width 1.6 Hz has ν_max ≈ 1.9 Hz.

### 6.3 AGN and ULX QPOs

**RE J1034+396.**
- 3730 ± 130 s, about 5.6σ, Q > 16 (Gierliński et al. 2008).
- Read as an analogue of the 67 Hz QPO in GRS 1915+105 (Middleton, Uttley & Done 2011), and non-stationary.
- Detected in 5 XMM observations (Alston et al. 2014).
- **It is an HFQPO analogue and does not calibrate type-C scaling.**

**ULXs.**
- NGC 5408 X-1:
  - Strohmayer et al. 2007: QPO and break frequencies 10–100 times lower than in Galactic black holes.
  - Strohmayer & Mushotzky 2009: 1000–9000 M☉ (likely 2000–5000) from type-C-like scaling of ν_QPO with photon index.
- M82 X-1:
  - Pasham, Strohmayer & Mushotzky 2014: twin 3.32/5.07 Hz QPOs (3:2), giving 428 ± 105 M☉ from inverse-mass scaling and 415 ± 63 M☉ from the RPM.
  - Pasham et al. 2013: the mHz QPO frequency tracks count rate (r = +0.97) but not hardness, which challenges the type-C analogy.
- **A uniform 88-ObsID ν_QPO track against photon index or hue, quoted with M = 8.48 M☉, is exactly the kind of calibrator these estimates rely on.**

### 6.4 Monitoring and alert pipelines (round-2 open item)

- **MAXI/GSC Nova-Alert** (Negoro et al. 2016, PASJ 68, S1): triggers at about 80 mCrab in a single scan; 15 novae, 177 ATels and 63 GCNs to March 2015. It discovered MAXI J1820.
- **OHMAN** (on-orbit MAXI→NICER link, since August 2022): response about 10 min instead of about 3 h; M15 burst followed 5.5 min after detection. Sources are the BAAS 2023 abstract and press releases; no refereed paper.
- **Swift/BAT Hard X-ray Transient Monitor** (Krimm et al. 2013, ApJS 209, 14): more than 900 sources; 245 detected, 17 new (2005–2013).
- **WATCHDOG** (Tetarenko et al. 2016, ApJS 222, 15): 132 outbursts from 1996–2015; about 40% hard-only.
- **No pipeline does automated Lorentzian or QPO tracking.** A modelingPS that can run on new NICER ToO data would be a methods contribution.

---

## 7. Updated recommendations (changes to the plan)

**Novelty and framing**

1. **Rewrite the novelty statement now.**
   - Drop any wording like "first per-ObsID Lorentzian catalogue".
   - Cite S&K 2020 and VizieR J/ApJ/889/142 as the prior catalogue.
   - Frame the work as an ObsID-keyed, current-calibration, 2–11 keV, statistically calibrated re-analysis and extension (Section 1.3).
   - Apply corrections C1–C10.
2. **Add a cross-match and validation step as a deliverable.**
   - Map S&K "Time" (days since MJD 58188) to each ObsID's start and mid MJD.
   - Compare ν0, Δ, rms and Q row by row.
   - Do the same against Bollemeijer's `qpofitsc.txt` (126–189) and Zhang Table A1 (138–197).
   - Report median offsets and scatter, attributed to band (2–11 vs 0–15 keV), reduction and Poisson treatment.
   - Confirm whether S&K reach the last hours of ObsID 197.
3. **Make fundamental-versus-harmonic a headline result.**
   - For each ObsID, report both Q > 2 peaks and their ratio.
   - Refit the RPM with both labellings, to show how a = 0.80 depends on the choice.
   - Compare with Fiori's a ≲ 0.15.

**Phase 1 reduction**

4. **Pin HEASoft 6.37.1 (≥ 6.35.2) with CALDB xti20240206 and `tasks=ALL`.**
   - Verify that 6.36 and 6.37.1 give identical GTIs.
   - Check that TIMEZERO (about −1 s) survives `ftextractor`.
5. **Run a screening audit before any PDS is made.** For each ObsID:
   - Run default, relaxed and NONE settings of `mpugtiscr`, `roundrobbinscr`, `lowmemscr` and `lowmemfrscr`.
   - Tabulate exposure, NUM_FPM and surviving MPUs.
   - Flag ObsIDs where autoscreening removes whole MPUs, ObsID 107 first.
6. **Treat the Poisson level for bright data properly.**
   - Never fix it at 2/rate, and do not rely on an analytic dead-time correction alone.
   - Fit it per ObsID above about 200 Hz.
   - Cross-check with MPU-disjoint cospectra.
   - Keep the FPM set fixed within each ObsID or segment and normalise to 52.
   - For the brightest ObsIDs (around 107–112), report rms against counts per FPM to show that the trends are not pile-up or dead time.
7. **GTI and segment rules.**
   - Use integer segments within GTIs. Use HEASARC's continuous segments of ≥ 16 s for ObsIDs 107–120.
   - As a secondary check, do Bollemeijer-style Poisson gap filling for gaps under 1 s per MPU.
   - Mask 40–70 Hz and 18–30 Hz in high-frequency work.
8. **Energy bands.**
   - Keep 2–11 keV as the primary, fixed band.
   - Also keep 0.5–0.9 / 2–3 / 5–10 keV (or similar) products from the same reprocessing, because noise frequencies depend on energy (Gao 2024/25).
   - Use ObsID 104 to reproduce Uttley & Malzac's Fig. 1 as a check of the lag pipeline.
9. **Screen the HIMS-end ObsIDs (about 189–197)** for ks-scale flares and dips (Buisson 2021; the Kajava dips) using flux and hardness light curves inside each GTI, before fitting the lowest Lorentzian.

**Fitting and software**

10. **Replace `scipy.stats.f` in `fitting.py` with SDDR or nested-sampling evidence.** Use priors following Gruberbauer 2009, with Malacaria-style prior-volume checks; FAMED's logic is a template.
11. **Replace the finite-difference Hessian** with posterior sampling for weak or bounded components.
12. **Engine.** Stingray's PSD likelihood, or a small CPU-JAX Whittle/Gamma(MW) model.
13. **GP cross-check.** Use Stingray `gpmodeling`, not new celerite2 code or QPOEstimation.
14. **Third code for validation.** Use ISIS+SITAR with a user-coded Whittle statistic in preference to XSPEC+flx2xsp.

**Templates and new products**

15. **Starting templates for early ObsIDs.**
    - Use Mudambi's five-Lorentzian layout.
    - Validate the lowest type-C fits against 47.7 mHz (Mudambi, MJD ≈ 58207) and 0.05 Hz with Q = 4.1 (Prabhakar, MJD 58208). Expect the rms to differ because of the band.
    - Test the 109 mHz peak (ratio 2.29) as a separate component.
16. **New products:**
    - Per-ObsID 3σ rms limits for HFQPOs over 50–500 Hz, with a gap-injection null test using the original GTI files, set against RPM predictions for a = 0.8 and a = 0.15.
    - The Körding variability-plane placement, fitting d log ν_l / d log Ṁ (1 expected; Cyg X-1 gives 1.3 ± 0.2) and θ(f) through the HIMS, with the caveats in Section 6.2.
    - An overlay of NS correlations (Altamirano 2008 VizieR) with a fitted shift factor.
    - A comparison with Kimura's H1–H6 phase boundaries, and of their 0.2–0.4 s shot timescales with 1/(2πν_max).
17. **Release format.**
    - ECSV on Zenodo with ObsID, MJD start/mid, exposure, band, component label (Klein-Wolt & van der Klis 2008), all BPK02 parameters with asymmetric errors, the test statistic and decision for every component, documented columns and units, and the software versions.
    - Fill the gaps found in S2 and S4: no ObsID or MJD column, no header, no errors.

---

## 8. Annotated bibliography

### 8.1 Verified

**Prior catalogues and per-ObsID tables (novelty)**
- **Stiele, H. & Kong, A. K. H. (2020)**, "A timing study of MAXI J1820+070 based on Swift/XRT and NICER monitoring in 2018/19", ApJ 889, 142. arXiv:1912.07625; doi:10.3847/1538-4357/ab64ef. The prior per-epoch NICER Lorentzian catalogue: Tables 3–6 hold 196/212/39/6 rows over days 1.58–116.9, in 0–15 keV, fitted with ISIS/GHATS (Section 1.1).
- **VizieR J/ApJ/889/142**, "Swift/XRT and NICER timing study of MAXI J1820+070 (Stiele+, 2020)". Six tables (88/101/196/212/39/6 rows), confirmed on the CfA mirror. No ObsID or MJD column. The baseline for the row-by-row comparison.
- **Bollemeijer, N., Uttley, P., Basak, A. et al. (2024)**, MNRAS 528, 558. arXiv:2312.09835; doi:10.1093/mnras/stad3912. Table D1: QPO fundamental for 63 ObsIDs (126–189), 0.077 → 0.51 → 0.22 Hz, no errors.
- **Bollemeijer, N. (2024)**, Zenodo doi:10.5281/zenodo.10391399. `qpofitsc.txt`: 63 × 16 Lorentzian parameters, no header or errors; HEAsoft 6.28.
- **Bollemeijer, N. & Uttley, P. (2025)**, MNRAS, doi:10.1093/mnras/staf1238; arXiv:2507.18437. Table 1: ν_QPO, hue and HR for 16 rise ObsIDs. Defers to S&K for PDS fits. The Zenodo DOI 10.5281/zenodo.14284350 returned 404.
- **Zhang, Y., Méndez, M., Steiner, J. F., García, F. et al. (2026)**, "A NICER view of the corona…", MNRAS 546, stag261; arXiv:2602.07247. Table A1: 27 ObsIDs (138–197) in 0.3–12 keV.
- **Bhargava, Y., Belloni, T., Bhattacharya, D., Motta, S., Ponti, G. (2021)**, "A timing-based estimate of the spin of the black hole in MAXI J1820+070", MNRAS 508, 3104; arXiv:2109.14371; doi:10.1093/mnras/stab2848. Author list checked on arXiv this round. Table 1: 21 ObsIDs (130–191) with 1σ errors; models only the higher QPO member; no HFQPO search.
- **Ma, X., Zhang, L., Tao, L. et al. (2023)**, ApJ 948, 116; arXiv:2303.00481. ν_QPO for 6 NICER epochs.
- **Li, Y., Yan, Z., Gao, C., Yu, W. (2024)**, MNRAS; arXiv:2407.08421. Table C1 covers ObsIDs 196/197 only; the "six > 3σ" remark refers to the decay (C1).
- **Basak, A., Uttley, P., Bollemeijer, N. et al. (2025)**, MNRAS; arXiv:2508.12053. Four ObsIDs in two hue bins; no Lorentzian table.
- **Fiori, M., Zampieri, L. et al. (2025)**, A&A 697, A222; arXiv:2504.06001. nicerl2 reprocessing of ObsIDs 101–278 (HEAsoft 6.29, CALDB 20210707, 1–12 keV); NICER QPO track only in Fig. 13.
- **Dziełak, M. A., De Marco, B., Zdziarski, A. A. (2021)**, MNRAS; arXiv:2102.11635. ObsIDs 103/104 only.
- **Kawamura, T., Done, C., Axelsson, M., Takahashi, T. (2023)**, MNRAS; arXiv:2209.14492. ObsIDs 106/108 only.

**New MAXI J1820 timing papers and context**
- **Kimura, M., Negoro, H., Yamada, S., Iwakiri, W., Sako, S., Ohsawa, R. (2024/25)**, "Evolution of X-ray and optical rapid variability during the low/hard state in the 2018 outburst of MAXI J1820+070 = ASASSN-18ey", PASJ, doi:10.1093/pasj/psae099; arXiv:2411.03602. Shots last about 0.2 s; hard shots are 0.31 s against 0.37 s soft; six phases H1–H6.
- **Uttley, P. & Malzac, J. (2025)**, "Large and complex X-ray time lags from black hole accretion disks with compact inner coronae", MNRAS 536, 3284; doi:10.1093/mnras/stae2514; arXiv:2312.08302. r_cor = 10 R_g; soft lag above about 1 Hz; τ ∝ ν^−0.7; Fig. 1 uses ObsID 104.
- **Davidson, E., Chauhan, J., Lohfink, A., Russell, T. D., Starling, R., Johnson, C. (2025)**, arXiv:2505.07962 (submitted to ApJL). Link between corona and jet; spectral.
- **Yang, P., Méndez, M., Rout, S. K. et al. (2026)**, "An atypical X-ray variability component in the black hole candidate AT2019wey", A&A accepted; arXiv:2606.17976. Paper verified; its use of MAXI J1820 is not.
- **Omama et al. (2023)**, ApJ 945, 92; arXiv:2302.11263. NICER cross-correlation lags; missing from rounds 1 and 2.
- **Fogantini, F. A., García, F., Méndez, M., König, O., Wilms, J. (2025)**, A&A 696, A237; arXiv:2503.03078. Analyses Cyg X-1 only (C3).

**Other missions**
- **Mudambi, S. P., Maqbool, B., Misra, R., Hebbar, S., Yadav, J. S., Gudennavar, S. B., Bubbly, S. G. (2020)**, ApJL 889, L17; arXiv:2001.00755; doi:10.3847/2041-8213/ab66bc. LAXPC20, 2018-03-30; five Lorentzians; QPO at 47.7 mHz.
- **Prabhakar, G., Mandal, S., Athulya, M. P., Nandi, A. (2022)**, MNRAS; arXiv:2204.13363; doi:10.1093/mnras/stac1176. LAXPC QPO at 0.05 Hz, Q = 4.1, 9.1σ, rms 9.1% (MJD 58208).
- **Banerjee, S., Dewangan, G. C., Knigge, C. et al. (2024)**, ApJ 964, 189; arXiv:2402.08237; doi:10.3847/1538-4357/ad24ef. Spectral; two coronae; R_in ≳ 51 R_g; a = 0.77 ± 0.21.
- **Chakraborty, S., Navale, N., Ratheesh, A., Bhattacharyya, S. (2020)**, MNRAS; arXiv:2009.02465. AstroSat plus NuSTAR spectra (author list from a search summary).
- **Buisson, D. J. K., Fabian, A. C., Gandhi, P. et al. (2021)**, "MAXI J1820+070 with NuSTAR II…", MNRAS 500, 3976; arXiv:2011.02472; doi:10.1093/mnras/staa3510. ks-scale flares and dips; soft lag of 20 s.
- **Gao, C., Yu, W., Yan, Z. (2024/25)**, "On the Broadening of the Characteristic Frequency Range towards Higher Photon Energies…", ApJ (to appear); arXiv:2409.06414. Noise frequencies depend on energy.
- **Fan, N. et al. (2024)**, ApJ 969, 61; arXiv:2404.12161. HXMT. **Peng et al. (2025)**, arXiv:2502.14508. HXMT.
- **Zhang, L. et al. (2018)**, ATel #11951 (55 Hz claim) and **ATel #11960**, "The 55 Hz signal we detected in MAXI J1820+070 is not a QPO" (retraction; GTI gaps). Seen in search-result text and the ATel index; the ATel pages are blocked to our fetcher.
- **NICER/ISS Science Nugget, 2018-03-22** (heasarc.gsfc.nasa.gov/docs/nicer/science_nuggets/20180322.html). "Marginally detected a 66 Hz QPO in the early data", "a tentative detection"; rate rising from about 900 to more than 20,000 c/s. Read directly this round.

**Software**
- **Stingray PR #739** (Joshi 2023), `stingray/modeling/gpmodeling.py`; MIT; tinygp + jaxns.
- **Foreman-Mackey, D. et al.**, tinygp, JOSS doi:10.21105/joss.04455.
- **Dupourqué, S., Barret, D., Diez, C. M., Guillot, S., Quintin, E. (2024)**, jaxspec, A&A 690, A317; doi:10.1051/0004-6361/202451736; arXiv:2409.05757.
- **Nowak, M.**, SITAR v1.2.8 (space.mit.edu/cxc/analysis/SITAR/); **Nowak, M. A. et al. (2008)**, "Beyond XSPEC…", arXiv:0806.0560.
- **Buchner, J.**, BXA (github.com/JohannesBuchner/BXA).
- **Corsaro, E. & De Ridder, J. (2014)**, DIAMONDS, A&A 571, A71; arXiv:1408.2515. **Corsaro, E., McKeever, J. M., Kuszlewicz, J. S. (2020)**, FAMED, A&A 640, A130; doi:10.1051/0004-6361/202037930; arXiv:2006.08245.
- **Gruberbauer, M., Kallinger, T., Weiss, W. W., Guenther, D. B. (2009)**, A&A 506, 1043; doi:10.1051/0004-6361/200811203; arXiv:0811.3345.
- **Hübner, M. T.**, QPOEstimation (GitHub; no licence).
- **Gúrpide, A. & Middleton, M. (2025)**, "Mind the gaps…", arXiv:2501.05602; MNRAS doi:10.1093/mnras/staf196 (per the repository); Zenodo 10.5281/zenodo.14600069.
- **Lefkir, M. et al. (2025)**, Pioran.jl, MNRAS 539, 1775; arXiv:2501.05886. Also arXiv:2510.14529, which notes that Pioran has no QPO component.
- **Yu, W. & Richards, G. T.**, EzTao, ascl:2201.001.
- **Marcel et al. (2020)**, JED-SAD V, arXiv:2005.10359. Recipe for automated QPO detection.
- **XSPEC/flx2xsp practice:** arXiv:2111.13403 (GX 339-4), 2402.14226 (RX J0440.9+4431), 2403.09341 (GRS 1915+105), 2410.07333 (Swift J1727.8-1613). Author lists not retrieved.

**NICER instrument and calibration**
- **HEASARC**, "NICER Data Analysis Tips & Caveats" (accessed 2026-11-05).
- **Markwardt, C. B. & Hare, J. (2024)**, CALDB release notes xti20240206. **Markwardt, C. B. & Steiner, J. (2020)**, release notes xti20200722. Source of the ~2 Crab pile-up threshold.
- **NICER Team (2018-08-14)**, NICERDAS Calibration Guide.
- **HEASARC analysis threads** for HEASoft 6.32/NICERDAS 11 and 6.33/NICERDAS 12; **HEASoft release notes** 6.33.0, 6.35.2, 6.36.0, 6.37.0, 6.37.1; **help pages** for nicerl2 and niautoscreen.
- **Webbe, R. & Young, A. J. (2024)**, RASTI; arXiv:2409.12574; doi:10.1093/rasti/rzae040.
- **Bollemeijer, N., Uttley, P., You, B. (2026)**, MNRAS; arXiv:2609.27608. Appendix C: bright-source recipe for ObsIDs 107–112.
- **Enoto, T. et al. (2021)**, Science 372, 187; doi:10.1126/science.abd4659; arXiv:2104.03492. Also the Enoto 2021 NICER workshop slides.
- **Tuo, Y. et al. (2022)**, ApJS; arXiv:2109.04709; doi:10.3847/1538-4365/ac4250.
- **Hamaguchi, K. (2019)**, IACHEC talk.
- **Deneva, J. S. et al. (2019)**, ApJ 874, 160 (already in round 2).
- **Bogdanov, S. et al. (2019)**, ApJL 887, L25; arXiv:1912.05706.
- **O'Connor, B. et al. (2023)**, GCN 34549 (already in round 2).
- **Prigozhin, G. et al. (2016)**, Proc. SPIE 9905, 99051I; doi:10.1117/12.2231718. DOI pairing confirmed; full text not retrieved.

**NS, AGN, ULX and monitoring**
- **McHardy, I. M., Koerding, E., Knigge, C., Uttley, P., Fender, R. P. (2006)**, Nature 444, 730; doi:10.1038/nature05389.
- **Körding, E. G., Migliari, S., Fender, R., Belloni, T., Knigge, C., McHardy, I. (2007)**, MNRAS 380, 301; arXiv:0706.2959; doi:10.1111/j.1365-2966.2007.12067.x.
- **Scaringi, S. et al. (2015)**, Sci. Adv. 1, e1500686; arXiv:1510.02471.
- **González-Martín, O. & Vaughan, S. (2012)**, A&A 544, A80; arXiv:1205.4255.
- **Altamirano, D. et al. (2008)**, ApJ 685, 436; arXiv:0806.0962; VizieR J/ApJ/685/436.
- **van Straaten, S., van der Klis, M., Wijnands, R. (2005)**, ApJ 619, 455; astro-ph/0410505.
- **van Straaten, S., van der Klis, M., Méndez, M. (2003)**, ApJ 596, 1155; astro-ph/0207646. Exact title not confirmed. Companion: astro-ph/0307041.
- **Boutelier, M., Barret, D., Miller, M. C. (2009)**, MNRAS 399; arXiv:0907.3223. **Barret, D. et al. (2011)**, arXiv:1011.4601.
- **Gierliński, M. et al. (2008)**, Nature 455, 369 (round 2). **Middleton, M., Uttley, P., Done, C. (2011)**, arXiv:1106.1294. **Alston, W. N. et al. (2014)**, arXiv:1407.7657.
- **Strohmayer, T. E. et al. (2007)**, astro-ph/0701390. **Strohmayer, T. E. & Mushotzky, R. F. (2009)**, arXiv:0911.1076. Volumes not retrieved for either.
- **Pasham, D. R., Strohmayer, T. E., Mushotzky, R. F. (2014)**, Nature 513, 74; arXiv:1501.03180. **Pasham, D. R. et al. (2013)**, arXiv:1308.1677.
- **Negoro, H. et al. (2016)**, PASJ 68, S1; doi:10.1093/pasj/psw016.
- **Krimm, H. A. et al. (2013)**, ApJS 209, 14; arXiv:1309.0755.
- **Tetarenko, B. E., Sivakoff, G. R., Heinke, C. O., Gladstone, J. C. (2016)**, ApJS 222, 15; arXiv:1512.00778.
- **OHMAN**: JAXA/ISAS topic 003208, RIKEN/Chuo press release (Oct 2022), BAAS 2023 abstract. No refereed paper.

### 8.2 Unverified or partly verified: check before citing

- **ATel #11423** (NICER team 2018). The title appeared in search results. Its claimed 66 Hz hint is unverified; cite the 2018-03-22 Science Nugget instead.
- **Remillard, R. A. et al. (2022)**, AJ 163, 130; arXiv:2105.09901. The Crab rate of 10,531 ± 60 c/s comes from a snippet only.
- **Tuo et al. 2022**: the NICER systematic of 5.5 µs comes from a snippet only.
- **arXiv:2010.08291**, a review of kHz QPOs in neutron stars. Authors not confirmed.
- **Strohmayer & Mushotzky 2003**, the M82 X-1 ~54 mHz QPO. From memory; not searched.
- **Sanna et al. 2010** (XTE J1701-462). Seen only through Barret 2011.
- **Barret, Olive & Miller 2006.** Seen only as a citation.
- **Yang P. 2026's use of MAXI J1820** as a comparison source (snippet).
- **The S&K counts of about 155 QPO rows and about 76 at σ ≥ 3** (day ≤ 117) are rough counts from the PDF text. Recount from the VizieR table.
- **The mapping of Mudambi's observation to ObsIDs 113–115** was inferred from the date. Check against the ObsID MJDs.
- **The Körding-plane prediction of ν_l ≈ 2–7 Hz** is our own estimate, not from a paper.
- **Pile-up moving soft photons into 2–11 keV and biasing rms** is our own inference.
- **The licences of BXA, DIAMONDS, FAMED and SITAR** were not seen.
- **The `lowmemscr` default** is inconsistent between the code default (200) and the help text (500).

---

## 9. Remaining open questions

1. **VizieR keyword search** for other MAXI J1820 per-exposure timing tables, especially HXMT papers in ApJ (Ma X. 2021 Nature Astronomy lists 3 tables and no supplement is mentioned). The CDS main site reset the connection.
2. **An ADS full-text search** for multi-source NICER samples that use MAXI J1820 without naming it in the title or abstract.
3. **Whether S&K's day-116.9 rows cover the type-C→B portion of ObsID 197.** Settle this during the cross-match.
4. **Quantitative pile-up and telemetry effects at 400–1000 c/s per FPM.** No refereed study and no per-MPU model exists. A small simulation of packet-limited gaps plus about 22 µs dead time per FPM would itself be a methods contribution.
5. **TIMEZERO handling under `ftextractor`** (HEASoft 6.37). Test it directly.
6. **The origin of the 2018 66 Hz hint.** Re-derive it from the early ObsIDs under the new pipeline, and include it in the HFQPO null tests.
7. **Fundamental-versus-harmonic identification at ObsIDs 130–155.** Ratios of 1.76–2.3 in published tables suggest a broad component blending with the pair. This needs the joint PDS plus cross-spectrum fit (Méndez 2024) to settle.
8. **The full text of Uttley & Malzac 2025 (second half)** and of Prigozhin 2016 were not read.