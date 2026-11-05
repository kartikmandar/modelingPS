# Literature review, round 2: extended directions

*Project: modelingPS. The pipeline fits sums of BPK02 Lorentzians (chi^2 + COBYLA, F-test for component addition, Hessian errors) to 88 NICER 2-11 keV (rms/mean)^2 power spectra of MAXI J1820+070. These are ObsIDs 1200120101-197, covering the 2018 rise from the hard state to the HIMS. It also tracks the type-C QPO. This round extends `literature/review_round1.md` and does not repeat it. Round-1 papers come up again only where a new specific point is added.*

> **Search-coverage caveat (read first).** Eight topic agents gathered this round's sources. The shared WebSearch budget (200 calls per turn) ran out partway through. Topic 7 (NICER reduction) was cut short. Topic 9 (reproducibility) ran no keyword searches and relied on fetching primary pages directly. The second completeness critic could not run at all. When this document was put together, one more search was tried and refused because of the budget. So **the novelty check in Section 10 is incomplete**, and every claim of novelty below is hedged ("to our knowledge, pending a full search"). Section 12.3 lists the open queries for a follow-up round.

---

## 0. Executive summary: what is new compared with round 1

Round 1 covered the MAXI J1820 NICER and HXMT timing papers, the BPK02 method, the Protassov F-test caveat, Whittle likelihoods, NICER dead time and 3C50, Stingray, vKompth and Lense-Thirring (LT) precession. Round 2 adds eight new directions:

1. **A population context (Section 1).** About 40 papers now catalogue type-C QPOs and broadband-noise Lorentzians in other black-hole transients: MAXI J1348-630, Swift J1727.8-1613, MAXI J1535-571, GX 339-4, EXO 1846-031, MAXI J1803-298, MAXI J1631-479, 4U 1630-47, GRS 1915+105, H1743-322 and Cyg X-1. They share three features:
   - A recurring **critical type-C frequency of about 2-3 Hz**, where lags change sign, the rms turns over or the Q behaviour breaks.
   - An **inclination and jet dependence of the QPO lag sign**.
   - **Corona size against nu_QPO** curves from vKompth.

   MAXI J1820 crosses 2-3 Hz during the rise, so it can be placed directly on these relations.
2. **A NICER hue benchmark that already includes MAXI J1820.** Zhang Y. et al. 2026 (arXiv:2602.07247) give a NICER power-colour convention and a nu_QPO-hue track for MAXI J1820, rising from about 0.2 Hz in the LHS to about 5 Hz in the HIMS. Basak et al. 2025 (arXiv:2508.12053) order the MAXI J1820 NICER PDS by hue. These papers are both a benchmark and a **partial novelty threat** (Section 10).
3. **Polarimetry (Section 2).** IXPE finds a hard-state corona extended in the disc plane, with PD about 4% and PA along the jet. The first QPO-phase-resolved polarimetry (Swift J1727) disagrees between papers. Zhao 2024 finds no PD or PA modulation. Ninoyu 2025 finds a marginal (2.9 sigma) PD modulation. A Fourier polarimetry-timing estimator has been verified on IXPE data (Ewing 2025). Round 1 lacked the MAXI J1820 optical, IR and radio timing papers (Paice 2021, Tetarenko 2021, Baglio 2025, Echiburu-Trujillo 2024); they are added here.
4. **Time-domain methods (Section 3).** Two papers on our own source and epoch were missing from round 1:
   - **Yu et al. 2023** (Hilbert-Huang transform, HHT): the QPO width comes mainly from frequency modulation.
   - **Zhou et al. 2022** (the QPO and the noise combine multiplicatively, which lowers the fundamental's rms by 2-3 times).

   Celerite Gaussian-process fits of black-hole QPOs now exist (Wang 2025, Dong 2026). So do wavelet S-factors and bicoherence studies of MAXI J1820 with HXMT (Ding 2023).
5. **Spectral-timing products from event data (Section 4).** Three QPO phase-resolved methods could be ported to NICER: Fourier waveform reconstruction, the cross-correlation function and HHT phase. Others are covariance spectra (Yang 2025 already has a MAXI J1820 two-Comptoniser result), rms spectra per Lorentzian, frequency-resolved reflection and **joint PDS plus cross-spectrum Lorentzian fitting**. That last one reveals "hidden" QPOs in Cyg X-1 (Fogantini 2025), and the paper reports the same features in MAXI J1820.
6. **Statistics (Section 5).**
   - Savage-Dickey Bayes factors are a direct replacement for the F-test (Malacaria 2024).
   - Nested-sampling evidences handle non-nested model comparisons.
   - Bayer and Seljak give the look-elsewhere correction.
   - Segment and parametric bootstraps check the Hessian errors.
   - Simulation-based inference handles dead-time bias.
   - **We found no hierarchical model across ObsIDs in the X-ray QPO literature.**
7. **Updates to NICER reduction (Section 7).** Reprocessing 2018 data with HEASoft 6.32-6.37 is **not equivalent** to the older archive products:
   - New `mingti` and `erodedilate` defaults change the GTIs.
   - Per-MPU lowmem screening was added.
   - A lowmem bug, fixed only in 6.35.2, screened MPU0 only.
   - Noise-ringer windows were added.
   - Older versions had a `.mkf` overwrite bug.

   Bright-source practice is to use MPU-disjoint cospectra (Bollemeijer 2025). Telemetry saturation produces a 55 Hz artefact.
8. **Theory to test (Section 8).** At least five nu_QPO mappings compete: rigid LT, LT slowed by coupling to the disc, JED-SAD (nu_K(r_J)/q), the accretion-ejection instability (AEI) and a transition layer. There are also new GRMHD results on tearing, MADs and magnetic truncation, and a disc-warp model of the type-C to type-B transition (Marcel 2026). JED-SAD fits of MAXI J1820 already exist (Marino 2021), so **the nu_QPO-r_J test can be done right away.**

Three findings change the plan most:

- **Plan (b) components need re-examination.** Lorentzians can hide in a power-only fit: the narrow QPO in Fogantini 2025 and the type-B shoulder in Jin 2025. The additive Lorentzian model itself has been challenged by Zhou 2022 and by bicoherence studies (Maccarone 2011). So the component counts the F-test chose, and the rms values of the fundamental, need to be checked again.
- **The screening version is a systematic.** Reprocessing defaults can change the GTIs and the active MPUs in the most telemetry-saturated ObsIDs. The screening version must be pinned and audited.
- **Novelty has to rest on method, not on coverage.** NICER type-C tracks of MAXI J1820 already exist (Stiele & Kong 2020 in round 1, Zhang Y. 2026, Basak 2025). What is new must be uniform, validated, statistically calibrated and publicly released Lorentzian decomposition, plus the new diagnostics.

---

## 1. Comparative black-hole transient catalogues and universal correlations

### 1.1 Same-band, same-method NICER analogues
- **MAXI J1348-630.**
  - Alabarta 2022 (arXiv:2205.11581, MNRAS 514, 2839) fits multiple Lorentzians to NICER PDS. The type-C rms rises up to about 2-3 keV and is flat out to 12 keV. Its lags are hard, and the QPO in the reflare matches the main-outburst QPO.
  - Alabarta et al. (arXiv:2409.14883, preprint submitted to ApJ, about 19 authors) fit vKompth to the whole outburst. They need two coupled coronae, about 1e3 km and about 1e4 km in size.
  - Wang X.-L. 2026 (arXiv:2603.22100, HXMT) follows the type-C from 0.24 to 10.3 Hz. They find nu-flux hysteresis loops that run in opposite directions in the main and mini outbursts. The rms spectrum hardens (peaking above 20 keV, above 10% amplitude) as the jet weakens.
  - Two caution papers: 2512.13124 (type-B and type-C coexisting) and 2406.12477 (an atypical hard-state LFQPO). Not every narrow hard-state Lorentzian is a type-C.
- **Swift J1727.8-1613.** This is the most direct comparison for the MAXI J1820 rise.
  - Rawat 2025 (arXiv:2504.06705, A&A 697, A229) uses 45 NICER observations with nu from 0.3 to 7 Hz. The 10 vs 0.5 keV lag is about +0.5 rad at 0.34 Hz, crosses zero near 1.2 Hz and reaches about -1.1 rad near 7 Hz. Two coronae of about 6.5e3 and 2e3 km shrink to about 2e3 km, and R_in falls from about 30-40 to about 10 Rg.
  - Liao 2024 (arXiv:2410.06574) fits vkompthdk through a flare, with nu rising from 0.33 to 2.63 Hz.
  - Jin P. 2025 (arXiv:2510.10353, A&A 706, A208) finds a **type-B in every observation, hidden in the HIMS as a broad shoulder on the type-C**. It reaches Q of 6 or more at the flare peak, which coincides with three radio ejections.
  - Bollemeijer 2025 (arXiv:2505.05060) finds the waveform pivots at about 15-20 keV.
  - 2405.09772 finds the frequency depends on energy.
  - 2605.18050 detects the type-C above 250 keV.
- **MAXI J1535-571.**
  - Rawat 2023 (arXiv:2301.04418) finds the corona shrinks from about 1e4 to 3e3 km as nu goes from 1.8 to 9 Hz. The relations break at **about 3.0 Hz, close to when the radio jet is quenched**.
  - Huang 2018 (arXiv:1808.05318) finds soft lags at the centroid and argues for high inclination.
  - Rawat 2023b (arXiv:2307.11460) tests the dynamic-frequency model.
- **GX 339-4 (2021).**
  - Stiele & Kong 2023 (arXiv:2304.06077) use the same authors and method as their MAXI J1820 paper (round 1). They find nu_QPO is anti-correlated with the size of the Comptonising region.
  - Jin Y.-J. 2025 (arXiv:2504.04763, A&A 697, A120) fits **three broadband-noise Lorentzians (L1, L2, L3)**. Their frequencies rise with nu_QPO, and each one's rms spectrum peaks in a different band. This is close to our setup.
  - Zhang Z. 2024 (arXiv:2305.18249, ApJ 971, 148) ties nu_QPO to R_in and mdot in GX 339-4 and EXO 1846-031.
- **Other NICER and HXMT catalogues.**
  - EXO 1846-031 (Liu 2021, arXiv:2009.10956): the rms-nu relation turns over at about 2 Hz.
  - MAXI J1803-298 (Zhu 2023, arXiv:2306.13993), a dipping source.
  - MAXI J1631-479 (Rout 2021, arXiv:2105.08426; Bu 2021, arXiv:2107.06333). The QPO rms is flat up to 100 keV while the flat-top noise drops below 0.1% above 30 keV.
  - 4U 1630-47 (arXiv:2501.13163): the Q behaviour changes at about 2.31 Hz. A 2025 outburst is in arXiv:2607.02228.
  - Swift J1658.2-4242 (arXiv:1910.04354).
  - A review of HXMT black-hole timing: Galaxies 13, 111 (2025).

### 1.2 Reference populations
- **GRS 1915+105.**
  - Zhang L. 2020 (arXiv:2003.08928) uses 620 RXTE observations. The lags change sign at about 2 Hz, and the QPO rms peaks near zero lag.
  - Karpouzas 2021 (arXiv:2103.09675) finds the rms-nu turnover at the same frequency.
  - Mendez 2022 (Nature Astronomy 6, 577, arXiv:2203.02963) finds hard lags below 2 Hz only when the radio is strong.
  - **Bellavita 2026 (arXiv:2609.15672)** uses NICER data from 2018-2020 at low radio flux (5 mJy or less). It finds only soft lags, which supports the claim that the hard lags are driven by the jet.
  - Garcia 2022 (arXiv:2204.13279) applies vKompth.
  - 2311.13467 finds the lags vary within single observations, a caution against averaging per ObsID.
  - 2510.10017 is an HXMT compilation.
- **H1743-322.**
  - Shui 2023 (arXiv:2211.16679) finds the type-C **rms scales with the non-thermal fraction**. At low nu, the rms in the decay is lower than in the rise.
  - The 2018 failed outburst (doi:10.1093/mnras/stac773) stays in the LHS, with nu between about 0.1 and 0.4 Hz.
- **Cyg X-1 with NICER.**
  - König 2024 (arXiv:2405.07754) finds two humps that split near 1 Hz, with a lag jump and a coherence drop that follow the Lorentzians.
  - **Fogantini 2025 (arXiv:2503.03078, A&A 696, A237)** fits the PDS and the real and imaginary cross spectrum together. A coherence dip and a lag "cliff" move from 1 to 6 Hz as Gamma goes from 1.8 to 2.4. Fitting them needs a low-rms, large-lag narrow Lorentzian, which they call a hidden type-C. **The paper says the same features appear in MAXI J1820.**
- **The AstroSat archive.** Ram & Chakraborty 2026 (arXiv:2605.17344) cover 14 sources and 29 QPOs. In high-inclination sources the lag changes sign near 2 Hz: hard lags reach about 36 ms below it and soft lags stay within about 10 ms above it. Low-inclination sources show no such change.

### 1.3 Correlations that MAXI J1820 can be placed on
| Relation | Population evidence | MAXI J1820 test |
|---|---|---|
| Critical nu of about 2-3 Hz (lag sign, rms turnover, Q break) | GRS 1915 at about 2 Hz; EXO 1846 at about 2 Hz; MAXI J1535 at about 3.0 Hz; 4U 1630 at about 2.31 Hz; AstroSat high-i sources at about 2 Hz. *Swift J1727 at about 1.2 Hz, which is lower, so the scale is not universal.* | Broken-line fits of rms, Q, FWHM and lag against nu_QPO |
| Lag sign against inclination and jet | High-i sources show soft lags or a sign change; low-i sources show none; in GRS 1915 hard lags appear only when the radio is strong | Two-band cross spectra plus the 2018 radio light curves |
| rms against non-thermal fraction, and rise/decay hysteresis | H1743 (Shui 2023); MAXI J1348 nu-flux loops | Needs reprocessed spectra and the decay ObsIDs |
| nu_QPO(R_in, mdot) | GX 339-4, EXO 1846, MAXI J1535, H1743 | nicerl2 spectral fits |
| Broadband-noise nu against nu_QPO (WK99, PBK99) | GX 339-4 L1, L2, L3; Cyg X-1 humps | Uses existing fits (cheap) |
| Corona size against nu_QPO (vKompth) | 1e4 km falling to about 2-3e3 km for nu about 1 to 9 Hz | Needs rms and lag spectra |

**Inclination correction (critic consolidation).** The topic agents quote MAXI J1820's inclination as "about 63-66 deg", "about 64 deg (jet)" and "about 63-75 deg". The paper should adopt one value with one cited source. The jet inclination from Atri 2020 (round 1) is the usual choice. If the orbital inclination differs, both should be stated, because the misalignment matters for LT models (Poutanen 2022, round 1).

---

## 2. Polarimetry and multiwavelength

### 2.1 IXPE results for the coronal geometry
- **Hard and HIMS states: PD about 4%, PA along the jet. The corona is extended in the disc plane, not a lamppost.**
  - Cyg X-1: 4.0 ± 0.2% (Krawczynski 2022, Science 378, 650, arXiv:2206.09972). Poutanen 2023 (arXiv:2302.11674) show this needs an outflowing medium of at least about 0.4c.
  - Swift J1727 hard state: 4.1 ± 0.2% (Veledina 2023, arXiv:2309.15928).
  - Swift J1727 HIMS: about 4% falling to about 3%, with constant PA aligned with the radio (Ingram 2024, arXiv:2311.05497). The same paper argues that the long NICER soft lags (above about 10 ms) are not light-travel delays.
  - Swift J1727 soft state: 1.2% or less (Svoboda 2024, arXiv:2403.04689).
  - Swift J1727 dim decaying hard state: 3.3 ± 0.4% (Podgorný 2024, arXiv:2404.19601). The geometry follows the state, not the luminosity.
- **High inclination.** IGR J17091-3624 has a type-C at about 0.2 Hz and PD = 9.1 ± 1.6% (Ewing 2025, arXiv:2503.22665). Debnath 2025 (arXiv:2504.16402) find PD rising to about 30% at 6-8 keV. **This is the best analogue for MAXI J1820's inclination.**
- **Other sources.**
  - GX 339-4 in the SIMS: 1.3%, aligned with the ejecta (Mastroserio 2025, arXiv:2408.06856). This is the closest analogue to MAXI J1820's type-C to type-B switch and ejection.
  - Cyg X-1 across 13 pointings: hard state about 4.0%, soft state about 2.2%, with no superorbital signal (Kravtsov 2025, arXiv:2505.03942).
  - 4U 1630-47: soft-state PD of 6-10% (Ratheesh 2024, doi:10.3847/1538-4357/ad226e). An IXPE plus NICER analysis is Rawat 2023 (arXiv:2303.02745).
  - LMC X-1: PD below 1.1% (Podgorný 2023, arXiv:2303.12034).
  - LMC X-3: PD about 3% and a spin of about 0.2 (Garg 2023, arXiv:2309.06830; Svoboda 2024, arXiv:2309.10813).
  - Review: Dovčiak 2024 (Galaxies 12, 54).

### 2.2 Polarimetry resolved by QPO phase, and the precession question
- **Predictions.** Ingram 2015 (arXiv:1505.00015) predict that PD and PA oscillate at the QPO frequency. eXTP forecasts give about 1.4% rms in PD and about 4 deg in PA.
- **Observations (Swift J1727).** These conflict:
  - **Zhao 2024 (arXiv:2401.08970, ApJL), using HHT:** Gamma is strongly modulated with QPO phase, but PD and PA are not. The authors call this inconsistent with simple LT precession.
  - **Ninoyu 2025 (ApJ 985, 69), using cycle stacking:** PD is 5.1% in high-flux phases and 3.8% in low-flux phases, at about 2.9 sigma, and the PA behaviour is inconsistent between observations.

  The predicted modulation is close to IXPE's sensitivity limit, so **precession has not been ruled out.** Because the PA stays constant through the HIMS, any precession has to be small in amplitude, or the region that modulates the QPO has to differ from the main polarised Comptoniser.
- **Methods and forecasts.**
  - The Fourier polarimetry-timing estimator of Ingram & Maccarone 2017 (arXiv:1707.06659) was implemented and verified on IXPE pulsars by Ewing 2025 (arXiv:2507.15461).
  - GRMHD plus Monte Carlo polarisation of tilted, truncated discs is in Fragile 2025 (arXiv:2505.11446).
  - eXTP 2025 (arXiv:2506.08105) forecasts a 3.6 sigma detection of the modulation in 350 ks at 500 mCrab, and notes IXPE's telemetry cap of about 15 million counts.
  - The original eXTP mission paper is Zhang S.N. 2019 (arXiv:1812.04020).
  - The STROBE-X concept is Ray 2019 (arXiv:1903.03035).
- **LT test with the fitted corona radius.** Ma, Done & Kubota 2025 (arXiv:2506.18857) find r_cor shrinking from 45 to 9 Rg in Swift J1727. The LT slope matches but the normalisation is offset. This is the same group as the MAXI J1820 dual-corona paper (Ma R. 2023, round 1), so the method is ready to apply to MAXI J1820.

### 2.3 Multiwavelength timing of MAXI J1820 missing from round 1
- **Paice 2021 (arXiv:2105.11769)** covers 80 days of ULTRACAM, HiPERCAM and NICER data, overlapping **our ObsID range**. They find an optical lag of about +0.2 s, a growing anticorrelation, and optical lags coherent with the QPO in at least two epochs.
- **Tetarenko 2021 (arXiv:2103.09318, MNRAS 504, 3862)** time 10 bands at once. Jet parameters: Gamma about 6.8, opening angle about 0.45 deg, and P_jet about 0.6 L(1-100 keV). Zdziarski 2022 (arXiv:2108.10929) caution that the jet Lorentz factor is not uniquely constrained.
- **Echiburu-Trujillo 2024 (arXiv:2311.11523):** the jet spectral break moves over about 3 decades in frequency during the outburst.
- **Baglio 2025 (arXiv:2510.16124, A&A 704, A198):** the optical rms peaks in the hard state and rises with wavelength, and winds are seen in the re-brightenings.
- **Optical polarimetry:** Veledina 2019 (arXiv:1808.09002) and Kosenkov 2020 (arXiv:2005.09699) find PD of 0.3-0.7% during the rise, linked to winds.
- **Kajava 2019 (arXiv:1906.06519):** XMM-Newton sees **X-ray dips** during the rise from a partially covering absorber, and a complex UV/X-ray cross-correlation function. If dips fall inside NICER GTIs they bias the (rms/mean)^2 and add low-frequency red noise, so they need to be screened (Section 6).
- **Fender, Belloni & Gallo 2004 (astro-ph/0409360)** give the hardness-intensity diagram (HID), hysteresis and jet-line framework needed to place the HIMS end of our sample relative to the July 2018 radio flare.

---

## 3. Time-domain and machine-learning methods

### 3.1 Papers on our source missing from round 1
- **Yu et al. 2023 (arXiv:2305.12317, ApJ 951, 130)** apply HHT to HXMT data over MJD 58191-58286. **The QPO peak is broadened mainly by frequency modulation**, so a Lorentzian FWHM measures frequency wandering more than a true coherence time. They also find negative QPO lags that grow in size with energy.
- **Zhou et al. 2022 (MNRAS 515, 1914, doi:10.1093/mnras/stac1789)** argue the QPO and the broadband noise combine by convolution (a multiplicative PDS), not addition. In the multiplicative model the fundamental's rms is **2-3 times lower** (about 3.8% against 9.1% at 1-2.6 keV), while the centroid and FWHM barely change. *This came only from an agent's takeaway, with no source entry, so it needs checking (Section 12.2).*
- **Ding et al. 2023 (RAA 23, 085024, doi:10.1088/1674-4527/acdc88; no arXiv id)** compute bicoherence with HXMT. They see a "cross" pattern at nu_QPO about 0.045 Hz and a "hypotenuse" pattern at about 0.23 Hz. They compare a type-B-like pattern with **NICER ObsID 1200120139, which is in our sample.**

### 3.2 Competing explanations for QPO incoherence
- Lachowicz & Done 2010 (arXiv:0908.2740), using wavelets and matching pursuit on XTE J1550: many independent wave packets, each at roughly constant frequency, with a range of lifetimes.
- Su 2015 (arXiv:1510.09122) and Hu 2014 (arXiv:1404.5108), using HHT and O-C analysis: an intermittent oscillation wandering over 3-5 Hz.
- Yu 2023 (above) supports frequency modulation for MAXI J1820.
- NICER's count rates make it possible to test which picture fits in 2-11 keV.

### 3.3 Wavelet S-factor
The S-factor is the fraction of time the QPO is "on", per GTI. It is about 0 for type-B and finite for type-C (Chen & Wang 2024, arXiv:2402.11293, NICER MAXI J1535). For type-C QPOs it **decreases with nu_QPO** in MAXI J1803 (Jin 2024, arXiv:2410.13367, MNRAS 535, 207). Chen 2022 (arXiv:2209.10408) separate QPO-on and QPO-off intervals and compare their spectra. The S-factor is cheap to compute for each ObsID.

### 3.4 Gaussian processes
celerite (Foreman-Mackey 2017, arXiv:1703.09710) works in O(N), and its SHO and exponential terms are Lorentzians in frequency. Applications so far:
- MAXI J1348-630 type-B with an SHO plus damped-random-walk (DRW) kernel: Q about 10 (Wang Y. 2025, arXiv:2510.05587, A&A 703, A134).
- **The full Swift J1727 type-C sweep from 0.07 to 5 Hz, with Q rising from 3 to 10** and viscosity alpha about 0.1 (Dong 2026, arXiv:2604.01901).
- Hübner 2022 (ApJ 936, 17), now built into Stingray.

Cautions: preprocessing and detrending can create false QPOs (Song & Mao 2024, arXiv:2405.04033), and a GP test rejected 5 of 6 claimed AGN QPOs (Zhang H. 2023, arXiv:2304.08044). *We found no GP analysis of MAXI J1820 (search incomplete).*

### 3.5 Bicoherence
- Arur & Maccarone 2019 (arXiv:1904.05895) and 2020 (arXiv:1910.01687): during the HIMS-to-SIMS transition the "web" pattern turns into a "cross" pattern in **high-inclination** sources and a "hypotenuse" pattern in low-inclination ones.
- Applications: MAXI J1535 (Zhu Z. 2024, arXiv:2403.13606) and Swift J1727 (arXiv:2408.17039).
- For details see Section 4.

### 3.6 Machine learning and simulation-based inference
- Simulation-based inference for periodograms biased by dead time (Huppenkothen & Bachetti 2022, arXiv:2104.03278).
- Unsupervised light-curve fingerprints (Orwat-Kapola 2022, arXiv:2110.10063).
- Autoencoder state maps (Ricketts 2023, arXiv:2301.10467).
- Supervised classification of GRS 1915 (Huppenkothen 2017, arXiv:1611.01332).
- K-Means state labels that match human labels 81-95% of the time (Sreehari & Nandi 2021, MNRAS 502, 1334).
- Koopman/EDMD (Miao 2026, arXiv:2609.01734): each eigenfunction contributes one Lorentzian to the PSD, which gives the BPK02 sum a dynamical-systems reading.
- **We found no published CNN or transformer that classifies QPO type from PDS** (search incomplete).

### 3.7 Shift-and-add and dynamic PDS
Dynamic PDS in 16-64 s segments can separate slow drift of the centroid within an ObsID from intrinsic coherence. Shifting and adding the PDS biases the QPO narrower when the per-segment frequency errors are large.

---

## 4. Spectral-timing techniques enabled by event data

### 4.1 QPO phase-resolved spectroscopy (three routes)
1. **Fourier waveform reconstruction.**
   - The method is Ingram & van der Klis 2015 (arXiv:1411.1967). First check that the fundamental-harmonic phase difference is stable.
   - A tomography follow-up is Ingram 2017 (arXiv:1610.00948). Two illumination patches fit best, and the reflection fraction varies with QPO phase at about 3.5 sigma.
   - **Nathan 2022 (arXiv:2201.01765)** is the closest NICER precedent: NICER plus NuSTAR on GRS 1915. The Fe-line centroid shifts with phase, but the inferred R_in is too large for LT. The paper discusses NICER systematics in detail.
2. **The cross-correlation function** (Stevens & Uttley 2016, round 1). Public code is at github.com/abigailStev; we could not confirm it is still maintained.
3. **HHT or VMD phase.**
   - Shui 2024 (arXiv:2407.18106, ApJ 973, 59), with NICER, NuSTAR and HXMT on Swift J1727: the disc temperature leads the other components by about 0.4-0.5 cycles. *This is the same paper that Topic 1 listed as "(authors not retrieved)". The entries are merged here.*
   - Shui 2024b (arXiv:2403.18272) uses a low-energy phase reference to recover the QPO above 170 keV.
   - For MAXI J1820 only HXMT HHT and VMD studies exist (2311.03635, round 1). **We found no NICER phase-resolved study of the soft band below 2 keV.**

### 4.2 Bicoherence and biphase
- The founding papers are Maccarone 2011 (arXiv:1101.1979), which shows the QPO is **phase-coupled to the broadband noise** and that the coupling depends on radio state, and Maccarone 2013 (arXiv:1308.3150), which explains how to read the biphase.
- Maccarone & Schnittman 2005 (astro-ph/0411266) show how to tell resonance broadening from phase jitter.
- Arur & Maccarone 2022 (arXiv:2205.12981) link the bispectrum to the radio state of GRS 1915.
- The biphase at (f_QPO, f_QPO) gives the fundamental-harmonic phase difference and so the waveform asymmetry. Bollemeijer 2025 publish a reproduction package (Zenodo record 14243337).
- **Taken together, these question the BPK02 assumption that the Lorentzians are independent and simply add.**

### 4.3 Covariance, rms and frequency-resolved spectra
- The covariance spectrum (Wilkinson & Uttley 2009, arXiv:0905.0587) shows the disc below 1 keV varying more than the Comptonised emission on timescales above 1 s. **Our 2-11 keV cut throws this regime away.**
- **Yang S.-K. 2025 (arXiv:2511.17285)** compute MAXI J1820 covariance spectra up to 150 keV. Coherence with 2-10 keV drops above 30 keV, and kTe differs between short and long timescales. They model this with two Comptonisation zones.
- QPO rms spectra lack a disc component and run opposite to the hardness of the time-averaged spectrum (Sobolewska & Życki 2006, astro-ph/0603024).
- The harmonic and the QPO have different spectra (Axelsson 2014, arXiv:1307.4396; Axelsson & Done 2016, MNRAS 458, 1778). A forward model from a precessing flow is Życki 2016 (arXiv:1610.07871).
- Frequency-resolved Fe-line spectroscopy is from Revnivtsev 1999 (astro-ph/9906198).
- Kawamura 2023 (arXiv:2209.14492) fit a propagating-fluctuation model with pivoting to MAXI J1820 NICER and HXMT data, and read the QPO as an external LT modulation.
- Basak 2025 (arXiv:2508.12053) compare MAXI J1820 and Cyg X-1. MAXI J1820 has soft lags at 1-10 Hz that Cyg X-1 lacks.

### 4.4 Joint PDS and cross-spectrum Lorentzian fitting
This follows Méndez 2024 (round 1) and is now used in Fogantini 2025 (hidden QPO) and Rawat 2025. Each Lorentzian gets a constant phase lag or a constant time lag. It extends our COBYLA and F-test framework directly.

### 4.5 Software
- Stingray 2 (Section 9).
- HENDRICS (ascl:1805.019) makes cross spectra and rms, lag and covariance-energy spectra from event files.
- Ingram 2022 (arXiv:2212.10238) is a tutorial on the error formulae.
- **No mainstream package implements bicoherence or IvdK15 waveform reconstruction**, so those need custom code or the Bollemeijer package.

---

## 5. Advanced statistics and hierarchical modelling

- **Savage-Dickey density ratio (SDDR), a replacement for the F-test.** Malacaria 2024 (arXiv:2310.16498) fit a Whittle likelihood with JAX and blackjax NUTS. They put a uniform prior on the added Lorentzian's amplitude that includes A = 0. They get log B = -30.3, which becomes -23.4 when the prior volume is 100 times larger. This is an explicit test of prior sensitivity.
- **Evidence from nested sampling.** dynesty (Speagle 2020, arXiv:1904.02180) and UltraNest (round 1) compute evidences valid for non-nested comparisons, such as a harmonic against a subharmonic against a broad hump.
- **Look-elsewhere effect.**
  - Bayer & Seljak 2020 (arXiv:2007.13821) show the evidence penalises it through the prior volume on the centroid, so **the prior range on the centroid must be stated and justified.**
  - Vaughan 2016 (arXiv:1606.02620) show red noise produces false periodicities and that multiplicity needs correcting. This applies to searching 88 ObsIDs.
  - Inglis 2015 (arXiv:1410.8162) show that a white-noise null inflates significance.
- **Frequentist tests for non-nested models.**
  - Algeri, Conrad & van Dyk 2016 (arXiv:1509.01010): a comprehensive model plus the Gross-Vitells correction.
  - A uniform threshold of ΔBIC ≥ 10 across 675 flares (Inglis 2016, arXiv:1610.07454).
  - A ΔAICc < 2 stopping rule (NGC 247 ULX-1, arXiv:2104.11163).
  - AIC for X-ray QPOs (arXiv:2501.01581).
- **Error checks.**
  - Bootstrap over segments with 1000 resamples, with phase wrapping handled (Nathan 2022).
  - Parametric bootstrap with 1000 refits and asymmetric intervals (M82 ULXs, arXiv:2505.16921).
  - Hessian errors are expected to fail for weak, broad or boundary-hugging Lorentzians.
- **Simulation-based inference (SBI)** for dead time, background and GTI gaps (Huppenkothen & Bachetti 2022).
- **Cross-validation across methods.** Report a QPO only when Fourier, wavelet and GP methods agree (Huppenkothen 2025, arXiv:2504.10153). A Bayesian X-ray workflow with posterior predictive checks and QQ plots is in Buchner & Boorman 2023 (arXiv:2309.05705).
- **Hierarchical modelling.** We found **no hierarchical or GP-prior model that pools Lorentzian parameters across the observations of one outburst.** Every evolution paper fits each PDS on its own and draws the track afterwards. This looks novel, pending a full search.
- **Label switching.** The general literature (Stephens 2000; ordering constraints) is relevant, but no X-ray paper deals with label switching of Lorentzian components.
- **Correction to round 1.** Huppenkothen et al. 2013 (ApJ 768, 87) is **arXiv:1212.1011**, not 1303.6648.

---

## 6. Component tracking and state classification

- **Power colours and hue.**
  - Heil, Uttley & Klein-Wolt 2015a (arXiv:1405.2024): PC1 = var(0.25-2 Hz)/var(0.0039-0.031 Hz) and PC2 = var(0.031-0.25 Hz)/var(2-16 Hz). These can be computed from our PDS without fitting.
  - Heil 2015b (arXiv:1405.2026): **the QPO shifts the hue by up to about 20 deg in high-inclination sources**, so report hue both with and without the QPO.
  - Gardenier & Uttley 2018 (arXiv:1809.06093): neutron stars follow the same track.
- **NICER convention and MAXI J1820 benchmark.** Zhang Y. et al. 2026 (arXiv:2602.07247, MNRAS 546, stag261, 15 authors; exact title begins "A NICER view of the corona…") use a PC origin of (4.5, 0.45) and a 45 deg reference axis. Hue ranges: LHS -20 to 140 deg, HIMS 140-220 deg, SIMS 220-300 deg. They give a nu_QPO-hue track for 9 sources, **including MAXI J1820**. Ram & Chakraborty 2026 track components ordered by hue, including harmonics and shoulders.
- **rms-intensity diagram and rms-flux relation.**
  - The rms-intensity diagram is from Muñoz-Darias 2011 (arXiv:1008.0558).
  - Heil 2012 (arXiv:1202.5877): the flux intercept of the rms-flux relation **becomes strongly positive from the hard state into the HIMS**.
  - The relation was discovered by Uttley & McHardy 2001 (astro-ph/0103367). Gleissner 2004 (astro-ph/0311039) extend it.
  - **The QPO has its own rms-flux slope, which changes sign with nu** (positive below about 4 Hz, negative above about 6 Hz; Heil 2011, arXiv:1011.6321). We found no test of this for MAXI J1820.
  - Scargle 2020 (arXiv:2001.08314) caution that a linear rms-flux relation does not uniquely imply multiplicative variability.
- **Component identity.**
  - **A new point from a round-1 paper:** in MAXI J1820 the lowest Lorentzian (L1) follows nu_QPO = 2 nu_L1, so it behaves as a subharmonic, not a break on the WK99 track (Yang Z.-X. 2022, 2204.00739). Automatic labelling should test the subharmonic hypothesis.
  - Jin P. 2025 (type-B shoulder) and Fogantini 2025 (hidden QPO) add more candidates.
  - Kawamura 2023b (arXiv:2304.12003) explain the long soft lags with propagating fluctuations whose sign flips. *The paper's link to MAXI J1820 NICER data was seen only in search snippets.*
- **Systematics: dips.** Kajava 2019 dips should be found from hardness light curves within each GTI, and the dipping segments flagged or removed.
- **Tracking methods.** The automated methods found are clustering and ML (Sreehari 2021, Huppenkothen 2017, Orwat-Kapola 2022, Ricketts 2023). **We found no published tracker ordered by hue or with neighbouring-observation priors.**

---

## 7. NICER reduction updates (2022-2026)

*Coverage was incomplete because of the search budget: pulse pile-up, NICER-specific Crab absolute offsets and the timing-thread author were not covered.*

| Version | Change relevant to 2018 MAXI J1820 data |
|---|---|
| **6.32 / NICERDAS 11** | New defaults: `mingti=5`, `erodedilate=5` (removes GTI pieces of 5 s or less), `keep_noisering=NO`, `thresh_range=-3-3`. niautoscreen may now return **zero good time** instead of keeping about 8 detectors. MPU_LOWMEM columns widened from 8 to 32 bits. The thread notes that at extreme rates (about 55 packets/s) MPU transmission limits cause many gaps per second. *Linking this to the 55 Hz feature in ObsID 107 is our inference, not a statement in the thread.* 6.32.0 had a bug in the per-MPU normalisation of nifpmsel; use 6.32.1 or later. |
| **6.33.2 / NICERDAS 12a / xti20240206** | 3C50 for the gain epoch "2022". SCORPEON v23 is the default (not interchangeable with v22). Per-MPU `lowmemscr` (200 lost events per GTI) **can drop single MPUs during a campaign**, which changes the FPM count, the dead time and the Poisson level. Before this version nicerl2 silently ignored `gtifiles`. |
| **6.34 / NICERDAS 13** | Bug fixed: with `tasks=SCREEN`, an archived `.mkf.gz` overwrote the new filter file, losing the NICERV5 columns. **With archive auxil/ files, run `tasks=ALL` on 6.34 or later.** New version keywords: NIBKGTYP, SCORPVER, BGYR3C50. |
| **6.35.2 / NICERDAS 14** | niautoscreen v1.6 **fixes lowmem screening, which had acted only on MPU0**, so results differ from 6.33-6.35.1. New: `lowmemfrscr` (lost-event fraction, suited to 20-50 kc/s), `autoscreen_tseg`, NICERV6 columns, and a `timerange` that removes TIME < 100 after MPU reboots. |
| **6.36** | barycorr handles FPM_SEL and filter files. A full reinstall is needed. |
| **6.37.1** | The newest seen. Official Dockerfiles (build 1.05), Section 9. *Correction: Topic 7 calls 6.36 "current", but Topic 9 found 6.37.1.* |

- **Noise ringers.** These are spurious events within about 110 µs after an undershoot. Screening them out adds up to about 5% effective dead time at 500 ct/s/FPM. This is not in the Webbe & Young model, so **fit the Poisson level rather than fixing it.**
- **Light leak.** It began in May 2023, so it does not affect 2018 data. But the tools have since been tuned for it, so run 2018 data with night or default threshold settings.
- **Bright-source practice** (Bollemeijer 2025; O'Connor 2023, GCN 34549):
  - Telemetry saturation is per MPU.
  - Use cospectra from **disjoint FPMs in disjoint MPUs**.
  - Measure the noise level at 200-250 Hz.
  - Treat any feature near 55 Hz as instrumental.
  - Expect many zero-GTI ObsIDs under the default `underonly_range`.
  - FAD (Bachetti & Huppenkothen 2018, round 1) can use MPU groups A and B as independent detectors.
  - An HXMT FAD application is 2407.07807, *seen as a title only*.
- **Timing.**
  - MJDREF = 56658 + 0.000777592592592593, TT. Add TIMEZERO (about -1 s). extractor 5.37 and earlier drops it, so use `niextract-events`. Accuracy is about 80 ns.
  - barycorr with `refframe=ICRS ephem=JPLEPH.430` and the `.orb` file.
  - Cross-mission absolute timing from the IACHEC Crab work: Bachetti 2026 (PASP 138, 014501, arXiv:2509.13865) and Terada 2025 (arXiv:2503.13807). Missions agree within 200 µs or better, and the ephemerides agree within about 5 µs or better.
  - NICER timing was stable in 2017-18 (Deneva 2019, arXiv:1902.07130; Vivekanand 2020, arXiv:1911.10028).
  - A 2021 workshop timing deck uses **MAXI J1820 ObsID 1200120107 (5912 GTIs)** as its window-function example.
- **SCORPEON** has no refereed paper. Cite the HEASARC threads and the version keyword. At our rates the background affects the rms by much less than 1%, but it should still be reported.
- **Coincidence to avoid.** Musoke 2023 GRMHD predicts an HFQPO near 55 Hz from tearing, while telemetry artefacts also sit near 55 Hz. Any HFQPO search (Section 8) must rule out the instrument first.

---

## 8. Theory to test

### 8.1 Five competing mappings for nu_QPO
1. **Rigid LT precession of a hot flow**, with the outer radius set by truncation.
   - Ingram & Done 2012a (arXiv:1108.0789): fits give r_o from about 68 to 13 Rg. The same paper ties the broadband noise (propagating fluctuations) and the QPO together.
   - Ingram & Done 2012b (arXiv:1208.0728): nu rises with flux and with softening.
   - Motta 2018 (arXiv:1709.02608): a flow wider than about 10-100 Rg aligns before it can precess coherently, and the frequencies set lower limits on spin.
2. **LT slowed by the outer disc** to about 5% of the isolated-torus rate (Bollimpalli 2023, arXiv:2210.02944). This biases any r_trunc derived from nu by about a factor of 20. Bollimpalli 2024 (arXiv:2312.14876) add epicyclic Lorentzians at the edge of the torus.
3. **JED-SAD and jet instability.**
   - nu_QPO = nu_K(r_J)/q, with q = 133 ± 4 for GX 339-4 (Marcel 2020, arXiv:2005.10359).
   - Ferreira 2022 (arXiv:2202.00438) propose the jet-instability origin.
   - Marcel 2022 (arXiv:2109.13592) give the radiative efficiency (15-40%).
   - **Marino 2021 (arXiv:2109.00218) already fit JED-SAD to 8 MAXI J1820 hard-state epochs**, and Barnier 2022 (arXiv:2109.02895) provide public XSPEC tables.
   - Limits: the Swift J1727 JED-SAD polarisation prediction fails (arXiv:2609.15301).
4. **The accretion-ejection instability (AEI).** omega is about 0.1-0.3 Omega(r_in), with a turnover near the ISCO (Varnière 2002, astro-ph/0203325).
5. **A transition layer.** The index-nu relation saturates near Gamma of about 2.7 (Titarchuk & Fiorito 2004, astro-ph/0403198).

### 8.2 GRMHD complications
- Jets follow the tilted disc axis (Liska 2018, arXiv:1707.06619).
- A disc-corona phase lag appears (Liska 2019, arXiv:1901.05970).
- **Disc tearing** gives an inner sub-disc that precesses for only 1-5 cycles at about 3 Hz, plus a 55 Hz HFQPO, so Q is limited and the QPO is intermittent (Musoke 2023, arXiv:2201.03085).
- Magnetic truncation (Liska 2022, arXiv:2201.03526) and flux-dependent outbursts (Liska 2023, arXiv:2309.15926).
- **MADs precess retrograde** under magnetic torque, with weak LT (Jiang 2025, arXiv:2507.13680).

### 8.3 Broadband noise
- Viscous damping puts a break at the lowest viscous frequency (Mushtukov 2018, arXiv:1707.07578).
- A break-frequency formula as a function of the fluctuation timescale, h/r and M (Turner & Reynolds 2021, arXiv:2103.09281).
- Analytic relativistic Green's functions (Mummery 2023, arXiv:2305.18314).
- GRMHD gives a linear rms-flux relation (Bollimpalli 2020, MNRAS 496, 3808).
- Three Compton zones in Mahmoud & Done 2018 (arXiv:1803.04811) correspond to 2-3 broadband Lorentzians.

### 8.4 The end of the HIMS
- **Disc-warp unification (Marcel 2026, arXiv:2511.10474, A&A 710, A387):** at the type-C to type-B transition the QPO amplitude drops sharply and the broadband noise disappears.
- Precessing-jet type-B QPOs (Kylafis 2020, arXiv:2008.02557).

### 8.5 HFQPOs and the relativistic precession model (critic addition)
- A systematic RXTE HFQPO search (Belloni, Sanna & Méndez 2012, arXiv:1207.2311) found only 11 detections among 22 sources.
- RPM mass and spin from the QPO triplet of GRO J1655-40 (Motta 2014, arXiv:1309.3652).
- For M of about 8.5 Msun, each observed nu_QPO together with an assumed spin predicts HFQPO frequencies. A NICER upper limit at those frequencies would test the RPM.

### 8.6 The dominant systematic: spin
Estimates range from a = 0.1-0.3 (continuum fitting) to a > 0.8 (reflection, the RPM and soft-state fits; arXiv:2603.12950). **Any r_trunc(t) derived from LT must be reported for both spin branches.**

### 8.7 Scaling to other masses (critic addition)
- The break timescale scales with mass and accretion rate (McHardy 2006, astro-ph/0612273).
- RE J1034+396 (Gierliński 2008, arXiv:0807.1899).
- QPOs in tidal disruption events: Pasham 2019 (arXiv:1810.10713), and Pasham 2024 (arXiv:2402.09689), a TDE precession reading with spin between 0.05 and 0.5.

---

## 9. Reproducibility and data release

- **Stingray 2** (Bachetti 2024, JOSS 9, 7389, doi:10.21105/joss.07389).
  - Products for a bright NICER observation in about 1 s.
  - Lag, covariance and coherence spectra, GP QPO models and power colours.
  - The modeling API: `PSDParEst.fit` (Whittle likelihood, with AIC and BIC), `sample()` (emcee), **`calibrate_lrt()` (a simulated likelihood-ratio test that is valid under Protassov)**, `calibrate_highest_outlier()` and `fit_lorentzians()`.
  - Cite both papers plus a DOI for the version used (concept DOI 10.5281/zenodo.1490116).
  - **Pitfall:** astropy `Lorentz1D` is parameterised by (amplitude, x_0, FWHM), while BPK02 uses (nu_max, Q, r). The conversion between them needs a unit test.
- **Three code paths for cross-validation:** our COBYLA code; Stingray or HENDRICS `HENmodel`; and XSPEC or ISIS through `HEN2xspec` or `ftflx2xsp`. `ftflx2xsp` makes a unit-diagonal response and accepts Hz as the x unit. XSPEC's chi^2 is Gaussian, so the PDS must be heavily rebinned.
- **Environment.**
  - HEASARC provides HEASoft 6.37.1 Dockerfiles but no registry image. The image defaults to the **remote** CALDB, so mount the local xti20240206 CALDB and GEOMAG_PATH.
  - heasoftpy provides `hsp.nicerl2(..., noprompt=True)` with `allow_failure`.
  - The nicerl2 thread says: run `tasks=ALL clobber=YES` once; do not use `incremental=YES` after a calibration change; use `cldir` to write outside the archive tree.
  - Do not run tools incrementally across HEASoft versions (NICER versions thread).
- **Workflow and practice.**
  - Snakemake (Mölder 2021, doi:10.12688/f1000research.29032.2).
  - Ten Simple Rules (Sandve 2013).
  - Best and Good-Enough practices (Wilson 2014, arXiv:1210.0530; Wilson 2017, arXiv:1609.00037).
  - Momcheva & Tollerud 2015 (arXiv:1507.03989): 90% of astronomers write their own code but only 8% have substantial training, which justifies cross-validating against other codes.
- **FAIR release** (Wilkinson 2016).
  - An ECSV or FITS catalogue with one row per ObsID and component, units, errors and test statistics.
  - The PDS and model arrays behind every figure.
  - Seeds, the container digest and the version keywords (NIBKGTYP, SCORPVER, ftversion, nicaldbver).
  - A Zenodo DOI obtained before submission.
- **Background model.** The choice (3C50, SCORPEON or Space Weather) changes the mean rate used in the normalisation, so record it as metadata and as a column for the systematic.
- **Open item.** We could not search for an existing public MAXI J1820 Lorentzian or QPO table on Zenodo, VizieR or GitHub.

---

## 10. Novelty check

**What already exists (verified this round or in round 1):**
- NICER type-C tracking of the MAXI J1820 rise with multiple Lorentzians: Stiele & Kong 2020 (round 1).
- **A NICER nu_QPO-hue track and vKompth fits for MAXI J1820 inside a 9-source sample** (Zhang Y. 2026, arXiv:2602.07247).
- PDS of MAXI J1820 ordered by hue and matched to Cyg X-1 (Basak 2025).
- Joint PDS and cross-spectrum features (Fogantini 2025 report them in MAXI J1820).
- HXMT studies over the whole outburst: Ma X. 2021/2023 (round 1), HHT (Yu 2023), bicoherence (Ding 2023), broadband noise with a subharmonic L1 (Yang 2022).
- A JED-SAD spectral evolution (Marino 2021).
- Covariance spectra (Yang S.-K. 2025).
- Optical and X-ray timing (Paice 2021).

**What we could not check (budget exhausted):**
- Later per-ObsID tables from HXMT, NICER or AstroSat. Possible AstroSat LAXPC papers, such as a recalled "Mudambi et al. 2020", are *unverified*.
- Whether Stiele & Kong 2020 or Bollemeijer (round 1) already publish a per-ObsID Lorentzian table for 1200120101-197.
- Recalled papers that might overlap: Axelsson & Veledina 2021 (frequency-resolved spectroscopy of MAXI J1820) and any VizieR or Zenodo tables. *Also unverified.*

**Verdict (hedged).** To our knowledge, and pending a full search, **a full Lorentzian catalogue of the 2018 rise has not been published as an open, uniform, statistically calibrated release.** That would mean all 88 ObsIDs with per-component tables, consistent rules for selecting components, validated errors, current-CALDB reprocessing with a screening audit, and released data. Coverage of the rise alone is **not** novel. What could be novel:
1. A uniform, pre-registered rule for selecting components (SDDR or calibrated LRT instead of the F-test), with a public table of decisions.
2. Validated errors (Hessian checked against segment bootstrap, parametric bootstrap and MCMC) and three-code agreement.
3. A **hierarchical model across ObsIDs** (no precedent found).
4. An audit of how the 2023-2025 screening changes affect the 2018 timing products (no precedent found).
5. NICER bicoherence across the sequence and the rms-flux relation of the QPO (not found for MAXI J1820).
6. A NICER GP time-domain cross-check (not found for MAXI J1820).
7. Testing the nu_QPO-r_J relation (JED-SAD) on a second source (not found).
8. Population placement against the 2-3 Hz critical frequency and the lag-jet relation, using MAXI J1820's radio coverage.
9. A FAIR release.

Anything that depends only on nu_QPO(t) or on hue tracking should be presented as a confirmation, not as a discovery.

---

## 11. Map of research directions (ranked)

Ranking weighs payoff against effort, with the existing PDS and fits taken as the starting point. Effort: L = days, M = weeks, H = months.

| # | Direction | What it needs | Effort | Payoff | Key refs |
|---|---|---|---|---|---|
| 1 | **Power colours and hue for all 88 ObsIDs, with and without the QPO; overlay the published MAXI J1820 track** | Existing PDS; the Heil bands; the Zhang 2026 convention | L | High: a fit-free state axis and a benchmark | 1405.2024, 1405.2026, 2602.07247, 2508.12053 |
| 2 | **Replace the F-test with SDDR, a calibrated LRT or evidence; publish the decision table** | JAX/NUTS or Stingray `calibrate_lrt`; Whittle likelihood | M | High: settles the Protassov problem; disagreements are themselves a result | 2310.16498, 2007.13821, 1606.02620, Stingray 2 |
| 3 | **Validate the Hessian errors** with segment bootstrap, parametric bootstrap and MCMC, and check coverage on injection tests | Simulations at NICER rates | L-M | High: underpins every number | 2201.01765, 2505.16921 |
| 4 | **Screening-sensitivity audit** (6.36/6.37 defaults vs relaxed; MPU survival; Poisson level) plus MPU-disjoint cospectra and FAD | Scripted nicerl2; Stingray | M | High: robustness and a new systematic | HEASoft 6.32-6.36 threads, 2505.05060, 1709.09700 |
| 5 | **Population placement:** broken-line fits of rms, Q and FWHM against nu_QPO; broadband-noise frequencies against nu_QPO; subharmonic test | Existing fits plus digitised literature tables | L-M | High: 2-3 Hz critical-frequency test | 2003.08928, 2301.04418, 2605.17344, 2504.04763, 2204.00739 |
| 6 | **Joint PDS and two-band cross-spectrum Lorentzian fits**; search for a hidden QPO and a type-B shoulder; lag against nu | Reprocessed events; about 30 more parameters | M | High: may change the component IDs; gives the lag-jet test | 2503.03078, 2510.10353, Méndez 2024 |
| 7 | **Additive vs multiplicative (convolution) QPO model**; bicoherence check of independence | New model function; bicoherence code | M | High: the fundamental's rms could change by 2-3 times | Zhou 2022 (check), 1101.1979 |
| 8 | **nu_QPO against r_J test** (JED-SAD), plus ranking the five theory mappings jointly against the broadband breaks | Marino 2021 r_J or Barnier tables; spectral fits; both spin branches | M | High: first test of the relation on a second source | 2005.10359, 2109.00218, 2109.02895, 2210.02944, 2603.12950 |
| 9 | **NICER bicoherence sequence** (web, cross, hypotenuse against nu and hue); biphase waveform | Custom code; noise-floor simulations | M | Medium-high: first with NICER; inclination test | 1910.01687, Ding 2023, 1308.3150 |
| 10 | **Hierarchical GP-prior model across ObsIDs** (log nu, log Q, rms against MJD or hue); label switching | numpyro/PyMC; Whittle likelihood | H | High, novel methodology | §5 (no precedent found) |
| 11 | **rms spectra and covariance spectra per Lorentzian**, 0.3-12 keV, then vKompth | Refits in 10-15 energy bands | M | Medium-high: corona size against nu; Yang 2025 two zones | 0905.0587, 2511.17285, 2504.06705, astro-ph/0603024 |
| 12 | **rms-intensity diagram, rms-flux intercept and QPO rms-flux slope against nu** | Light curves in segments | L-M | Medium: test of the Heil 2011 sign change | 1008.0558, 1202.5877, 1011.6321 |
| 13 | **Wavelet S-factor, HHT intermittency and dynamic PDS** per ObsID | PyEMD or wavelets; red-noise simulations | M | Medium: tests FM against wave packets; Q bias | 2305.12317, 0908.2740, 2402.11293, 2410.13367 |
| 14 | **GP time-domain cross-check** (celerite2 SHO plus DRW) on a subset | Binned 2-11 keV light curves | M | Medium: first for MAXI J1820 | 2510.05587, 2604.01901, 1703.09710 |
| 15 | **Energy-dependent centroid** (2-4, 4-7, 7-11 keV) | Refits | L | Medium: constrains differential precession | 2405.09772 |
| 16 | **Dip screening** (Kajava 2019) and telemetry-quality covariates | HK and .mkf files; hardness light curves | L | Medium: removes systematics | 1906.06519, HEASoft threads |
| 17 | **Phase-resolved spectroscopy of the soft band** (IvdK15, CCF or HHT) in bright HIMS ObsIDs | Phase method; spectral fits | M-H | Medium-high: first study of the disc below 2 keV | 1411.1967, 2201.01765, 2407.18106 |
| 18 | **Multiwavelength alignment** (optical rms and lags, jet break against nu_QPO; radio against lag sign) | Digitised tables | L-M | Medium | 2105.11769, 2311.11523, 2510.16124, 2203.02963 |
| 19 | **Extend to the decay and reflare ObsIDs**: hysteresis at matched hue; rms against non-thermal fraction | New downloads; same pipeline | M | Medium-high | astro-ph/0409360, 2211.16679 |
| 20 | **HFQPO upper limits** (50-500 Hz, FAD or cospectra) against RPM predictions | Bright HIMS ObsIDs | L-M | Medium: spin dispute; must exclude the 55 Hz artefact | 1207.2311, 1309.3652 |
| 21 | **Reproducibility stack:** three-code validation, Snakemake, a container, pytest injection tests, a FAIR Zenodo release | Engineering | M | High for credibility and reuse | Stingray 2, heasoftpy, FAIR, Sandve 2013 |
| 22 | **Polarisation-modulation prediction** for MAXI J1820, a polarimetry-timing simulator, and eXTP/STROBE-X forecasts | Fitted rms, a geometry model | M | Medium (forward-looking) | 1505.00015, 2507.15461, 2506.08105, 1812.04020, 1903.03035 |
| 23 | **SBI for the dead-time bias** at peak rates; a CNN or flow classifier of QPO type | Simulator; sbi | H | Medium (novel ML) | 2104.03278 |
| 24 | **Corona radius against nu with SSsed** as an LT test | Broadband spectra | M-H | Medium | 2506.18857 |
| 25 | **Mass-scaling and TDE placement** for the discussion | Literature values | L | Low-medium | astro-ph/0612273, 0807.1899, 2402.09689 |

---

## 12. Annotated bibliography

### 12.1 Verified (seen on a search result page or a primary page this round)

**§1 Catalogues and correlations**
- 2205.11581: Alabarta+2022, MNRAS 514, 2839. MAXI J1348 multi-Lorentzian fits with NICER; shape of rms against energy.
- 2409.14883: Alabarta+ (preprint submitted to ApJ, about 19 authors; cite as "et al."). Two coronae from vKompth.
- 2603.22100: Wang X.-L.+2026, A&A. MAXI J1348 type-C from 0.24 to 10.3 Hz; nu-flux hysteresis.
- 2512.13124: 2025 (authors not retrieved). Type-B and type-C coexisting in MAXI J1348.
- 2406.12477: 2024 (authors not retrieved). Atypical hard-state LFQPO in MAXI J1348.
- 2504.06705: Rawat+2025, A&A 697, A229. Swift J1727 with NICER; lag crosses zero at 1.2 Hz; corona shrinks.
- 2410.06574: Liao+2024/25, ApJ. Swift J1727 flare, vkompthdk.
- 2510.10353: Jin P.+2025, A&A 706, A208. Type-B hidden as a shoulder in the HIMS.
- 2505.05060: Bollemeijer, Uttley & You 2025, MNRAS (doi:10.1093/mnras/staf750). Waveform pivot; reduction recipe for bright sources.
- 2405.09772: 2024 (authors not retrieved). Energy-dependent frequency of the Swift J1727 type-C.
- 2605.18050: 2026 (authors not retrieved). Swift J1727 type-C above 250 keV.
- 2301.04418: Rawat+2023, MNRAS. MAXI J1535 break at 3.0 Hz.
- 1808.05318: Huang+2018, ApJ. MAXI J1535 soft lags.
- 2307.11460: Rawat+2023. Test of the dynamic-frequency model.
- 2304.06077: Stiele & Kong 2023. GX 339-4 with the same method as their MAXI J1820 paper.
- 2504.04763: Jin Y.-J.+2025, A&A 697, A120. Three broadband-noise Lorentzians in GX 339-4.
- 2305.18249: Zhang Z.+2024, ApJ 971, 148. nu_QPO against R_in and mdot.
- 2009.10956: Liu H.-X.+2021, RAA 21, 70. EXO 1846, turnover at 2 Hz.
- 2306.13993: Zhu+2023. MAXI J1803.
- 2105.08426: Rout+2021, MNRAS. MAXI J1631 with NICER.
- 2107.06333: Bu+2021, ApJ. MAXI J1631 QPO compared with broadband noise above 30 keV.
- 2501.13163: 2025 (authors not retrieved). 4U 1630, Q break at 2.31 Hz.
- 2607.02228: 2026 (authors not retrieved). 4U 1630, 2025 outburst.
- 2003.08928: Zhang L.+2020, MNRAS. GRS 1915 lags against nu, 620 observations.
- 2203.02963: Méndez+2022, Nature Astronomy 6, 577 (doi:10.1038/s41550-022-01617-y). Corona-jet coupling.
- 2609.15672: Bellavita+2026, A&A. NICER GRS 1915, only soft lags at low radio.
- 2103.09675: Karpouzas+2021. rms-nu turnover in GRS 1915.
- 2204.13279: García+2022, MNRAS. vKompth applied to GRS 1915.
- 2311.13467: 2023 (authors not retrieved). Fast-varying QPO lags.
- 2510.10017: 2025 (authors not retrieved). HXMT GRS 1915 compilation.
- 2211.16679: Shui+2023, ApJ 943, 165. H1743 rms against non-thermal fraction.
- doi:10.1093/mnras/stac773: 2022 (authors not retrieved). H1743 2018 failed outburst.
- 2405.07754: König+2024, A&A 687, A284. Cyg X-1 with NICER, two humps.
- 2503.03078: Fogantini+2025, A&A 696, A237. Hidden QPO in Cyg X-1.
- 2605.17344: Ram & Chakraborty 2026, MNRAS. AstroSat lag-nu population.
- doi:10.3390/galaxies13050111: 2025 (authors not retrieved). Review of HXMT black-hole timing.
- 1910.04354: 2019/20 (authors not retrieved). Swift J1658.2-4242.

**§2 Polarimetry and multiwavelength**
- 2309.15928: Veledina+2023, ApJL 958, L16.
- 2311.05497: Ingram+2024, ApJ 968, 76.
- 2401.08970: Zhao Q.-C.+2024, ApJL. No PD or PA modulation with QPO phase.
- doi:10.3847/1538-4357/adcabd: Ninoyu+2025, ApJ 985, 69. Marginal PD modulation.
- 2403.04689: Svoboda+2024. 2404.19601: Podgorný+2024. 2506.18857: Ma, Done & Kubota 2025.
- 2206.09972: Krawczynski+2022, Science 378, 650. 2302.11674: Poutanen+2023. 2505.03942: Kravtsov+2025.
- 2408.06856: Mastroserio+2025. doi:10.3847/1538-4357/ad226e: Ratheesh+2024. 2303.02745: Rawat, Garg & Méndez 2023.
- 2303.12034: Podgorný+2023. 2309.06830: Garg+2023. 2309.10813: Svoboda+2024.
- 2503.22665: Ewing+2025. 2504.16402: Debnath+2025. Galaxies 12, 54: Dovčiak+2024.
- 1505.00015: Ingram+2015. 1707.06659: Ingram & Maccarone 2017. 2507.15461: Ewing+2025 (method).
- 2505.11446: Fragile+2025. 2506.08105: eXTP 2025. 2504.10981: Brigitte & Svoboda 2025.
- 1812.04020: Zhang S.N.+2019. 1903.03035: Ray+2019.
- 2105.11769: Paice+2021. 2103.09318: Tetarenko+2021. 1808.09002: Veledina+2019. 2005.09699: Kosenkov+2020.
- 2311.11523: Echiburu-Trujillo+2024. 2108.10929: Zdziarski+2022. 2510.16124: Baglio+2025.
- 1906.06519: Kajava+2019. astro-ph/0409360: Fender, Belloni & Gallo 2004.

**§3 Time domain and ML**
- 2305.12317: Yu+2023, ApJ 951, 130 (on our source).
- 1510.09122: Su+2015. 1404.5108: Hu+2014. 0908.2740: Lachowicz & Done 2010.
- 2403.18272: Shui+2024 (MAXI J1535 above 170 keV).
- 2209.10408: Chen+2022. 2402.11293: Chen & Wang 2024. 2410.13367: Jin+2024, MNRAS 535, 207.
- 2510.05587: Wang Y.+2025, A&A 703, A134. 2604.01901: Dong+2026. 1703.09710: Foreman-Mackey+2017.
- 2405.04033: Song & Mao 2024. 2304.08044: Zhang H.+2023. 2104.03278: Huppenkothen & Bachetti 2022.
- 2110.10063: Orwat-Kapola+2022. 2301.10467: Ricketts+2023. 1611.01332: Huppenkothen+2017. Sreehari & Nandi 2021, MNRAS 502, 1334.
- 2609.01734: Miao+2026 (Koopman). 2403.13606: Zhu Z.+2024. 2408.17039: Swift J1727 bicoherence, ApJ.
- doi:10.1088/1674-4527/acdc88: **Ding Q., Ji L., Bu Q.-C., Dong T., Chang J. 2023, RAA 23, 085024, "Nonlinear Variability Observed with Insight-HXMT in MAXI J1820+070 and MAXI J1535-571"**. Topic 4 gave the full author list; it replaces the partial entry from Topic 3.

**§4 Spectral timing**
- 1411.1967: Ingram & van der Klis 2015. 1610.00948: Ingram+2017. 2201.01765: Nathan+2022, MNRAS 511, 255.
- 2407.18106: **Shui Q.-C.+2024, ApJ 973, 59** (merged entry).
- 1904.05895, 1910.01687, 2205.12981: Arur & Maccarone 2019, 2020, 2022.
- 1101.1979: Maccarone+2011. 1308.3150: Maccarone 2013. astro-ph/0411266: Maccarone & Schnittman 2005.
- 0905.0587: Wilkinson & Uttley 2009. 2511.17285: Yang S.-K.+2025 (MAXI J1820 covariance spectra).
- astro-ph/9906198: Revnivtsev+1999. astro-ph/0603024: Sobolewska & Życki 2006.
- 1307.4396: Axelsson+2014. MNRAS 458, 1778: Axelsson & Done 2016. 1610.07871: Życki+2016.
- 2209.14492: Kawamura+2023. 2508.12053: Basak+2025 (14 authors, MNRAS accepted). 2212.10238: Ingram 2022.
- HENDRICS: ascl:1805.019. Stevens cross_correlation code (GitHub; maintenance not confirmed).

**§5 Statistics**
- 2310.16498: Malacaria+2024, A&A. 1904.02180: Speagle 2020. 1606.02620: Vaughan+2016.
- 2007.13821: Bayer & Seljak 2020. 1509.01010: Algeri+2016. 1610.07454: Inglis+2016. 1410.8162: Inglis+2015.
- 2104.11163: NGC 247 ULX-1. 2501.01581: mHz oscillations in an SMBH. 2505.16921: M82 ULX QPOs with NuSTAR.
- 2504.10153: Huppenkothen+2025. 2309.05705: Buchner & Boorman 2023.
- 1212.1011: Huppenkothen+2013 (round-1 id corrected).

**§6 Tracking and states**
- 1405.2024, 1405.2026: Heil+2015a,b. 1809.06093: Gardenier & Uttley 2018.
- 2602.07247: **Zhang Y.+2026, "A NICER view of the corona…", MNRAS 546, stag261, 15 authors**.
- 1008.0558: Muñoz-Darias+2011. 1202.5877: Heil+2012. 1011.6321: Heil+2011.
- astro-ph/0103367: Uttley & McHardy 2001. astro-ph/0311039: Gleissner+2004. 2001.08314: Scargle 2020.
- 2304.12003: Kawamura+2023b. 2204.00739: Yang Z.-X.+2022 (round 1; new point on the subharmonic L1).

**§7 NICER reduction**
- HEASARC threads: heasoft632, heasoft633, heasoft634, RelNotes 6.35.2, RelNotes 6.36, nicer-versions, SCORPEON overview, noise ringers, light leak, TIME definition, nicerl2, and the 2021 timing workshop slides.
- GCN 34549: O'Connor+2023.
- 2509.13865: Bachetti+2026, PASP 138, 014501. 2503.13807: Terada+2025.
- 1902.07130: Deneva+2019, ApJ 874, 160. 1911.10028: Vivekanand 2020, A&A 633, A57.

**§8 Theory**
- 1108.0789, 1208.0728: Ingram & Done 2012a,b. 1709.02608: Motta+2018.
- 1707.06619: Liska+2018. 1901.05970: Liska+2019. 2201.03085: Musoke+2023. 2201.03526: Liska+2022. 2309.15926: Liska+2023.
- 2210.02944: Bollimpalli+2023. 2312.14876: Bollimpalli+2024. 2507.13680: Jiang+2025.
- 2005.10359: Marcel+2020. 2202.00438: Ferreira+2022. 2109.00218: Marino+2021. 2109.02895: Barnier+2022. 2109.13592: Marcel+2022.
- 2511.10474: Marcel+2026, A&A 710, A387. 2609.15301: JED-SAD study of Swift J1727.
- 1707.07578: Mushtukov+2018. 2103.09281: Turner & Reynolds 2021. 2305.18314: Mummery 2023. Bollimpalli+2020, MNRAS 496, 3808.
- 1803.04811: Mahmoud & Done 2018. 2008.02557: Kylafis+2020. astro-ph/0203325: Varnière+2002. astro-ph/0403198: Titarchuk & Fiorito 2004.
- 2603.12950: MAXI J1820 soft-state spin (authors not captured).
- 1207.2311: Belloni, Sanna & Méndez 2012. 1309.3652: Motta+2014.
- astro-ph/0612273: McHardy+2006. 0807.1899: Gierliński+2008. 1810.10713: Pasham+2019. 2402.09689: Pasham+2024.

**§9 Reproducibility**
- doi:10.21105/joss.07389: Stingray 2 (Bachetti+2024). Stingray modeling API documentation. heasoftpy. HEASoft Docker page (6.37.1). ftflx2xsp help page.
- doi:10.1093/mnras/stab3437 (SBI, as above). Hübner+2022, ApJ 936, 17. Heil+2015, doi:10.1093/mnras/stv191.
- doi:10.1007/978-981-16-4544-0_137-1: Bachetti & Huppenkothen 2022, "Fourier methods".
- 2206.14220: Astropy 2022. doi:10.1038/sdata.2016.18: FAIR. doi:10.12688/f1000research.29032.2: Snakemake.
- doi:10.1371/journal.pcbi.1003285: Sandve+2013. 1210.0530: Wilson+2014. 1609.00037: Wilson+2017. 1507.03989: Momcheva & Tollerud 2015.

### 12.2 Unverified or partly verified: check before citing
- **Zhou et al. 2022, MNRAS 515, 1914 (doi:10.1093/mnras/stac1789).** The multiplicative QPO model for MAXI J1820 appears only in an agent's takeaway, with no source entry. Check the authors, title and the numbers (3.8% vs 9.1%).
- **Springford, Eadie & Thomson 2020, AJ 159, 205** (multitaper periodogram). Seen only in other papers' reference lists (verified = false).
- **Hübner+2022 arXiv id.** Not checked; the DOI was seen in the Stingray 2 reference list.
- **Ratheesh+2024 (4U 1630-47) arXiv id.** Not confirmed; the DOI is verified.
- **arXiv:2407.07807** (HXMT/ME FAD application). Title only.
- **Su 2018 NCU thesis** (mean high-amplitude QPO duration of about 1.45 s). Seen only in the thesis abstract.
- **Kawamura+2023b link to MAXI J1820 NICER data.** From search snippets only.
- **Topic 7's attribution of the 55 Hz feature in ObsID 1200120107 to MPU packet limits.** An inference, not stated in the HEASARC thread.
- **Topic 6's Ding 2023 nu_QPO range** ("0.02 to 0.51 Hz to 0.22 Hz") is garbled. Take the values from the paper.
- **Unconfirmed recollections; do not cite.** Mudambi+2020 (AstroSat MAXI J1820 QPOs), Axelsson & Veledina 2021 (frequency-resolved spectroscopy of MAXI J1820), and any Cavecchi PDS or QPO paper (his known work is on neutron-star bursts and burst oscillations). Two guessed arXiv ids, 2001.03424 and 2106.08385, resolved to **unrelated** papers and have been discarded.
- **Many entries with "authors not retrieved"** (2512.13124, 2406.12477, 2405.09772, 2605.18050, 2501.13163, 2607.02228, 2311.13467, 2510.10017, 1910.04354, 2104.11163, 2501.01581, 2505.16921, 2603.12950, stac773, Galaxies 13, 111). The ids or DOIs were seen, but the author lists need filling in.

### 12.3 Corrections from the critics and open queries for round 3
**Corrections applied**
- Huppenkothen 2013 is arXiv:1212.1011.
- 2407.18106 is Shui+2024, ApJ 973, 59 (entries merged).
- Duplicates across topics removed: 2504.04763, 2410.13367, 2402.11293, Ding 2023, 2104.03278, 2510.05587, 2604.01901, 2201.01765, 2508.12053, 2511.17285 and 1910.01687.
- Zhang Y. 2026: exact title begins "A NICER view…", 15 authors.
- Alabarta 2409.14883: a preprint, cited with "et al.".
- Méndez 2022 DOI added.
- HEASoft "current" changed from 6.36 to 6.37.1.
- The source of the inclination value must be made consistent.
- The 2-3 Hz critical frequency is not universal (Swift J1727 at about 1.2 Hz).

**Open queries, blocked by the search budget**
- NICER MAXI J1820 timing papers from 2024-2026.
- Per-ObsID Lorentzian or QPO tables for 1200120101-197 in Stiele & Kong 2020, Bollemeijer, HXMT and AstroSat papers, and on VizieR or Zenodo.
- AstroSat LAXPC papers on MAXI J1820.
- NICER pulse pile-up and NICER-specific Crab absolute offsets.
- Automated outburst-monitoring pipelines.
- Comparisons with neutron-star kHz QPOs beyond Klein-Wolt & van der Klis 2008.
- Other open PDS-fitting codes.

Sources: the one search tried while assembling this round (`MAXI J1820+070 NICER power spectra Lorentzian evolution type-C QPO 2018 outburst`) was refused because the turn's search budget was spent, and it returned no URLs. All entries above come from the round-2 topic-agent and critic inputs, which were verified as stated in each entry.