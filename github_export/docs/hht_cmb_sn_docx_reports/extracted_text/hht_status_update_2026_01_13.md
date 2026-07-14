# HHT_Status_Update_2026-01-13

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\HHT_Status_Update_2026-01-13.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2026-01-15T08:30:00Z`

## Extracted Text

HHT / Seasonal Robustness Investigation – Status Update

Update date: 2026-01-15Scope: CEEMDAN/HHT locking and robustness tests on SN residual fields; seed-grid stability; DR16×CMB cross-correlation workplan.

Prepared for internal research tracking. This document summarises what has been observed so far, what is statistically robust (and what is not), and what tests are currently in progress.

Executive summary

Analysis of DR16 Quaia quasar dataset Across the runs captured to date, two patterns stand out:

(1) Rayleigh test p-values are frequently small (often < 0.01) in both raw and detrended series for several seeds, suggesting non-uniformity in the relevant phase distribution under the specific pipeline settings;

(2) the primary IMF2↔λ phase-locking (PLV permutation p-values) is not consistently significant across seeds, and significance can move between raw and detrended channels depending on seed and TRIALS depth.

Several headline metrics are effectively constant across seeds (e.g., Wald p ≈ 0.264, Running-L best p ≈ 3.576×10⁻⁴, Running-λ spin-null global p ≈ 0.357).

This indicates these outputs are either (a) deterministic given the underlying data and not materially affected by the Monte Carlo components, or

(b) being read from a fixed section of the log/summary.

The current TRIALS=200 seed generator run (40 seeds, odd integers 3…81) is in progress; early results (seeds 3–15) show Rayleigh p-values remain frequently small, while PLV p-values remain variable and often non-significant.

Headline interpretation

At present, the evidence for a stable, seed-invariant ‘locking’ signature tied to λ is mixed: the phase-uniformity (Rayleigh) signal is more repeatable than the PLV-based locking against λ.

Accordingly, the near-term objective being investigated… (i) separate deterministic effects from Monte Carlo variability, (ii) establish which metrics survive rotation/null tests, and (iii) promote any promising signals to an explicitly pre-registered hypothesis test to control multiple-comparisons risk.

1) Metrics currently being applied

1.1 CEEMDAN / IMF extraction (signal decomposition)

What it is:Decomposition of the primary time series (or ordered series) into a set of Intrinsic Mode Functions (IMFs) using CEEMDAN (Complete Ensemble Empirical Mode Decomposition with Adaptive Noise). This is a data-driven, non-parametric method designed to separate oscillatory components without assuming stationarity.

Why it matters:

It is well suited to nonlinear / nonstationary structure (e.g., “chirp-like” frequency drifts).

It supports phase-based tests (locking) on specific IMFs (commonly IMF2/IMF3 in your pipeline).

Key outputs used downstream:

IMF2 (and sometimes IMF3) time series

Instantaneous phase and frequency (via Hilbert transform after IMF extraction)

1.2 PLV (Phase-Locking Value) tests

What it is:PLV quantifies how consistently the phase difference between two signals stays constant over the series. In the current analysis, the most notable results include:

PLV(IMF2 vs λ): does the IMF2 phase align with a periodic “calendar” variable (often λ, interpreted as a solar/ecliptic phase proxy depending on the build of the calendar)?

PLV phase stability (raw↔det): whether the phase structure persists after detrending / preprocessing.

Interpretation:

PLV near 0: weak or no consistent phase relation.

Higher PLV: stronger phase consistency, but significance must come from null tests (permutation, circular shift, sky-rotation, block resampling).

Associated p-values you report:

PLV perm p (raw/det): permutation null where phase pairing is randomized.

PLV circ-shift p: circularly shifting one series (preserves autocorrelation but breaks alignment).

PLV block-λ p: blockwise re-sampling / block-shuffle null (more conservative when there is long memory).

1.3 Rayleigh test (circular uniformity)

What it is:A classical test for whether phases are uniformly distributed on . Currently applied  to phase angles derived from CEEMDAN IMFs and/or a phase residual.

Interpretation:

Small p-value: phases are not uniform; indicates clustering around a preferred direction on the circle.

Important caution: Rayleigh significance can be inflated by autocorrelation or unmodelled structure; it should be interpreted alongside more robust nulls… (IAAFT/block/AR1/spin tests).

1.4 Energy / “E” metrics (E2 and E(IMF2+IMF3))

Current running tests… a set of “energy-like” scalar summaries:

E2 surrogate test: compares observed E2 to a surrogate distribution; in your logs it often returns p=1 (observed smaller than typical surrogates), suggesting no excess structure under that definition of E2.

E(imf2+imf3) with IAAFT surrogates: frequently returns low p (~0.0012) in your prints because N_IAAFT is 800 and you are hitting the resolution floor. This indicates “observed looks more extreme than most IAAFT surrogates” under that particular statistic.

E(imf2+imf3) with BLOCK / AR1+block surrogates: typically returns very high p (~0.99), meaning those conservative nulls consider the observed statistic unexceptional.

How to interpret the apparent contradiction:

IAAFT preserves the amplitude distribution and approximately preserves the power spectrum, but it can still allow structures that block-based nulls treat differently.

Block and AR1+block nulls can be much more conservative for slow drift / long-memory processes.

When IAAFT shows “significance” but block/AR1 + block does not, the safest reading is: your statistic is sensitive to low-frequency structure or dependence that the conservative nulls treat as typical.

1.5 Spin tests (sky-rotation null)

What it is:A null test where the sky field (or mapping between sky coordinate and the ordered index) is rotated / scrambled in a way that preserves some geometric features but breaks the hypothesised physical alignment.

In current investigation output results:

Outputting spin p, plus median/q95 of the null distribution.

In the runs analysed to date, spin p is generally not small (e.g., ~0.25), which argues against a strong sky-locked effect for that metric.

1.6 Extended Regression / Dipole fits

Current Operating fit…  linear model:

Where λ, β are coordinate proxies (often ecliptic longitude/latitude or similar depending on your calendar definition), and RA is right ascension.

dipole direction (λ°) with uncertainty and amplitude

Wald p (joint significance of the harmonic terms, or a subset)

In your logs, Wald p is frequently not small (~0.26), suggesting the global harmonic regression is not strongly significant under that specification.

1.7 Running-L / scanning statistic

Currently running-window amplitude and p-value (“Running-L best L_center … p=…”) and also compare against a spin null.

In current runs, runL_best_p is consistently ~3.6e-4 across seeds, which is a red flag for a non-seed-dependent deterministic scan feature. Often indicates:

the scan is operating on a fixed input (or fixed ordering) and the “seed” is not perturbing the relevant stage, or

the scan is using a statistic dominated by large-scale structure that is stable across CEEMDAN seeds.

This is an area being audited explicitly: aim being to confirm which stages are truly stochastic vs deterministic.

2) DR16 dataset description (what it is and what it contains)

“DR16” usually refers to SDSS Data Release 16, which includes the eBOSS final spectroscopic sample (and associated pipelines/catalogues). In your context, when you say “DR16 dataset” you appear to be referring to large quasar/galaxy redshift catalogues used for anisotropy/dipole analyses and/or BAO/LSS work.

Typical DR16 spectroscopic catalogue content includes (varies by specific table/product):

Sky position: Right Ascension (RA), Declination (Dec)

Redshift: (with measurement uncertainties and flags)

Object class / target class: QSO, LRG, ELG, etc.

Photometry / magnitudes: multi-band magnitudes and errors (depending on product)

Quality flags: pipeline flags, redshift confidence, selection bits

Derived quantities: absolute magnitudes, K-corrections (in some value-added products), completeness weights, systematics weights (in clustering catalogues)

If you are using the eBOSS DR16 clustering catalogues, they additionally include:

Weights: FKP weights, systematics weights, redshift failure weights

Random catalogues: for selection function / window corrections

Survey mask / completeness fields

If you are using Quaia (you referenced Quaia previously) that is a separate quasar astrometric catalogue derived from Gaia crossmatches; it is not “DR16” per se, but it can be used alongside DR16/eBOSS.

Best practice for the report: explicitly name the exact file/product you are using (e.g., “eBOSS DR16 QSO clustering catalogue” vs “SDSS DR16 specObj QSO table” vs “Quaia_G20.fits”), because “DR16” alone is ambiguous.

3) Programs used in this investigation (what each does)

HHT_Surrogates_and_Locking_Compare_V4.py

Primary analysis engine. Responsibilities include:

Load raw and detrended residual series (e.g., sn_residuals_enriched.csv, sn_residuals_detrended.csv)

Build the “calendar” mapping (solar mode in your runs)

Perform CEEMDAN decomposition and compute IMF phases

Compute PLV metrics and null p-values (perm, shift, block, spin)

Compute Rayleigh p-values

Compute energy statistics vs multiple surrogate families (IAAFT, block, AR1+block)

Run regression/dipole fits and running-window scans

Write summary CSV (hht_seasonal_robustness_summary_v4.csv) with schema versioning

run_hht_seedgrid_v2.py / v3.py

Orchestration / batch driver. Responsibilities include:

Set env vars (TRIALS, N_PERM, etc.)

Run the analysis script across a seed list (or generated seeds)

Capture stdout/stderr to per-seed logs

Robustly parse metrics from CSV tail and/or from log text (to avoid pandas ParserError)

Write an aggregate CSV for quick comparisons across seeds

run_hht_seedgrid_redo_NA_v1.py (and schema-fix snippets)

Utility scripts to:

Re-run failed seeds or “NA” results

Repair mixed-schema summary files (schema3 vs schema4) and preserve provenance

4) What has been done to remove local galaxy velocities (and what remains)

This is a critical point, because “local velocities” can mean two distinct things:

4.1 Observer motion (Solar System / Local Group) — the kinematic dipole

What it is:Our motion with respect to the CMB rest frame induces a dipole anisotropy and can bias redshift distributions and angular number counts.

What your current pipeline appears to do:

It includes harmonic regression terms (sin/cos of angles) and reports dipole directions.

It performs “spin” nulls and running-window scans.

What is not clearly evidenced in current output results:

An explicit correction of observed redshifts for our peculiar velocity relative to the CMB (i.e., transforming to the CMB frame)

Or explicit aberration/boost corrections when constructing sky maps

Next Stage Working Recommendation (for DR16 work):

Transform observed redshifts and/or inferred velocities into a consistent rest frame (CMB frame or Local Group frame) before building anisotropy fields.

Apply consistent masks and completeness weights.

Compare dipole direction to the known CMB dipole direction as a sanity check.

4.2 Galaxy peculiar velocities (non-Hubble flows)

What it is:Peculiar velocities are significant at low redshift; they broaden redshift-space distortions and can imprint correlations at large angular scales.

What’s already done (based on current working output logs):

Detrending: loaded sn_residuals_detrended.csv (so the SN analysis is explicitly attempting to remove smooth trends/structure).

Conservative nulls: block and AR1+block surrogates treat long-memory drift as “typical”, reducing false positives from correlated structures.

What needs to be added when DR16 is involved:

Redshift cuts: e.g., exclude very low-z objects where peculiar velocities dominate (common practice: z > 0.1 or z > 0.2 depending on tracer).

RSD-aware modelling: if using clustering catalogues, apply weights and treat line-of-sight effects carefully.

Systematics regression: regress out known observing-condition systematics (seeing, depth, stellar density, dust) when building sky fields.

4.3 For SN Ia specifically

SN cosmology analyses usually already incorporate:

Host-galaxy velocity dispersion terms and peculiar-velocity covariance contributions (often as an additional error term at low z)

Standard corrections (light-curve standardization)

Optional CMB-frame redshift corrections for the heliocentric → CMB transformation

The current HHT workflow is working off residual series rather than raw distance moduli, so you should document:

Which upstream pipeline produced sn_residuals_enriched.csv and what frame/corrections it used

What “detrended” means (polynomial? smoothing? model subtraction? windowing?)

History band summary

Datasets and inputs

Primary inputs used by the current HHT script runs (as printed in logs):

sn_residuals_enriched.csv (raw residual series)

sn_residuals_detrended.csv (detrended residual series)

These files are treated as time/sequence-ordered residuals produced upstream from the FR/log1pz pipeline. They are the only inputs explicitly confirmed in the HHT run logs included in this chat.

DR16 (eBOSS/SDSS) and CMB data have been discussed as next-step integrations. At present, these are not shown as direct inputs to the HHT_Surrogates_and_Locking_Compare_V4 pipeline outputs captured here.

Programs and scripts in use

Core analysis script:

HHT_Surrogates_and_Locking_Compare_V4.py (CEEMDAN decomposition; PLV/Rayleigh/rotation/null tests; surrogate tests; regression diagnostics; summary CSV writing).

Orchestration / robustness runners:

run_hht_seedgrid_v2.py (seed grid runner with tail-read parsing + log-text fallback).

run_hht_seedgrid_redo_NA_v1.py (rerun seeds where parsing or file writing produced NA).

run_hht_seedgrid_v3.py (seed generator run: N_SEEDS, SEED_START, SEED_STEP; TRIALS bump to 200; larger grid).

Data hygiene utilities:

schema-fix snippet to separate schema-3 (59 col) vs schema-4 (65 col) rows and restore a clean summary CSV for downstream parsing.

Method and test suite

The HHT_Surrogates_and_Locking_Compare pipeline applies CEEMDAN to extract intrinsic mode functions (IMFs) and then evaluates multiple hypotheses about phase structure and sky-coordinate dependence. Key tests include:

E2 surrogate test: Energy-based statistic for IMF2 (or related component), compared against surrogate distributions; reports E_obs and p.

PLV locking: Phase locking value between IMF2 phase and λ (ecliptic longitude) or analogous coordinate; permutation p-values and shift/block nulls.

Spin test on sky rotations: Rotation-based null to test whether the apparent locking survives sky coordinate rotations.

Rayleigh test: Circular uniformity test on phases; small p indicates departure from uniformity (phase clustering).

Extended regression: Residual ~ const + sinλ + cosλ + sinβ + z + sinRA + cosRA; dipole direction estimates and Wald test.

Running-L (windowed) regression: Scanning for best ‘L_center’ producing maximal amplitude; reports best center, amplitude, and p.

IAAFT / BLOCK / AR1+BLOCK surrogates (E(imf2+imf3)): Multiple null models to quantify whether IMF2+IMF3 energy/fraction is consistent with noise preserving different autocorrelation structures.

Results to date

Seed-grid results (TRIALS=50)

The table below reflects the TRIALS=50 seedgrid output captured in this chat. Note: for seed=47, the TRIALS=50 p-values were not printed in the captured snippet; fields are shown as NA.

TRIALS=50 – key p-values by seed

Seed-grid results (TRIALS=100)

When TRIALS was increased to 100, Rayleigh p-values generally decreased for several seeds, while PLV p-values continued to vary widely by seed and by channel (raw vs detrended).

TRIALS=100 – key p-values by seed

TRIALS=100 vs TRIALS=50 – comparison

Direct comparisons (same seed) show that increasing TRIALS can materially shift the reported p-values for both PLV and Rayleigh. This is expected because TRIALS affects Monte Carlo stability and the effective discretisation of p-values. However, the qualitative pattern remains: Rayleigh p-values are frequently small, whereas PLV-based ‘locking’ is not consistently significant.

Examples (same seed, different TRIALS):

Seed 3: PLV_p_raw: 0.286 → 0.0347 (more significant); Rayleigh_p_det: 0.00109 → 0.000832 (similar order).

Seed 7: PLV_p_raw: 0.00310 → 0.0166 (less significant); Rayleigh_p_raw: 0.000290 → 0.00176 (less extreme but still small).

Seed 11: PLV_p_det: 0.00245 → 0.0140 (less significant); Rayleigh_p_det: 0.000270 → 0.00116 (still small).

Current run in progress (TRIALS=200; 40 seeds)

The seed-generator run (run_hht_seedgrid_v3.py) is currently executing with TRIALS=200, N_SEEDS=40, SEED_START=3, SEED_STEP=2 (odd seeds 3…81). Interim results provided so far are listed below.

TRIALS=200 (interim, seeds 3–15) – key p-values

What the current results may imply

If the Rayleigh signal remains consistently small across the full TRIALS=200 grid, it would support the presence of systematic phase clustering under this pipeline (CEEMDAN settings + coordinate transform + residual construction). However, absent a correspondingly robust PLV locking signal that survives rotations/shifts/blocks, the most conservative interpretation is that the pipeline detects some structured phase behaviour, but its alignment with λ as a specific physical coordinate remains unproven.

The fact that certain outputs are constant across seeds (Wald p, Running-L best p/center/amp, and spin-null global p) is a priority diagnostic: these should be confirmed as intentionally deterministic metrics (computed once from the data), and, if so, relocated in reporting as ‘global’ rather than ‘seed-dependent’ results.

Supplementary figures (from pixel-shuffle MC)

The following plots were supplied in the chat. They illustrate p-values versus redshift for multiple Galactic latitude cuts. As presented, the curves show strong bin-to-bin variability, and any interpretation should be grounded in the exact definition of p(A) and p(|A_par|), the binning strategy, masking, and the null construction (pixel-shuffle Monte Carlo).

Figure 1. Pixel-shuffle MC: p(A) versus redshift (multiple |b| cuts).

Figure 2. Pixel-shuffle MC: p(|A_par|) versus redshift (multiple |b| cuts).

Tests in progress and recommended next targets

In progress

Complete TRIALS=200 seedgrid (40 seeds) and generate aggregate statistics (fraction significant at α=0.05; stability by channel).

Post-run: compute false discovery rate (FDR) adjusted q-values for PLV and Rayleigh across the grid.

Verify deterministic metrics: confirm why Wald/Running-L/spin-global p are constant; ensure parsing is reading the correct lines.

Recommended next targets

Robustness across analysis settings (not only RNG seed):

Sweep BLOCK_LEN (e.g., 16/32/64) and CEEMDAN parameters (noise width, ensemble size) while holding seed fixed.

Confirm results are stable to coordinate choice (ecliptic vs galactic vs equatorial), and to masking strategy.

Null-model strengthening:

Increase N_PERM and N_SHIFT on a subset of seeds to ensure p-value resolution is not limiting inference (e.g., N_PERM≥100k).

Add ‘phase randomisation’ surrogates matched to the power spectrum as an independent null check.

Promote promising signals to pre-registered tests:

If Rayleigh remains consistently significant, define one primary endpoint (e.g., Rayleigh p on detrended IMF2 phase) and one secondary endpoint (PLV vs λ), then hold analysis settings fixed for confirmatory runs.

DR16 × CMB integration (planned): cross-correlation approach

Direct subtraction of DR16 and CMB maps is generally not physically meaningful because they are different observables at different epochs. A more defensible approach is to construct comparable sky fields (e.g., a DR16 dipole residual map and a CMB temperature map), apply consistent masks, and evaluate the cross-correlation (or harmonic-space cross-power) against null simulations.

Implementation outline (code-level tasks):

Build HEALPix maps from DR16 quantities (e.g., number-count residuals or redshift residuals per pixel).

Load a CMB map (e.g., Planck SMICA) at compatible NSIDE, downgrade/upgrade as required.

Construct joint mask (Galactic plane + survey footprint + point-source masks).

Compute harmonic coefficients (a_lm) and cross-spectrum C_l^{DR16×CMB}; also compute real-space correlation w(θ).

Generate null distributions via: (i) DR16 pixel-shuffles, (ii) phase randomisation, and (iii) map rotations.

Where a signal is suggested, localise by multipole ranges (dipole-only; low-l; etc.) and by redshift slices.

## Extracted Tables

### Table 1

Date | Change / action | Outcome / note
2026-01-11 | Seed-grid runner initially failed with pandas ParserError due to mixed-schema CSV (59 vs 65 columns). | Implemented robust tail-read parsing and schema-fix workflow; re-runs possible.
2026-01-11 | TRIALS=50 seedgrid (10 seeds) executed; some NA fields due to parsing/resume logic. | Observed variable PLV significance; Rayleigh often small; deterministic-looking Running-L/Wald outputs.
2026-01-11 | Mixed-schema summary CSV repaired (schema-4 fixed; schema-3 extracted for provenance). | Seedgrid aggregation restored; NA mostly eliminated.
2026-01-12/13 | TRIALS increased to 100; seedgrid rerun (10 seeds). | Rayleigh p-values generally reduced vs TRIALS=50; PLV significance shifted among seeds; still non-robust overall.
2026-01-13 | Seedgrid v3 started (TRIALS=200; N_SEEDS=40; generator seeds 3..81 step 2). | In progress; early seeds maintain small Rayleigh p-values; PLV p-values remain variable.

### Table 2

Seed | PLV_p_raw | PLV_p_det | Rayleigh_p_raw | Rayleigh_p_det | wald_p_det | runL_best_p | runningL_spin_global_p_det
3 | 0.285536 | 0.0119994 | 0.0478379 | 0.0010938 | 0.264339 | 0.0003576 | NA
7 | 0.00309985 | 0.0982951 | 0.000290059 | 0.0142065 | 0.264339 | 0.0003576 | NA
11 | 0.519074 | 0.00244988 | 0.109719 | 0.000269839 | 0.264339 | 0.0003576 | 0.356822
19 | 0.116644 | 0.828259 | 0.0175462 | 0.290568 | 0.264339 | 0.0003576 | 0.356822
23 | 0.0105495 | 0.00419979 | 0.000765376 | 0.000278613 | 0.264339 | 0.0003576 | 0.356822
29 | 0.020599 | 0.0515974 | 0.00192495 | 0.0100894 | 0.264339 | 0.0003576 | 0.356822
31 | 0.00359982 | 0.0039498 | 0.000370325 | 0.000480415 | 0.264339 | 0.0003576 | 0.356822
37 | 0.528624 | 0.0476976 | 0.107546 | 0.00733468 | 0.264339 | 0.0003576 | 0.356822
41 | 0.151042 | 0.112444 | 0.0181682 | 0.0125939 | 0.264339 | 0.0003576 | 0.356822
47 | NA | NA | NA | NA | 0.264339 | 0.0003576 | 0.356822

### Table 3

seed | PLV_p_raw | PLV_p_det | Rayleigh_p_raw | Rayleigh_p_det | wald_p_det | runL_best_p | runningL_spin_global_p_det
3 | 0.0347483 | 0.00779961 | 0.00395593 | 0.000832433 | 0.264339 | 0.0003576 | 0.356822
7 | 0.0166492 | 0.873256 | 0.0017627 | 0.279458 | 0.264339 | 0.0003576 | 0.356822
11 | 0.109195 | 0.0140493 | 0.0175311 | 0.00116166 | 0.264339 | 0.0003576 | 0.356822
19 | 0.080246 | 0.169192 | 0.0108615 | 0.0240435 | 0.264339 | 0.0003576 | 0.356822
23 | 0.0924954 | 0.0272986 | 0.0124545 | 0.00287748 | 0.264339 | 0.0003576 | 0.356822
29 | 0.247088 | 0.0849958 | 0.0441774 | 0.00949909 | 0.264339 | 0.0003576 | 0.356822
31 | 0.736613 | 0.0437978 | 0.158204 | 0.00534982 | 0.264339 | 0.0003576 | 0.356822
37 | 0.464427 | 0.0228489 | 0.0823575 | 0.00230467 | 0.264339 | 0.0003576 | 0.356822
41 | 0.0495475 | 0.0386481 | 0.00594272 | 0.00626715 | 0.264339 | 0.0003576 | 0.356822
47 | 0.0542473 | 0.273686 | 0.0065827 | 0.0451994 | 0.264339 | 0.0003576 | 0.356822

### Table 4

seed | PLV_p_raw | PLV_p_det | Rayleigh_p_raw | Rayleigh_p_det
3 | 0.0368982 | 0.0837958 | 0.00454048 | 0.0110751
5 | 0.153492 | 0.0443978 | 0.0203636 | 0.00623008
7 | 0.0896455 | 0.0574471 | 0.012129 | 0.0072404
9 | 0.0719964 | 0.433778 | 0.0092839 | 0.0949933
11 | 0.040698 | 0.154792 | 0.00429267 | 0.0257259
13 | 0.0706965 | 0.125744 | 0.010038 | 0.0149429
15 | 0.150842 | 0.0626469 | 0.0255274 | 0.0084387
