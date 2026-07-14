# HHT_Status_Update_2026-01-13 (1)

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\HHT_Status_Update_2026-01-13 (1).docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2013-12-23T23:15:00Z`

## Extracted Text

HHT / Seasonal Robustness Investigation – Status Update

Update date: 2026-01-15Scope: CEEMDAN/HHT locking and robustness tests on SN residual fields; seed-grid stability; DR16×CMB cross-correlation workplan.

Prepared for internal research tracking. This document summarises what has been observed so far, what is statistically robust (and what is not), and what tests are currently in progress.

Executive summary

Across the runs captured to date, two patterns stand out: (1) Rayleigh test p-values are frequently small (often < 0.01) in both raw and detrended series for several seeds, suggesting non-uniformity in the relevant phase distribution under the specific pipeline settings; (2) the primary IMF2↔λ phase-locking (PLV permutation p-values) is not consistently significant across seeds, and significance can move between raw and detrended channels depending on seed and TRIALS depth.

Several headline metrics are effectively constant across seeds (e.g., Wald p ≈ 0.264, Running-L best p ≈ 3.576×10⁻⁴, Running-λ spin-null global p ≈ 0.357). This indicates these outputs are either (a) deterministic given the underlying data and not materially affected by the Monte Carlo components, or (b) being read from a fixed section of the log/summary.

The current TRIALS=200 seed generator run (40 seeds, odd integers 3…81) is in progress; early results (seeds 3–15) show Rayleigh p-values remain frequently small, while PLV p-values remain variable and often non-significant.

Headline interpretation

At present, the evidence for a stable, seed-invariant ‘locking’ signature tied to λ is mixed: the phase-uniformity (Rayleigh) signal is more repeatable than the PLV-based locking against λ. Accordingly, the near-term objective should be to (i) separate deterministic effects from Monte Carlo variability, (ii) establish which metrics survive rotation/null tests, and (iii) promote any promising signals to an explicitly pre-registered hypothesis test to control multiple-comparisons risk.

History band summary

Datasets and inputs

Primary inputs used by the current HHT script runs (as printed in logs):

• sn_residuals_enriched.csv (raw residual series)

• sn_residuals_detrended.csv (detrended residual series)

These files are treated as time/sequence-ordered residuals produced upstream from the FR/log1pz pipeline. They are the only inputs explicitly confirmed in the HHT run logs included in this chat.

DR16 (eBOSS/SDSS) and CMB data have been discussed as next-step integrations. At present, these are not shown as direct inputs to the HHT_Surrogates_and_Locking_Compare_V4 pipeline outputs captured here.

Programs and scripts in use

Core analysis script:

• HHT_Surrogates_and_Locking_Compare_V4.py (CEEMDAN decomposition; PLV/Rayleigh/rotation/null tests; surrogate tests; regression diagnostics; summary CSV writing).

Orchestration / robustness runners:

• run_hht_seedgrid_v2.py (seed grid runner with tail-read parsing + log-text fallback).

• run_hht_seedgrid_redo_NA_v1.py (rerun seeds where parsing or file writing produced NA).

• run_hht_seedgrid_v3.py (seed generator run: N_SEEDS, SEED_START, SEED_STEP; TRIALS bump to 200; larger grid).

Data hygiene utilities:

• schema-fix snippet to separate schema-3 (59 col) vs schema-4 (65 col) rows and restore a clean summary CSV for downstream parsing.

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

• Complete TRIALS=200 seedgrid (40 seeds) and generate aggregate statistics (fraction significant at α=0.05; stability by channel).

• Post-run: compute false discovery rate (FDR) adjusted q-values for PLV and Rayleigh across the grid.

• Verify deterministic metrics: confirm why Wald/Running-L/spin-global p are constant; ensure parsing is reading the correct lines.

Recommended next targets

1) Robustness across analysis settings (not only RNG seed):

• Sweep BLOCK_LEN (e.g., 16/32/64) and CEEMDAN parameters (noise width, ensemble size) while holding seed fixed.

• Confirm results are stable to coordinate choice (ecliptic vs galactic vs equatorial), and to masking strategy.

2) Null-model strengthening:

• Increase N_PERM and N_SHIFT on a subset of seeds to ensure p-value resolution is not limiting inference (e.g., N_PERM≥100k).

• Add ‘phase randomisation’ surrogates matched to the power spectrum as an independent null check.

3) Promote promising signals to pre-registered tests:

• If Rayleigh remains consistently significant, define one primary endpoint (e.g., Rayleigh p on detrended IMF2 phase) and one secondary endpoint (PLV vs λ), then hold analysis settings fixed for confirmatory runs.

DR16 × CMB integration (planned): cross-correlation approach

Direct subtraction of DR16 and CMB maps is generally not physically meaningful because they are different observables at different epochs. A more defensible approach is to construct comparable sky fields (e.g., a DR16 dipole residual map and a CMB temperature map), apply consistent masks, and evaluate the cross-correlation (or harmonic-space cross-power) against null simulations.

Implementation outline (code-level tasks):

• Build HEALPix maps from DR16 quantities (e.g., number-count residuals or redshift residuals per pixel).

• Load a CMB map (e.g., Planck SMICA) at compatible NSIDE, downgrade/upgrade as required.

• Construct joint mask (Galactic plane + survey footprint + point-source masks).

• Compute harmonic coefficients (a_lm) and cross-spectrum C_l^{DR16×CMB}; also compute real-space correlation w(θ).

• Generate null distributions via: (i) DR16 pixel-shuffles, (ii) phase randomisation, and (iii) map rotations.

• Where a signal is suggested, localise by multipole ranges (dipole-only; low-l; etc.) and by redshift slices.

## Extracted Tables

### Table 1

Date | Change / action | Outcome / note
2026-01-11 | Seed-grid runner initially failed with pandas ParserError due to mixed-schema CSV (59 vs 65 columns). | Implemented robust tail-read parsing and schema-fix workflow; re-runs possible.
2026-01-11 | TRIALS=50 seedgrid (10 seeds) executed; some NA fields due to parsing/resume logic. | Observed variable PLV significance; Rayleigh often small; deterministic-looking Running-L/Wald outputs.
2026-01-11 | Mixed-schema summary CSV repaired (schema-4 fixed; schema-3 extracted for provenance). | Seedgrid aggregation restored; NA mostly eliminated.
2026-01-12/13 | TRIALS increased to 100; seedgrid rerun (10 seeds). | Rayleigh p-values generally reduced vs TRIALS=50; PLV significance shifted among seeds; still non-robust overall.
2026-01-13 | Seedgrid v3 started (TRIALS=200; N_SEEDS=40; generator seeds 3..81 step 2). | In progress; early seeds maintain small Rayleigh p-values; PLV p-values remain variable.

### Table 2

seed | PLV_p_raw | PLV_p_det | Rayleigh_p_raw | Rayleigh_p_det | wald_p_det | runL_best_p | runningL_spin_global_p_det
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
