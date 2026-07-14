# HHT_BAO_Cosmology_Overview

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\HHT_BAO_Cosmology_Overview.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2025-09-29T16:43:00Z`

## Extracted Text

HHT-Based Cosmological Analysis Overview

1. Introduction

This report provides an overview of investigations into the use of the Hilbert–Huang Transform (HHT) for analysing cosmological data, with a focus on supernova (SN) residuals, baryon acoustic oscillation (BAO) scales, and symmetry-breaking tests in the Full Relativity (FR) model.

Unlike traditional Fourier or wavelet methods, HHT offers a local, data-driven approach to extracting intrinsic mode functions (IMFs) and instantaneous frequencies, allowing the detection of non-stationary, quasi-periodic structures or scale couplings that might otherwise be hidden.

2. Background, Variables, and Acronyms

Key equations and definitions:

• SN residuals: Δμ(z) = μ_obs(z) − μ_model(z)

• BAO observables: DM/rd and DH/rd, where rd is the sound horizon at drag epoch.

• Pearson correlation used for locking tests: r = Cov(ω̄, 1/k_eff) / (σ_ω̄ σ_1/k).

• IMFs: Intrinsic Mode Functions derived from HHT decomposition. IMF2 and IMF3 represent the 2nd and 3rd oscillatory modes in SN residuals.

• Noether strength: Control parameter (0–12 tested) governing conservation symmetry constraints (combo_soft, conserve_time, etc.).

• Acronyms: SN = Supernova, BAO = Baryon Acoustic Oscillation, FR = Full Relativity model, PBH = Primordial Black Hole, CMB = Cosmic Microwave Background.

3. Methodology

The HHT method combines Empirical Mode Decomposition (EMD) with Hilbert transforms to extract instantaneous frequencies ω(z) from SN residuals as a function of log(1+z).

These are compared against effective BAO scales (1/k_eff) derived from clustering measurements.

Pearson correlation coefficients (r) are computed for IMF2 and IMF3 to test alignment with DM/rd and DH/rd. Significance is assessed using phase-randomized surrogate datasets.

4. Results to Date

4.1 Harmonic Ratios

Harmonic tests compared the IMF3/IMF2 ratio across multiple packs (baseline, conserve_time, lambda_on, lambda_grid).

Results show clear deviations from both the 0.5 harmonic ratio and the 1.0 independence line, with lambda_on and lambda_grid exhibiting stronger alignment towards harmonic coupling.

See file: harmonic_test_harmonic_test.csv and figure: harmonic_test_ratio_summary.png.

4.2 BAO Locking

BAO correlation tests demonstrated non-random alignment between IMF2/IMF3 instantaneous frequencies and BAO effective scales.

Strength sweeps (0–12) showed IMF3 becoming significantly aligned (p < 0.01) at higher Noether strengths (s≈8).

See: bao_compare_dm_dh_freq_scale_fits.csv and bao_compare_dm_dh_summary.md.

4.3 Directional Locking (DM vs DH)

Directionality tests computed Δr = r_DM − r_DH.

Results indicate IMF2 tends to favor DM/rd locking, while IMF3 locks more strongly with DH/rd at higher strengths.

See: noether_grid_IMF2_IMF3_vs_strength.csv.

4.4 Redshift Windows

When HHT analysis was restricted to specific redshift windows, IMF2 showed weak correlations in z=[0.05,1.5], while IMF3 exhibited more robust alignment over z=[0.02,1.9]. See: hht_sn_FRPBH_windows_IMF2_IMF3_by_window.csv.

5. Interpretation and Implications

These findings suggest that residual structure in SN datasets contains oscillatory modes which may couple to BAO scales in a way not predicted by ΛCDM.

If validated, this could indicate that dark energy phenomenology is not simply a cosmological constant, but linked to underlying dynamical processes such as symmetry breaking (tested via Noether strength), variable-c scaling in FR, or PBH entropy growth models.

The harmonic coupling observed between IMF2 and IMF3 provides additional evidence for non-random, physically meaningful structure in SN residuals, potentially supporting models with fractal or scale-locked evolution.

6. HHT vs Standard Cosmology Approaches

Advantages of HHT analysis:

• Local, adaptive decomposition that does not assume stationarity.

• Ability to detect scale locking and harmonic structure.

• Complements Fourier/ΛCDM residual analyses by exposing new patterns.

Limitations and caveats:

• Sensitive to noise and decomposition choices.

• Requires large surrogate testing for robust significance.

• Interpretation depends on mapping IMF frequencies to cosmological scales.

7. Next Steps

• Extend surrogate testing to >10,000 permutations for robust significance.

• Expand redshift-window analyses to isolate low-z vs high-z contributions.

• Cross-validate with eBOSS DR16 BAO datasets and JWST high-z supernovae.

• Integrate with FR + PBH cosmological model fits for joint constraints.

• Prepare publication-ready plots and results tables for peer review.
