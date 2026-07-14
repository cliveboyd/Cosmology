# HHT_BAO_Locking_Report_1V4

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\HHT_BAO_Locking_Report_1V4.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2025-09-29T16:03:00Z`

## Extracted Text

HHT–BAO Locking Analysis Report

Version: 1V4

Date: September 2025

Executive Summary

This update builds on versions 1V0–1V3 of the HHT–BAO locking analysis. The most recent analyses extend the Noether soft-breaking grid, introduce harmonic ratio tests, and explores redshift-window dependencies. Key results show that IMF2 remains the primary carrier of BAO locking at low Noether strengths, but decoheres at higher values. IMF3, in contrast, demonstrates conditional behavior: independent in most contexts, but harmonic (1:2) under specific λ-driven scenarios.

Results Update

1. Extended Noether Grid (s = 0–12)

- IMF2: Stable locking at s = 0–4; decoheres at s = 6 and 10; partial recovery at s = 8.- IMF3: Weak baseline signal; flips correlation sign at s = 6–8; acts independent at higher strengths.Reference files:- noether_grid_IMF2_IMF3_vs_strength.csv- noether_grid_summary_direction_scores.csv

2. Harmonic Ratio Tests (IMF3 / IMF2)

- λ_on: Mean ratio ≈ 0.5 → harmonic locking.- λ_grid: Mean ratio ≈ 1.0 → independent oscillators.- Baseline: Scatter, wide variance around −0.5.- Conserve_time: ~0.06, decoupled.Reference files:- harmonic_test_harmonic_test.csv- harmonic_test_ratio_summary.png

3. Redshift Window Dependencies

- Full window (z = 0.02–1.9): IMF2 modest lock (p ≈ 0.06), IMF3 insignificant.- Restricted (z = 0.05–1.5): IMF2 weak (p ≈ 0.91), IMF3 non-significant.- Window (z = 0.10–1.2): signals nearly vanish.Interpretation: Locking strongest only when full redshift range is retained.Reference file:- hht_sn_FRPBH_windows_IMF2_IMF3_by_window.csv

4. DM vs DH Fits

- IMF2: Consistent DM/DH correlations at low strengths.- IMF3: Opposing DM vs DH tendencies.- Δr highlights systematic IMF3 tension.Reference files:- bao_compare_dm_dh_freq_scale_fits.csv- bao_compare_dm_dh_summary.md- noether_grid_dm_dh_soft0_freq_scale_fits.csv- noether_grid_dm_dh_soft10_freq_scale_fits.csv

Interpretation

- IMF2 is the most stable BAO-lock carrier but fragile under Noether perturbations or redshift cuts.- IMF3 is conditionally harmonic but often independent or oppositional.- Harmonic ratio results reinforce distinct dynamical regimes driven by Λ(z) or coupling assumptions.

Pros and Cons of HHT Approach (Update)

Pros:- Sensitive to quasi-periodic microstructures in SN residuals.- Detects conditional harmonic relations (IMF3/IMF2 ≈ 0.5).- Flexible across redshift windows and parameter sweeps.Cons:- Locking signals fragile under data cuts.- Dependent on smoothing and permutation count.- IMF3 often noisy, challenging interpretation.

Next Steps

1. Extend harmonic ratio analysis using raw instantaneous frequencies.2. Introduce Planck priors and BAO covariance into HHT runs.3. Apply wavelet decomposition as cross-check.4. Rank FR, PBH, and ΛCDM models by AIC/BIC with HHT signals.5. Investigate IMF3’s oppositional DM vs DH behavior.

Reference Dataset: Final Curated Drop

Curated results: /Users/boyde/Desktop/HHT_BAO_FINAL_DROP_CLEAN/Key files:- hht_bao_lock_real_baseline_s9_bao_lock_heatmap_table_with_p.csv- noether_grid_IMF2_IMF3_vs_strength.csv- noether_grid_summary_direction_scores.csv- bao_compare_dm_dh_freq_scale_fits.csv- bao_compare_dm_dh_summary.md- harmonic_test_harmonic_test.csv- harmonic_test_ratio_summary.png- hht_sn_FRPBH_windows_IMF2_IMF3_by_window.csv
