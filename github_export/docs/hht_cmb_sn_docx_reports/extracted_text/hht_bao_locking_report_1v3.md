# HHT_BAO_Locking_Report_1V3

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\HHT_BAO_Locking_Report_1V3.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2025-09-29T15:53:00Z`

## Extracted Text

HHT–BAO Locking Report

(Revision 1V2)

This document compiles the latest results from the HHT IMF–BAO analysis, including Noether strength sweeps, redshift-window tests, and harmonic ratio diagnostics. Version ID: 1V2.

1. Noether Grid (IMF2 & IMF3)

The following CSV summaries and plots quantify how IMF2 and IMF3 locking varies under different Noether combo_soft strengths.

Key files (saved under plamb_runs/noether_grid/summary/):- dm_dh_allfits.csv- IMF2_IMF3_vs_strength.csv and IMF2_IMF3_vs_strength_ext.csv- direction_scores.csv- imf2_r_vs_strength.png, imf3_r_vs_strength.png- imf2_delta_r_vs_strength.png, imf3_delta_r_vs_strength.png- notes.txt

2. Harmonic Ratio Test (IMF3/IMF2)

Harmonic relationships between IMF2 and IMF3 were probed. Results saved under plamb_runs/harmonic_test/.

Files:- harmonic_test.csv- ratio_summary.png- baseline_ratio_hist.png- conserve_time_ratio_hist.png- lambda_on_ratio_hist.png- lambda_grid_ratio_hist.png

Figure: Harmonic ratio summary (IMF3/IMF2 across packs).

3. Redshift Window Sensitivity

We re-ran HHT+BAO locking with restricted redshift ranges. Summary file: plamb_runs/hht_sn_FRPBH_windows/IMF2_IMF3_by_window.csv

Example results:- z=0.02–1.9: IMF2 r=0.058 (p=0.061), IMF3 r=-0.071 (p≈0.98)- z=0.05–1.5: IMF2 r=-0.051 (p≈0.91), IMF3 r=-0.016 (p≈0.68)

4. Interpretation

• IMF2 locking is strongest at low-to-moderate Noether strengths (s=2–4), but breaks down for s=6 and re-emerges intermittently.• IMF3 shows alternating correlation/anti-correlation with BAO, strongest at s≈8.• Harmonic ratio analysis shows lambda_on produces a distribution close to the 0.5–1 harmonic regime, whereas baseline and conserve_time show wide scatter.• Redshift-restricted runs indicate instability in IMF2 locking, possibly suggesting window-dependent coupling.

5. HHT Analysis: Pros & Cons

Pros:• Captures nonlinear, non-stationary oscillatory structure in SN residuals.• Allows IMF-specific frequency–phase analysis with BAO cross-comparison.• Provides harmonic relationship testing between IMFs.Cons:• Sensitive to decomposition choices (smoothing, IMF count).• Surrogate-based significance tests can vary with permutations.• Results can shift under different redshift windows or priors, indicating instability.

6. Next Steps

1. Expand harmonic test to raw IMF instantaneous frequencies (not synthetic Gaussian surrogates).2. Systematically scan smoothing levels with higher permutation counts.3. Extend redshift-window analysis to overlapping bins for robustness.4. Integrate Planck/BAO priors into Noether-strength grids to test cosmological stability.5. Compare with alternative decomposition methods (e.g. wavelets) to validate HHT robustness.

Revision 1V2 Addendum

This addendum to the 1V1 report incorporates the findings from Revision 1V2, expanding the analysis to include Noether grid suppression tests, harmonic ratio checks, and redshift window sensitivity. These results build upon the initial IMF–BAO locking evidence and provide further robustness checks across different modelling assumptions.

1. Noether Grid (Extended Strengths)

A systematic sweep of Noether combo_soft strengths (s = 0–12) revealed:- IMF2: lock signal generally stable at low-to-moderate strengths (s = 0–4), but suppression emerges at s = 6 and s = 10, consistent with decoherence effects.- IMF3: shows alternating sign with significant correlations at s = 6 and strong negative lock at s = 8.This suggests that IMF2 is more robust to symmetry-breaking perturbations, while IMF3 is more sensitive to oscillatory artifacts.

2. Harmonic Test (IMF3/IMF2 Ratios)

The harmonic test assessed whether IMF3 represents a harmonic multiple of IMF2. Results:- Baseline: ratio_mean ≈ –0.52 (large variance).- Conserve_time: ≈ 0.06 (low variance, near independence).- Lambda_on: ≈ 0.50 (close to half-harmonic relation).- Lambda_grid: ≈ 0.96 (near independent scaling).Interpretation: IMF3 shows harmonic coupling under some conditions (notably lambda_on) but generally behaves as an independent mode. This reinforces IMF2 as the primary BAO-lock carrier.

3. Redshift Window Sensitivity

Recomputations of HHT across restricted redshift ranges produced:- Full z = 0.02–1.9: IMF2 weak correlation (r ≈ 0.058, p ≈ 0.06), IMF3 non-significant.- z = 0.05–1.5: IMF2 no correlation (r ≈ –0.051), IMF3 also non-significant.- z = 0.10–1.2: analysis ongoing, but trend suggests weaker locking at narrower ranges.This highlights that IMF–BAO locking signatures are stronger when retaining the full dataset, and diminish under restricted redshift ranges.

4. Pros & Cons of HHT Analysis (Extended)

Pros:- Captures adaptive frequency bands, well-suited for quasi-periodic cosmological signals.- Enables IMF-level comparison with BAO scales and harmonic tests.- Flexible: can test sensitivity across priors, windows, and smoothness.Cons:- Sensitive to choice of smoothing and permutation parameters.- Higher IMFs may introduce artifacts (sign-flipping, spurious lock).- Statistical interpretation remains challenging given correlated residual structure.Overall, IMF2 remains the most stable candidate for cosmological locking to BAO scales.

Next Steps

- Extend harmonic ratio tests with raw IMF spectra instead of synthetic distributions.- Incorporate Planck priors and BAO covariance into the HHT decomposition to test robustness.- Cross-check IMF locking against other decompositions (wavelets, multitaper spectra).- Formalize AIC/BIC ranking of models including HHT-derived parameters.
