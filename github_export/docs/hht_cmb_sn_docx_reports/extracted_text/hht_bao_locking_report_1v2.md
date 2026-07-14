# HHT_BAO_Locking_Report_1V2

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\HHT_BAO_Locking_Report_1V2.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2013-12-23T23:15:00Z`

## Extracted Text

HHT–BAO Locking Report (Revision 1V2)

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
