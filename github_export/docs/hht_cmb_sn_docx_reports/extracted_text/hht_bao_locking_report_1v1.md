# HHT_BAO_Locking_Report_1V1

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\HHT_BAO_Locking_Report_1V1.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2025-09-27T21:21:00Z`
- DOCX modified: `2025-09-27T21:21:00Z`

## Extracted Text

Report 1V0: Report: HHT–BAO Locking Analysis

Executive Summary

We investigated Hilbert–Huang Transform (HHT) decomposition of supernova (SN) Ia residuals to test for possible scale-locking with baryon acoustic oscillation (BAO) effective scales.

By extracting intrinsic mode functions (IMFs) and comparing their instantaneous-frequency ridges with BAO-derived wavenumbers, we tested baseline,

Noether-constrained, and Λ(z)-coupled models.

Results demonstrate a robust IMF 2 correlation across all models and smoothing choices, and an emergent IMF 3 lock only under Λ(z) coupling.

Background

This investigation is embedded in the PLamb cosmological framework, which explores variable-mass and variable-c dynamics relative to standard Friedmann–Robertson–Walker (FRW) cosmology.

The analysis leverages supernova (Pantheon+SH0ES) data and compares to BAO distance scales.

Models tested:

- FRPBH baseline (variable mass, variable c, PBH-driven acceleration).

- Noether symmetry (conserve_time: freezes accretion terms).

- Λ(z) coupled OU process with gamma coupling (fluctuation-induced).

Acronyms

• HHT – Hilbert–Huang Transform.

• IMF – Intrinsic Mode Function (oscillatory component from EMD/EEMD).

• BAO – Baryon Acoustic Oscillation.

• k_eff – effective BAO wavenumber from DM/rd measurements.

• SN – Supernova Ia distance modulus data (Pantheon+SH0ES).

Key Equations

HHT decomposition: Residual(z) = Σ IMF_i(z) + trend(z).

IMF energy:  E_i = Σ IMF_i(z)^2.

BAO effective scale: k_eff(z) = 2π / (D_M(z)/r_d).

Method

SN residuals were detrended and placed on a uniform log(1+z) grid with tapering. Empirical Mode Decomposition (EMD/EEMD) extracted IMFs, and Hilbert transforms provided instantaneous frequencies.

IMF ridges were compared to BAO scales using correlation (Pearson, Spearman) and permutation-based significance testing.

FDR correction identified robust IMFs.

Smoothings of IMFs 5, 7, 9, 11 were tested to confirm stability.

Results

Findings by model:

• Baseline FRPBH: IMF 2 significant across all smoothings (r≈0.14, p≈5e-4).

• Noether conserve_time: IMF 2 strengthened (r≈0.22, p≈5e-4).

• Λ(z) coupling: IMF 2 weaker but still significant (r≈0.09, p≈8e-3).

– Additionally, IMF 3 becomes robustly significant (r≈0.14, p≈5e-4).

• Higher IMFs (4–6) showed no consistent correlations (p>0.1).

• Robustness analysis confirmed IMF 2’s persistence across all smoothing kernels.

Interpretation

IMF 2 represents a universal locking to BAO transverse scale.

The symmetry constraint (conserve_time) suppresses extraneous structure and enhances the IMF 2 lock.

Λ(z) coupling shifts coherence to IMF 3, suggesting an additional frequency band induced by fluctuations.

Conclusions and Next Steps

Evidence indicates HHT-extracted IMFs track BAO scales.

IMF 2 is a stable universal feature; IMF 3 emerges under Λ(z) coupling.

Next steps include mapping IMF median frequencies to physical lengths, testing radial (DH/rd) versus transverse (DM/rd) BAO dependence, and jackknife resampling of SN data to confirm robustness.

HHT–BAO Locking Report (Revision 1V1)

New Results: Noether Strength Scan

We extended the BAO–IMF locking analysis to include varying soft Noether constraint strengths (s = 0, 2, 4, 6, 8, 10). This tested whether increasing the degree of symmetry breaking/smoothing systematically changes the correlation structure.

Key findings:

IMF 2 · DM/rd: max |r|=0.122 at s=2 (p=0.00122,  trend ↑)

IMF 2 · DH/rd: max |r|=-0.114 at s=2 (p=0.00252,  trend ↓)

IMF 3 · DM/rd: max |r|=-0.253 at s=8 (p=1.1e-11,   trend ↓)

IMF 3 · DH/rd: max |r|=0.246 at s=8 (p=4.77e-11, trend ↑)

Interpretation

• IMF 2 continues to capture a robust BAO lock, with the strongest signal around low to moderate soft Noether strength (s=2).

This suggests the BAO-scale correlation is most stable under partial conservation constraints.• IMF 3 emerges as a strong harmonic under higher soft strength (s=8), suggesting a secondary fluctuation mode consistent with Λ(z)-like variability.• DM (transverse) and DH (radial) show opposing correlation directions in IMF 2 and IMF 3, confirming a directional split in how BAO locking manifests.

Pros and Cons of HHT Analysis

Advantages:

Adaptive, data-driven decomposition without requiring basis functions.

Captures local oscillatory modes (IMFs) aligned with cosmological scales.

Provides instantaneous frequency analysis that can be directly compared with BAO scales.

Reveals hidden coupling signals (IMF2 with BAO, IMF3 under Λ(z) fluctuations).

Limitations:

Mode mixing: IMF content may split or merge depending on data length/smoothing.

Statistical significance requires heavy permutation testing (computationally intensive).

Interpretation of higher-order IMFs can be ambiguous.

Sensitivity to boundary conditions and redshift coverage gaps.

Next Steps

Extend frequency–BAO regression to include other fluctuation models (OU, piecewise, spline Λ(z)).

Cross-check IMF2 and IMF3 against Planck distance priors for consistency.

Incorporate directionality tests into model selection (does IMF2 prefer DM while IMF3 prefers DH?).

Automate grid scans over Noether and Λ(z) fluctuation parameters, generating summary plots.
