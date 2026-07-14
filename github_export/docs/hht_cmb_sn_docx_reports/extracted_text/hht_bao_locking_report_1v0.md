# HHT_BAO_Locking_Report_1V0

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\HHT_BAO_Locking_Report_1V0.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2025-09-27T06:11:00Z`

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
