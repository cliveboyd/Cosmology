# BAO–HHT Integration: 2-bin ρ(z) Update

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\BAO_HHT_2bin_update.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2025-10-10T05:05:00Z`

## Extracted Text

BAO–HHT Integration: 2‑bin ρ(z) Update

Generated: 2025-10-10 04:58

Overview

You tested a two‑bin BAO covariance ρ(z) with a split at z=0.7, ρ_lo=+0.374 (z<0.7) and ρ_hi=−0.117 (z≥0.7).

Relative to the constant‑ρ baseline (+0.174), the 2‑bin model reduces the BAO χ² by ≈10.4 (237.17 → 226.77) and the total χ² by ≈13.17, while keeping (H0, Ωm) nearly unchanged.

Model comparison (constant‑ρ vs 2‑bin ρ(z))

Note: For a fair comparison, AIC/BIC for the 2‑bin case include a +1 hyper‑parameter penalty for {ρ_lo, ρ_hi} compared to constant‑ρ. PLamb’s printed AIC/BIC count only (H0, Ωm).

Per‑z χ² contributions (2‑bin − constant‑ρ)

HHT interpretation & projection notes

Mapping residuals to log(1+z) for HHT: IMF “frequency” is oscillations per unit log(1+z); bandwidth reflects temporal coherence across cosmic time. IMF energy (amplitude²) vs log(1+z) localizes which redshift ranges drive χ².

In these fits, the two‑bin covariance reduces χ² mainly near z≈0.38 and z≈0.70, consistent with LOOCV diagnostics.

This points to a mild redshift dependence in the DM–DH correlation.

Key equations / definitions (flat ΛCDM)

E(z) = H(z)/H0 = √[ Ωm(1+z)^3 + (1−Ωm) ]

D_M(z) = (c/H0) ∫₀ᶻ dz′/E(z′)

D_H(z) = c/H(z) = c/(H0 E(z))

Predictions: D_M/rd and D_H/rd; χ² uses 2×2 DM–DH blocks per‑z.

Two‑bin covariance: ρ(z) = ρ_lo for z<z_split, ρ_hi for z≥z_split (off‑diagonal term).

Artifacts & paths

• Constant‑ρ HHT export: plamb_runs/tools/bao_hht_const_rho0174/

• 2‑bin HHT export (long‑chain best): plamb_runs/tools/bao_hht_2bin_z0p7_longbest/

• Per‑z Δχ² CSV: plamb_runs/tools/bao_chi2_perz_const_vs_2bin.csv

• LOOCV (2‑bin): plamb_runs/tools/bao_loocv_2bin_z0p7/bao_loocv.csv

• ρ‑surface scan (z=0.7): plamb_runs/tools/bao_2bin_rho_surface_z0p7.csv

Recommended next checks

• Verify stability when adding/altering low‑z BAO points.

• Consider a limited 3‑bin test if justified by data volume.

• Cross‑validate HHT energy localization vs. the Δχ² map.

• Re‑evaluate IC penalties if promoting ρ(z) to a parametric form.

## Extracted Tables

### Table 1

Model | H0 | Ωm | χ²_BAO | χ²_total | AIC (eff.) | BIC (eff.) | Δχ²_total
Constant ρ | 73.135540 | 0.236538 | 237.172 | 2334.681 | 2338.681 | 2349.576 | —
2‑bin ρ(z)  (z<0.7:+0.374, z≥0.7:−0.117) | 73.155974 | 0.236852 | 226.774 | 2321.514 | 2327.514 | 2343.856 | -13.167

### Table 2

z | Δχ² (2‑bin − const) | Interpretation
0.380 | -9.805 | Improves fit at this z‑pair (notably z≈0.38, 0.70).
0.510 | +3.858 | Slightly worse fit at this z‑pair.
0.700 | -9.735 | Improves fit at this z‑pair (notably z≈0.38, 0.70).
0.845 | +2.520 | Slightly worse fit at this z‑pair.
0.850 | -0.387 | Improves fit at this z‑pair (notably z≈0.38, 0.70).
1.480 | -0.707 | Improves fit at this z‑pair (notably z≈0.38, 0.70).
2.330 | +1.017 | Slightly worse fit at this z‑pair.
