# RBAO dipole prelim report

- Source: `RBAO_dipole_prelim_report.docx`
- Source size: 98762 bytes
- Source modified: 2025-10-25T07:45:40
- Extracted: 2026-07-14
- Word count estimate: 700

## Extracted Text
Preliminary R_BAO Dipole Analysis

Generated: 2025-10-23 06:45 UTC

Overview

We map the BAO ratio R_BAO = (D_H / r_d) / (D_M / r_d) = D_H / D_M from the unified catalog onto a HEALPix sky and fit a dipole.

The map shown uses NSIDE=64, Gaussian smoothing FWHM=8°, redshift range z∈[0.1, 1.2], and a Galactic mask |b|>20°.

R_BAO values are constructed from columns DH_over_rd and DM_over_rd.

Inputs & Settings

• Catalog columns used: ra_deg, dec_deg, z, DH_over_rd, DM_over_rd, kind (BAO).

• Pixelization: HEALPix NSIDE=64 (RING), smoothed with FWHM=8°.

• Mask: Galactic latitude |b|>20°.

• R_BAO mode: DH_over_rd/DM_over_rd.

• FR axis (for reference): (l,b)=(264°, 48°).

• KdS axis (for reference): (l,b)=(300°, 30°).

Method (Summary)

1) Construct R_BAO per object and bin to HEALPix pixels with optional weights (if σ columns provided). 

2) Smooth the partial-sky map to mitigate sparsity. 

3) Fit a dipole on pixels via a centered linear model: R(n) = R0 + a·n. 

4) Quantify hemispheric contrast Δ across the fitted axis.

5) Evaluate null significance via permutation with re-estimation of the dipole axis per shuffle.

6) Report angles between the fitted dipole and supplied FR/KdS axes.

Equations (plain text)

• R_BAO = (D_H/r_d) / (D_M/r_d) = D_H / D_M
• Dipole model: R( n̂ ) = R₀ + a · n̂, with amplitude |a|
• Hemisphere contrast across axis d̂: Δ = ⟨R⟩_{hemisph(d̂)}^+ − ⟨R⟩_{hemisph(d̂)}^−
• 1D projection along fitted axis: s = cosψ = n̂ · d̂; regression R ≈ R₀ + β (s − ⟨s⟩)
• Practical amplitudes from β:
    – A_sigma = |β| · σ_s  (typical variation over footprint)
    – A_range = |β| · (max(s) − min(s))/2  (extreme hemispherical swing)

Latest Results

Quantity | Value

R₀ (mean level) | 0.320486

Dipole amplitude A_eff | 0.006847  (~2.136% of mean)

Hemispheric contrast Δ | 0.038987  (~12.165% of mean)

Permutation p-value | 0.904 (re-fit axis per shuffle)

Dipole pole (l, b) | (160.81°, -49.49°)

Angle to FR axis | 131.62°  (axis-min 48.38°)

Angle to KdS axis | 143.70°  (axis-min 36.30°)

1D regression R² | 0.533

corr(R, cosψ) | 0.730

A_sigma | 0.026119  (~8.150% of mean)

A_range | 0.046312  (~14.450% of mean)

Interpretation

• The smoothed R_BAO map exhibits a shallow gradient: A_eff ≈ 0.00685 (~2.1% of mean), 
with hemispheric contrast Δ ≈ 0.039 (~12% of mean).

• Under a permutation test that re-estimates the dipole each shuffle, p ≈ 0.90—i.e., the observed amplitude is consistent with chance given the current partial-sky footprint and smoothing.

• The fitted dipole axis is well separated from both reference axes (FR and KdS), so no notable alignment is present in this configuration.

• 1D diagnostics over the footprint indicate moderate correlation with the projected cosine (corr ≈ 0.73), but this is not compelling given the high p-value under permutation.

Are the marked directions near known features?

The fitted dipole pole (l≈161°, b≈−49°) lies in the southern Galactic cap and does not coincide with commonly-cited large-scale anomalies (e.g., the CMB kinematic dipole direction, 
Great Attractor/Shapley region, or the well-known CMB Cold Spot).

A targeted cross-match against catalogs (CMB features, LSS overdensities, and previous anisotropy claims) would be the right next step if we want to probe this further.

Outstanding Issues & Next Steps

1) Incorporate per-pixel inverse-variance weights using σ(DH_over_rd) and σ(DM_over_rd) when available.

2) Test stability versus NSIDE, FWHM smoothing, and the Galactic mask threshold.

3) Slice the map in narrower redshift bins to check for redshift evolution of the dipole.

4) Quantify and mitigate partial-sky biases (e.g., fit in harmonic space; include survey window/mask explicitly).

5) Try jackknife by survey footprint or hemisphere to assess robustness.

6) Compare against independent datasets (DESI, eBOSS subsets) and synthetic mocks.

Reproduction (command used)

%run '/Users/boyde/.spyder-py3/SM_SkyMap_FR_KdS.py' \
  '/Users/boyde/.spyder-py3/plamb_runs/SkyMap/data/unified_catalog.csv' \
  --out '/Users/boyde/.spyder-py3/plamb_runs/SkyMap/figs/SM_RBAO_z0p1_1p2.png' \
  --nside 64 --zmin 0.1 --zmax 1.2 --fwhm-smooth-deg 8.0 \
  --coord-in icrs --rbao-mode DH_over_rd/DM_over_rd --gal-mask-abs-b-deg 20.0 \
  --fr-vec-lb-deg "264,48" --kds-vec-lb-deg "300,30" \
  --dipole-search-nside 32 --cmap coolwarm

Figure: Latest sky map
