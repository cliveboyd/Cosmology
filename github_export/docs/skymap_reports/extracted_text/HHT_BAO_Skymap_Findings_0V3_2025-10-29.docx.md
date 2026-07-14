# HHT BAO Skymap Findings 0V3 2025-10-29

- Source: `HHT_BAO_Skymap_Findings_0V3_2025-10-29.docx`
- Source size: 133580 bytes
- Source modified: 2025-10-30T00:27:30
- Extracted: 2026-07-14
- Word count estimate: 1134

## Extracted Text
Executive Summary & Technical Brief

Prepared for: PLamb / Clive S. Boyd   •   Date: 2025-10-29

Executive Summary

This document consolidates progress to date on the Hilbert–Huang Transform (HHT) residual analysis of supernova datasets and the construction of a BAO-derived anisotropy sky map.

The analysis focuses on: 

(i) intrinsic mode functions (IMFs) extracted via EMD/EEMD/CEEMDAN and evidence for IMF2/IMF3 “locking”; 

(ii) the interplay between BAO distance measures — the line-of-sight Hubble distance D_H and the transverse comoving distance D_M — and their ratio R_BAO ≡ D_H/D_M; and 

(iii) a new sky‑mapping pipeline (SM_SkyMap_FR_KdS.py) that constructs a HEALPix map of R_BAO and fits its dipole while comparing to FR/KdS model axes.

• HHT: IMF2 (and sometimes IMF3) shows localized energy ridges consistent with quasi‑stationary or gently‑chirped content; phase‑randomized surrogate tests are used to assess significance.

• BAO ratio: R_BAO = (D_H/r_d) / (D_M/r_d) = D_H / D_M reduces absolute calibration sensitivity and emphasizes geometry.

• Sky map: Mask‑aware HEALPix binning, Gaussian smoothing, dipole fitting with permutation p‑values, and redshift tomography implemented.

• Reproducibility: PNG/FITS/JSON/TXT outputs are written; tomography CSVs summarize per‑bin axes and amplitudes.

In parallel, we study two complementary cosmological distance measures: D_H (line‑of‑sight) and D_M (transverse). 

Their ratio R_BAO highlights geometry and cancels a common calibration scale. 

Projecting R_BAO across the sky tests for directional differences (a dipole), and repeating this in redshift slices tracks evolution.

HHT IMF Findings to Date

Methods: EMD, EEMD, and CEEMDAN decompose residuals into IMFs (IMF1, IMF2, IMF3, …). 
IMF2 often captures repeatable structure; IMF3 can reflect broader‑scale trends.

‘Locking’ refers to alignment between IMF energy ridges and external reference scales.

Robustness: Phase‑randomized surrogates (typically 500–1000) and parameter sweeps (pad, grids, seeds) are used to assess stability; outputs are aggregated to CSV/JSON for review.

Observations: Concentrated IMF2 energy appears within specific redshift windows, with hints of gentle chirp‑like evolution.

These are candidates for further validation and physical interpretation.

Caveats: HHT is sensitive to endpoints, noise realisations, and method settings; convergence checks across methods/seeds are essential.

BAO Distance Measures and Ratio R_BAO

Key quantities:

D_H(z) 	= c / H(z)

D_M(z) 	= (1+z) · D_A(z)

R_BAO(z) = (D_H/r_d) / (D_M/r_d) = D_H / D_M

Interpretation: R_BAO compares line‑of‑sight and transverse scalings at the same redshift.

Because r_d cancels, R_BAO is less sensitive to absolute calibration and more sensitive to geometry.

Coherent deviations can be summarized by a dipole.

R_BAO SkyMap Pipeline (SM_SkyMap_FR_KdS.py)

Input: Unified catalog with RA/Dec or Galactic l/b, redshift z, and BAO columns (DH_over_rd, DM_over_rd) or (DH, DM, rd).

Optional HHT metrics (imf2_E, imf3_E, hht_lock) for overlay.

Map: HEALPix binning with mask‑aware smoothing (numerator/denominator normalisation) and optional Galactic |b| mask in analysis/display.

Dipole: Centered linear regression for axis/amplitude with ridge stabilisation; significance from permutations over valid pixels.

Tomography: Repeat mapping in user‑defined z‑bins; bins with too few valid pixels are skipped.

CSV sidecars capture per‑bin results.

Model axes: Overlay and compare measured dipole to FR and KdS axes; report minimum axis separation in degrees.

Skymap HEALPix Map Output

Practical Guidance & Reproducibility

Example Spyder command:

%run '/Users/boyde/.spyder-py3/SM_SkyMap_FR_KdS.py' \

'/Users/boyde/.spyder-py3/plamb_runs/SkyMap/data/unified_catalog.csv' \

--out '/Users/boyde/.spyder-py3/plamb_runs/SkyMap/figs/SM_RBAO_z0p1_1p2.png' \

--nside 128 \
--zmin 0.1 \
--zmax 1.2 \
--fwhm-smooth-deg 8.0 \
--coord-in icrs \
--rbao-mode auto \
--gal-mask-abs-b-deg 20.0 \
--highlight-quantile 0.90 \
--lock-min 0.60 \
--fr-vec-lb-deg '264,48' \
--kds-vec-lb-deg '300,30' \
--dipole-search-nside 128 \
--cmap coolwarm \
--inset-loc 'lower right' \
--save-stats \
--save-json

Console Out:

[analysis] global: valid_pix = 4443/196608 (2.3%)

[dipole-fit ⟂ 1] R0 = 2.17777 A_eff = 0.0685007 (frac = 3.145% )  cond  ≈ 7.58e+01

[hemi-contrast] Δ = 0.100039  (frac = 4.594% )

[null|permute] n=500  p(A_null ≥ A_obs) ≈ 0.0020

[angle ] dipole–FR   =   48.27 deg   (axis-min =   48.27 deg)

[angle ] dipole–KdS  =   35.12 deg   (axis-min =   35.12 deg)

[dipole-dir] l =  339.53   b =   48.22

[fit-quality] R^2 ≈ 0.485 (≈ corr^2), corr(y,cosθ_axis) ≈ 0.697

[amplitude]  A_sigma ≈ 0.060355  (2.771% of mean), A_range ≈ 0.122090  (5.606% of mean)

[Y1-check] |a_l=1| ≈ 0.295202  (mask-biased sanity check)

[X-power] no template_map provided; skipping cross-power.

[analysis] z∈[0.10,0.40]: valid_pix = 4434/196608 (2.3%)

[tomography] z∈[0.10,0.40]: l=158.77 b_deg=-48.36  R0=2.42322  A_eff=5.05034e-14 (frac=0.00%)

[analysis] z∈[0.40,0.70]: valid_pix = 4001/196608 (2.0%)

[tomography] z∈[0.40,0.70]: l=48.41 b_deg=13.73  R0=1.62465  A_eff=0.140957 (frac=8.68%)

[analysis] z∈[0.70,1.00]: valid_pix = 3592/196608 (1.8%)

[tomography] z∈[0.70,1.00]: l=52.88 b_deg=27.20  R0=0.981772  A_eff=0.0509945 (frac=5.19%)

[analysis] z∈[1.00,1.20]: valid_pix = 2740/196608 (1.4%)

[tomography] z∈[1.00,1.20]: l=12.44 b_deg=58.87  R0=0.823463  A_eff=1.81724e-15 (frac=0.00%)

[tomography] saved → /Users/boyde/.spyder-py3/plamb_runs/SkyMap/figs/SM_RBAO_z0p1_1p2_tomography.csv

[X-power] C1(map, FR)  ≈ -0.131126

[X-power] C1(map, KdS) ≈ -0.164876

setting the output map dtype to [dtype('float64')]

[ok] saved stats → /Users/boyde/.spyder-py3/plamb_runs/SkyMap/figs/SM_RBAO_z0p1_1p2_stats.txt

<Figure size 864x576 with 0 Axes>

Outputs:

PNG sky map,

FITS map,

JSON summary (dipole, tomography),

and optional TXT stats;

tomography CSVs saved alongside the PNG.

Acronyms

Acronym | Meaning

BAO | Baryon Acoustic Oscillation

HHT | Hilbert–Huang Transform

IMF | Intrinsic Mode Function

EMD/EEMD/CEEMDAN | Empirical Mode Decomposition / Ensemble EMD / Complete Ensemble EMD with Adaptive Noise

HEALPix | Hierarchical Equal Area isoLatitude Pixelization

FR | Full Relativity (model axis)

KdS | Kerr–de Sitter (model axis)

R_BAO | Distance ratio D_H/D_M

NSIDE | HEALPix resolution parameter

Variables & Symbols

Symbol | Description | Units

C | Speed of light | km s⁻¹

H(z) | Hubble expansion rate at redshift z | km s⁻¹ Mpc⁻¹

D_H(z) | Hubble distance c/H(z) | Mpc

D_A(z) | Angular diameter distance | Mpc

D_M(z) | Transverse comoving distance = (1+z)·D_A | Mpc

r_d | Sound horizon at the drag epoch | Mpc

R_BAO(z) | BAO distance ratio D_H/D_M | Dimensionless

Z | Redshift | —

l, b | Galactic longitude/latitude | Deg

NSIDE | HEALPix map resolution parameter | —

imf2_E, imf3_E | IMF2/IMF3 energy metric | arb. Units

hht_lock | HHT locking score or filter | 0–1 (typ.)

Limitations & Next Steps

• Validate stability across HHT methods and seeds; record effect sizes and false‑positive rates.

• Use inverse‑variance weights whenever available (σ columns).

• Check systematics: Galactic masking thresholds, smoothing scales, and catalog sub‑selections.

• Extend tomography with adaptive binning to maintain pixel counts.

• Cross‑validate with independent datasets or mocks.
