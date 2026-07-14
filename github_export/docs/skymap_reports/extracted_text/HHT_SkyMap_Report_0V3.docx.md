# HHT SkyMap Report 0V3

- Source: `HHT_SkyMap_Report_0V3.docx`
- Source size: 5787 bytes
- Source modified: 2025-10-30T00:33:42
- Extracted: 2026-07-14
- Word count estimate: 1362

## Extracted Text
HHT–BAO Anisotropy & Sky Map Report (0V3)

Generated 2025-10-29 13:32 (Australia/Perth)

EXECUTIVE SUMMARY

- We constructed all-sky HEALPix maps of the BAO ratio R_BAO ≡ D_H / D_M from a unified catalog, applying optional Galactic latitude masks and Gaussian smoothing. A dipole component was fitted with a centered linear model and evaluated with permutation nulls.

- Across a quick grid of configurations (NSIDE, smoothing, masking, and redshift binning), we consistently find a significant dipole amplitude (permutation p≈0.005 with n=200 permutations for each run).

- However, the recovered dipole direction is not stable across all configurations, clustering into a few regions on the sky (e.g., near (l,b)≈(302°,-12°); (42°,55°); and (5°,15°)). This sensitivity suggests either mask/coverage effects, redshift slicing dependencies, or residual systematics.

- Angular separations to the reference FR and KdS axes are generally large (Δ_FR ≈ 69–86°, Δ_KdS ≈ 40–72° in the sampled runs), so there is no compelling alignment with either axis at this stage.

- IMF2/IMF3 overlays and locking filters are implemented in the plotting pipeline, but interpretation should be reserved until we finalize the map configuration and stability checks.

BACKGROUND (Plain-English for a Technical Reader)

BAO distances and the ratio R_BAO:

• D_H(z) = c / H(z) is the ‘Hubble distance’ that converts a small change in redshift to a radial comoving distance.

• D_M(z) is the transverse comoving distance that sets the angular scale for a fixed comoving ruler at redshift z.

• r_d is the sound horizon at the drag epoch (used to non-dimensionalize BAO measurements).

• In many catalogs, the reported quantities are D_H/r_d and D_M/r_d. Their ratio cancels r_d:

R_BAO(z) = (D_H/r_d) / (D_M/r_d) = D_H / D_M.

If the Universe is perfectly isotropic on large scales, sky maps of R_BAO (after appropriate masking and smoothing) should not show a strong dipole component. A reproducible dipole can indicate anisotropy or systematics (e.g., survey footprint, calibration, or redshift-dependent selection).

HHT IMF FINDINGS (Brief Recap)

We apply the Hilbert–Huang Transform (HHT) to supernova residuals to decompose signals into Intrinsic Mode Functions (IMFs).

Key notions:

• IMF2/IMF3: The 2nd and 3rd IMFs often carry mid-frequency structure; their instantaneous energy and phase are tracked vs. z.

• Locking: We quantify phase/energy ‘locking’ between IMF components and external scales (e.g., BAO-related 1/k features).

• Robustness: We evaluate significance via phase-randomized surrogates and grid sweeps (EMD/EEMD/CEEMDAN; seeds; pads).

To date, IMF2 ‘locking’ has shown hints of structure near BAO scales, but the effect size depends on method, seeds, and masking. We therefore treat HHT-based signals as suggestive pending full robustness and cross-validation with independent subsamples.

METHODS (Sky Map Pipeline)

Inputs: unified catalog with columns (suggested): ra_deg, dec_deg, z, DH_over_rd, DM_over_rd (or DH, DM, rd), and optionally imf2_E, imf3_E, hht_lock.

Coordinates: RA/Dec are converted to Galactic (l,b).

Map-making (NSIDE=N): Bin points to HEALPix pixels; average R_BAO per pixel with optional inverse-variance weights; apply Gaussian smoothing (FWHM). A Galactic |b|≥b_min mask can be applied consistently to prevent leakage.

Dipole fit: Solve a centered linear model Y ≈ R0 + a·n with ridge regularization; report effective amplitude A_eff and axis (l,b).

Null test: Permute pixel values (n permutations), refit dipole each time, and estimate p(A_null ≥ A_obs).

Tomography: Repeat map-making and dipole fit in redshift bins to track evolution and directional stability.

Cross-power C1: Correlate the map with cosθ templates of candidate axes (FR, KdS) to gauge directional affinity.

KEY EQUATIONS

D_H(z) = c / H(z)

D_M(z) =

⎧  c ∫₀ᶻ dz′ / H(z′)                                 (k=0)

⎨  (c/√|k|) sin( √|k| ∫₀ᶻ dz′ / H(z′) )              (k<0)  [open]

⎩  (c/√k)  sinh( √k  ∫₀ᶻ dz′ / H(z′) )               (k>0)  [closed]

R_BAO(z) = D_H(z) / D_M(z)

Dipole (centered fit on valid pixels):

Let y be the map values, y0 = mean(y), and s_i = n_i·d where d is the unit dipole direction, n_i the pixel unit vector.

Solve (V_cᵀ V_c + λI) a = V_cᵀ (y − y0), where V rows are pixel unit vectors and V_c is column-centered.

The effective dipole amplitude A_eff ≈ ½ [max(y_dip) − min(y_dip)], y_dip = V_c a.

QUICK GRID: RECENT SKYMAP RUNS (Oct 29, 2025)

Permutation count per run: n=200 → p_min ≈ 1/(n+1) ≈ 0.004975

FR axis = (l,b) = (264°, 48°);  KdS axis = (300°, 30°)

Summary (angles are axis-min ambiguities, 0–90°):

tag                        l_deg        b_deg        ΔFR_min   ΔKdS_min   R0        A_eff      frac      p_perm

ns32_s8_mask20_noaxes     306.0221    -9.7816       68.6781   40.2014    2.258415  0.316877   0.140310  0.004975

ns32_s6_nomask_auto       302.4955   -12.0909       69.1195   42.1595    2.337374  0.359368   0.153749  0.004975

ns64_s8_mask20_axes        42.4437    55.3468       70.9370   72.2301    2.226457  0.208888   0.093821  0.004975

ns64_s8_mask20_kds_only    42.4437    55.3468       70.9370   72.2301    2.226457  0.208888   0.093821  0.004975

ns64_s6_mask20_bins2_fr     5.0758    14.9948       86.0947   61.1906    2.275749  0.201339   0.088472  0.004975

ns64_s6_mask20_bins5        5.0758    14.9948       86.0947   61.1906    2.275749  0.201339   0.088472  0.004975

INTERPRETATION

• Significance: All sampled configurations yield p≈0.005 (with n=200 perms), indicating the observed dipole amplitudes are unlikely under pixel-value permutations. Increase permutations (e.g., n≥2000) to firm up p-values and control for multiple tests.

• Directional stability: The recovered dipole direction moves notably with masking, smoothing, and redshift binning. This can arise from (i) sky coverage and selection anisotropies, (ii) redshift-evolving systematics, or (iii) genuine large-scale anisotropy whose projection mixes with mask boundaries. Stability improves confidence; instability urges caution.

• FR/KdS alignment: No tight alignment (≲30°) with FR or KdS axes in these quick runs. Δ_FR is ~69–86°, Δ_KdS ~40–72°.

• Amplitude scale: A_eff/mean ≈ 9–15%. This is large enough to track but small enough to be susceptible to subtle systematics.

RECOMMENDED NEXT STEPS

1) Robust statistics:

– Increase permutations to n≥2000 (per run) and apply FDR/Bonferroni across the grid.

– Add ‘spin’ tests (random rigid rotations of the map relative to the mask) to better model mask-induced couplings.

2) Coverage control:

– Enforce identical sky masks across runs (analysis *and* display) and consider apodized masks.

– Try NSIDE up/down-sampling to a common resolution before fitting; validate with Y₁ amplitude checks.

3) Weighting and errors:

– Use inverse-variance weights derived from σ(D_H/r_d) and σ(D_M/r_d); propagate to R_BAO uncertainty carefully.

– Jackknife by sky sectors (e.g., 12 or 48 HEALPix superpixels) to assess regional leverage.

4) Redshift tomography:

– Require a minimum number of valid pixels per bin (already in code); compare directions across adjacent bins for continuity; test alternative bin edges.

5) Cross-power with templates:

– Compute C₁(map,cosθ_FR) and C₁(map,cosθ_KdS) for each run; bootstrap errors over pixels to get uncertainties.

6) HHT overlays:

– Freeze a specific locking threshold and highlight quantile, then test stability under binning and masking.

7) Systematics audits:

– Re-run maps on subcatalogs (by survey, footprint, or calibration version); inspect ecliptic/Equatorial masks;

– Validate with null maps (scramble z within sky patches; or randomize positions within patches).

ACRONYMS & SYMBOLS

BAO     — Baryon Acoustic Oscillation

SN      — Supernova (Type Ia)

HHT     — Hilbert–Huang Transform

IMF     — Intrinsic Mode Function (IMF2, IMF3)

NSIDE   — HEALPix resolution parameter (N_pix = 12·NSIDE²)

R_BAO   — D_H / D_M; the ratio of radial to transverse BAO distances

D_H     — c/H(z), ‘Hubble distance’

D_M     — Transverse comoving distance

r_d     — Sound horizon at the drag epoch

Y₁      — Spherical-harmonic dipole (ℓ=1)

A_eff   — Effective half-range of the fitted dipole signal

p_perm  — Permutation-based p-value for dipole amplitude

FR axis — (l,b)=(264°,48°)

KdS axis— (l,b)=(300°,30°)

VARIABLE NOTES

ra_deg, dec_deg : Right ascension and declination (ICRS), degrees

l_deg,  b_deg   : Galactic longitude and latitude, degrees

z               : Redshift

DH_over_rd      : D_H / r_d  (if provided in the catalog)

DM_over_rd      : D_M / r_d  (if provided in the catalog)

imf2_E, imf3_E  : IMF energies (from HHT decomposition)

hht_lock        : Locking scalar (0–1) indicating selection strength

REPRODUCIBILITY (CLI)

Example quick grid command that produced the summary:

%run '/Users/boyde/.spyder-py3/sm_quickgrid.py' \

--catalog /Users/boyde/.spyder-py3/plamb_runs/SkyMap/data/unified_catalog.csv \

--script  /Users/boyde/.spyder-py3/SM_SkyMap_FR_KdS.py \

--outdir  /Users/boyde/.spyder-py3/plamb_runs/SkyMap/figs/quickgrid
