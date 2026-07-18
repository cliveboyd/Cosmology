# SDSS DR16Q v4 Independent Cross-Catalogue Validation Preregistration

Date locked: 2026-07-18 19:44:55 AEST  
Analysis status at lock: no SDSS DR16Q v4 catalogue rows had been downloaded, opened or analysed in this branch.

## Purpose and Claim Boundary

This analysis tests whether the previously identified Quaia redshift-anisotropy candidate has a like-for-like counterpart in the independent SDSS DR16Q v4 spectroscopic quasar catalogue. It is a validation test, not a new search. A failure of any promotion gate retains the result as a catalogue-specific or unresolved survey-selection feature; it does not constitute evidence against the underlying SU(2) framework in general.

## Locked External Target

The cross-catalogue direction and amplitude were specified before accessing SDSS outcomes. The locked comparator is the **raw** Quaia dipole in the primary row `1.00 <= z < 1.50`, `|b| >= 35 deg`:

- ICRS direction: `RA = 145.3283284733 deg`, `Dec = -26.7809296987 deg`.
- Redshift-gradient amplitude: `A = 0.003049663343`.
- Source row: `su2_quaia_scanlaw_colour_residual_global_observed_windows_2026-07-15.csv`.

The template-controlled Quaia direction is not substituted for this comparator because the independent SDSS catalogue does not share the Quaia selection/template model.

## Catalogue and Integrity Rules

- Catalogue: [SDSS DR16Q v4](https://data.sdss.org/sas/dr16/eboss/qso/DR16Q/DR16Q_v4.fits).
- Expected transport size: `1,718,484,480` bytes.
- Expected binary-table row count: `750,414`.
- Required columns: `SDSS_NAME`, `RA`, `DEC`, `Z`, `IS_QSO_FINAL`, `Z_CONF`, `ZWARNING`, and `SN_MEDIAN_ALL`.
- The downloaded file must have the expected byte count, a readable FITS binary table, the expected row count and all required columns. Its SHA-256 digest is recorded after download; no unverified expected digest is assumed.
- Primary quality selection: `IS_QSO_FINAL == 1`, finite `RA`, `DEC` and `Z`, `0 <= RA < 360`, and `-90 <= DEC <= 90`.
- Duplicate `SDSS_NAME` values, if any, are resolved before redshift/window cuts by retaining the row with the highest `Z_CONF`, then highest finite `SN_MEDIAN_ALL`, then lowest original row index.
- `ZWARNING` and `Z_CONF` are not primary exclusion cuts: the catalogue's final visual classification is the official primary quasar decision. Two declared quality controls are nevertheless reported: `ZWARNING == 0` and `Z_CONF >= 2`.

## Fixed Windows and Latitude Cuts

The upper redshift edge is exclusive throughout.

| Role | Tag | Redshift interval |
| --- | --- | --- |
| primary | `z1p00_1p50` | `1.00 <= z < 1.50` |
| perturbation | `z0p95_1p45` | `0.95 <= z < 1.45` |
| perturbation | `z1p05_1p55` | `1.05 <= z < 1.55` |
| perturbation | `z1p00_1p40` | `1.00 <= z < 1.40` |
| perturbation | `z1p10_1p50` | `1.10 <= z < 1.50` |

The fixed Galactic latitude cuts are `|b| >= 20, 25, 30, 35 deg`; `35 deg` is the primary cut. The primary global family is the Cartesian product of five windows, four latitude cuts and four statistics, giving 80 tests.

## Mask and Coordinate Treatment

The analysis is conditional on the observed DR16Q angular positions after the primary quality selection. No full-sky completion or post-outcome footprint trimming is permitted. Each row's design matrix is evaluated only at its observed positions, so the SDSS footprint and its harmonic mode coupling enter the fitted Gram matrix directly. Galactic coordinates use the fixed ICRS-to-Galactic rotation implemented by `astropy.coordinates.SkyCoord`.

A row is invalid, rather than silently modified, if it has fewer than 1,000 objects, a rank-deficient full design matrix, or a full-design condition number above `1e8`.

## Harmonic Model and Statistics

For each fixed window/cut row, redshift is modelled jointly with real spherical harmonics through `l = 3`:

\[
\begin{aligned}
z_i &= a_0 + \sum_{l=1}^{3}\sum_{m=-l}^{l} a_{lm}Y^{\mathrm{R}}_{lm}(l_i,b_i) + \epsilon_i, \\
F_l &= \frac{(\mathrm{RSS}_{-l}-\mathrm{RSS}_{1:3})/(2l+1)}{\mathrm{RSS}_{1:3}/(N-16)}, \\
F_{1:3} &= \frac{(\mathrm{RSS}_{0}-\mathrm{RSS}_{1:3})/15}{\mathrm{RSS}_{1:3}/(N-16)}.
\end{aligned}
\]

Here `RSS_-l` is the residual sum of squares after omitting the entire `l` group while retaining the other declared multipoles. The four registered statistics are `F_1`, `F_2`, `F_3` and `F_1:3`. For comparable family maximisation they are transformed to `S = -log10(P_F)`, using the nominal F survival probability only as a monotonic scale; inferential probabilities come exclusively from the spatial-block null.

The dipole amplitude and direction are obtained from the `l = 1` Cartesian Galactic basis fitted jointly with `l = 2,3`; the fitted direction is also transformed to ICRS for the locked cross-catalogue comparison. Multipole powers use the squared Euclidean norm of each real-harmonic coefficient group divided by `2l+1` and are descriptive under the partial-sky mask.

## Spatial-Block Null and Global Family

The null is a fixed-design, restricted-residual Rademacher wild bootstrap. Within each row, residuals from the intercept-only null are multiplied by a common `+1/-1` sign for every object in the same equal-area Galactic longitude by `sin(b)` cell. The same cell-sign realisation is shared across all overlapping windows and latitude cuts, preserving their covariance. Positions and redshifts are never reassigned between footprint regions.

Three fixed nominal cell scales are used independently:

| Scale | Mocks | Seed |
| --- | ---: | ---: |
| `8 deg` | 20,000 | 20,268,790 |
| `12 deg` | 20,000 | 20,272,826 |
| `16 deg` | 20,000 | 20,276,862 |

At each scale, the global statistic is the maximum `S` over the complete 80-test family. The Monte Carlo probability is

\[
\begin{aligned}
p_{\mathrm{global}} &= \frac{1 + \#\{S^{(j)}_{\max} \ge S_{\max,\mathrm{obs}}\}}{1 + N_{\mathrm{mock}}}, \\
p_{\mathrm{point}} &= \frac{1 + \#\{S^{(j)}_{r,k} \ge S_{r,k,\mathrm{obs}}\}}{1 + N_{\mathrm{mock}}}.
\end{aligned}
\]

No block scale may be selected after outcomes are inspected. The conservative cross-scale readout is the largest of the three global probabilities.

## Promotion Gates

Promotion requires **all** gates below. No weighted score or compensating success is allowed.

1. **Integrity:** the catalogue and every registered row pass the locked file, sample-size, rank and conditioning checks.
2. **Global rarity:** `p_global <= 0.01` independently at 8, 12 and 16 degrees.
3. **Mask-calibrated joint coherence:** the primary row's joint `F_1:3` has empirical `p_point <= 0.05` independently at 8, 12 and 16 degrees. Individual `l=1,2,3` probabilities are retained so any single-multipole domination remains visible.
4. **Cross-catalogue direction:** the primary SDSS dipole lies within 30 degrees of the locked raw Quaia ICRS direction.
5. **Window direction stability:** at `|b| >= 35 deg`, every perturbation-window dipole lies within 30 degrees of the primary-window dipole.
6. **Latitude direction stability:** in the primary window, every latitude-cut dipole lies within 30 degrees of the primary-cut dipole.
7. **Amplitude stability:** every primary/perturbation analysis-row dipole amplitude divided by the primary-row amplitude lies in `[0.5, 2.0]`. This includes all fixed window perturbations at `|b| >= 35 deg` and all fixed latitude cuts in the primary window.

The declared `ZWARNING == 0` and `Z_CONF >= 2` controls remain diagnostics and are not additional promotion gates or members of the 80-test family.

If a direction is numerically unstable because its dipole amplitude is effectively zero, the corresponding stability gate fails rather than being waived.

## Reporting Rules

All observed rows, per-scale pointwise probabilities, global maxima, null quantiles, condition numbers, direction separations, amplitude ratios and gate decisions are retained. Every failed gate and failed subcomponent must be listed explicitly. The report must call the result `not promoted` unless every gate passes. Any later exploratory changes require a new dated registration and must not be merged into this confirmatory family.
