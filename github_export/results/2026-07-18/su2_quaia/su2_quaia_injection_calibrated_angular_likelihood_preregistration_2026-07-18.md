# SU2 / Quaia Injection-Calibrated Angular-Likelihood Preregistration

Date locked: 18 July 2026 (Australia/Sydney)

Status: locked before the new observed held-out statistic or injection-calibrated threshold was calculated.

## Scientific Question

Does a redshift-dependent angular term shaped by the independently fitted SU2R-minus-LCDM scalar distance-modulus residual improve prediction of Quaia redshifts on spatially held-out sky regions after dust, catalogue depth, Gaia scanning geometry, WISE colour and catalogue-quality templates have been fitted jointly?

This is a targeted follow-up to the previously selected `z ~ 1-1.5` angular candidate. It is not an independent discovery test and cannot replace cross-catalogue validation.

## Locked Sample and Geometry

- Catalogue: local Quaia `G < 20.0` FITS catalogue.
- Primary redshift window: `1.00 <= z < 1.50`.
- Primary latitude cut: `|b| >= 15 deg`.
- Angular axis: the pre-existing raw Quaia `bcut=15 deg` direction, `(l, b) = (138.9906218, -36.6041026) deg`.
- Reference redshift-dipole amplitude: `A_ref = 0.002548691334890594`.
- Held-out folds: twelve fixed Galactic-longitude sectors. Equal-area 12-degree longitude-by-sine-latitude cells are assigned to the sector containing their cell centre, so a bootstrap cell is never split between train and test data.
- Spatial null blocks: the same 12-degree equal-area Galactic cells, with one Rademacher sign per cell.

## Locked Scalar Shape

The scalar curve is

```text
s(z) = mu_SU2R(z) - mu_LCDM(z),
```

using the existing DESI DR2 plus Pantheon plus Planck best-fit files:

- `plamb_runs/updated_datasets_20260710/runs/desi_SU2R_Pantheon_SN_BAO_Planck/SU2R_emcee_bestfit.txt`
- `plamb_runs/updated_datasets_20260710/runs/desi_LCDM_Pantheon_SN_BAO_Planck/LCDM_emcee_bestfit.txt`

The fixed scalar-modulated angular basis is

```text
w_i          = standardise[s(z_i)],
x_SU2,i      = w_i n_i,
d_locked     = A_ref u_locked.
```

Here `n_i` is the Galactic unit vector and `u_locked` is the locked angular axis. The three components of `x_SU2,i` allow the angular vector to be recovered rather than imposed. The standardisation of `w_i` is unsupervised and is calculated once over the eligible sample. The scalar shape is not refitted in a hold-out fold.

This interaction avoids treating `s(z)` alone as a predictor of `z`, which would leak a deterministic transform of the response into the model.

## Joint Template Model

All continuous nuisance columns are median-centred and divided by their standard deviation over the eligible sample. The fixed template stack is:

- dust: SFD `E(B-V)`, `log(1+E(B-V))` and squared `E(B-V)`;
- depth and footprint: Quaia selection function, random-catalogue density, AllWISE W1/W2 median and low-coverage fields, and Gaia source density when available;
- Gaia scanning geometry: DR3 scan-count, mean cosine and sine of twice the scan angle, and scan-angle resultant;
- WISE and Gaia colour: W1-W2, BP-RP and their standardised interaction;
- catalogue quality: Quaia redshift uncertainty, Gaia G magnitude, W1 and W2 magnitudes, and proper-motion uncertainties.

Columns with zero variance or non-finite values are removed by a documented deterministic rule. Rows must be finite for every retained locked column. No outcome-driven template selection is allowed.

## Candidate Models and Held-Out Score

For response `z_i`, the models are

```text
M_T:  z_i = alpha + T_i beta + epsilon_i,
M_S:  z_i = alpha + T_i beta + d . x_SU2,i + epsilon_i.
```

Each model is fitted by ordinary least squares on eleven sectors and scored on the untouched twelfth sector. The Gaussian residual variance is estimated from the training portion only. The primary statistic is

```text
Delta_CVLPD = sum_f [LPD_test,f(M_S) - LPD_test,f(M_T)],
```

where positive values favour the locked SU2-shaped angular term. Fold-level log-predictive-density differences and squared-error changes will also be reported. They are diagnostics, not additional selection statistics.

## Injection and Power Calibration

The template-only full-sample fit supplies fitted values and residuals. In calibration replicate `r`, one common sign is drawn for each 12-degree spatial cell:

```text
z_i^(r,a) = zhat_T,i + xi_cell(i),r e_T,i + a d_locked . x_SU2,i,
```

with `xi` in `{-1,+1}` and locked amplitude scales `a = 0, 0.5, 1.0, 1.5, 2.0`. The same spatial-sign draw is reused across amplitude scales within a replicate. The primary run uses 2,000 paired replicates and seed `180718`.

The first 1,000 null replicates are a threshold-calibration partition. The detection threshold is their empirical 99th percentile of `Delta_CVLPD`. The second 1,000 null replicates are an untouched validation partition used to estimate the realised false-positive rate. Power and physical recovery are also evaluated on this validation partition. This replaces the former fixed `Delta BIC < -10` requirement.

A non-zero injection is counted as recovered only when all of the following hold in the same replicate:

- `Delta_CVLPD` is positive and strictly above the frozen threshold;
- the pooled cross-fold recovered vector has positive orientation relative to the injected vector;
- recovered-to-injected amplitude lies in `[0.5, 1.5]`;
- recovered direction is within `30 deg` of the injected direction.

## Locked Decision Rule

The angular likelihood passes only if all conditions hold conjunctively:

1. the independent validation-null false-positive rate is at most `0.01` under the frozen threshold;
2. joint held-out, direction and amplitude recovery at `1.0 A_ref` is at least `0.80` under that same threshold;
3. observed pooled `Delta_CVLPD` is positive and strictly greater than the frozen threshold;
4. the observed pooled cross-fold angular vector has positive orientation, an amplitude ratio in `[0.5, 1.5]` relative to `A_ref`, and direction separation at most `30 deg` from the locked axis;
5. at least eight of the twelve held-out sectors have positive `Delta LPD`.

If condition 2 fails, the result is labelled **underpowered** as an analysis outcome. The threshold remains frozen and is not relaxed, moved or re-estimated. If false-positive control and power are adequate but one or more observed-data conditions fail, the result is a **gate failure**. Passing this gate would justify follow-up only; promotion still requires the separately preregistered independent-catalogue test.

In-sample AIC, BIC and likelihood-ratio values may be reported as diagnostics. No in-sample BIC result can substitute for positive pooled held-out predictive gain or any other conjunctive condition.

## Fixed Outputs

The run will write full calibration draws, fold scores, power summaries, configuration, decision JSON, a technical report and a compact readout under `plamb_runs/diagnostics/su2_quaia_injection_calibrated_angular_likelihood_20260718/`. Compact, hash-audited outputs will be exported to `github_export/results/2026-07-18/su2_quaia/`.
