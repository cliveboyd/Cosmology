# SU2 / Quaia Matched Catalogue-Quality Gate

Date: July 15, 2026

## Purpose

This gate tests whether the locked Quaia `0.95 < z < 1.45` angular mode survives when catalogue-quality variables are explicitly balanced between the two hemispheres defined by the fitted dipole axis.

The analysis uses propensity-overlap weights. For each latitude cut, a logistic model predicts the positive dipole-axis hemisphere from catalogue-quality variables. Objects are then overlap-weighted so the two hemispheres have matched quality distributions before the redshift dipole is refit.

## Model

```text
z_i       = c0 + d . n_i + epsilon_i
A_dipole  = ||d||
SNR_d     = sqrt(d^T Cov(d)^(-1) d)
```

## Locked Configuration

- redshift window: `0.95 <= z < 1.45`
- latitude cuts: `10, 15, 20, 25, 30, 35`
- matched-quality shuffle mocks per group: `500`

## Matched-Quality Summary

| fit_group | amp_ratio_vs_baseline | max_direction_sep_vs_baseline_deg | snr_max | max_abs_smd_after | max_weighted_ks_after | readout |
| --- | --- | --- | --- | --- | --- | --- |
| zerr_only | 1.01004 | 0.475513 | 4.8135 | 0.000326864 | 0.00544453 | survives matched-quality gate |
| gaia_colour_quality | 0.62955 | 10.7436 | 3.12893 | 0.000273522 | 0.00706574 | does not pass matched-quality gate |
| wise_colour_quality | 0.724522 | 3.86062 | 3.77471 | 0.000234986 | 0.00529918 | does not pass matched-quality gate |
| pm_error_quality | 1.02777 | 2.59327 | 4.85898 | 0.000323531 | 0.00887279 | survives matched-quality gate |
| all_catalogue_quality | 0.309709 | 121.39 | 1.87638 | 0.000365446 | 0.00760149 | does not pass matched-quality gate |

## Matched-Quality Shuffle P-Values

| fit_group | n_mocks | p_snr_max_ge_observed | p_pair_sep_le_observed | p_coherence_ge_observed | p_joint_snr_and_pair_sep |
| --- | --- | --- | --- | --- | --- |
| zerr_only | 500 | 0.00199601 | 0.00199601 | 0.00199601 | 0.00199601 |
| gaia_colour_quality | 500 | 1 | 0.00199601 | 0.756487 | 0.00199601 |
| wise_colour_quality | 500 | 0.113772 | 0.00199601 | 0.00199601 | 0.00199601 |
| pm_error_quality | 500 | 0.00199601 | 0.00199601 | 0.00199601 | 0.00199601 |
| all_catalogue_quality | 500 | 1 | 0.986028 | 1 | 0.986028 |

## Interpretation Gate

- Amp/direction gate: amplitude ratio >= `0.8` and maximum direction rotation <= `30.0` degrees.
- Balance gate: post-weight max absolute SMD <= `0.1` and post-weight weighted KS <= `0.1`.
- If all-catalogue quality passes both gates, the catalogue-quality explanation is weakened.
- If all-catalogue quality fails the amp/direction gate after good balance, the catalogue-quality explanation remains the leading caveat.

## Outputs

- summary: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_matched_quality_gate_20260715\su2_quaia_matched_quality_summary.csv`
- fit rows: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_matched_quality_gate_20260715\su2_quaia_matched_quality_fit_rows.csv`
- balance diagnostics: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_matched_quality_gate_20260715\su2_quaia_matched_quality_balance.csv`
- p-values: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_matched_quality_gate_20260715\su2_quaia_matched_quality_pvalues.csv`
