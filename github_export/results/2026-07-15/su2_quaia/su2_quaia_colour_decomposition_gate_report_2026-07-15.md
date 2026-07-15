# SU2 / Quaia Colour-Decomposition Gate

Date: July 15, 2026

## Purpose

This follow-up decomposes the catalogue-quality suppression found in the matched-quality gate. The locked target remains the Quaia `0.95 <= z < 1.45` angular mode, with the same latitude cuts and redshift-dipole model. Only the matched feature groups change.

The aim is to distinguish four explanations:

- Gaia optical magnitude or colour.
- WISE infrared magnitude or colour.
- Quaia redshift-uncertainty structure.
- Interactions between colour space and the Quaia redshift estimator.

## Model

```text
z_i       = c0 + d . n_i + epsilon_i
A_dipole  = ||d||
SNR_d     = sqrt(d^T Cov(d)^(-1) d)
p_i       = P(H_i = 1 | x_i)
w_i       = 1 - p_i,  H_i = 1
w_i       = p_i,      H_i = 0
```

Here `H_i` is the hemisphere defined by the baseline dipole axis for the relevant latitude cut, and `x_i` is the tested catalogue-quality feature block. Polynomial degree-2 propensity features are used, so multi-variable blocks include squares and pairwise interactions.

## Locked Configuration

- redshift window: `0.95 <= z < 1.45`
- latitude cuts: `10, 15, 20, 25, 30, 35`
- matched-feature shuffle mocks per group: `300`
- propensity strata for redshift shuffles: `10`

## Ranked Suppressors

| rank_by_amp_suppression | fit_group | amp_ratio_vs_baseline | max_direction_sep_vs_baseline_deg | snr_max | p_joint_snr_and_pair_sep | readout |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | all_catalogue_quality | 0.309709 | 121.39 | 1.87638 | 0.990033 | does not pass matched-quality gate |
| 2 | colour_redshift_estimator_compact | 0.41956 | 27.5163 | 2.45191 | 0.013289 | does not pass matched-quality gate |
| 3 | gaia_wise_colour_cross | 0.508139 | 18.5926 | 2.78799 | 0.00332226 | does not pass matched-quality gate |
| 4 | zerr_x_gaia_wise_colour | 0.527275 | 17.9618 | 2.90093 | 0.00332226 | does not pass matched-quality gate |
| 5 | gaia_bp_rp_only | 0.626609 | 10.9648 | 3.11152 | 0.00332226 | does not pass matched-quality gate |
| 6 | gaia_g_plus_bp_rp | 0.62955 | 10.7436 | 3.12893 | 0.00332226 | does not pass matched-quality gate |
| 7 | zerr_x_gaia_colour | 0.639289 | 10.8356 | 3.16836 | 0.00332226 | does not pass matched-quality gate |
| 8 | wise_w1_plus_w2 | 0.717794 | 3.18802 | 3.85253 | 0.00332226 | does not pass matched-quality gate |
| 9 | wise_w1_w2_plus_mags | 0.724522 | 3.86062 | 3.77471 | 0.00332226 | does not pass matched-quality gate |
| 10 | zerr_x_wise_colour | 0.809327 | 2.77263 | 4.14339 | 0.00332226 | survives matched-quality gate |
| 11 | wise_w1_w2_only | 0.819913 | 2.34925 | 4.2001 | 0.00332226 | survives matched-quality gate |
| 12 | wise_w1_only | 0.963743 | 0.496725 | 4.64187 | 0.00332226 | survives matched-quality gate |

## Full Summary

| fit_group | amp_ratio_vs_baseline | max_direction_sep_vs_baseline_deg | snr_max | max_abs_smd_after | max_weighted_ks_after | readout |
| --- | --- | --- | --- | --- | --- | --- |
| zerr_only | 1.01004 | 0.475513 | 4.8135 | 0.000326864 | 0.00544453 | survives matched-quality gate |
| gaia_g_only | 1.00063 | 0.416871 | 4.77949 | 0.000388115 | 0.00705779 | survives matched-quality gate |
| gaia_bp_rp_only | 0.626609 | 10.9648 | 3.11152 | 3.60446e-05 | 0.00446648 | does not pass matched-quality gate |
| gaia_g_plus_bp_rp | 0.62955 | 10.7436 | 3.12893 | 0.000273522 | 0.00706574 | does not pass matched-quality gate |
| wise_w1_only | 0.963743 | 0.496725 | 4.64187 | 7.49129e-05 | 0.00533162 | survives matched-quality gate |
| wise_w2_only | 1.00693 | 0.398399 | 4.79055 | 0.000301626 | 0.00534204 | survives matched-quality gate |
| wise_w1_w2_only | 0.819913 | 2.34925 | 4.2001 | 0.000116845 | 0.00457839 | survives matched-quality gate |
| wise_w1_plus_w2 | 0.717794 | 3.18802 | 3.85253 | 0.000292097 | 0.0054702 | does not pass matched-quality gate |
| wise_w1_w2_plus_mags | 0.724522 | 3.86062 | 3.77471 | 0.000234986 | 0.00529918 | does not pass matched-quality gate |
| gaia_wise_colour_cross | 0.508139 | 18.5926 | 2.78799 | 0.000338307 | 0.00518034 | does not pass matched-quality gate |
| zerr_x_gaia_colour | 0.639289 | 10.8356 | 3.16836 | 0.000305935 | 0.00553668 | does not pass matched-quality gate |
| zerr_x_wise_colour | 0.809327 | 2.77263 | 4.14339 | 0.000266252 | 0.00580446 | survives matched-quality gate |
| zerr_x_gaia_wise_colour | 0.527275 | 17.9618 | 2.90093 | 0.000336554 | 0.00568575 | does not pass matched-quality gate |
| colour_redshift_estimator_compact | 0.41956 | 27.5163 | 2.45191 | 0.000333639 | 0.00629225 | does not pass matched-quality gate |
| pm_error_quality | 1.02777 | 2.59327 | 4.85898 | 0.000323531 | 0.00887279 | survives matched-quality gate |
| all_catalogue_quality | 0.309709 | 121.39 | 1.87638 | 0.000365446 | 0.00760149 | does not pass matched-quality gate |

## Matched-Feature Shuffle P-Values

| fit_group | n_mocks | p_snr_max_ge_observed | p_pair_sep_le_observed | p_coherence_ge_observed | p_joint_snr_and_pair_sep |
| --- | --- | --- | --- | --- | --- |
| zerr_only | 300 | 0.00332226 | 0.00332226 | 0.00332226 | 0.00332226 |
| gaia_g_only | 300 | 0.00332226 | 0.00332226 | 0.00332226 | 0.00332226 |
| gaia_bp_rp_only | 300 | 1 | 0.00332226 | 0.820598 | 0.00332226 |
| gaia_g_plus_bp_rp | 300 | 1 | 0.00332226 | 0.76412 | 0.00332226 |
| wise_w1_only | 300 | 0.0498339 | 0.00332226 | 0.00332226 | 0.00332226 |
| wise_w2_only | 300 | 0.00332226 | 0.00332226 | 0.00332226 | 0.00332226 |
| wise_w1_w2_only | 300 | 0.0199336 | 0.00332226 | 0.00332226 | 0.00332226 |
| wise_w1_plus_w2 | 300 | 0.0996678 | 0.00332226 | 0.00332226 | 0.00332226 |
| wise_w1_w2_plus_mags | 300 | 0.13289 | 0.00332226 | 0.00332226 | 0.00332226 |
| gaia_wise_colour_cross | 300 | 1 | 0.00332226 | 1 | 0.00332226 |
| zerr_x_gaia_colour | 300 | 1 | 0.00332226 | 0.581395 | 0.00332226 |
| zerr_x_wise_colour | 300 | 0.0365449 | 0.00332226 | 0.00332226 | 0.00332226 |
| zerr_x_gaia_wise_colour | 300 | 1 | 0.00332226 | 0.996678 | 0.00332226 |
| colour_redshift_estimator_compact | 300 | 1 | 0.013289 | 1 | 0.013289 |
| pm_error_quality | 300 | 0.00332226 | 0.00332226 | 0.00332226 | 0.00332226 |
| all_catalogue_quality | 300 | 1 | 0.990033 | 1 | 0.990033 |

## Interpretation Rule

A feature block is treated as a plausible catalogue explanation when it both balances successfully and strongly suppresses the mode. In practice, the most informative rows are those with low amplitude ratio, low SNR, and large direction rotation after good post-weight balance.

## Outputs

- ranked suppressors: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_colour_decomposition_gate_20260715\su2_quaia_colour_decomposition_ranked_suppressors.csv`
- summary: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_colour_decomposition_gate_20260715\su2_quaia_colour_decomposition_summary.csv`
- fit rows: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_colour_decomposition_gate_20260715\su2_quaia_colour_decomposition_fit_rows.csv`
- balance diagnostics: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_colour_decomposition_gate_20260715\su2_quaia_colour_decomposition_balance.csv`
- p-values: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_colour_decomposition_gate_20260715\su2_quaia_colour_decomposition_pvalues.csv`
