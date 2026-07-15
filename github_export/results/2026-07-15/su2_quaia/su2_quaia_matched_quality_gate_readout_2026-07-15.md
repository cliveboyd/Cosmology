# SU2 / Quaia Matched Catalogue-Quality Gate Readout

Date: July 15, 2026

## Bottom Line

The locked mode does not pass the all-catalogue-quality matching gate; catalogue-quality remains the leading promotion-blocking caveat.

## Gate Summary

| fit_group | amp_ratio_vs_baseline | max_direction_sep_vs_baseline_deg | snr_max | max_abs_smd_after | max_weighted_ks_after | readout |
| --- | --- | --- | --- | --- | --- | --- |
| zerr_only | 1.01004 | 0.475513 | 4.8135 | 0.000326864 | 0.00544453 | survives matched-quality gate |
| gaia_colour_quality | 0.62955 | 10.7436 | 3.12893 | 0.000273522 | 0.00706574 | does not pass matched-quality gate |
| wise_colour_quality | 0.724522 | 3.86062 | 3.77471 | 0.000234986 | 0.00529918 | does not pass matched-quality gate |
| pm_error_quality | 1.02777 | 2.59327 | 4.85898 | 0.000323531 | 0.00887279 | survives matched-quality gate |
| all_catalogue_quality | 0.309709 | 121.39 | 1.87638 | 0.000365446 | 0.00760149 | does not pass matched-quality gate |

## Shuffle P-Values

| fit_group | n_mocks | p_snr_max_ge_observed | p_pair_sep_le_observed | p_coherence_ge_observed | p_joint_snr_and_pair_sep |
| --- | --- | --- | --- | --- | --- |
| zerr_only | 500 | 0.00199601 | 0.00199601 | 0.00199601 | 0.00199601 |
| gaia_colour_quality | 500 | 1 | 0.00199601 | 0.756487 | 0.00199601 |
| wise_colour_quality | 500 | 0.113772 | 0.00199601 | 0.00199601 | 0.00199601 |
| pm_error_quality | 500 | 0.00199601 | 0.00199601 | 0.00199601 | 0.00199601 |
| all_catalogue_quality | 500 | 1 | 0.986028 | 1 | 0.986028 |

## All-Catalogue Quality Detail

- amplitude ratio vs baseline: `0.309709`
- maximum direction rotation vs baseline: `121.39` deg
- max post-weight absolute SMD: `0.000365446`
- max post-weight weighted KS: `0.00760149`
- matched-quality joint p-value: `0.986028`
