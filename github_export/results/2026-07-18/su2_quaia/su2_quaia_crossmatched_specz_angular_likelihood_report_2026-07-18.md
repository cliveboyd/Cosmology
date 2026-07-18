# Leakage-Free Quaia x DR16Q Angular Likelihood

Date: 18 July 2026

## Bottom Line

The preregistered outcome is **underpowered**. The leakage-free design did not recover the locked one-times signal with at least 80% power; the threshold remains frozen.

## Independent Response

The response is `z_Quaia - z_spec`, where `z_spec` is the independently measured DR16Q spectroscopic redshift. The SU2R-minus-Lambda-CDM scalar curve is evaluated only at `z_spec`.

- reciprocal one-to-one matches after all masks: `46687`
- maximum accepted separation: `0.524992` arcsec
- retained nuisance templates: `40`
- candidate design rank: `44`
- candidate design condition number: `3164.58`

## Leakage Audit

- exact predictor equality after permuting the Quaia response: `True`
- maximum absolute predictor change: `0`
- response changed under permutation: `True`
- predictor SHA-256: `e6f04c7a9d01f73b6ca0b2c5f27705399633a7fe58b81ce63a0e8a8ca65cfae9`

## Calibration

- frozen 99th-percentile threshold: `3.38511`
- untouched validation-null false-positive rate: `0.009`
- joint recovery at the locked 1x signal: `0`

## Observed Held-Out Result

- pooled Delta CVLPD: `-4.0084`
- empirical validation-null p-value: `0.727273`
- positive sectors: `4/12`
- recovered amplitude ratio: `1.08372`
- recovered direction: `(l,b)=(165.123, 14.9059) deg`
- separation from locked direction: `57.1066` deg

## Conjunctive Rule

| condition | pass |
| --- | --- |
| predictor_permutation_leakage_audit | True |
| validation_null_false_positive_rate_le_0p01 | True |
| joint_1x_recovery_ge_0p80 | False |
| observed_pooled_held_out_gain_positive_and_above_threshold | False |
| observed_direction_amplitude_orientation | False |

## Power Calibration

| amplitude_scale | partition | n_replicates | predictive_detection_rate | joint_recovery_rate | delta_cvlpd_median | amplitude_ratio_median | direction_separation_median_deg |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | threshold_calibration | 1000 | 0.009 |  | -2.58972 |  |  |
| 0 | untouched_validation | 1000 | 0.009 |  | -2.7264 |  |  |
| 0.5 | fresh_injection | 1000 | 0.032 | 0 | -2.20856 | 2.9986 | 64.1497 |
| 1 | fresh_injection | 1000 | 0.08 | 0 | -1.66563 | 1.70051 | 45.4385 |
| 1.5 | fresh_injection | 1000 | 0.203 | 0.016 | 0.340012 | 1.34243 | 36.3315 |
| 2 | fresh_injection | 1000 | 0.481 | 0.154 | 3.00948 | 1.20021 | 29.8867 |

## Held-Out Sectors

| fold | N_test | delta_lpd | delta_sse | recovered_amp | recovered_l_deg | recovered_b_deg |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 1443 | 0.0548888 | -0.00442428 | 0.000784344 | 216.514 | 43.4395 |
| 1 | 3990 | 0.0384755 | -0.00419932 | 0.00276384 | 160.567 | 6.69908 |
| 2 | 5294 | -1.00521 | 0.0693851 | 0.00336975 | 189.025 | 1.54792 |
| 3 | 8500 | -0.591623 | 0.0405093 | 0.0024807 | 204.476 | 13.7664 |
| 4 | 5751 | -1.75458 | 0.121993 | 0.00497216 | 149.083 | 1.60751 |
| 5 | 8357 | -0.253833 | 0.019343 | 0.0029834 | 153.718 | 31.3695 |
| 6 | 5953 | 0.150929 | -0.00760237 | 0.00262567 | 155.288 | 22.8743 |
| 7 | 3990 | -0.0935903 | 0.00617939 | 0.00297989 | 133.068 | 18.8457 |
| 8 | 885 | -0.198514 | 0.0143987 | 0.00220915 | 196.609 | 9.59707 |
| 9 | 877 | 0.027403 | -0.00230331 | 0.00245382 | 166.055 | 16.7102 |
| 10 | 476 | -0.0920262 | 0.00625132 | 0.00296756 | 159.668 | 13.035 |
| 11 | 1171 | -0.290725 | 0.019882 | 0.00339222 | 159.949 | 9.36102 |

## Claim Boundary

This replacement repairs the response leakage in the first injection branch. A pass would establish only a conditional cross-catalogue predictive association. A calibration or power failure is underpowered; an observed failure after adequate calibration rejects only the locked candidate in this sample and nuisance model.
