# SU2 / Quaia Injection-Calibrated Angular Likelihood

Date: 18 July 2026

## Bottom Line

The preregistered conjunctive outcome is **calibration failure and underpowered**. The untouched null false-positive rate exceeded 1%, and the locked 1x effect was not jointly recovered in at least 80% of validation injections. Both are calibration outcomes: the threshold remains unchanged and the observed candidate is not promotable.

## Locked Models

```text
M_T: z_i = alpha + T_i beta + epsilon_i
M_S: z_i = alpha + T_i beta + d . [standardise(Delta mu_SU2R-LCDM(z_i)) n_i] + epsilon_i
```

`T_i` jointly contains SFD dust, Quaia and WISE depth, Gaia source-density and scanning-geometry fields, WISE/Gaia colour and catalogue-quality variables. Twelve fixed Galactic-longitude sky sectors are held out in turn.

## Frozen Calibration

- threshold-calibration null replicates: `1000`
- untouched validation null replicates: `1000`
- frozen 99th-percentile Delta CVLPD threshold: `413.749`
- validation-null false-positive rate: `0.014`
- joint held-out/direction/amplitude recovery at 1x: `0`

## Observed Held-Out Result

- pooled Delta CVLPD: `-11653.4`
- empirical validation-null p-value: `1`
- positive held-out sectors: `0/12`
- pooled recovered amplitude: `0.00236282`
- amplitude ratio to locked A_ref: `0.92707`
- direction separation from locked axis: `113.615 deg`

## Conjunctive Rule

| condition | pass |
| --- | --- |
| validation_null_false_positive_rate_le_0p01 | False |
| joint_1x_recovery_ge_0p80 | False |
| observed_pooled_held_out_gain_positive_and_above_threshold | False |
| observed_direction_within_30_deg | False |
| observed_amplitude_ratio_in_0p5_to_1p5 | True |
| observed_orientation_positive | False |
| at_least_8_of_12_positive_held_out_sectors | False |

No in-sample BIC result is used in this rule. The reported BIC difference is diagnostic only and cannot substitute for positive held-out gain.

## Power Calibration

| amplitude_scale | partition | n_replicates | predictive_detection_rate | joint_recovery_rate | delta_cvlpd_median | amplitude_ratio_median | direction_separation_median_deg |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | threshold_calibration | 1000 | 0.009 |  | -320.984 |  |  |
| 0 | untouched_validation | 1000 | 0.014 |  | -320.053 |  |  |
| 0.5 | threshold_calibration | 1000 | 0.014 | 0 | -320.39 | 7.14994 | 79.6717 |
| 0.5 | untouched_validation | 1000 | 0.016 | 0 | -318.257 | 7.28787 | 80.4684 |
| 1 | threshold_calibration | 1000 | 0.021 | 0 | -310.142 | 3.60803 | 70.8498 |
| 1 | untouched_validation | 1000 | 0.016 | 0 | -313.267 | 3.70675 | 71.733 |
| 1.5 | threshold_calibration | 1000 | 0.033 | 0 | -293.026 | 2.46664 | 62.5409 |
| 1.5 | untouched_validation | 1000 | 0.023 | 0 | -300.124 | 2.57817 | 62.854 |
| 2 | threshold_calibration | 1000 | 0.042 | 0 | -273.435 | 1.92662 | 54.6093 |
| 2 | untouched_validation | 1000 | 0.03 | 0 | -277.322 | 2.03721 | 55.8081 |

A non-zero replicate counts as recovered only when its held-out statistic exceeds the frozen threshold, its recovered direction is within 30 degrees, its amplitude ratio is in [0.5, 1.5], and its orientation is positive.

## Held-Out Sectors

| fold | N_test | delta_lpd | delta_sse | recovered_amp | recovered_l_deg | recovered_b_deg |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 11070 | -477.41 | 11.0649 | 0.0121176 | 195.609 | 8.50855 |
| 1 | 18786 | -1614.83 | 37.105 | 0.0215839 | 220.448 | 4.44433 |
| 2 | 12903 | -714.385 | 16.5428 | 0.0143557 | 248.893 | 5.71908 |
| 3 | 18853 | -1493.84 | 34.5154 | 0.0203953 | 279.295 | 5.07261 |
| 4 | 12516 | -611.656 | 14.1899 | 0.0130753 | 306.454 | 8.65941 |
| 5 | 16967 | -1014.98 | 23.553 | 0.0165516 | 335.557 | 3.73217 |
| 6 | 12083 | -438.827 | 10.194 | 0.0103251 | 3.97876 | 8.60761 |
| 7 | 20359 | -1558.01 | 35.991 | 0.0194244 | 39.5101 | 6.30146 |
| 8 | 13347 | -508.093 | 11.8237 | 0.0105905 | 76.0778 | 11.801 |
| 9 | 18811 | -1234.48 | 28.6437 | 0.0175418 | 108.844 | 4.49849 |
| 10 | 12021 | -542.637 | 12.6058 | 0.0122284 | 143.305 | 6.39372 |
| 11 | 18023 | -1444.25 | 33.2763 | 0.0205902 | 166.284 | 6.09123 |

## Template Inventory

The joint fit retained `23` nuisance templates and `185739` catalogue rows. The machine-readable inventory records deterministic omissions and standardisation constants.

## Scope and Limitations

This is a targeted, conditional power audit on the already selected Quaia candidate. The scalar shape comes from existing late-time best fits, and the injection null conditions on the observed catalogue, mask and nuisance fields. A pass would not constitute an independent detection. A failure does not exclude SU2 generally. An underpowered result does not permit the threshold to be moved.

The SU2 scalar curve is evaluated at catalogue redshift and enters only through a centred angular interaction. This avoids adding the scalar curve alone as a circular predictor of redshift, but the result remains a phenomenological association test rather than a first-principles quasar-selection likelihood.
