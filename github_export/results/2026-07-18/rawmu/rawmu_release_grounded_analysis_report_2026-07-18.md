# Release-Grounded Raw-MU FR versus Lambda-CDM Analysis

Generated: 2026-07-18 20:08:33

## Primary Result

The primary comparison uses released `STAT+SYS` covariance, `zHD` for Pantheon+ and DES, and one flat profiled release intercept. FR and Lambda-CDM use the identical nuisance parameterisation, priors and rows in every cell.

- N = 3422
- LCDM Omega_m = 0.330680
- Delta BIC (FR - LCDM) = 94.344983
- Same-nuisance verified = True

Negative Delta BIC favours FR; positive Delta BIC favours Lambda-CDM.

## Covariance Sensitivities

| analysis_variant                                     | evidence_status                       |    N |   LCDM_Omega_m |   delta_BIC_FR_minus_LCDM | same_nuisance_verified   |
|:-----------------------------------------------------|:--------------------------------------|-----:|---------------:|--------------------------:|:-------------------------|
| released_total_primary                               | primary                               | 3422 |       0.33068  |                   94.345  | True                     |
| pantheon_calib_statonly_grouped                      | covariance_reconstruction_sensitivity | 3422 |       0.336221 |                   98.1247 | True                     |
| pantheon_calib_csp_recal_hstcalspec_statonly_grouped | covariance_reconstruction_sensitivity | 3422 |       0.3357   |                   96.6008 | True                     |
| des_allsys_statonly_grouped                          | covariance_reconstruction_sensitivity | 3422 |       0.32655  |                   86.3631 | True                     |
| combined_calib_des_statonly_grouped                  | covariance_reconstruction_sensitivity | 3422 |       0.332351 |                   89.702  | True                     |

The stat-only grouped-offset rows are covariance-reconstruction sensitivities. They are correlated with the primary release likelihood and are not independent evidence.

### Decisive Reconstruction Limitation

Dataset/survey common modes leave 94.0% of the Pantheon+ `CALIB` component Frobenius norm and 99.9% of the DES all-systematics Frobenius norm unexplained. They therefore cannot replace the full release covariance. All headline comparisons and locked hold-outs use released total covariance; grouped modes remain sensitivity diagnostics only.

## Dataset Hold-Outs

| holdout           |   N_train |   N_test |   LCDM_Omega_m_train |   delta_BIC_train_FR_minus_LCDM |   delta_chi2_test_FR_minus_LCDM | strong_opposite_test   |
|:------------------|----------:|---------:|---------------------:|--------------------------------:|--------------------------------:|:-----------------------|
| PantheonPlusSH0ES |      1842 |     1580 |             0.330199 |                         45.4027 |                         49.5596 | False                  |
| DES_Dovekie_raw   |      1602 |     1820 |             0.33246  |                         59.6923 |                         35.3832 | False                  |
| Union3p1_UNITY1p8 |      3400 |       22 |             0.329929 |                         76.8622 |                         17.4842 | False                  |

## High-Redshift Hold-Outs

| holdout   |   N_train |   N_test |   LCDM_Omega_m_train |   delta_chi2_test_FR_minus_LCDM | test_sign_matches_full   |
|:----------|----------:|---------:|---------------------:|--------------------------------:|:-------------------------|
| z>=0.50   |      2337 |     1085 |             0.345514 |                         33.0018 | True                     |
| z>=0.80   |      3255 |      167 |             0.33517  |                         17.1754 | True                     |
| z>=1.00   |      3376 |       46 |             0.333856 |                         12.7732 | True                     |

## Survey Hold-Out Gate

Eligible survey hold-outs: 17.
Sign-preserving fraction: 0.588235.
Strong opposite survey tests: 0.

## External-Budget Modes

| dataset           | group        |   N |   budget_sigma_mag |   delta_mode_FR_minus_LCDM_mag |   abs_delta_over_budget | within_budget   |
|:------------------|:-------------|----:|-------------------:|-------------------------------:|------------------------:|:----------------|
| PantheonPlusSH0ES | IDSURVEY_64  |  53 |         0.00321523 |                     -0.03681   |                11.4486  | False           |
| PantheonPlusSH0ES | IDSURVEY_57  |  85 |         0.00383506 |                     -0.0407742 |                10.632   | False           |
| PantheonPlusSH0ES | IDSURVEY_63  |  28 |         0.00560802 |                     -0.038889  |                 6.93453 | False           |
| PantheonPlusSH0ES | IDSURVEY_65  |  35 |         0.00597239 |                     -0.0367685 |                 6.15642 | False           |
| PantheonPlusSH0ES | IDSURVEY_10  | 203 |         0.00693693 |                      0.0365831 |                 5.27368 | False           |
| PantheonPlusSH0ES | IDSURVEY_15  | 269 |         0.00637122 |                      0.0298566 |                 4.68617 | False           |
| DES_Dovekie_raw   | IDSURVEY_64  |  31 |         0.0124574  |                     -0.0557632 |                 4.47631 | False           |
| PantheonPlusSH0ES | IDSURVEY_51  |  38 |         0.00960541 |                     -0.0395037 |                 4.11266 | False           |
| DES_Dovekie_raw   | IDSURVEY_150 | 117 |         0.0196764  |                     -0.0548801 |                 2.78913 | False           |
| DES_Dovekie_raw   | IDSURVEY_65  |  22 |         0.0230639  |                     -0.0575033 |                 2.49322 | False           |

The budget gate applies to the FR-minus-Lambda-CDM change in survey residual mode, not to the arbitrary release normalisation intercept.

## Robustness Gate

- `full_and_primary_high_z_same_sign`: **PASS**
- `no_strong_opposite_dataset_holdout`: **PASS**
- `survey_sign_fraction_at_least_0p8`: **FAIL**
- `no_strong_opposite_survey_holdout`: **PASS**
- `all_model_difference_modes_within_budget`: **FAIL**
- `same_nuisance_all_cells`: **PASS**
- `no_covariance_prior_double_counting`: **PASS**

Overall robustness gate: **FAIL**.

A failed gate blocks promotion of a full-sample model preference.

## Source Metadata

- [DES_SN5YR_release](https://github.com/des-science/DES-SN5YR)
- [PantheonPlusSH0ES_CALIB_raw](https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/sytematic_groupings/Pantheon%2BSH0ES_122221_CALIB.cov)
- [PantheonPlusSH0ES_CSP_RECAL_raw](https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/sytematic_groupings/Pantheon%2BSH0ES_122221_CSP_RECAL.cov)
- [PantheonPlusSH0ES_HSTCALSPEC_raw](https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/sytematic_groupings/Pantheon%2BSH0ES_122221_HSTCALSPEC.cov)
- [PantheonPlusSH0ES_STATONLY_raw](https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES_STATONLY.cov)
- [PantheonPlusSH0ES_release](https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR)
- [PantheonPlusSH0ES_systematic_groupings](https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/sytematic_groupings)
- [PantheonPlus_calibration_paper](https://arxiv.org/abs/2110.03486)
- [Union3_paper](https://doi.org/10.3847/1538-4357/adc0a5)
