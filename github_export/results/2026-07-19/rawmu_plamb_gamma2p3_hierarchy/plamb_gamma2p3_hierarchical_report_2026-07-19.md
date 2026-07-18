# Hierarchical raw-MU sensitivity: fixed gamma = 2.3

Generated: 2026-07-19 08:58:10

## Claim boundary

This is an outcome-aware sensitivity analysis. `gamma_c=2.3` was selected after inspection of the 18 July supernova profiles, so this run is not independent evidence and cannot promote PLAMB.

## Model

```text
c(z) / c0                    = 1 + 2.3 z
|dz/dT| / H0                = (1 + z)^p
D_L(z)                      = (c0/H0) integral[0,z] (1+2.3u)/(1+u)^p du
alpha                       = 0
```

Pantheon+ and DES use `zHD`; Union3.1 uses its released `z`. Every model has one unpenalised release intercept and one shape parameter, so delta BIC equals delta chi2.

## Calibration treatment

The primary cell uses each release's published `STAT+SYS` covariance. Calibration budgets already represented there are not added again. Sensitivity cells replace the relevant total covariance by `STATONLY + X Lambda X^T`, which is the covariance obtained after analytically marginalising zero-mean Gaussian dataset/survey offsets. Because the grouped reconstructions omit most of the original covariance structure, they are correlated diagnostics rather than alternative headline likelihoods.

## Results

| analysis_variant                                     | evidence_status                                  |   PATH_p |   LCDM_Omega_m |   delta_BIC_PATH_minus_LCDM | same_nuisance_verified   |
|:-----------------------------------------------------|:-------------------------------------------------|---------:|---------------:|----------------------------:|:-------------------------|
| released_total_primary                               | primary_released_total                           | 0.797429 |       0.33068  |                   0.0119585 | True                     |
| pantheon_calib_statonly_grouped                      | correlated_covariance_reconstruction_sensitivity | 0.804876 |       0.336221 |                  -0.816882  | True                     |
| pantheon_calib_csp_recal_hstcalspec_statonly_grouped | correlated_covariance_reconstruction_sensitivity | 0.804414 |       0.3357   |                  -0.740532  | True                     |
| des_allsys_statonly_grouped                          | correlated_covariance_reconstruction_sensitivity | 0.794539 |       0.32655  |                   0.565488  | True                     |
| combined_calib_des_statonly_grouped                  | correlated_covariance_reconstruction_sensitivity | 0.801971 |       0.332351 |                  -0.480082  | True                     |

Primary delta BIC (PATH - Lambda-CDM) = 0.011958.
Sensitivity delta-BIC range = [-0.816882, 0.565488].

## External calibration budgets

| dataset           | group        |   N |   budget_sigma_mag |   delta_mode_PATH_minus_LCDM_mag |   abs_delta_over_budget | within_budget   |
|:------------------|:-------------|----:|-------------------:|---------------------------------:|------------------------:|:----------------|
| PantheonPlusSH0ES | IDSURVEY_64  |  53 |         0.00321523 |                      -0.00220137 |                0.68467  | True            |
| PantheonPlusSH0ES | IDSURVEY_57  |  85 |         0.00383506 |                      -0.00230863 |                0.60198  | True            |
| PantheonPlusSH0ES | IDSURVEY_10  | 203 |         0.00693693 |                       0.00313801 |                0.452362 | True            |
| PantheonPlusSH0ES | IDSURVEY_63  |  28 |         0.00560802 |                      -0.00226327 |                0.403578 | True            |
| PantheonPlusSH0ES | IDSURVEY_15  | 269 |         0.00637122 |                       0.00256693 |                0.402895 | True            |
| PantheonPlusSH0ES | IDSURVEY_65  |  35 |         0.00597239 |                      -0.00220661 |                0.369468 | True            |
| DES_Dovekie_raw   | IDSURVEY_64  |  31 |         0.0124574  |                      -0.00371702 |                0.298379 | True            |
| PantheonPlusSH0ES | IDSURVEY_51  |  38 |         0.00960541 |                      -0.00226223 |                0.235517 | True            |
| DES_Dovekie_raw   | IDSURVEY_150 | 117 |         0.0196764  |                      -0.00369284 |                0.187679 | True            |
| DES_Dovekie_raw   | IDSURVEY_65  |  22 |         0.0230639  |                      -0.00378853 |                0.164262 | True            |

## Gates

- `same_nuisance_all_cells`: **PASS**
- `optimisation_success_all_cells`: **PASS**
- `primary_absolute_delta_BIC_at_most_2`: **PASS**
- `all_covariance_sensitivities_absolute_delta_BIC_at_most_2`: **PASS**
- `all_model_difference_modes_within_external_budget`: **PASS**
- `no_calibration_prior_double_counting`: **PASS**
- `independent_gamma_derivation`: **FAIL**

These gates describe calibration robustness of the post-hoc cell only. They do not repair the failed independent-gamma and high-redshift predictive gates.

## Reconstruction limitations

- `pantheon_calib_statonly_grouped` leaves 94.0% of the target Frobenius norm unexplained.
- `pantheon_calib_csp_recal_hstcalspec_statonly_grouped` leaves 97.9% of the target Frobenius norm unexplained.
- `des_allsys_statonly_grouped` leaves 99.9% of the target Frobenius norm unexplained.
