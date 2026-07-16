# FR c(z)=c0(1+z) Alpha=0 Calibration Ladder vs LCDM

Date: 2026-07-16 20:19:15

## Test

FR no-loss model:

`DL_FR = (c/H0) z (1 + z/2)`

LCDM comparison models:

- `LCDM_Om03_H0fixed675`: flat LCDM with Om fixed at 0.3.
- `LCDM_Omfree_H0fixed675`: flat LCDM with Om fitted.

H0 fixed at `67.5` km/s/Mpc for all rows. Calibration offsets and priors are identical within each comparison cell.

## Loaded Blocks

| frame             | dataset           | z_col   |   N_used | cov_note                                 | subset_note                |
|:------------------|:------------------|:--------|---------:|:-----------------------------------------|:---------------------------|
| HD                | PantheonPlusSH0ES | zHD     |     1580 | Pantheon-style covariance with leading N | covariance subset inverted |
| HD                | DES_Dovekie_raw   | zHD     |     1820 | DES packed inverse covariance            | precision matrix           |
| HD                | Union3p1_UNITY1p8 | z       |       22 | Union compressed inverse covariance      | precision matrix           |
| HEL               | PantheonPlusSH0ES | zHEL    |     1576 | Pantheon-style covariance with leading N | covariance subset inverted |
| HEL               | DES_Dovekie_raw   | zHEL    |     1820 | DES packed inverse covariance            | precision matrix           |
| HEL               | Union3p1_UNITY1p8 | z       |       22 | Union compressed inverse covariance      | precision matrix           |
| CMB_PANTHEON_ONLY | PantheonPlusSH0ES | zCMB    |     1578 | Pantheon-style covariance with leading N | covariance subset inverted |
| CMB_PANTHEON_ONLY | DES_Dovekie_raw   | zHD     |     1820 | DES packed inverse covariance            | precision matrix           |
| CMB_PANTHEON_ONLY | Union3p1_UNITY1p8 | z       |       22 | Union compressed inverse covariance      | precision matrix           |

## FR Alpha=0 Calibration-Budget Ladder

| prior_config                |   dataset_sigma_mag |   survey_sigma_mag |   data_chi2_map |   prior_chi2 |   posterior_objective |   marginal_score |   BIC_base |   profiled_rms_mag |   max_abs_offset_pull_sigma |
|:----------------------------|--------------------:|-------------------:|----------------:|-------------:|----------------------:|-----------------:|-----------:|-------------------:|----------------------------:|
| budget_010mmag_ds020        |               0.02  |              0.01  |         3512.64 |     243.061  |               3755.7  |          3765.78 |    3765.78 |           0.229474 |                     7.7551  |
| budget_016mmag_ds032        |               0.032 |              0.016 |         3300.16 |     228.28   |               3528.44 |          3545.7  |    3545.7  |           0.226778 |                     6.69001 |
| budget_025mmag_ds050        |               0.05  |              0.025 |         3180.98 |     168.723  |               3349.71 |          3376.9  |    3376.9  |           0.224997 |                     5.45004 |
| budget_050mmag_ds100_stress |               0.1   |              0.05  |         3108.5  |      75.5105 |               3184.01 |          3232.23 |    3232.23 |           0.223504 |                     3.13713 |

## Direct HD Model Comparison

| prior_config                |    N |   FR_BIC_base |   LCDMfree_Om |   LCDMfree_BIC_base |   delta_BIC_FR_minus_LCDMfree |   delta_posterior_FR_minus_LCDMfree |   delta_BIC_FR_minus_LCDM03 |
|:----------------------------|-----:|--------------:|--------------:|--------------------:|------------------------------:|------------------------------------:|----------------------------:|
| budget_010mmag_ds020        | 3422 |       3765.78 |      0.529011 |             3887.99 |                    -122.207   |                          -114.069   |                    -763.496 |
| budget_016mmag_ds032        | 3422 |       3545.7  |      0.491427 |             3650.47 |                    -104.768   |                           -96.6301  |                    -415.767 |
| budget_025mmag_ds050        | 3422 |       3376.9  |      0.432642 |             3448.05 |                     -71.1529  |                           -63.015   |                    -185.514 |
| budget_050mmag_ds100_stress | 3422 |       3232.23 |      0.358069 |             3233.99 |                      -1.76276 |                             6.37522 |                     -14.153 |

Negative `delta_BIC_FR_minus_LCDMfree` favors FR; positive favors LCDM.

## Best Model By Cell

| frame             | offset_mode      | prior_config                | model                     |   param_value |   BIC_base |   posterior_objective |   marginal_score |
|:------------------|:-----------------|:----------------------------|:--------------------------|--------------:|-----------:|----------------------:|-----------------:|
| CMB_PANTHEON_ONLY | dataset+idsurvey | budget_010mmag_ds020        | FR_C1PZ_ALPHA0_H0fixed675 |    nan        |    3704.27 |               3694.2  |          3704.27 |
| CMB_PANTHEON_ONLY | dataset+idsurvey | budget_016mmag_ds032        | FR_C1PZ_ALPHA0_H0fixed675 |    nan        |    3512.1  |               3494.84 |          3512.1  |
| CMB_PANTHEON_ONLY | dataset+idsurvey | budget_025mmag_ds050        | FR_C1PZ_ALPHA0_H0fixed675 |    nan        |    3367.95 |               3340.77 |          3367.95 |
| CMB_PANTHEON_ONLY | dataset+idsurvey | budget_050mmag_ds100_stress | LCDM_Omfree_H0fixed675    |      0.357747 |    3248.67 |               3192.33 |          3240.53 |
| CMB_PANTHEON_ONLY | none             | none                        | FR_C1PZ_ALPHA0_H0fixed675 |    nan        |    4115.91 |               4115.91 |          4115.91 |
| HD                | dataset+idsurvey | budget_010mmag_ds020        | FR_C1PZ_ALPHA0_H0fixed675 |    nan        |    3765.78 |               3755.7  |          3765.78 |
| HD                | dataset+idsurvey | budget_016mmag_ds032        | FR_C1PZ_ALPHA0_H0fixed675 |    nan        |    3545.7  |               3528.44 |          3545.7  |
| HD                | dataset+idsurvey | budget_025mmag_ds050        | FR_C1PZ_ALPHA0_H0fixed675 |    nan        |    3376.9  |               3349.71 |          3376.9  |
| HD                | dataset+idsurvey | budget_050mmag_ds100_stress | FR_C1PZ_ALPHA0_H0fixed675 |    nan        |    3232.23 |               3184.01 |          3232.23 |
| HD                | none             | none                        | FR_C1PZ_ALPHA0_H0fixed675 |    nan        |    4227.9  |               4227.9  |          4227.9  |
| HEL               | dataset+idsurvey | budget_010mmag_ds020        | FR_C1PZ_ALPHA0_H0fixed675 |    nan        |    3733.5  |               3723.42 |          3733.5  |
| HEL               | dataset+idsurvey | budget_016mmag_ds032        | FR_C1PZ_ALPHA0_H0fixed675 |    nan        |    3528.72 |               3511.48 |          3528.72 |
| HEL               | dataset+idsurvey | budget_025mmag_ds050        | FR_C1PZ_ALPHA0_H0fixed675 |    nan        |    3377.43 |               3350.26 |          3377.43 |
| HEL               | dataset+idsurvey | budget_050mmag_ds100_stress | FR_C1PZ_ALPHA0_H0fixed675 |    nan        |    3255.43 |               3207.25 |          3255.43 |
| HEL               | none             | none                        | FR_C1PZ_ALPHA0_H0fixed675 |    nan        |    4173.71 |               4173.71 |          4173.71 |

## Largest HD Offset Pulls

| prior_config         | offset_label                     | offset_type   |   prior_sigma_mag |   offset_mag |   pull_sigma |
|:---------------------|:---------------------------------|:--------------|------------------:|-------------:|-------------:|
| budget_010mmag_ds020 | PantheonPlusSH0ES:IDSURVEY_1     | survey        |             0.01  |   -0.077551  |     -7.7551  |
| budget_010mmag_ds020 | PantheonPlusSH0ES:IDSURVEY_150   | survey        |             0.01  |   -0.068513  |     -6.8513  |
| budget_016mmag_ds032 | PantheonPlusSH0ES:IDSURVEY_150   | survey        |             0.016 |   -0.10704   |     -6.69001 |
| budget_016mmag_ds032 | PantheonPlusSH0ES:IDSURVEY_1     | survey        |             0.016 |   -0.103768  |     -6.48549 |
| budget_010mmag_ds020 | PantheonPlusSH0ES:IDSURVEY_10    | survey        |             0.01  |   -0.0603811 |     -6.03811 |
| budget_010mmag_ds020 | PantheonPlusSH0ES:IDSURVEY_15    | survey        |             0.01  |   -0.0571432 |     -5.71432 |
| budget_025mmag_ds050 | PantheonPlusSH0ES:IDSURVEY_150   | survey        |             0.025 |   -0.136251  |     -5.45004 |
| budget_016mmag_ds032 | PantheonPlusSH0ES:IDSURVEY_10    | survey        |             0.016 |   -0.0869401 |     -5.43376 |
| budget_016mmag_ds032 | PantheonPlusSH0ES:IDSURVEY_15    | survey        |             0.016 |   -0.0774479 |     -4.84049 |
| budget_025mmag_ds050 | PantheonPlusSH0ES:IDSURVEY_1     | survey        |             0.025 |   -0.116069  |     -4.64278 |
| budget_016mmag_ds032 | PantheonPlusSH0ES:IDSURVEY_OTHER | survey        |             0.016 |   -0.068905  |     -4.30656 |
| budget_025mmag_ds050 | PantheonPlusSH0ES:IDSURVEY_10    | survey        |             0.025 |   -0.101736  |     -4.06943 |
| budget_025mmag_ds050 | PantheonPlusSH0ES:IDSURVEY_OTHER | survey        |             0.025 |   -0.101549  |     -4.06197 |
| budget_025mmag_ds050 | PantheonPlusSH0ES:IDSURVEY_64    | survey        |             0.025 |   -0.0972943 |     -3.89177 |
| budget_010mmag_ds020 | PantheonPlusSH0ES:IDSURVEY_OTHER | survey        |             0.01  |   -0.0384975 |     -3.84975 |
| budget_016mmag_ds032 | PantheonPlusSH0ES:IDSURVEY_64    | survey        |             0.016 |   -0.0613233 |     -3.8327  |
| budget_016mmag_ds032 | PantheonPlusSH0ES:IDSURVEY_57    | survey        |             0.016 |   -0.0579055 |     -3.61909 |
| budget_016mmag_ds032 | PantheonPlusSH0ES:IDSURVEY_5     | survey        |             0.016 |   -0.0565446 |     -3.53404 |
| budget_025mmag_ds050 | PantheonPlusSH0ES:IDSURVEY_15    | survey        |             0.025 |   -0.0873392 |     -3.49357 |
| budget_025mmag_ds050 | PantheonPlusSH0ES:IDSURVEY_57    | survey        |             0.025 |   -0.0870332 |     -3.48133 |

## Readout Boundary

This is still a promotion-gate diagnostic. A stable FR alpha=0 preference over added redshift-dimming factors does not by itself establish that FR beats LCDM. The decisive comparison here is the identical-offset FR-vs-LCDM BIC/base-score table, interpreted together with the calibration-budget dependence and offset pulls.

## Output Files

- `fr_c1pz_noloss_lcdm_ladder_summary.csv`
- `fr_c1pz_noloss_lcdm_ladder_comparison.csv`
- `fr_c1pz_noloss_lcdm_ladder_offsets.csv`
- `fr_c1pz_noloss_lcdm_ladder_blocks.csv`
- `fr_c1pz_noloss_lcdm_ladder_config.json`
