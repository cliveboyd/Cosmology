# Joint Raw-MU Hierarchical Offset Likelihood

Date: July 14, 2026

Datasets: DES-Dovekie raw, Pantheon+SH0ES, and Union3.1 UNITY1.8 compressed MU nodes.

Covariance treatment: each dataset uses its available full covariance or inverse covariance; dataset blocks are combined as independent block-diagonal precision blocks.

Model: `MU = 5 log10((c/H0) z (1 + beta z)) + 25 + dataset_offset + survey_offset`.

Calibration priors are zero-mean Gaussian random effects. Dataset offsets use `sigma_dataset`; nested IDSURVEY offsets use `sigma_survey`. Union3.1 has no IDSURVEY labels, so it receives only a dataset-level offset in `dataset+survey` mode.

Ranking columns:

- `posterior_objective = data_chi2_at_MAP + prior_chi2`.
- `marginal_score` adds the Gaussian random-effect Occam/log-determinant term, up to a common covariance constant.
- `BIC_random_marginal_base` additionally counts only the base H0/beta parameters, because calibration offsets have been marginalized with explicit priors.

## Configuration

- frames: `HD,HEL`
- random-effect modes: `none,dataset,dataset+survey`
- prior configs: `tight:0.02:0.01,moderate:0.05:0.02,loose:0.10:0.05,very_loose:0.20:0.10`
- z_min: `0.01`
- min IDSURVEY/group size before OTHER pooling: `10`
- H0/beta bounds: H0 `45.0` to `90.0`, beta `-0.25` to `1.25`

## Best Rows by Marginal Random-Effect BIC-Like Score

| frame   | random_mode    | prior_config   | model                  |    N |   k_base |   k_random |   dataset_sigma_mag |   survey_sigma_mag |      H0 |     beta |   data_chi2_map |   prior_chi2 |   marginal_score |   BIC_random_marginal_base |
|:--------|:---------------|:---------------|:-----------------------|-----:|---------:|-----------:|--------------------:|-------------------:|--------:|---------:|----------------:|-------------:|-----------------:|---------------------------:|
| HD      | dataset+survey | moderate       | FR_BETA05_H0fixed675   | 3422 |        0 |         27 |                0.05 |               0.02 | 67.5    | 0.5      |         3097.75 |    19.9073   |          3147.06 |                    3147.06 |
| HD      | dataset+survey | loose          | FR_BETA05_H0fixed675   | 3422 |        0 |         27 |                0.1  |               0.05 | 67.5    | 0.5      |         3082.88 |     9.59394  |          3150.49 |                    3150.49 |
| HD      | dataset+survey | moderate       | FR_BETAfree_H0fixed675 | 3422 |        1 |         27 |                0.05 |               0.02 | 67.5    | 0.521314 |         3095.49 |    18.0114   |          3142.9  |                    3151.04 |
| HD      | dataset+survey | moderate       | FR_BETA05_H0free       | 3422 |        1 |         27 |                0.05 |               0.02 | 69.2264 | 0.5      |         3097.26 |    16.8908   |          3143.55 |                    3151.69 |
| HD      | dataset+survey | moderate       | FR_BETAfree_H0free     | 3422 |        2 |         27 |                0.05 |               0.02 | 69.7141 | 0.526515 |         3095.49 |    12.5296   |          3137.41 |                    3153.69 |
| HD      | dataset+survey | loose          | FR_BETAfree_H0fixed675 | 3422 |        1 |         27 |                0.1  |               0.05 | 67.5    | 0.514357 |         3081.96 |     8.98427  |          3148.96 |                    3157.09 |
| HD      | dataset+survey | loose          | FR_BETA05_H0free       | 3422 |        1 |         27 |                0.1  |               0.05 | 69.3859 | 0.5      |         3082.82 |     8.60836  |          3149.44 |                    3157.58 |
| HD      | dataset        | loose          | FR_BETAfree_H0fixed675 | 3422 |        1 |          3 |                0.1  |               0    | 67.5    | 0.541273 |         3128.86 |     2.26102  |          3150.2  |                    3158.34 |
| HD      | dataset        | very_loose     | FR_BETAfree_H0fixed675 | 3422 |        1 |          3 |                0.2  |               0    | 67.5    | 0.54204  |         3128.85 |     0.5711   |          3152.65 |                    3160.79 |
| HD      | dataset        | moderate       | FR_BETAfree_H0fixed675 | 3422 |        1 |          3 |                0.05 |               0    | 67.5    | 0.53837  |         3129.08 |     8.69616  |          3152.71 |                    3160.85 |
| HD      | dataset        | moderate       | FR_BETAfree_H0free     | 3422 |        2 |          3 |                0.05 |               0    | 69.8057 | 0.543189 |         3128.88 |     2.84341  |          3146.65 |                    3162.93 |
| HD      | dataset+survey | loose          | FR_BETAfree_H0free     | 3422 |        2 |         27 |                0.1  |               0.05 | 69.657  | 0.515809 |         3081.95 |     7.63514  |          3147.6  |                    3163.88 |
| HD      | dataset+survey | tight          | FR_BETAfree_H0free     | 3422 |        2 |         27 |                0.02 |               0.01 | 69.8779 | 0.539498 |         3111.22 |    22.5273   |          3148.58 |                    3164.86 |
| HD      | dataset        | loose          | FR_BETAfree_H0free     | 3422 |        2 |          3 |                0.1  |               0    | 69.7897 | 0.542524 |         3128.85 |     0.721643 |          3148.65 |                    3164.92 |
| HD      | dataset        | very_loose     | FR_BETAfree_H0free     | 3422 |        2 |          3 |                0.2  |               0    | 69.7853 | 0.542356 |         3128.85 |     0.181094 |          3152.26 |                    3168.54 |
| HD      | dataset        | loose          | FR_BETA05_H0fixed675   | 3422 |        0 |          3 |                0.1  |               0    | 67.5    | 0.5      |         3150.15 |     1.53427  |          3170.76 |                    3170.76 |
| HD      | dataset        | moderate       | FR_BETA05_H0fixed675   | 3422 |        0 |          3 |                0.05 |               0    | 67.5    | 0.5      |         3150.18 |     6.08487  |          3171.2  |                    3171.2  |
| HD      | dataset        | tight          | FR_BETAfree_H0free     | 3422 |        2 |          3 |                0.02 |               0    | 69.9106 | 0.54756  |         3129.86 |    16.045    |          3155.45 |                    3171.72 |
| HD      | dataset+survey | very_loose     | FR_BETA05_H0fixed675   | 3422 |        0 |         27 |                0.2  |               0.1  | 67.5    | 0.5      |         3080.23 |     3.56303  |          3172.1  |                    3172.1  |
| HD      | dataset+survey | tight          | FR_BETA05_H0free       | 3422 |        1 |         27 |                0.02 |               0.01 | 69.0874 | 0.5      |         3117.39 |    32.7185   |          3164.94 |                    3173.08 |
| HD      | dataset        | very_loose     | FR_BETA05_H0fixed675   | 3422 |        0 |          3 |                0.2  |               0    | 67.5    | 0.5      |         3150.15 |     0.38439  |          3173.76 |                    3173.76 |
| HEL     | dataset+survey | moderate       | FR_BETA05_H0fixed675   | 3418 |        0 |         27 |                0.05 |               0.02 | 67.5    | 0.5      |         3131.97 |    14.2088   |          3175.56 |                    3175.56 |
| HD      | dataset        | moderate       | FR_BETA05_H0free       | 3422 |        1 |          3 |                0.05 |               0    | 68.861  | 0.5      |         3150.17 |     3.85505  |          3168.96 |                    3177.1  |
| HD      | dataset        | loose          | FR_BETA05_H0free       | 3422 |        1 |          3 |                0.1  |               0    | 68.8594 | 0.5      |         3150.15 |     0.973409 |          3170.2  |                    3178.34 |

## H0=67.5 Beta-Free Focus

| frame   | random_mode    | prior_config   |   k_random |   dataset_sigma_mag |   survey_sigma_mag |     beta |   beta_minus_0p5 |   data_chi2_map |   prior_chi2 |   marginal_score |   BIC_random_marginal_base |
|:--------|:---------------|:---------------|-----------:|--------------------:|-------------------:|---------:|-----------------:|----------------:|-------------:|-----------------:|---------------------------:|
| HD      | dataset+survey | moderate       |         27 |                0.05 |               0.02 | 0.521314 |       0.0213139  |         3095.49 |    18.0114   |          3142.9  |                    3151.04 |
| HD      | dataset+survey | loose          |         27 |                0.1  |               0.05 | 0.514357 |       0.0143566  |         3081.96 |     8.98427  |          3148.96 |                    3157.09 |
| HD      | dataset        | loose          |          3 |                0.1  |               0    | 0.541273 |       0.0412727  |         3128.86 |     2.26102  |          3150.2  |                    3158.34 |
| HD      | dataset        | very_loose     |          3 |                0.2  |               0    | 0.54204  |       0.0420398  |         3128.85 |     0.5711   |          3152.65 |                    3160.79 |
| HD      | dataset        | moderate       |          3 |                0.05 |               0    | 0.53837  |       0.0383696  |         3129.08 |     8.69616  |          3152.71 |                    3160.85 |
| HD      | dataset+survey | very_loose     |         27 |                0.2  |               0.1  | 0.512272 |       0.0122725  |         3079.37 |     3.38404  |          3171.06 |                    3179.2  |
| HD      | dataset+survey | tight          |         27 |                0.02 |               0.01 | 0.514866 |       0.0148657  |         3114.78 |    49.8767   |          3179.49 |                    3187.63 |
| HD      | dataset        | tight          |          3 |                0.02 |               0    | 0.523403 |       0.0234029  |         3134.62 |    43.982    |          3188.15 |                    3196.29 |
| HD      | none           | none           |          0 |                0    |               0    | 0.466532 |      -0.0334679  |         4109.87 |     0        |          4109.87 |                    4118    |
| HEL     | dataset+survey | moderate       |         27 |                0.05 |               0.02 | 0.514769 |       0.0147689  |         3130.4  |    13.7751   |          3173.55 |                    3181.68 |
| HEL     | dataset        | loose          |          3 |                0.1  |               0    | 0.525911 |       0.0259106  |         3155.53 |     1.83924  |          3176.44 |                    3184.58 |
| HEL     | dataset        | moderate       |          3 |                0.05 |               0    | 0.523622 |       0.0236222  |         3155.67 |     7.12216  |          3177.72 |                    3185.86 |
| HEL     | dataset        | very_loose     |          3 |                0.2  |               0    | 0.526515 |       0.0265148  |         3155.52 |     0.463736 |          3179.21 |                    3187.35 |
| HEL     | dataset+survey | loose          |         27 |                0.1  |               0.05 | 0.512536 |       0.0125358  |         3118.94 |     7.55044  |          3184.47 |                    3192.61 |
| HEL     | dataset+survey | tight          |         27 |                0.02 |               0.01 | 0.506924 |       0.00692352 |         3143.95 |    41.675    |          3200.45 |                    3208.59 |
| HEL     | dataset        | tight          |          3 |                0.02 |               0    | 0.511779 |       0.0117788  |         3159.54 |    37.3315   |          3206.41 |                    3214.55 |
| HEL     | dataset+survey | very_loose     |         27 |                0.2  |               0.1  | 0.512091 |       0.012091   |         3116.3  |     3.06174  |          3207.62 |                    3215.76 |
| HEL     | none           | none           |          0 |                0    |               0    | 0.465369 |      -0.034631   |         4047.16 |     0        |          4047.16 |                    4055.3  |

## Dataset Blocks

| frame   | dataset           | z_col   |   N_source |   N_used | cov_note                                 | subset_note                | source_path                                                                                                                                       | cov_path                                                                                                                                          |
|:--------|:------------------|:--------|-----------:|---------:|:-----------------------------------------|:---------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------|
| HD      | PantheonPlusSH0ES | zHD     |       1701 |     1580 | Pantheon-style covariance with leading N | covariance subset inverted | C:\Users\clive\Documents\Cosmology\Pantheon+SH0ES.dat                                                                                             | C:\Users\clive\Documents\Cosmology\Pantheon+SH0ES_STAT+SYS.cov                                                                                    |
| HD      | DES_Dovekie_raw   | zHD     |       1820 |     1820 | DES packed inverse covariance            | precision matrix           | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\DES_SN5YR_Dovekie\DES-Dovekie_HD.csv                               | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\DES_SN5YR_Dovekie\STAT+SYS.npz                                     |
| HD      | Union3p1_UNITY1p8 | z       |         22 |       22 | Union compressed inverse covariance      | precision matrix           | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\Union3_release\mu_mat_union3.1_UNITY1.8_template_cosmo=2_0_mu.fits | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\Union3_release\mu_mat_union3.1_UNITY1.8_template_cosmo=2_0_mu.fits |
| HEL     | PantheonPlusSH0ES | zHEL    |       1701 |     1576 | Pantheon-style covariance with leading N | covariance subset inverted | C:\Users\clive\Documents\Cosmology\Pantheon+SH0ES.dat                                                                                             | C:\Users\clive\Documents\Cosmology\Pantheon+SH0ES_STAT+SYS.cov                                                                                    |
| HEL     | DES_Dovekie_raw   | zHEL    |       1820 |     1820 | DES packed inverse covariance            | precision matrix           | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\DES_SN5YR_Dovekie\DES-Dovekie_HD.csv                               | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\DES_SN5YR_Dovekie\STAT+SYS.npz                                     |
| HEL     | Union3p1_UNITY1p8 | z       |         22 |       22 | Union compressed inverse covariance      | precision matrix           | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\Union3_release\mu_mat_union3.1_UNITY1.8_template_cosmo=2_0_mu.fits | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\Union3_release\mu_mat_union3.1_UNITY1.8_template_cosmo=2_0_mu.fits |

## Largest Offset Pulls

| frame   | random_mode    | prior_config   | model                  | offset_group                     | component_type   |   prior_sigma_mag |   residual_offset_mag |   pull_sigma |
|:--------|:---------------|:---------------|:-----------------------|:---------------------------------|:-----------------|------------------:|----------------------:|-------------:|
| HD      | dataset        | tight          | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET        | dataset          |              0.02 |            -0.127502  |     -6.3751  |
| HD      | dataset+survey | tight          | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET        | dataset          |              0.02 |            -0.124072  |     -6.20359 |
| HEL     | dataset        | tight          | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET        | dataset          |              0.02 |            -0.119597  |     -5.97984 |
| HD      | dataset+survey | tight          | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET        | dataset          |              0.02 |            -0.118788  |     -5.93941 |
| HD      | dataset        | tight          | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET        | dataset          |              0.02 |            -0.118682  |     -5.93408 |
| HEL     | dataset+survey | tight          | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET        | dataset          |              0.02 |            -0.116169  |     -5.80847 |
| HEL     | dataset        | tight          | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET        | dataset          |              0.02 |            -0.115147  |     -5.75734 |
| HEL     | dataset+survey | tight          | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET        | dataset          |              0.02 |            -0.113704  |     -5.68521 |
| HD      | dataset        | tight          | FR_BETA05_H0free       | PantheonPlusSH0ES:DATASET        | dataset          |              0.02 |            -0.0765105 |     -3.82553 |
| HEL     | dataset        | tight          | FR_BETA05_H0free       | PantheonPlusSH0ES:DATASET        | dataset          |              0.02 |            -0.0736176 |     -3.68088 |
| HD      | dataset+survey | tight          | FR_BETA05_H0free       | PantheonPlusSH0ES:DATASET        | dataset          |              0.02 |            -0.0712793 |     -3.56396 |
| HEL     | dataset+survey | tight          | FR_BETA05_H0free       | PantheonPlusSH0ES:DATASET        | dataset          |              0.02 |            -0.0696229 |     -3.48114 |
| HEL     | dataset        | tight          | FR_BETAfree_H0free     | PantheonPlusSH0ES:DATASET        | dataset          |              0.02 |            -0.0644057 |     -3.22028 |
| HEL     | dataset+survey | tight          | FR_BETAfree_H0free     | PantheonPlusSH0ES:DATASET        | dataset          |              0.02 |            -0.0630783 |     -3.15391 |
| HD      | dataset        | tight          | FR_BETAfree_H0free     | PantheonPlusSH0ES:DATASET        | dataset          |              0.02 |            -0.062958  |     -3.1479  |
| HD      | dataset+survey | tight          | FR_BETAfree_H0free     | PantheonPlusSH0ES:DATASET        | dataset          |              0.02 |            -0.0620036 |     -3.10018 |
| HD      | dataset        | moderate       | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET        | dataset          |              0.05 |            -0.137064  |     -2.74127 |
| HD      | dataset+survey | moderate       | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET        | dataset          |              0.05 |            -0.135361  |     -2.70721 |
| HD      | dataset+survey | moderate       | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET        | dataset          |              0.05 |            -0.128002  |     -2.56003 |
| HEL     | dataset        | moderate       | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET        | dataset          |              0.05 |            -0.12776   |     -2.5552  |
| HEL     | dataset+survey | moderate       | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET        | dataset          |              0.05 |            -0.125679  |     -2.51359 |
| HD      | dataset        | moderate       | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET        | dataset          |              0.05 |            -0.122209  |     -2.44417 |
| HD      | dataset        | tight          | FR_BETA05_H0free       | DES_Dovekie_raw:DATASET          | dataset          |              0.02 |             0.0487627 |      2.43814 |
| HEL     | dataset+survey | moderate       | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET        | dataset          |              0.05 |            -0.120573  |     -2.41145 |
| HEL     | dataset        | moderate       | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET        | dataset          |              0.05 |            -0.118586  |     -2.37173 |
| HD      | dataset        | tight          | FR_BETAfree_H0free     | DES_Dovekie_raw:DATASET          | dataset          |              0.02 |             0.0468484 |      2.34242 |
| HD      | dataset+survey | tight          | FR_BETA05_H0free       | DES_Dovekie_raw:IDSURVEY_10      | survey           |              0.01 |             0.0232886 |      2.32886 |
| HEL     | dataset        | tight          | FR_BETA05_H0free       | DES_Dovekie_raw:DATASET          | dataset          |              0.02 |             0.0465016 |      2.32508 |
| HEL     | dataset        | tight          | FR_BETAfree_H0free     | DES_Dovekie_raw:DATASET          | dataset          |              0.02 |             0.0451952 |      2.25976 |
| HEL     | dataset+survey | tight          | FR_BETAfree_H0free     | DES_Dovekie_raw:DATASET          | dataset          |              0.02 |             0.040722  |      2.0361  |
| HD      | dataset+survey | tight          | FR_BETAfree_H0free     | DES_Dovekie_raw:DATASET          | dataset          |              0.02 |             0.0396794 |      1.98397 |
| HEL     | dataset+survey | tight          | FR_BETA05_H0free       | DES_Dovekie_raw:DATASET          | dataset          |              0.02 |             0.0388275 |      1.94137 |
| HD      | dataset+survey | tight          | FR_BETA05_H0free       | DES_Dovekie_raw:DATASET          | dataset          |              0.02 |             0.0369918 |      1.84959 |
| HD      | dataset        | tight          | FR_BETAfree_H0fixed675 | Union3p1_UNITY1p8:DATASET        | dataset          |              0.02 |            -0.0354955 |     -1.77477 |
| HEL     | dataset+survey | loose          | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:IDSURVEY_OTHER | survey           |              0.05 |            -0.0887351 |     -1.7747  |
| HEL     | dataset+survey | loose          | FR_BETAfree_H0free     | PantheonPlusSH0ES:IDSURVEY_OTHER | survey           |              0.05 |            -0.0883723 |     -1.76745 |
| HD      | dataset+survey | loose          | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:IDSURVEY_OTHER | survey           |              0.05 |            -0.0866727 |     -1.73345 |
| HD      | dataset+survey | tight          | FR_BETA05_H0fixed675   | DES_Dovekie_raw:IDSURVEY_10      | survey           |              0.01 |             0.0173174 |      1.73174 |
| HD      | dataset+survey | loose          | FR_BETAfree_H0free     | PantheonPlusSH0ES:IDSURVEY_OTHER | survey           |              0.05 |            -0.0862528 |     -1.72506 |
| HEL     | dataset+survey | loose          | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:IDSURVEY_OTHER | survey           |              0.05 |            -0.0859156 |     -1.71831 |

## Interpretation Notes

- If beta near 0.5 survives only for loose calibration priors, the raw-MU signal remains calibration-dependent.
- If dataset+survey priors improve data chi2 but are disfavored by the marginal random-effect score, the free-offset improvement is mostly Occam-penalized survey structure.
- Large offset pulls identify which calibration blocks need direct audit before any model-promotion claim.