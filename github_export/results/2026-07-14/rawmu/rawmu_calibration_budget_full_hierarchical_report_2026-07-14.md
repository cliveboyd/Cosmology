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

- frames: `HD,HEL,CMB_PANTHEON_ONLY`
- random-effect modes: `none,dataset,dataset+survey`
- prior configs: `budget_010mmag_ds020:0.020:0.010,budget_016mmag_ds032:0.032:0.016,budget_025mmag_ds050:0.050:0.025,budget_050mmag_ds100_stress:0.100:0.050`
- z_min: `0.01`
- min IDSURVEY/group size before OTHER pooling: `10`
- H0/beta bounds: H0 `45.0` to `90.0`, beta `-0.25` to `1.25`

## Best Rows by Marginal Random-Effect BIC-Like Score

| frame             | random_mode    | prior_config                | model                  |    N |   k_base |   k_random |   dataset_sigma_mag |   survey_sigma_mag |      H0 |     beta |   data_chi2_map |   prior_chi2 |   marginal_score |   BIC_random_marginal_base |
|:------------------|:---------------|:----------------------------|:-----------------------|-----:|---------:|-----------:|--------------------:|-------------------:|--------:|---------:|----------------:|-------------:|-----------------:|---------------------------:|
| HD                | dataset+survey | budget_025mmag_ds050        | FR_BETA05_H0fixed675   | 3422 |        0 |         27 |               0.05  |              0.025 | 67.5    | 0.5      |         3092.82 |    19.1536   |          3146.03 |                    3146.03 |
| HD                | dataset+survey | budget_025mmag_ds050        | FR_BETA05_H0free       | 3422 |        1 |         27 |               0.05  |              0.025 | 69.2777 | 0.5      |         3092.36 |    15.929    |          3142.34 |                    3150.48 |
| HD                | dataset+survey | budget_050mmag_ds100_stress | FR_BETA05_H0fixed675   | 3422 |        0 |         27 |               0.1   |              0.05  | 67.5    | 0.5      |         3082.88 |     9.59394  |          3150.49 |                    3150.49 |
| HD                | dataset+survey | budget_025mmag_ds050        | FR_BETAfree_H0fixed675 | 3422 |        1 |         27 |               0.05  |              0.025 | 67.5    | 0.518263 |         3091.21 |    17.8551   |          3143.12 |                    3151.26 |
| HD                | dataset+survey | budget_025mmag_ds050        | FR_BETAfree_H0free     | 3422 |        2 |         27 |               0.05  |              0.025 | 69.702  | 0.523532 |         3091.17 |    12.5071   |          3137.74 |                    3154.01 |
| HD                | dataset+survey | budget_016mmag_ds032        | FR_BETAfree_H0free     | 3422 |        2 |         27 |               0.032 |              0.016 | 69.7669 | 0.530845 |         3100.44 |    15.7566   |          3139.21 |                    3155.49 |
| HD                | dataset+survey | budget_016mmag_ds032        | FR_BETA05_H0fixed675   | 3422 |        0 |         27 |               0.032 |              0.016 | 67.5    | 0.5      |         3104.33 |    28.6384   |          3155.98 |                    3155.98 |
| HD                | dataset+survey | budget_016mmag_ds032        | FR_BETA05_H0free       | 3422 |        1 |         27 |               0.032 |              0.016 | 69.1825 | 0.5      |         3103.11 |    21.9175   |          3148.04 |                    3156.18 |
| HD                | dataset+survey | budget_050mmag_ds100_stress | FR_BETAfree_H0fixed675 | 3422 |        1 |         27 |               0.1   |              0.05  | 67.5    | 0.514357 |         3081.96 |     8.98427  |          3148.96 |                    3157.09 |
| HD                | dataset+survey | budget_050mmag_ds100_stress | FR_BETA05_H0free       | 3422 |        1 |         27 |               0.1   |              0.05  | 69.3859 | 0.5      |         3082.82 |     8.60836  |          3149.44 |                    3157.58 |
| HD                | dataset        | budget_050mmag_ds100_stress | FR_BETAfree_H0fixed675 | 3422 |        1 |          3 |               0.1   |              0     | 67.5    | 0.541273 |         3128.86 |     2.26102  |          3150.2  |                    3158.34 |
| HD                | dataset+survey | budget_016mmag_ds032        | FR_BETAfree_H0fixed675 | 3422 |        1 |         27 |               0.032 |              0.016 | 67.5    | 0.519264 |         3100.92 |    28.2004   |          3152.14 |                    3160.28 |
| HD                | dataset        | budget_025mmag_ds050        | FR_BETAfree_H0fixed675 | 3422 |        1 |          3 |               0.05  |              0     | 67.5    | 0.53837  |         3129.08 |     8.69616  |          3152.71 |                    3160.85 |
| HD                | dataset        | budget_025mmag_ds050        | FR_BETAfree_H0free     | 3422 |        2 |          3 |               0.05  |              0     | 69.8057 | 0.543189 |         3128.88 |     2.84341  |          3146.65 |                    3162.93 |
| HD                | dataset+survey | budget_050mmag_ds100_stress | FR_BETAfree_H0free     | 3422 |        2 |         27 |               0.1   |              0.05  | 69.657  | 0.515809 |         3081.95 |     7.63514  |          3147.6  |                    3163.88 |
| HD                | dataset        | budget_016mmag_ds032        | FR_BETAfree_H0free     | 3422 |        2 |          3 |               0.032 |              0     | 69.8358 | 0.544438 |         3129.01 |     6.74608  |          3148.05 |                    3164.32 |
| HD                | dataset+survey | budget_010mmag_ds020        | FR_BETAfree_H0free     | 3422 |        2 |         27 |               0.02  |              0.01  | 69.8779 | 0.539498 |         3111.22 |    22.5273   |          3148.58 |                    3164.86 |
| HD                | dataset        | budget_050mmag_ds100_stress | FR_BETAfree_H0free     | 3422 |        2 |          3 |               0.1   |              0     | 69.7897 | 0.542524 |         3128.85 |     0.721643 |          3148.65 |                    3164.92 |
| HD                | dataset        | budget_016mmag_ds032        | FR_BETAfree_H0fixed675 | 3422 |        1 |          3 |               0.032 |              0     | 67.5    | 0.533437 |         3130.04 |    19.8383   |          3162.17 |                    3170.31 |
| HD                | dataset        | budget_050mmag_ds100_stress | FR_BETA05_H0fixed675   | 3422 |        0 |          3 |               0.1   |              0     | 67.5    | 0.5      |         3150.15 |     1.53427  |          3170.76 |                    3170.76 |
| HD                | dataset        | budget_025mmag_ds050        | FR_BETA05_H0fixed675   | 3422 |        0 |          3 |               0.05  |              0     | 67.5    | 0.5      |         3150.18 |     6.08487  |          3171.2  |                    3171.2  |
| HD                | dataset        | budget_010mmag_ds020        | FR_BETAfree_H0free     | 3422 |        2 |          3 |               0.02  |              0     | 69.9106 | 0.54756  |         3129.86 |    16.045    |          3155.45 |                    3171.72 |
| CMB_PANTHEON_ONLY | dataset+survey | budget_025mmag_ds050        | FR_BETA05_H0fixed675   | 3420 |        0 |         27 |               0.05  |              0.025 | 67.5    | 0.5      |         3122.24 |    16.5497   |          3172.83 |                    3172.83 |
| HD                | dataset+survey | budget_010mmag_ds020        | FR_BETA05_H0free       | 3422 |        1 |         27 |               0.02  |              0.01  | 69.0874 | 0.5      |         3117.39 |    32.7185   |          3164.94 |                    3173.08 |

## H0=67.5 Beta-Free Focus

| frame             | random_mode    | prior_config                |   k_random |   dataset_sigma_mag |   survey_sigma_mag |     beta |   beta_minus_0p5 |   data_chi2_map |   prior_chi2 |   marginal_score |   BIC_random_marginal_base |
|:------------------|:---------------|:----------------------------|-----------:|--------------------:|-------------------:|---------:|-----------------:|----------------:|-------------:|-----------------:|---------------------------:|
| CMB_PANTHEON_ONLY | dataset+survey | budget_025mmag_ds050        |         27 |               0.05  |              0.025 | 0.515864 |       0.0158641  |         3120.89 |     15.7033  |          3170.63 |                    3178.77 |
| CMB_PANTHEON_ONLY | dataset        | budget_050mmag_ds100_stress |          3 |               0.1   |              0     | 0.536845 |       0.0368446  |         3154.85 |      1.94957 |          3175.87 |                    3184.01 |
| CMB_PANTHEON_ONLY | dataset+survey | budget_016mmag_ds032        |         27 |               0.032 |              0.016 | 0.51638  |       0.0163802  |         3129.65 |     24.3952  |          3177.04 |                    3185.18 |
| CMB_PANTHEON_ONLY | dataset        | budget_025mmag_ds050        |          3 |               0.05  |              0     | 0.534194 |       0.0341944  |         3155.03 |      7.50449 |          3177.47 |                    3185.6  |
| CMB_PANTHEON_ONLY | dataset+survey | budget_050mmag_ds100_stress |         27 |               0.1   |              0.05  | 0.513231 |       0.0132313  |         3111.85 |      8.38815 |          3178.22 |                    3186.36 |
| CMB_PANTHEON_ONLY | dataset        | budget_016mmag_ds032        |          3 |               0.032 |              0     | 0.529689 |       0.0296888  |         3155.85 |     17.1439  |          3185.28 |                    3193.41 |
| CMB_PANTHEON_ONLY | dataset+survey | budget_010mmag_ds020        |         27 |               0.02  |              0.01  | 0.512293 |       0.0122934  |         3141.58 |     43.2007  |          3199.6  |                    3207.74 |
| CMB_PANTHEON_ONLY | dataset        | budget_010mmag_ds020        |          3 |               0.02  |              0     | 0.520512 |       0.0205124  |         3159.73 |     38.1165  |          3207.39 |                    3215.53 |
| CMB_PANTHEON_ONLY | none           | none                        |          0 |               0     |              0     | 0.468172 |      -0.0318277  |         4009.33 |      0       |          4009.33 |                    4017.47 |
| HD                | dataset+survey | budget_025mmag_ds050        |         27 |               0.05  |              0.025 | 0.518263 |       0.0182626  |         3091.21 |     17.8551  |          3143.12 |                    3151.26 |
| HD                | dataset+survey | budget_050mmag_ds100_stress |         27 |               0.1   |              0.05  | 0.514357 |       0.0143566  |         3081.96 |      8.98427 |          3148.96 |                    3157.09 |
| HD                | dataset        | budget_050mmag_ds100_stress |          3 |               0.1   |              0     | 0.541273 |       0.0412727  |         3128.86 |      2.26102 |          3150.2  |                    3158.34 |
| HD                | dataset+survey | budget_016mmag_ds032        |         27 |               0.032 |              0.016 | 0.519264 |       0.0192641  |         3100.92 |     28.2004  |          3152.14 |                    3160.28 |
| HD                | dataset        | budget_025mmag_ds050        |          3 |               0.05  |              0     | 0.53837  |       0.0383696  |         3129.08 |      8.69616 |          3152.71 |                    3160.85 |
| HD                | dataset        | budget_016mmag_ds032        |          3 |               0.032 |              0     | 0.533437 |       0.0334372  |         3130.04 |     19.8383  |          3162.17 |                    3170.31 |
| HD                | dataset+survey | budget_010mmag_ds020        |         27 |               0.02  |              0.01  | 0.514866 |       0.0148657  |         3114.78 |     49.8767  |          3179.49 |                    3187.63 |
| HD                | dataset        | budget_010mmag_ds020        |          3 |               0.02  |              0     | 0.523403 |       0.0234029  |         3134.62 |     43.982   |          3188.15 |                    3196.29 |
| HD                | none           | none                        |          0 |               0     |              0     | 0.466532 |      -0.0334679  |         4109.87 |      0       |          4109.87 |                    4118    |
| HEL               | dataset+survey | budget_025mmag_ds050        |         27 |               0.05  |              0.025 | 0.513209 |       0.0132088  |         3127.11 |     13.8644  |          3175.01 |                    3183.14 |
| HEL               | dataset        | budget_050mmag_ds100_stress |          3 |               0.1   |              0     | 0.525911 |       0.0259106  |         3155.53 |      1.83924 |          3176.44 |                    3184.58 |
| HEL               | dataset        | budget_025mmag_ds050        |          3 |               0.05  |              0     | 0.523622 |       0.0236222  |         3155.67 |      7.12216 |          3177.72 |                    3185.86 |
| HEL               | dataset+survey | budget_016mmag_ds032        |         27 |               0.032 |              0.016 | 0.512121 |       0.0121206  |         3134.3  |     22.3464  |          3179.64 |                    3187.78 |
| HEL               | dataset+survey | budget_050mmag_ds100_stress |         27 |               0.1   |              0.05  | 0.512536 |       0.0125358  |         3118.94 |      7.55044 |          3184.47 |                    3192.61 |
| HEL               | dataset        | budget_016mmag_ds032        |          3 |               0.032 |              0     | 0.519728 |       0.0197276  |         3156.33 |     16.4386  |          3185.05 |                    3193.19 |
| HEL               | dataset+survey | budget_010mmag_ds020        |         27 |               0.02  |              0.01  | 0.506924 |       0.00692352 |         3143.95 |     41.675   |          3200.45 |                    3208.59 |
| HEL               | dataset        | budget_010mmag_ds020        |          3 |               0.02  |              0     | 0.511779 |       0.0117788  |         3159.54 |     37.3315  |          3206.41 |                    3214.55 |
| HEL               | none           | none                        |          0 |               0     |              0     | 0.465369 |      -0.034631   |         4047.16 |      0       |          4047.16 |                    4055.3  |

## Dataset Blocks

| frame             | dataset           | z_col   |   N_source |   N_used | cov_note                                 | subset_note                | source_path                                                                                                                                       | cov_path                                                                                                                                          |
|:------------------|:------------------|:--------|-----------:|---------:|:-----------------------------------------|:---------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------|
| HD                | PantheonPlusSH0ES | zHD     |       1701 |     1580 | Pantheon-style covariance with leading N | covariance subset inverted | C:\Users\clive\Documents\Cosmology\Pantheon+SH0ES.dat                                                                                             | C:\Users\clive\Documents\Cosmology\Pantheon+SH0ES_STAT+SYS.cov                                                                                    |
| HD                | DES_Dovekie_raw   | zHD     |       1820 |     1820 | DES packed inverse covariance            | precision matrix           | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\DES_SN5YR_Dovekie\DES-Dovekie_HD.csv                               | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\DES_SN5YR_Dovekie\STAT+SYS.npz                                     |
| HD                | Union3p1_UNITY1p8 | z       |         22 |       22 | Union compressed inverse covariance      | precision matrix           | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\Union3_release\mu_mat_union3.1_UNITY1.8_template_cosmo=2_0_mu.fits | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\Union3_release\mu_mat_union3.1_UNITY1.8_template_cosmo=2_0_mu.fits |
| HEL               | PantheonPlusSH0ES | zHEL    |       1701 |     1576 | Pantheon-style covariance with leading N | covariance subset inverted | C:\Users\clive\Documents\Cosmology\Pantheon+SH0ES.dat                                                                                             | C:\Users\clive\Documents\Cosmology\Pantheon+SH0ES_STAT+SYS.cov                                                                                    |
| HEL               | DES_Dovekie_raw   | zHEL    |       1820 |     1820 | DES packed inverse covariance            | precision matrix           | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\DES_SN5YR_Dovekie\DES-Dovekie_HD.csv                               | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\DES_SN5YR_Dovekie\STAT+SYS.npz                                     |
| HEL               | Union3p1_UNITY1p8 | z       |         22 |       22 | Union compressed inverse covariance      | precision matrix           | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\Union3_release\mu_mat_union3.1_UNITY1.8_template_cosmo=2_0_mu.fits | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\Union3_release\mu_mat_union3.1_UNITY1.8_template_cosmo=2_0_mu.fits |
| CMB_PANTHEON_ONLY | PantheonPlusSH0ES | zCMB    |       1701 |     1578 | Pantheon-style covariance with leading N | covariance subset inverted | C:\Users\clive\Documents\Cosmology\Pantheon+SH0ES.dat                                                                                             | C:\Users\clive\Documents\Cosmology\Pantheon+SH0ES_STAT+SYS.cov                                                                                    |
| CMB_PANTHEON_ONLY | DES_Dovekie_raw   | zHD     |       1820 |     1820 | DES packed inverse covariance            | precision matrix           | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\DES_SN5YR_Dovekie\DES-Dovekie_HD.csv                               | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\DES_SN5YR_Dovekie\STAT+SYS.npz                                     |
| CMB_PANTHEON_ONLY | Union3p1_UNITY1p8 | z       |         22 |       22 | Union compressed inverse covariance      | precision matrix           | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\Union3_release\mu_mat_union3.1_UNITY1.8_template_cosmo=2_0_mu.fits | C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\Union3_release\mu_mat_union3.1_UNITY1.8_template_cosmo=2_0_mu.fits |

## Largest Offset Pulls

| frame             | random_mode    | prior_config         | model                  | offset_group              | component_type   |   prior_sigma_mag |   residual_offset_mag |   pull_sigma |
|:------------------|:---------------|:---------------------|:-----------------------|:--------------------------|:-----------------|------------------:|----------------------:|-------------:|
| HD                | dataset        | budget_010mmag_ds020 | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.127502  |     -6.3751  |
| HD                | dataset+survey | budget_010mmag_ds020 | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.124072  |     -6.20359 |
| HEL               | dataset        | budget_010mmag_ds020 | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.119597  |     -5.97984 |
| CMB_PANTHEON_ONLY | dataset        | budget_010mmag_ds020 | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.118827  |     -5.94133 |
| HD                | dataset+survey | budget_010mmag_ds020 | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.118788  |     -5.93941 |
| HD                | dataset        | budget_010mmag_ds020 | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.118682  |     -5.93408 |
| HEL               | dataset+survey | budget_010mmag_ds020 | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.116169  |     -5.80847 |
| HEL               | dataset        | budget_010mmag_ds020 | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.115147  |     -5.75734 |
| CMB_PANTHEON_ONLY | dataset+survey | budget_010mmag_ds020 | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.114682  |     -5.73412 |
| HEL               | dataset+survey | budget_010mmag_ds020 | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.113704  |     -5.68521 |
| CMB_PANTHEON_ONLY | dataset        | budget_010mmag_ds020 | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.111101  |     -5.55507 |
| CMB_PANTHEON_ONLY | dataset+survey | budget_010mmag_ds020 | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.110317  |     -5.51584 |
| HD                | dataset        | budget_016mmag_ds032 | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET | dataset          |             0.032 |            -0.13407   |     -4.1897  |
| HD                | dataset+survey | budget_016mmag_ds032 | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET | dataset          |             0.032 |            -0.130986  |     -4.09333 |
| HEL               | dataset        | budget_016mmag_ds032 | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET | dataset          |             0.032 |            -0.125225  |     -3.91329 |
| CMB_PANTHEON_ONLY | dataset        | budget_016mmag_ds032 | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET | dataset          |             0.032 |            -0.124888  |     -3.90275 |
| HD                | dataset+survey | budget_016mmag_ds032 | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET | dataset          |             0.032 |            -0.124267  |     -3.88335 |
| HD                | dataset        | budget_010mmag_ds020 | FR_BETA05_H0free       | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.0765105 |     -3.82553 |
| HEL               | dataset+survey | budget_016mmag_ds032 | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET | dataset          |             0.032 |            -0.121933  |     -3.81041 |
| HD                | dataset        | budget_016mmag_ds032 | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET | dataset          |             0.032 |            -0.12122   |     -3.78811 |
| CMB_PANTHEON_ONLY | dataset+survey | budget_016mmag_ds032 | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET | dataset          |             0.032 |            -0.120367  |     -3.76148 |
| HEL               | dataset        | budget_010mmag_ds020 | FR_BETA05_H0free       | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.0736176 |     -3.68088 |
| HEL               | dataset+survey | budget_016mmag_ds032 | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET | dataset          |             0.032 |            -0.117699  |     -3.67809 |
| HEL               | dataset        | budget_016mmag_ds032 | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET | dataset          |             0.032 |            -0.117622  |     -3.67568 |
| CMB_PANTHEON_ONLY | dataset+survey | budget_016mmag_ds032 | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET | dataset          |             0.032 |            -0.114659  |     -3.58308 |
| CMB_PANTHEON_ONLY | dataset        | budget_010mmag_ds020 | FR_BETA05_H0free       | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.0714805 |     -3.57403 |
| HD                | dataset+survey | budget_010mmag_ds020 | FR_BETA05_H0free       | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.0712793 |     -3.56396 |
| CMB_PANTHEON_ONLY | dataset        | budget_016mmag_ds032 | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET | dataset          |             0.032 |            -0.113484  |     -3.54637 |
| HEL               | dataset+survey | budget_010mmag_ds020 | FR_BETA05_H0free       | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.0696229 |     -3.48114 |
| CMB_PANTHEON_ONLY | dataset+survey | budget_010mmag_ds020 | FR_BETA05_H0free       | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.065711  |     -3.28555 |
| HEL               | dataset        | budget_010mmag_ds020 | FR_BETAfree_H0free     | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.0644057 |     -3.22028 |
| HEL               | dataset+survey | budget_010mmag_ds020 | FR_BETAfree_H0free     | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.0630783 |     -3.15391 |
| HD                | dataset        | budget_010mmag_ds020 | FR_BETAfree_H0free     | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.062958  |     -3.1479  |
| HD                | dataset+survey | budget_010mmag_ds020 | FR_BETAfree_H0free     | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.0620036 |     -3.10018 |
| CMB_PANTHEON_ONLY | dataset        | budget_010mmag_ds020 | FR_BETAfree_H0free     | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.0592824 |     -2.96412 |
| CMB_PANTHEON_ONLY | dataset+survey | budget_010mmag_ds020 | FR_BETAfree_H0free     | PantheonPlusSH0ES:DATASET | dataset          |             0.02  |            -0.0574951 |     -2.87476 |
| HD                | dataset        | budget_025mmag_ds050 | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET | dataset          |             0.05  |            -0.137064  |     -2.74127 |
| HD                | dataset+survey | budget_025mmag_ds050 | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET | dataset          |             0.05  |            -0.134865  |     -2.69729 |
| HD                | dataset+survey | budget_025mmag_ds050 | FR_BETA05_H0fixed675   | PantheonPlusSH0ES:DATASET | dataset          |             0.05  |            -0.128778  |     -2.57555 |
| HEL               | dataset        | budget_025mmag_ds050 | FR_BETAfree_H0fixed675 | PantheonPlusSH0ES:DATASET | dataset          |             0.05  |            -0.12776   |     -2.5552  |

## Interpretation Notes

- If beta near 0.5 survives only for loose calibration priors, the raw-MU signal remains calibration-dependent.
- If dataset+survey priors improve data chi2 but are disfavored by the marginal random-effect score, the free-offset improvement is mostly Occam-penalized survey structure.
- Large offset pulls identify which calibration blocks need direct audit before any model-promotion claim.