# FR c(z)=c0(1+z) No-Loss Dimming Sequence

Date: 2026-07-16 20:06:34

## Model Under Test

Flat/no-expansion path with no post-emission photon-energy loss:

`r(z) = (c/H0) * integral_0^z (1+z') dz' = (c/H0) * z * (1+z/2)`

The sequence tests:

`DL_alpha(z) = r(z) * (1+z)^alpha`

- `alpha=0`: no extra redshift dimming beyond geometric distance.
- `alpha=0.5`: one flux/rate redshift factor.
- `alpha=1`: two redshift flux factors on the same flat path.

H0 fixed at `67.5` km/s/Mpc. Datasets and redshift frames follow the July 14 raw-MU audit policy.

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

## Primary HD Reference

| prior_config         |   alpha |    N |   data_chi2 |   prior_chi2 |   posterior_objective |   chi2_dof_data |   profiled_rms_mag |   max_abs_offset_pull_sigma |
|:---------------------|--------:|-----:|------------:|-------------:|----------------------:|----------------:|-------------------:|----------------------------:|
| budget_025mmag_ds050 |     0   | 3422 |     3180.98 |      168.723 |               3349.71 |        0.929569 |           0.224997 |                     5.45004 |
| budget_025mmag_ds050 |     0.5 | 3422 |     5421.15 |     1460.56  |               6881.71 |        1.58421  |           0.272233 |                    15.61    |
| budget_025mmag_ds050 |     1   | 3422 |    11994.7  |     4473.47  |              16468.2  |        3.50517  |           0.371298 |                    31.7022  |

## No-Offset Control

| frame             |   alpha |    N |   data_chi2 |   posterior_objective |   chi2_dof_data |   profiled_rms_mag |
|:------------------|--------:|-----:|------------:|----------------------:|----------------:|-------------------:|
| CMB_PANTHEON_ONLY |     0   | 3420 |     4115.91 |               4115.91 |         1.20348 |           0.237608 |
| CMB_PANTHEON_ONLY |     0.5 | 3420 |    34837.1  |              34837.1  |        10.1863  |           0.441196 |
| CMB_PANTHEON_ONLY |     1   | 3420 |   120199    |             120199    |        35.1459  |           0.78738  |
| HD                |     0   | 3422 |     4227.9  |               4227.9  |         1.23551 |           0.238853 |
| HD                |     0.5 | 3422 |    35139.7  |              35139.7  |        10.2688  |           0.442149 |
| HD                |     1   | 3422 |   120718    |             120718    |        35.2769  |           0.78798  |
| HEL               |     0   | 3418 |     4173.71 |               4173.71 |         1.2211  |           0.237897 |
| HEL               |     0.5 | 3418 |    35224.1  |              35224.1  |        10.3055  |           0.444971 |
| HEL               |     1   | 3418 |   120966    |             120966    |        35.3908  |           0.79209  |

## Best Alpha By Frame And Offset Mode

| frame             | offset_mode      | prior_config         |   alpha |   posterior_objective |   data_chi2 |   prior_chi2 |   profiled_rms_mag |
|:------------------|:-----------------|:---------------------|--------:|----------------------:|------------:|-------------:|-------------------:|
| CMB_PANTHEON_ONLY | dataset+idsurvey | budget_025mmag_ds050 |       0 |               3340.77 |     3196.78 |      143.989 |           0.225504 |
| CMB_PANTHEON_ONLY | none             | none                 |       0 |               4115.91 |     4115.91 |        0     |           0.237608 |
| HD                | dataset+idsurvey | budget_025mmag_ds050 |       0 |               3349.71 |     3180.98 |      168.723 |           0.224997 |
| HD                | none             | none                 |       0 |               4227.9  |     4227.9  |        0     |           0.238853 |
| HEL               | dataset+idsurvey | budget_025mmag_ds050 |       0 |               3350.26 |     3201.17 |      149.09  |           0.225622 |
| HEL               | none             | none                 |       0 |               4173.71 |     4173.71 |        0     |           0.237897 |

## Largest HD Offset Pulls

| offset_mode      | prior_config         |   alpha | offset_label                     | offset_type   |   prior_sigma_mag |   offset_mag |   pull_sigma |
|:-----------------|:---------------------|--------:|:---------------------------------|:--------------|------------------:|-------------:|-------------:|
| dataset+idsurvey | budget_025mmag_ds050 |     1   | DES_Dovekie_raw:IDSURVEY_10      | survey        |             0.025 |    -0.792555 |    -31.7022  |
| dataset+idsurvey | budget_025mmag_ds050 |     1   | PantheonPlusSH0ES:IDSURVEY_4     | survey        |             0.025 |    -0.717235 |    -28.6894  |
| dataset+idsurvey | budget_025mmag_ds050 |     1   | PantheonPlusSH0ES:IDSURVEY_10    | survey        |             0.025 |    -0.65262  |    -26.1048  |
| dataset+idsurvey | budget_025mmag_ds050 |     1   | PantheonPlusSH0ES:IDSURVEY_15    | survey        |             0.025 |    -0.562494 |    -22.4997  |
| dataset+idsurvey | budget_025mmag_ds050 |     1   | Union3p1_UNITY1p8                | dataset       |             0.05  |    -0.995512 |    -19.9102  |
| dataset+idsurvey | budget_025mmag_ds050 |     1   | PantheonPlusSH0ES:IDSURVEY_1     | survey        |             0.025 |    -0.434862 |    -17.3945  |
| dataset+idsurvey | budget_025mmag_ds050 |     0.5 | PantheonPlusSH0ES:IDSURVEY_4     | survey        |             0.025 |    -0.39025  |    -15.61    |
| dataset+idsurvey | budget_025mmag_ds050 |     0.5 | DES_Dovekie_raw:IDSURVEY_10      | survey        |             0.025 |    -0.388682 |    -15.5473  |
| dataset+idsurvey | budget_025mmag_ds050 |     0.5 | PantheonPlusSH0ES:IDSURVEY_10    | survey        |             0.025 |    -0.377178 |    -15.0871  |
| dataset+idsurvey | budget_025mmag_ds050 |     1   | PantheonPlusSH0ES:IDSURVEY_OTHER | survey        |             0.025 |    -0.374455 |    -14.9782  |
| dataset+idsurvey | budget_025mmag_ds050 |     1   | PantheonPlusSH0ES:IDSURVEY_150   | survey        |             0.025 |    -0.343942 |    -13.7577  |
| dataset+idsurvey | budget_025mmag_ds050 |     0.5 | PantheonPlusSH0ES:IDSURVEY_15    | survey        |             0.025 |    -0.324916 |    -12.9967  |
| dataset+idsurvey | budget_025mmag_ds050 |     0.5 | PantheonPlusSH0ES:IDSURVEY_1     | survey        |             0.025 |    -0.275466 |    -11.0186  |
| dataset+idsurvey | budget_025mmag_ds050 |     0.5 | Union3p1_UNITY1p8                | dataset       |             0.05  |    -0.505113 |    -10.1023  |
| dataset+idsurvey | budget_025mmag_ds050 |     0.5 | PantheonPlusSH0ES:IDSURVEY_150   | survey        |             0.025 |    -0.240097 |     -9.60387 |
| dataset+idsurvey | budget_025mmag_ds050 |     0.5 | PantheonPlusSH0ES:IDSURVEY_OTHER | survey        |             0.025 |    -0.238002 |     -9.52008 |
| dataset+idsurvey | budget_025mmag_ds050 |     1   | PantheonPlusSH0ES:IDSURVEY_57    | survey        |             0.025 |    -0.182302 |     -7.29208 |
| dataset+idsurvey | budget_025mmag_ds050 |     1   | PantheonPlusSH0ES:IDSURVEY_5     | survey        |             0.025 |    -0.176904 |     -7.07617 |
| dataset+idsurvey | budget_025mmag_ds050 |     1   | PantheonPlusSH0ES:IDSURVEY_64    | survey        |             0.025 |    -0.171655 |     -6.86621 |
| dataset+idsurvey | budget_025mmag_ds050 |     1   | PantheonPlusSH0ES:IDSURVEY_51    | survey        |             0.025 |    -0.158445 |     -6.33782 |

## Readout Boundary

This is a model-consistency diagnostic. It does not promote FR/PLAMB as a cosmological detection by itself. The relevant question is whether the clarified no-loss path prefers alpha=0, 0.5, or 1 under the same calibration and redshift-frame controls used in the existing raw-MU audit.

## Output Files

- `fr_c1pz_noloss_dimming_summary.csv`
- `fr_c1pz_noloss_dimming_offsets.csv`
- `fr_c1pz_noloss_dimming_blocks.csv`
- `fr_c1pz_noloss_dimming_config.json`
