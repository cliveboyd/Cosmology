# FR c(z)=c0(1+z) No-Loss Follow-Up Gates

Date: 2026-07-16 20:49:46

## Scope

This follow-up reuses the identical-offset FR alpha=0 versus LCDM ladder programme and asks whether the apparent FR advantage is broad across catalogues, redshift bands, and calibration-pull structure.

Negative `delta_BIC_FR_minus_LCDMfree` favours FR. Negative `delta_posterior_FR_minus_LCDMfree` favours FR before the base-parameter BIC penalty.

## Executed Cases

| gate          | case_label       | datasets                 | frames                   |   z_min | max_z   | status   |
|:--------------|:-----------------|:-------------------------|:-------------------------|--------:|:--------|:---------|
| dataset_combo | pantheon_only    | pantheon                 | HD,HEL,CMB_PANTHEON_ONLY |    0.01 |         | ran      |
| dataset_combo | des_only         | des_raw                  | HD,HEL                   |    0.01 |         | ran      |
| dataset_combo | union31_only     | union31                  | HD                       |    0.01 |         | ran      |
| dataset_combo | pantheon_des     | pantheon,des_raw         | HD,HEL,CMB_PANTHEON_ONLY |    0.01 |         | ran      |
| dataset_combo | pantheon_union31 | pantheon,union31         | HD,HEL,CMB_PANTHEON_ONLY |    0.01 |         | ran      |
| dataset_combo | des_union31      | des_raw,union31          | HD,HEL                   |    0.01 |         | ran      |
| dataset_combo | all_three        | pantheon,des_raw,union31 | HD,HEL,CMB_PANTHEON_ONLY |    0.01 |         | ran      |
| redshift_band | all_z            | pantheon,des_raw,union31 | HD,HEL,CMB_PANTHEON_ONLY |    0.01 |         | ran      |
| redshift_band | low_z_0p01_0p10  | pantheon,des_raw,union31 | HD,HEL,CMB_PANTHEON_ONLY |    0.01 | 0.1     | ran      |
| redshift_band | mid_z_0p10_0p50  | pantheon,des_raw,union31 | HD,HEL,CMB_PANTHEON_ONLY |    0.1  | 0.5     | ran      |
| redshift_band | high_z_0p50_plus | pantheon,des_raw,union31 | HD,HEL,CMB_PANTHEON_ONLY |    0.5  |         | ran      |

## Leave-One-Dataset-Out: HD Dataset+Survey

| case_label       | prior_config                |    N |   FR_BIC_base |   LCDMfree_Om |   LCDMfree_BIC_base |   delta_BIC_FR_minus_LCDMfree |   delta_posterior_FR_minus_LCDMfree | fr_beats_lcdm_bic   | fr_beats_lcdm_posterior   |
|:-----------------|:----------------------------|-----:|--------------:|--------------:|--------------------:|------------------------------:|------------------------------------:|:--------------------|:--------------------------|
| all_three        | budget_010mmag_ds020        | 3422 |     3765.78   |      0.529011 |           3887.99   |                   -122.207    |                         -114.069    | True                | True                      |
| all_three        | budget_016mmag_ds032        | 3422 |     3545.7    |      0.491427 |           3650.47   |                   -104.768    |                          -96.6301   | True                | True                      |
| all_three        | budget_025mmag_ds050        | 3422 |     3376.9    |      0.432642 |           3448.05   |                    -71.1529   |                          -63.015    | True                | True                      |
| all_three        | budget_050mmag_ds100_stress | 3422 |     3232.23   |      0.358069 |           3233.99   |                     -1.76276  |                            6.37522  | True                | False                     |
| des_only         | budget_010mmag_ds020        | 1820 |     1662.52   |      0.394263 |           1676.82   |                    -14.2993   |                           -6.79269  | True                | True                      |
| des_only         | budget_016mmag_ds032        | 1820 |     1661.62   |      0.367879 |           1669.77   |                     -8.15738  |                           -0.650785 | True                | True                      |
| des_only         | budget_025mmag_ds050        | 1820 |     1661.6    |      0.341613 |           1661.97   |                     -0.372315 |                            7.13428  | True                | False                     |
| des_only         | budget_050mmag_ds100_stress | 1820 |     1663.78   |      0.314781 |           1654.37   |                      9.40546  |                           16.9121   | False               | False                     |
| des_union31      | budget_010mmag_ds020        | 1842 |     1712.57   |      0.384792 |           1713.77   |                     -1.19803  |                            6.32058  | True                | False                     |
| des_union31      | budget_016mmag_ds032        | 1842 |     1712.27   |      0.362114 |           1705.03   |                      7.24226  |                           14.7609   | False               | False                     |
| des_union31      | budget_025mmag_ds050        | 1842 |     1713.01   |      0.342023 |           1696.76   |                     16.2481   |                           23.7667   | False               | False                     |
| des_union31      | budget_050mmag_ds100_stress | 1842 |     1716.5    |      0.323422 |           1690.34   |                     26.1666   |                           33.6852   | False               | False                     |
| pantheon_des     | budget_010mmag_ds020        | 3400 |     3715.73   |      0.573229 |           3776.41   |                    -60.6809   |                          -52.5494   | True                | True                      |
| pantheon_des     | budget_016mmag_ds032        | 3400 |     3495.05   |      0.53309  |           3572.14   |                    -77.0898   |                          -68.9583   | True                | True                      |
| pantheon_des     | budget_025mmag_ds050        | 3400 |     3325.49   |      0.464704 |           3395.79   |                    -70.3003   |                          -62.1688   | True                | True                      |
| pantheon_des     | budget_050mmag_ds100_stress | 3400 |     3179.5    |      0.367126 |           3197.43   |                    -17.925    |                           -9.79349  | True                | True                      |
| pantheon_only    | budget_010mmag_ds020        | 1580 |     2053.21   |      0.6      |           1935.52   |                    117.691    |                          125.056    | False               | False                     |
| pantheon_only    | budget_016mmag_ds032        | 1580 |     1833.43   |      0.6      |           1783.92   |                     49.513    |                           56.8782   | False               | False                     |
| pantheon_only    | budget_025mmag_ds050        | 1580 |     1663.89   |      0.6      |           1665.37   |                     -1.47959  |                            5.88559  | True                | False                     |
| pantheon_only    | budget_050mmag_ds100_stress | 1580 |     1515.73   |      0.448457 |           1533.55   |                    -17.8273   |                          -10.4621   | True                | True                      |
| pantheon_union31 | budget_010mmag_ds020        | 1602 |     2103.26   |      0.6      |           2087.83   |                     15.4235   |                           22.8025   | False               | False                     |
| pantheon_union31 | budget_016mmag_ds032        | 1602 |     1884.08   |      0.575461 |           1911.94   |                    -27.86     |                          -20.481    | True                | True                      |
| pantheon_union31 | budget_025mmag_ds050        | 1602 |     1715.3    |      0.496644 |           1757.68   |                    -42.3808   |                          -35.0018   | True                | True                      |
| pantheon_union31 | budget_050mmag_ds100_stress | 1602 |     1568.45   |      0.39111  |           1578.31   |                     -9.85725  |                           -2.47824  | True                | True                      |
| union31_only     | budget_010mmag_ds020        |   22 |       50.0505 |      0.366527 |             38.6275 |                     11.423    |                           14.514    | False               | False                     |
| union31_only     | budget_016mmag_ds032        |   22 |       50.6505 |      0.351942 |             37.9912 |                     12.6593   |                           15.7504   | False               | False                     |
| union31_only     | budget_025mmag_ds050        |   22 |       51.4102 |      0.342681 |             37.8678 |                     13.5424   |                           16.6335   | False               | False                     |
| union31_only     | budget_050mmag_ds100_stress |   22 |       52.7264 |      0.336368 |             38.5527 |                     14.1738   |                           17.2648   | False               | False                     |

### Dataset-Combination Tally

| case_label       |   bic_wins |   posterior_wins |   n_rows |
|:-----------------|-----------:|-----------------:|---------:|
| all_three        |          4 |                3 |        4 |
| des_only         |          3 |                2 |        4 |
| des_union31      |          1 |                0 |        4 |
| pantheon_des     |          4 |                4 |        4 |
| pantheon_only    |          2 |                1 |        4 |
| pantheon_union31 |          3 |                3 |        4 |
| union31_only     |          0 |                0 |        4 |

## Redshift-Band Check: HD Dataset+Survey

| case_label       | prior_config                |    N |   FR_BIC_base |   LCDMfree_Om |   LCDMfree_BIC_base |   delta_BIC_FR_minus_LCDMfree |   delta_posterior_FR_minus_LCDMfree | fr_beats_lcdm_bic   | fr_beats_lcdm_posterior   |
|:-----------------|:----------------------------|-----:|--------------:|--------------:|--------------------:|------------------------------:|------------------------------------:|:--------------------|:--------------------------|
| all_z            | budget_010mmag_ds020        | 3422 |      3765.78  |      0.529011 |            3887.99  |                    -122.207   |                          -114.069   | True                | True                      |
| all_z            | budget_016mmag_ds032        | 3422 |      3545.7   |      0.491427 |            3650.47  |                    -104.768   |                           -96.6301  | True                | True                      |
| all_z            | budget_025mmag_ds050        | 3422 |      3376.9   |      0.432642 |            3448.05  |                     -71.1529  |                           -63.015   | True                | True                      |
| all_z            | budget_050mmag_ds100_stress | 3422 |      3232.23  |      0.358069 |            3233.99  |                      -1.76276 |                             6.37522 | True                | False                     |
| high_z_0p50_plus | budget_010mmag_ds020        | 1084 |      1018.85  |      0.443065 |             993.177 |                      25.6688  |                            32.6572  | False               | False                     |
| high_z_0p50_plus | budget_016mmag_ds032        | 1084 |      1011.77  |      0.457836 |             983.59  |                      28.1819  |                            35.1703  | False               | False                     |
| high_z_0p50_plus | budget_025mmag_ds050        | 1084 |      1001.83  |      0.465785 |             974.873 |                      26.9549  |                            33.9433  | False               | False                     |
| high_z_0p50_plus | budget_050mmag_ds100_stress | 1084 |       987.288 |      0.442182 |             966.498 |                      20.7906  |                            27.779   | False               | False                     |
| low_z_0p01_0p10  | budget_010mmag_ds020        |  823 |      1086.98  |      0.6      |            1106.25  |                     -19.2645  |                           -12.5516  | True                | True                      |
| low_z_0p01_0p10  | budget_016mmag_ds032        |  823 |      1022.13  |      0.6      |            1038.73  |                     -16.6018  |                            -9.88886 | True                | True                      |
| low_z_0p01_0p10  | budget_025mmag_ds050        |  823 |       947.671 |      0.6      |             961.12  |                     -13.4496  |                            -6.73661 | True                | True                      |
| low_z_0p01_0p10  | budget_050mmag_ds100_stress |  823 |       857.501 |      0.6      |             866.823 |                      -9.32266 |                            -2.6097  | True                | True                      |
| mid_z_0p10_0p50  | budget_010mmag_ds020        | 1515 |      1472.41  |      0.559147 |            1537.65  |                     -65.2398  |                           -57.9166  | True                | True                      |
| mid_z_0p10_0p50  | budget_016mmag_ds032        | 1515 |      1433.78  |      0.529381 |            1495.71  |                     -61.9287  |                           -54.6056  | True                | True                      |
| mid_z_0p10_0p50  | budget_025mmag_ds050        | 1515 |      1403.69  |      0.473264 |            1454.22  |                     -50.5301  |                           -43.2069  | True                | True                      |
| mid_z_0p10_0p50  | budget_050mmag_ds100_stress | 1515 |      1382.4   |      0.370809 |            1399.5   |                     -17.1051  |                            -9.78191 | True                | True                      |

### Redshift-Band Tally

| case_label       |   bic_wins |   posterior_wins |   n_rows |
|:-----------------|-----------:|-----------------:|---------:|
| all_z            |          4 |                3 |        4 |
| high_z_0p50_plus |          0 |                0 |        4 |
| low_z_0p01_0p10  |          4 |                4 |        4 |
| mid_z_0p10_0p50  |          4 |                4 |        4 |

## Calibration-Pull Audit

| gate          | case_label       | prior_config                | offset_dataset    | offset_type   |   n_offsets |   max_abs_pull_sigma |   n_abs_pull_ge_3 |   n_abs_pull_ge_5 |   median_abs_pull_sigma |
|:--------------|:-----------------|:----------------------------|:------------------|:--------------|------------:|---------------------:|------------------:|------------------:|------------------------:|
| dataset_combo | pantheon_des     | budget_010mmag_ds020        | PantheonPlusSH0ES | survey        |          14 |             7.7551   |                 8 |                 4 |                3.16711  |
| dataset_combo | pantheon_union31 | budget_010mmag_ds020        | PantheonPlusSH0ES | survey        |          14 |             7.7551   |                 8 |                 4 |                3.16711  |
| dataset_combo | all_three        | budget_010mmag_ds020        | PantheonPlusSH0ES | survey        |          14 |             7.7551   |                 8 |                 4 |                3.16711  |
| dataset_combo | pantheon_only    | budget_010mmag_ds020        | PantheonPlusSH0ES | survey        |          14 |             7.7551   |                 8 |                 4 |                3.16711  |
| redshift_band | all_z            | budget_010mmag_ds020        | PantheonPlusSH0ES | survey        |          14 |             7.7551   |                 8 |                 4 |                3.16711  |
| dataset_combo | pantheon_union31 | budget_016mmag_ds032        | PantheonPlusSH0ES | survey        |          14 |             6.69001  |                 8 |                 3 |                3.57657  |
| dataset_combo | pantheon_des     | budget_016mmag_ds032        | PantheonPlusSH0ES | survey        |          14 |             6.69001  |                 8 |                 3 |                3.57657  |
| redshift_band | all_z            | budget_016mmag_ds032        | PantheonPlusSH0ES | survey        |          14 |             6.69001  |                 8 |                 3 |                3.57657  |
| dataset_combo | pantheon_only    | budget_016mmag_ds032        | PantheonPlusSH0ES | survey        |          14 |             6.69001  |                 8 |                 3 |                3.57657  |
| dataset_combo | all_three        | budget_016mmag_ds032        | PantheonPlusSH0ES | survey        |          14 |             6.69001  |                 8 |                 3 |                3.57657  |
| redshift_band | all_z            | budget_025mmag_ds050        | PantheonPlusSH0ES | survey        |          14 |             5.45004  |                 8 |                 1 |                3.44522  |
| dataset_combo | all_three        | budget_025mmag_ds050        | PantheonPlusSH0ES | survey        |          14 |             5.45004  |                 8 |                 1 |                3.44522  |
| dataset_combo | pantheon_only    | budget_025mmag_ds050        | PantheonPlusSH0ES | survey        |          14 |             5.45004  |                 8 |                 1 |                3.44522  |
| dataset_combo | pantheon_union31 | budget_025mmag_ds050        | PantheonPlusSH0ES | survey        |          14 |             5.45004  |                 8 |                 1 |                3.44522  |
| dataset_combo | pantheon_des     | budget_025mmag_ds050        | PantheonPlusSH0ES | survey        |          14 |             5.45004  |                 8 |                 1 |                3.44522  |
| redshift_band | low_z_0p01_0p10  | budget_016mmag_ds032        | PantheonPlusSH0ES | survey        |          10 |             5.3529   |                 4 |                 1 |                2.43248  |
| redshift_band | low_z_0p01_0p10  | budget_010mmag_ds020        | PantheonPlusSH0ES | survey        |          10 |             4.94475  |                 2 |                 0 |                1.78645  |
| redshift_band | low_z_0p01_0p10  | budget_025mmag_ds050        | PantheonPlusSH0ES | survey        |          10 |             4.74514  |                 4 |                 0 |                2.79018  |
| redshift_band | mid_z_0p10_0p50  | budget_010mmag_ds020        | PantheonPlusSH0ES | survey        |           4 |             4.22382  |                 3 |                 0 |                3.30992  |
| redshift_band | mid_z_0p10_0p50  | budget_016mmag_ds032        | PantheonPlusSH0ES | survey        |           4 |             4.19389  |                 3 |                 0 |                3.40627  |
| redshift_band | mid_z_0p10_0p50  | budget_025mmag_ds050        | PantheonPlusSH0ES | survey        |           4 |             3.45499  |                 2 |                 0 |                2.895    |
| dataset_combo | all_three        | budget_050mmag_ds100_stress | PantheonPlusSH0ES | survey        |          14 |             3.13713  |                 1 |                 0 |                2.25906  |
| redshift_band | all_z            | budget_050mmag_ds100_stress | PantheonPlusSH0ES | survey        |          14 |             3.13713  |                 1 |                 0 |                2.25906  |
| dataset_combo | pantheon_union31 | budget_050mmag_ds100_stress | PantheonPlusSH0ES | survey        |          14 |             3.13713  |                 1 |                 0 |                2.25906  |
| dataset_combo | pantheon_des     | budget_050mmag_ds100_stress | PantheonPlusSH0ES | survey        |          14 |             3.13713  |                 1 |                 0 |                2.25906  |
| dataset_combo | pantheon_only    | budget_050mmag_ds100_stress | PantheonPlusSH0ES | survey        |          14 |             3.13713  |                 1 |                 0 |                2.25906  |
| redshift_band | low_z_0p01_0p10  | budget_050mmag_ds100_stress | PantheonPlusSH0ES | survey        |          10 |             2.94183  |                 0 |                 0 |                2.35825  |
| redshift_band | high_z_0p50_plus | budget_025mmag_ds050        | PantheonPlusSH0ES | survey        |           3 |             2.61327  |                 0 |                 0 |                2.14737  |
| redshift_band | high_z_0p50_plus | budget_016mmag_ds032        | PantheonPlusSH0ES | survey        |           3 |             2.29201  |                 0 |                 0 |                1.88715  |
| redshift_band | high_z_0p50_plus | budget_050mmag_ds100_stress | PantheonPlusSH0ES | survey        |           3 |             2.21337  |                 0 |                 0 |                1.79532  |
| redshift_band | mid_z_0p10_0p50  | budget_050mmag_ds100_stress | PantheonPlusSH0ES | survey        |           4 |             2.03501  |                 0 |                 0 |                1.74691  |
| redshift_band | high_z_0p50_plus | budget_010mmag_ds020        | PantheonPlusSH0ES | survey        |           3 |             1.70085  |                 0 |                 0 |                1.39871  |
| dataset_combo | des_only         | budget_010mmag_ds020        | DES_Dovekie_raw   | survey        |           5 |             1.2865   |                 0 |                 0 |                0.601941 |
| dataset_combo | des_union31      | budget_010mmag_ds020        | DES_Dovekie_raw   | survey        |           5 |             1.2865   |                 0 |                 0 |                0.601941 |
| dataset_combo | pantheon_des     | budget_010mmag_ds020        | DES_Dovekie_raw   | survey        |           5 |             1.2865   |                 0 |                 0 |                0.601941 |
| redshift_band | all_z            | budget_010mmag_ds020        | DES_Dovekie_raw   | survey        |           5 |             1.2865   |                 0 |                 0 |                0.601941 |
| dataset_combo | all_three        | budget_010mmag_ds020        | DES_Dovekie_raw   | survey        |           5 |             1.2865   |                 0 |                 0 |                0.601941 |
| redshift_band | high_z_0p50_plus | budget_010mmag_ds020        | DES_Dovekie_raw   | survey        |           1 |             1.03077  |                 0 |                 0 |                1.03077  |
| dataset_combo | des_only         | budget_025mmag_ds050        | DES_Dovekie_raw   | survey        |           5 |             0.963326 |                 0 |                 0 |                0.765917 |
| dataset_combo | all_three        | budget_025mmag_ds050        | DES_Dovekie_raw   | survey        |           5 |             0.963326 |                 0 |                 0 |                0.765917 |

## Readout Boundary

This is still a promotion-gate diagnostic, not a discovery claim. A robust FR result should not be carried by a single catalogue, a single redshift band, or implausibly large calibration pulls. The output tables are intended to identify exactly which gate is currently strongest or weakest.

## Output Files

- `fr_c1pz_followup_dataset_comparison.csv`
- `fr_c1pz_followup_redshift_comparison.csv`
- `fr_c1pz_followup_all_summary.csv`
- `fr_c1pz_followup_all_offsets.csv`
- `fr_c1pz_followup_pull_summary.csv`
- `fr_c1pz_followup_run_manifest.csv`
