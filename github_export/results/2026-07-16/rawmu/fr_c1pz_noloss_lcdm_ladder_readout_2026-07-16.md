# FR Alpha=0 Calibration Ladder vs LCDM Readout

Date: 2026-07-16

## Executive readout

FR alpha=0 beats LCDM Om-free by base BIC in `4` of `4` HD calibration-budget rows and by raw posterior objective in `3` of `4` rows.
At the loosest 50/100 mmag stress budget, the distinction matters: FR is better by base BIC (`-1.763`), but LCDM has the slightly better raw posterior objective (`FR-LCDM = 6.375`).

| prior_config                |   FR_BIC_base |   LCDMfree_Om |   LCDMfree_BIC_base |   delta_BIC_FR_minus_LCDMfree |   delta_posterior_FR_minus_LCDMfree |
|:----------------------------|--------------:|--------------:|--------------------:|------------------------------:|------------------------------------:|
| budget_010mmag_ds020        |       3765.78 |      0.529011 |             3887.99 |                    -122.207   |                          -114.069   |
| budget_016mmag_ds032        |       3545.7  |      0.491427 |             3650.47 |                    -104.768   |                           -96.6301  |
| budget_025mmag_ds050        |       3376.9  |      0.432642 |             3448.05 |                     -71.1529  |                           -63.015   |
| budget_050mmag_ds100_stress |       3232.23 |      0.358069 |             3233.99 |                      -1.76276 |                             6.37522 |

## FR calibration ladder

| prior_config                |   data_chi2_map |   prior_chi2 |   BIC_base |   profiled_rms_mag |   max_abs_offset_pull_sigma |
|:----------------------------|----------------:|-------------:|-----------:|-------------------:|----------------------------:|
| budget_010mmag_ds020        |         3512.64 |     243.061  |    3765.78 |           0.229474 |                     7.7551  |
| budget_016mmag_ds032        |         3300.16 |     228.28   |    3545.7  |           0.226778 |                     6.69001 |
| budget_025mmag_ds050        |         3180.98 |     168.723  |    3376.9  |           0.224997 |                     5.45004 |
| budget_050mmag_ds100_stress |         3108.5  |      75.5105 |    3232.23 |           0.223504 |                     3.13713 |

## Claim boundary

A negative FR-minus-LCDM delta favors the clarified FR no-loss distance law under identical offsets. The result should still be treated as calibration-gated, because the tightest budgets expose large offset pulls.
