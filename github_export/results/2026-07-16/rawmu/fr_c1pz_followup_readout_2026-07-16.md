# FR No-Loss Follow-Up Gates Readout

Date: 2026-07-16

## Executive readout

Leave-one-dataset-out gate, HD dataset+survey rows:

| case_label       |   bic_wins |   posterior_wins |   n_rows |
|:-----------------|-----------:|-----------------:|---------:|
| all_three        |          4 |                3 |        4 |
| des_only         |          3 |                2 |        4 |
| des_union31      |          1 |                0 |        4 |
| pantheon_des     |          4 |                4 |        4 |
| pantheon_only    |          2 |                1 |        4 |
| pantheon_union31 |          3 |                3 |        4 |
| union31_only     |          0 |                0 |        4 |

Redshift-band gate, HD dataset+survey rows:

| case_label       |   bic_wins |   posterior_wins |   n_rows |
|:-----------------|-----------:|-----------------:|---------:|
| all_z            |          4 |                3 |        4 |
| high_z_0p50_plus |          0 |                0 |        4 |
| low_z_0p01_0p10  |          4 |                4 |        4 |
| mid_z_0p10_0p50  |          4 |                4 |        4 |

## Largest calibration-pull structures

| gate          | case_label       | prior_config         | offset_dataset    | offset_type   |   max_abs_pull_sigma |   n_abs_pull_ge_3 |   n_abs_pull_ge_5 |
|:--------------|:-----------------|:---------------------|:------------------|:--------------|---------------------:|------------------:|------------------:|
| dataset_combo | pantheon_des     | budget_010mmag_ds020 | PantheonPlusSH0ES | survey        |              7.7551  |                 8 |                 4 |
| dataset_combo | pantheon_union31 | budget_010mmag_ds020 | PantheonPlusSH0ES | survey        |              7.7551  |                 8 |                 4 |
| dataset_combo | all_three        | budget_010mmag_ds020 | PantheonPlusSH0ES | survey        |              7.7551  |                 8 |                 4 |
| dataset_combo | pantheon_only    | budget_010mmag_ds020 | PantheonPlusSH0ES | survey        |              7.7551  |                 8 |                 4 |
| redshift_band | all_z            | budget_010mmag_ds020 | PantheonPlusSH0ES | survey        |              7.7551  |                 8 |                 4 |
| dataset_combo | pantheon_union31 | budget_016mmag_ds032 | PantheonPlusSH0ES | survey        |              6.69001 |                 8 |                 3 |
| dataset_combo | pantheon_des     | budget_016mmag_ds032 | PantheonPlusSH0ES | survey        |              6.69001 |                 8 |                 3 |
| redshift_band | all_z            | budget_016mmag_ds032 | PantheonPlusSH0ES | survey        |              6.69001 |                 8 |                 3 |
| dataset_combo | pantheon_only    | budget_016mmag_ds032 | PantheonPlusSH0ES | survey        |              6.69001 |                 8 |                 3 |
| dataset_combo | all_three        | budget_016mmag_ds032 | PantheonPlusSH0ES | survey        |              6.69001 |                 8 |                 3 |
| redshift_band | all_z            | budget_025mmag_ds050 | PantheonPlusSH0ES | survey        |              5.45004 |                 8 |                 1 |
| dataset_combo | all_three        | budget_025mmag_ds050 | PantheonPlusSH0ES | survey        |              5.45004 |                 8 |                 1 |
| dataset_combo | pantheon_only    | budget_025mmag_ds050 | PantheonPlusSH0ES | survey        |              5.45004 |                 8 |                 1 |
| dataset_combo | pantheon_union31 | budget_025mmag_ds050 | PantheonPlusSH0ES | survey        |              5.45004 |                 8 |                 1 |
| dataset_combo | pantheon_des     | budget_025mmag_ds050 | PantheonPlusSH0ES | survey        |              5.45004 |                 8 |                 1 |

## Claim boundary

The right reading is gate-based: a favourable FR-vs-LCDM delta is only persuasive if it persists across catalogue combinations and redshift bands without relying on implausibly large calibration pulls.
