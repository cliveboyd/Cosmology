# Release-Grounded Raw-MU Calibration Prior Registry

Generated: 2026-07-18 20:11:39

## Interpretation

Projected widths are calibration-mode budgets. They are not automatically additive priors. Every row in this registry is already represented in the corresponding release `STAT+SYS` product, so `sigma_add_mag = 0` and the primary likelihood does not add it again.

Pantheon+ uses its calibration-only covariance. DES uses `C_total-C_stat` only as an all-systematics upper bound. Union3.1 records its documented 0.01 mag important-dataset zero-point scale, already propagated in the compressed release covariance.

## Dataset Modes

| dataset           |    N | component_kind                                        |   sigma_projected_mag |   sigma_documented_mag |   sigma_hierarchy_mag |   sigma_add_mag | prior_action         |
|:------------------|-----:|:------------------------------------------------------|----------------------:|-----------------------:|----------------------:|----------------:|:---------------------|
| DES_Dovekie_raw   | 1820 | all_systematics_covariance_upper_bound                |           0.0104247   |                 nan    |                     0 |               0 | covariance_contained |
| PantheonPlusSH0ES | 1580 | CALIB                                                 |           0.000331554 |                 nan    |                     0 |               0 | covariance_contained |
| PantheonPlusSH0ES | 1580 | CSP_RECAL                                             |           0           |                 nan    |                     0 |               0 | covariance_contained |
| PantheonPlusSH0ES | 1580 | HSTCALSPEC                                            |           0           |                 nan    |                     0 |               0 | covariance_contained |
| PantheonPlusSH0ES | 1580 | CALIB+CSP_RECAL+HSTCALSPEC                            |           0           |                 nan    |                     0 |               0 | covariance_contained |
| Union3p1_UNITY1p8 |   22 | documented_zero_point_budget_in_compressed_covariance |           0.01        |                   0.01 |                   nan |               0 | covariance_contained |

## Survey Modes

| dataset           | component_kind                         | group          |    N |   sigma_projected_mag |   sigma_hierarchy_mag | explicit_base   |
|:------------------|:---------------------------------------|:---------------|-----:|----------------------:|----------------------:|:----------------|
| DES_Dovekie_raw   | all_systematics_covariance_upper_bound | IDSURVEY_10    | 1623 |           0.0117667   |           0.0117667   | STATONLY        |
| DES_Dovekie_raw   | all_systematics_covariance_upper_bound | IDSURVEY_150   |  117 |           0.0196764   |           0.0196764   | STATONLY        |
| DES_Dovekie_raw   | all_systematics_covariance_upper_bound | IDSURVEY_64    |   31 |           0.0124574   |           0.0124574   | STATONLY        |
| DES_Dovekie_raw   | all_systematics_covariance_upper_bound | IDSURVEY_65    |   22 |           0.0230639   |           0.0230639   | STATONLY        |
| DES_Dovekie_raw   | all_systematics_covariance_upper_bound | IDSURVEY_OTHER |   27 |           0.0231729   |           0.0231729   | STATONLY        |
| PantheonPlusSH0ES | CALIB                                  | IDSURVEY_1     |  321 |           0.00807598  |           0.00807598  | STATONLY        |
| PantheonPlusSH0ES | CSP_RECAL                              | IDSURVEY_1     |  321 |           0.00105765  |           0.00105765  | STATONLY        |
| PantheonPlusSH0ES | HSTCALSPEC                             | IDSURVEY_1     |  321 |           0.000873832 |           0.000873832 | STATONLY        |
| PantheonPlusSH0ES | CALIB+CSP_RECAL+HSTCALSPEC             | IDSURVEY_1     |  321 |           0.00819168  |           0.00819168  | STATONLY        |
| PantheonPlusSH0ES | CALIB                                  | IDSURVEY_10    |  203 |           0.00693693  |           0.00693693  | STATONLY        |
| PantheonPlusSH0ES | CSP_RECAL                              | IDSURVEY_10    |  203 |           0.00100712  |           0.00100712  | STATONLY        |
| PantheonPlusSH0ES | HSTCALSPEC                             | IDSURVEY_10    |  203 |           0.00354828  |           0.00354828  | STATONLY        |
| PantheonPlusSH0ES | CALIB+CSP_RECAL+HSTCALSPEC             | IDSURVEY_10    |  203 |           0.00785656  |           0.00785656  | STATONLY        |
| PantheonPlusSH0ES | CALIB                                  | IDSURVEY_15    |  269 |           0.00637122  |           0.00637122  | STATONLY        |
| PantheonPlusSH0ES | CSP_RECAL                              | IDSURVEY_15    |  269 |           0.00100071  |           0.00100071  | STATONLY        |
| PantheonPlusSH0ES | HSTCALSPEC                             | IDSURVEY_15    |  269 |           0.00245353  |           0.00245353  | STATONLY        |
| PantheonPlusSH0ES | CALIB+CSP_RECAL+HSTCALSPEC             | IDSURVEY_15    |  269 |           0.00690027  |           0.00690027  | STATONLY        |
| PantheonPlusSH0ES | CALIB                                  | IDSURVEY_150   |  173 |           0.0150257   |           0.0150257   | STATONLY        |
| PantheonPlusSH0ES | CSP_RECAL                              | IDSURVEY_150   |  173 |           0.00106797  |           0.00106797  | STATONLY        |
| PantheonPlusSH0ES | HSTCALSPEC                             | IDSURVEY_150   |  173 |           0.000298266 |           0.000298266 | STATONLY        |
| PantheonPlusSH0ES | CALIB+CSP_RECAL+HSTCALSPEC             | IDSURVEY_150   |  173 |           0.0150666   |           0.0150666   | STATONLY        |
| PantheonPlusSH0ES | CALIB                                  | IDSURVEY_4     |  160 |           0.0189364   |           0.0189364   | STATONLY        |
| PantheonPlusSH0ES | CSP_RECAL                              | IDSURVEY_4     |  160 |           0.00107597  |           0.00107597  | STATONLY        |
| PantheonPlusSH0ES | HSTCALSPEC                             | IDSURVEY_4     |  160 |           0.0173944   |           0.0173944   | STATONLY        |
| PantheonPlusSH0ES | CALIB+CSP_RECAL+HSTCALSPEC             | IDSURVEY_4     |  160 |           0.0257354   |           0.0257354   | STATONLY        |
| PantheonPlusSH0ES | CALIB                                  | IDSURVEY_5     |   74 |           0.0193121   |           0.0193121   | STATONLY        |
| PantheonPlusSH0ES | CSP_RECAL                              | IDSURVEY_5     |   74 |           0.0135243   |           0.0135243   | STATONLY        |
| PantheonPlusSH0ES | HSTCALSPEC                             | IDSURVEY_5     |   74 |           0.00117973  |           0.00117973  | STATONLY        |
| PantheonPlusSH0ES | CALIB+CSP_RECAL+HSTCALSPEC             | IDSURVEY_5     |   74 |           0.0236063   |           0.0236063   | STATONLY        |
| PantheonPlusSH0ES | CALIB                                  | IDSURVEY_50    |   27 |           0.027088    |           0.027088    | STATONLY        |
| PantheonPlusSH0ES | CSP_RECAL                              | IDSURVEY_50    |   27 |           0.00127582  |           0.00127582  | STATONLY        |
| PantheonPlusSH0ES | HSTCALSPEC                             | IDSURVEY_50    |   27 |           0.000355556 |           0.000355556 | STATONLY        |
| PantheonPlusSH0ES | CALIB+CSP_RECAL+HSTCALSPEC             | IDSURVEY_50    |   27 |           0.0271204   |           0.0271204   | STATONLY        |
| PantheonPlusSH0ES | CALIB                                  | IDSURVEY_51    |   38 |           0.00960541  |           0.00960541  | STATONLY        |
| PantheonPlusSH0ES | CSP_RECAL                              | IDSURVEY_51    |   38 |           4.95129e-05 |           4.95129e-05 | STATONLY        |
| PantheonPlusSH0ES | HSTCALSPEC                             | IDSURVEY_51    |   38 |           0.000607895 |           0.000607895 | STATONLY        |
| PantheonPlusSH0ES | CALIB+CSP_RECAL+HSTCALSPEC             | IDSURVEY_51    |   38 |           0.00962475  |           0.00962475  | STATONLY        |
| PantheonPlusSH0ES | CALIB                                  | IDSURVEY_56    |   36 |           0.0445636   |           0.0445636   | STATONLY        |
| PantheonPlusSH0ES | CSP_RECAL                              | IDSURVEY_56    |   36 |           0.00136671  |           0.00136671  | STATONLY        |
| PantheonPlusSH0ES | HSTCALSPEC                             | IDSURVEY_56    |   36 |           0.001875    |           0.001875    | STATONLY        |
| PantheonPlusSH0ES | CALIB+CSP_RECAL+HSTCALSPEC             | IDSURVEY_56    |   36 |           0.044624    |           0.044624    | STATONLY        |
| PantheonPlusSH0ES | CALIB                                  | IDSURVEY_57    |   85 |           0.00383506  |           0.00383506  | STATONLY        |
| PantheonPlusSH0ES | CSP_RECAL                              | IDSURVEY_57    |   85 |           0.00097925  |           0.00097925  | STATONLY        |
| PantheonPlusSH0ES | HSTCALSPEC                             | IDSURVEY_57    |   85 |           0.000998824 |           0.000998824 | STATONLY        |
| PantheonPlusSH0ES | CALIB+CSP_RECAL+HSTCALSPEC             | IDSURVEY_57    |   85 |           0.00408218  |           0.00408218  | STATONLY        |
| PantheonPlusSH0ES | CALIB                                  | IDSURVEY_63    |   28 |           0.00560802  |           0.00560802  | STATONLY        |
| PantheonPlusSH0ES | CSP_RECAL                              | IDSURVEY_63    |   28 |           0.00127154  |           0.00127154  | STATONLY        |
| PantheonPlusSH0ES | HSTCALSPEC                             | IDSURVEY_63    |   28 |           0.000192857 |           0.000192857 | STATONLY        |
| PantheonPlusSH0ES | CALIB+CSP_RECAL+HSTCALSPEC             | IDSURVEY_63    |   28 |           0.0057536   |           0.0057536   | STATONLY        |
| PantheonPlusSH0ES | CALIB                                  | IDSURVEY_64    |   53 |           0.00321523  |           0.00321523  | STATONLY        |
| PantheonPlusSH0ES | CSP_RECAL                              | IDSURVEY_64    |   53 |           3.16846e-05 |           3.16846e-05 | STATONLY        |
| PantheonPlusSH0ES | HSTCALSPEC                             | IDSURVEY_64    |   53 |           0.00158491  |           0.00158491  | STATONLY        |
| PantheonPlusSH0ES | CALIB+CSP_RECAL+HSTCALSPEC             | IDSURVEY_64    |   53 |           0.00358478  |           0.00358478  | STATONLY        |
| PantheonPlusSH0ES | CALIB                                  | IDSURVEY_65    |   35 |           0.00597239  |           0.00597239  | STATONLY        |
| PantheonPlusSH0ES | CSP_RECAL                              | IDSURVEY_65    |   35 |           0.000919965 |           0.000919965 | STATONLY        |
| PantheonPlusSH0ES | HSTCALSPEC                             | IDSURVEY_65    |   35 |           0.00481714  |           0.00481714  | STATONLY        |
| PantheonPlusSH0ES | CALIB+CSP_RECAL+HSTCALSPEC             | IDSURVEY_65    |   35 |           0.00772791  |           0.00772791  | STATONLY        |
| PantheonPlusSH0ES | CALIB                                  | IDSURVEY_OTHER |   78 |           0.0157944   |           0.0157944   | STATONLY        |
| PantheonPlusSH0ES | CSP_RECAL                              | IDSURVEY_OTHER |   78 |           0.00107893  |           0.00107893  | STATONLY        |
| PantheonPlusSH0ES | HSTCALSPEC                             | IDSURVEY_OTHER |   78 |           0.00521538  |           0.00521538  | STATONLY        |
| PantheonPlusSH0ES | CALIB+CSP_RECAL+HSTCALSPEC             | IDSURVEY_OTHER |   78 |           0.0166682   |           0.0166682   | STATONLY        |

## Explicit-Hierarchy Reconstruction

| dataset           | component                  | base     | status                                        |   relative_frobenius_error |   captured_frobenius_power_fraction |   grouping_covariance_relative_frobenius_error |   released_total_relative_frobenius_error |   diagonal_rmse_mag2 |
|:------------------|:---------------------------|:---------|:----------------------------------------------|---------------------------:|------------------------------------:|-----------------------------------------------:|------------------------------------------:|---------------------:|
| PantheonPlusSH0ES | CALIB                      | STATONLY | stat_only_grouping_reconstruction_sensitivity |                   0.940071 |                          0.116267   |                                    0.227823    |                               0.864597    |          0.000617995 |
| PantheonPlusSH0ES | CSP_RECAL                  | STATONLY | stat_only_grouping_reconstruction_sensitivity |                   0.998049 |                          0.00389751 |                                    0.167467    |                               0.867304    |          3.55514e-05 |
| PantheonPlusSH0ES | HSTCALSPEC                 | STATONLY | stat_only_grouping_reconstruction_sensitivity |                   0.978065 |                          0.0433896  |                                    0.175634    |                               0.866721    |          0.000266136 |
| PantheonPlusSH0ES | CALIB+CSP_RECAL+HSTCALSPEC | STATONLY | stat_only_grouping_reconstruction_sensitivity |                   0.978713 |                          0.0421204  |                                    0.492078    |                               0.864192    |          0.00080548  |
| DES_Dovekie_raw   | ALL_SYSTEMATICS            | STATONLY | all_systematics_common_mode_sensitivity       |                   0.999158 |                          0.0016824  |                                    1.08848e-05 |                               1.08848e-05 |          0.0129146   |

Pantheon+ grouping files include `STATONLY`; the program first forms each grouping increment as `C_grouping - C_stat`, then reconstructs `C_stat + C_grouped_offsets`. DES starts from `STATONLY` and approximates all released systematics with common modes. The released total-covariance likelihood remains primary; reconstruction rows are sensitivities.

## Covariance Audit

```json
{
  "DES_Dovekie_raw": {
    "N_selected": 1820,
    "N_source": 1820,
    "negative_systematics_diagonal_count": 0,
    "stat_loader_note": "DES packed inverse covariance",
    "systematics_diagonal_max_mag2": 0.15905692802018573,
    "systematics_diagonal_min_mag2": 0.00011370281466080781,
    "total_loader_note": "DES packed inverse covariance"
  },
  "PantheonPlusSH0ES": {
    "CSP_RECAL_to_all_systematics_frobenius_ratio": 0.10496443467993903,
    "FRAGILISTIC": {
      "N_labels": 102,
      "eigenvalue_min_mag2": 1.1067452027441792e-06,
      "likelihood_use": "source calibration-parameter covariance only; not mapped directly without the release SN-distance Jacobian",
      "shape": [
        102,
        102
      ],
      "sigma_max_mag": 0.02658604818973166,
      "sigma_median_mag": 0.010442028358164674,
      "sigma_min_mag": 0.0032376764755204083
    },
    "HSTCALSPEC_to_all_systematics_frobenius_ratio": 0.10950166292297617,
    "N_selected": 1580,
    "N_source": 1701,
    "calibration_increment_diagonal_max_mag2_selected": 0.01106273,
    "calibration_increment_diagonal_min_mag2_selected": 3.57000000000135e-06,
    "calibration_symmetry_max_abs": 0.0,
    "calibration_to_all_systematics_frobenius_ratio": 0.1502654494867115,
    "grouping_files_include_STATONLY": true,
    "stat_symmetry_max_abs": 0.0,
    "total_symmetry_max_abs": 0.0
  },
  "Union3p1_UNITY1p8": {
    "N_compressed_nodes": 22,
    "matrix_shape": [
      23,
      23
    ],
    "total_compressed_common_mode_sigma_mag": 0.004196273933288975,
    "total_compressed_common_mode_variance_mag2": 1.7608714923200526e-05
  }
}
```

## Source Metadata

| name                          | role                                             |    bytes | sha256                                                           | source_url                                                                                                                                                             |
|:------------------------------|:-------------------------------------------------|---------:|:-----------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| DES_data                      | distance table                                   |   148002 | 2f57019d783eaa976df80a41b0054171a2d994ee9808d715ce850c2df5720aaf | https://github.com/des-science/DES-SN5YR                                                                                                                               |
| DES_stat                      | STATONLY precision                               |    19366 | ba596d4f035530e53bf5bf2c24359e2dc00ef6449661a30c83bdbe50ae26f6b2 | https://github.com/des-science/DES-SN5YR                                                                                                                               |
| DES_total                     | STAT+SYS precision                               |  6244951 | ffd3124b32148b1372bd95fda9299269f0352a9f8eee02d416c610e38495463b | https://github.com/des-science/DES-SN5YR                                                                                                                               |
| PantheonPlusSH0ES_CSP_RECAL   | CSP recalibration covariance                     | 32157100 | 820169cb4b7424b2ca0ba1a4144848e63fff03831ee7103c4668e5e16b4ca6bd | https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/sytematic_groupings/Pantheon%2BSH0ES_122221_CSP_RECAL.cov  |
| PantheonPlusSH0ES_FRAGILISTIC | 102-parameter Fragilistic calibration covariance |    90672 | c58704e80061f1aa1b74d97861c5b115ee412bea685e4fef88545525d16dfc70 | https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/2_CALIBRATION/FRAGILISTIC_COVARIANCE.npz                                         |
| PantheonPlusSH0ES_HSTCALSPEC  | HST CALSPEC covariance                           | 33147626 | 9f0ad4e15207f216adabd7900c1fa31aed6749e9d989984026ddab55e4c72b94 | https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/sytematic_groupings/Pantheon%2BSH0ES_122221_HSTCALSPEC.cov |
| PantheonPlusSH0ES_calibration | calibration-only covariance                      | 33316040 | 76672f971a32a2966dde9629796b38cd665f784341a4e3fbc28c44fad6f330f3 | https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/sytematic_groupings/Pantheon%2BSH0ES_122221_CALIB.cov      |
| PantheonPlusSH0ES_data        | distance table                                   |   579283 | 1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8 | https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR                                                                      |
| PantheonPlusSH0ES_stat        | STATONLY covariance                              | 31827416 | 9f177129a332735d3637affd20054080d5260815f3ca0809120c05b2c902297f | https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES_STATONLY.cov                              |
| PantheonPlusSH0ES_total       | STAT+SYS covariance                              | 33284960 | abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc | https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR                                                                      |
| Union3p1                      | compressed distance/precision                    |     8640 | bfebb15f64f2c39eb33fca91a6b370cfaebbe8946fcdd58c8d1019ea302482b9 | https://doi.org/10.3847/1538-4357/adc0a5                                                                                                                               |

The Fragilistic file is a 102 x 102 calibration-parameter covariance. It is archived and hashed here, but is not mapped directly into SN-distance offsets without the corresponding release Jacobian; the released distance-covariance grouping remains the applicable product.

- [PantheonPlusSH0ES_release](https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR)
- [PantheonPlusSH0ES_systematic_groupings](https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/sytematic_groupings)
- [PantheonPlus_calibration_paper](https://arxiv.org/abs/2110.03486)
- [PantheonPlusSH0ES_STATONLY_raw](https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES_STATONLY.cov)
- [PantheonPlusSH0ES_CALIB_raw](https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/sytematic_groupings/Pantheon%2BSH0ES_122221_CALIB.cov)
- [PantheonPlusSH0ES_CSP_RECAL_raw](https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/sytematic_groupings/Pantheon%2BSH0ES_122221_CSP_RECAL.cov)
- [PantheonPlusSH0ES_HSTCALSPEC_raw](https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/sytematic_groupings/Pantheon%2BSH0ES_122221_HSTCALSPEC.cov)
- [PantheonPlusSH0ES_FRAGILISTIC_raw](https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/2_CALIBRATION/FRAGILISTIC_COVARIANCE.npz)
- [DES_SN5YR_release](https://github.com/des-science/DES-SN5YR)
- [Union3_paper](https://doi.org/10.3847/1538-4357/adc0a5)

## Double-Counting Guard

A residual-offset sensitivity may be run only for a demonstrably external component absent from the release covariance. It must use the identical design and prior for FR and Lambda-CDM and must not be presented as independent evidence.
