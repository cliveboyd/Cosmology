# PLAMB boundary and high-redshift audit

**Analysis date:** 2026-07-18  
**Status:** post hoc diagnostic; cannot promote a PLAMB claim  
**Locked data:** 3,422 supernovae with released total covariance

## Headline result

Under the widened cap gamma_c<=8, the full-sample PATH fit gives gamma_c=2.276517, p=0.784634 and chi-squared=3047.635738. Its boundary fraction is 0.284565, and the fixed-grid profile rise at gamma_c=8 is 2197.651774.

The widened PATH curve is almost identical to Lambda-CDM in raw chi-squared, but still has Delta BIC=+8.111201 because it uses one additional shape parameter.

The corresponding z<0.5 fit gives gamma_c=1.524995, p=0.246775. Its conditional z>=0.5 prediction has Delta chi-squared=+321.928234 relative to Lambda-CDM.

A narrow outcome-visible grid interval, gamma_c=2.30-2.40, would beat Lambda-CDM on the held-out quadratic score. Its best grid point is gamma_c=2.30, but that point is worse than the independently selected low-redshift optimum by Delta chi-squared=4.597673. It is therefore an exploratory tension diagnostic, not a valid retuning of the predictor.

**Diagnostic decision: FAIL.** This decision concerns boundedness and transfer only; a pass would still not be confirmatory evidence because the wider bound was outcome-motivated.

## Cap fits

| Model | Sample | gamma cap | gamma_c | p | alpha | chi-squared | high-z Delta chi-squared | Bound fraction |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `PATH_FREE` | full | 2 | 2.000000 | 0.628192 | 0.000000 | 3052.869589 | - | 0.000000 |
| `PATH_FREE` | full | 4 | 2.276516 | 0.784633 | 0.000000 | 3047.635738 | - | 0.430871 |
| `PATH_FREE` | full | 8 | 2.276517 | 0.784634 | 0.000000 | 3047.635738 | - | 0.284565 |
| `PATH_FREE` | low_z | 2 | 1.524998 | 0.246777 | 0.000000 | 2101.391879 | +321.926337 | 0.237501 |
| `PATH_FREE` | low_z | 4 | 1.524998 | 0.246777 | 0.000000 | 2101.391879 | +321.926550 | 0.381249 |
| `PATH_FREE` | low_z | 8 | 1.524995 | 0.246775 | 0.000000 | 2101.391879 | +321.928234 | 0.190624 |
| `GENERAL_FREE` | full | 2 | 1.831979 | 1.753326 | 0.701593 | 3047.174295 | - | 0.084010 |
| `GENERAL_FREE` | full | 4 | 1.831979 | 1.753333 | 0.701597 | 3047.174295 | - | 0.457995 |
| `GENERAL_FREE` | full | 8 | 1.831980 | 1.753334 | 0.701597 | 3047.174295 | - | 0.228997 |
| `GENERAL_FREE` | low_z | 2 | 2.000000 | -0.180063 | -0.452072 | 2101.386163 | +318.911013 | 0.000000 |
| `GENERAL_FREE` | low_z | 4 | 2.050233 | -0.225670 | -0.500000 | 2101.385593 | +318.003376 | 0.487442 |
| `GENERAL_FREE` | low_z | 8 | 2.050233 | -0.225670 | -0.500000 | 2101.385593 | +318.003494 | 0.256279 |

## Additive conditional decomposition

The rows below are exact Cholesky-innovation contributions for the cap-8 PATH low-redshift fit relative to the independently fitted low-redshift Lambda-CDM control.

| Release | Redshift band | PATH chi-squared | LCDM chi-squared | Delta chi-squared |
|---|---|---:|---:|---:|
| DES_Dovekie_raw | `z0p50_0p65` | 415.351816 | 418.040031 | -2.688215 |
| DES_Dovekie_raw | `z0p65_0p80` | 248.729243 | 245.125477 | +3.603766 |
| DES_Dovekie_raw | `z0p80_1p00` | 130.555172 | 115.637015 | +14.918158 |
| DES_Dovekie_raw | `z1p00_plus` | 31.798426 | 16.642999 | +15.155426 |
| PantheonPlusSH0ES | `z0p50_0p65` | 83.739911 | 75.901578 | +7.838334 |
| PantheonPlusSH0ES | `z0p65_0p80` | 44.874408 | 35.605421 | +9.268987 |
| PantheonPlusSH0ES | `z0p80_1p00` | 8.860610 | 6.748007 | +2.112603 |
| PantheonPlusSH0ES | `z1p00_plus` | 31.643318 | 21.777131 | +9.866186 |
| Union3p1_UNITY1p8 | `z0p50_0p65` | 0.559193 | 1.678671 | -1.119478 |
| Union3p1_UNITY1p8 | `z0p65_0p80` | 12.477909 | 2.458456 | +10.019453 |
| Union3p1_UNITY1p8 | `z0p80_1p00` | 3.379072 | 4.242308 | -0.863236 |
| Union3p1_UNITY1p8 | `z1p00_plus` | 259.925520 | 6.109270 | +253.816250 |

## Largest survey contributions

| Release | Survey | Delta chi-squared |
|---|---|---:|
| Union3p1_UNITY1p8 | `COMPRESSED` | +261.852989 |
| DES_Dovekie_raw | `10` | +30.989135 |
| PantheonPlusSH0ES | `106` | +14.279424 |
| PantheonPlusSH0ES | `10` | +8.692538 |
| PantheonPlusSH0ES | `4` | +6.632806 |
| PantheonPlusSH0ES | `101` | +2.880763 |
| PantheonPlusSH0ES | `15` | +1.781977 |
| PantheonPlusSH0ES | `100` | -5.181397 |

## Largest ordered-innovation contributions

| Release | Source index | z | Survey | Redshift band | Delta chi-squared |
|---|---:|---:|---|---|---:|
| Union3p1_UNITY1p8 | 21 | 2.262260 | `COMPRESSED` | `z1p00_plus` | +246.280504 |
| Union3p1_UNITY1p8 | 14 | 0.750000 | `COMPRESSED` | `z0p65_0p80` | +4.408860 |
| DES_Dovekie_raw | 1813 | 1.046800 | `10` | `z1p00_plus` | +4.382451 |
| PantheonPlusSH0ES | 1699 | 1.911650 | `106` | `z1p00_plus` | +4.337333 |
| Union3p1_UNITY1p8 | 18 | 1.094415 | `COMPRESSED` | `z1p00_plus` | +4.210380 |
| DES_Dovekie_raw | 1817 | 1.121320 | `10` | `z1p00_plus` | +4.196545 |
| Union3p1_UNITY1p8 | 13 | 0.700000 | `COMPRESSED` | `z0p65_0p80` | +4.120775 |
| Union3p1_UNITY1p8 | 20 | 1.390961 | `COMPRESSED` | `z1p00_plus` | +3.661509 |
| PantheonPlusSH0ES | 1700 | 2.261370 | `106` | `z1p00_plus` | +3.195348 |
| PantheonPlusSH0ES | 1698 | 1.801190 | `106` | `z1p00_plus` | +3.088341 |
| DES_Dovekie_raw | 1805 | 1.017000 | `10` | `z1p00_plus` | +3.028060 |
| DES_Dovekie_raw | 1756 | 0.875650 | `10` | `z0p80_1p00` | +2.970844 |

Individual values are ordered Gaussian innovations after sorting by redshift and source index. Their allocation is order-dependent, while their sum is the exact, order-invariant conditional quadratic form.

Union3.1 is a 22-node compressed likelihood. Its source indices are compressed redshift nodes, not individual supernovae. The dominant z=2.262260 entry therefore identifies endpoint leverage in the compressed release rather than one aberrant event.

## Covariance sensitivities

| Readout | Model | High-z chi-squared | Delta chi-squared vs LCDM |
|---|---|---:|---:|
| `released_conditional` | `PATH_FREE` | 1271.894599 | +321.928234 |
| `released_marginal` | `PATH_FREE` | 1117.511127 | +183.038487 |
| `released_conditional` | `GENERAL_FREE` | 1267.969858 | +318.003494 |
| `released_marginal` | `GENERAL_FREE` | 1115.604726 | +181.132086 |
| `released_conditional` | `LCDM` | 949.966364 | +0.000000 |
| `released_marginal` | `LCDM` | 934.472640 | +0.000000 |
| `diagonal_refit` | `PATH_FREE` | 984.936838 | +90.627560 |
| `diagonal_refit` | `GENERAL_FREE` | 984.900153 | +90.590875 |
| `diagonal_refit` | `LCDM` | 894.309278 | +0.000000 |

## Exploratory release-exclusion sensitivity

This refit was triggered by the dominant Union3.1 compressed endpoint and was not a registered gate.

| Excluded release | PATH gamma_c | PATH p | N test | Delta high-z chi-squared |
|---|---:|---:|---:|---:|
| `NONE` | 1.524995 | 0.246775 | 1085 | +321.928235 |
| `PantheonPlusSH0ES` | 0.809860 | -0.323832 | 875 | +678.868786 |
| `DES_Dovekie_raw` | 1.636172 | 0.320817 | 223 | +250.965900 |
| `Union3p1_UNITY1p8` | 1.729055 | 0.392636 | 1072 | +38.606432 |

The conditional score is primary. Marginal and diagonal rows diagnose whether released cross-covariance or individual variances control the sign.
The general low-redshift control reaches alpha=-0.5, its registered lower bound, so it does not provide an identified alternative explanation of the transfer failure.

## Registered diagnostic gates

| Gate | Pass |
|---|---|
| `cap8_full_optimum_at_least_5_percent_from_boundary` | **True** |
| `profile_rise_at_gamma8_at_least_3p84` | **True** |
| `cap8_best_two_optimisers_agree_within_1e_4` | **True** |
| `cap8_PATH_high_z_delta_chi2_not_positive` | **False** |
| `no_release_redshift_bin_loss_above_10` | **False** |
| `no_single_innovation_above_25_percent_of_positive_excess` | **False** |
| `no_survey_above_50_percent_of_positive_excess` | **True** |
| `covariance_sensitivity_signs_match` | **True** |

## Interpretation boundary

This audit can establish whether the earlier gamma boundary was artificial and whether the associated curve transfers beyond z=0.5. It cannot convert a post hoc bound extension into evidence for clock-rate or light-speed evolution. A physical claim still requires a pre-data derivation of gamma_c, p and the flux exponent, followed by a fresh external-data test.

## Reproduction

```powershell
python github_export/code/rawmu/audit_plamb_boundary_highz_2026_07_18.py
```
