# Peter clock-law v1: physical and supernova consistency audit

Generated: 2026-07-19 17:43:52

## Decision

**THE EXACT `1+z/2` PATH IS IMPLEMENTED CORRECTLY BUT FAILS THE SUPERNOVA AND PHYSICAL-CLOSURE GATES. DO NOT ESCALATE TO A LIGHT-CURVE-LEVEL FIT.**

The exact no-loss, constant-luminosity branch gives `Delta BIC=+94.344983` relative to flat Lambda-CDM for `N=3422`. Positive values favour Lambda-CDM. This reproduces the already-known fixed-law result and is not new evidence by itself.

More importantly, the proposed atomic-clock statement is inconsistent with the registered spectroscopic-redshift relation. With `a=1` and `C=1+z`, the required atomic-frequency ratio is `Q=1`, not `Q=(1+z)^-1`. With fixed `h`, the stated mass and light-speed laws also make the Compton frequency grow as `1+z`, rather than fall as `(1+z)^-1`.

## Frozen status

Protocol SHA-256: `3f595a34f908771475fd631ba593917dd8f5fefc76dee7d570967fa6ea0fdf6c`.

The earlier exact-path Delta BIC was disclosed in the protocol. The newly frozen content was the physical closure, deterministic source and photon branches, covariance sensitivities and escalation rule.

## Registered equations

$$
\begin{aligned}
a(z)                          &= 1, \\
C(z)                          &= 1+z, \\
R(z)                          &= H_0, \\
\chi(z)                       &= \frac{c_0}{H_0}z\left(1+\frac{z}{2}\right), \\
1+z                           &= \frac{C(z)}{a(z)Q(z)}, \\
Q_{\rm required}(z)           &= 1, \\
D_{\rm eff}(z)                &= \chi(z)(1+z)^s, \\
s                             &= \tau-\frac{\ell}{2}, \\
L_{\rm Ia}(z)/L_{\rm Ia,0}    &= (1+z)^\ell.
\end{aligned}
$$

The supernova likelihood identifies only the composite `s`. It cannot separate photon-distance transfer `tau` from intrinsic luminosity evolution `ell`.

## Primary likelihood

| model                      | parameter_name   |   parameter_value |   parameter_sigma |     chi2 |      BIC |   delta_BIC_vs_LCDM |
|:---------------------------|:-----------------|------------------:|------------------:|---------:|---------:|--------------------:|
| PETER_NOLOSS_LCONST        | s_fixed          |         0         |      nan          |  3150.15 |  3174.56 |             94.345  |
| PETER_NOLOSS_PARTICLECOUNT | s_fixed          |        -0.5       |      nan          |  6957.83 |  6982.24 |           3902.03   |
| PETER_NOLOSS_CHANDRA       | s_fixed          |        -1.75      |      nan          | 45438    | 45462.4  |          42382.2    |
| STATIC_DUALITY_LCONST      | s_fixed          |         2         |      nan          | 54115.8  | 54140.2  |          51060      |
| STATIC_DUALITY_CHANDRA     | s_fixed          |         0.25      |      nan          |  3728.67 |  3753.08 |            672.869  |
| COMPOSITE_S_FREE           | s                |         0.0376049 |        0.00869096 |  3131.42 |  3163.98 |             83.7609 |
| FLAT_LCDM                  | Omega_m          |         0.33068   |      nan          |  3047.66 |  3080.21 |              0      |

The diagnostic composite fit gives `s=0.037605 +/- 0.008691`. It is equivalent to `ell=-0.075210 +/- 0.017382` under the no-loss interpretation, or `ell=3.924790 +/- 0.017382` under the static distance-duality interpretation. These are exactly degenerate descriptions of the same magnitude curve.

## Deterministic branches

| model                      |   tau_photon_distance_power |   ell_luminosity_power |   s_composite_power | description                                                            |
|:---------------------------|----------------------------:|-----------------------:|--------------------:|:-----------------------------------------------------------------------|
| PETER_NOLOSS_LCONST        |                           0 |                    0   |                0    | Peter no-loss, constant standardised luminosity                        |
| PETER_NOLOSS_PARTICLECOUNT |                           0 |                    1   |               -0.5  | Peter no-loss, illustrative particle-count luminosity                  |
| PETER_NOLOSS_CHANDRA       |                           0 |                    3.5 |               -1.75 | Peter no-loss, fixed-constant Chandrasekhar luminosity proxy           |
| STATIC_DUALITY_LCONST      |                           2 |                    0   |                2    | Static distance-duality, constant standardised luminosity              |
| STATIC_DUALITY_CHANDRA     |                           2 |                    3.5 |                0.25 | Static distance-duality, fixed-constant Chandrasekhar luminosity proxy |

With fixed `hbar`, `G` and particle masses proportional to `(1+z)^-1`, the Chandrasekhar proxy predicts `M_Ch/M_Ch,0=(1+z)^(7/2)`. The corresponding no-loss branch gives `Delta BIC=+42382.152295`.

## Physical consistency readout

|   z_observed |   C_c_over_c0 |   Q_required_by_spectroscopic_closure |   Q_asserted_clock_control |   z_predicted_if_asserted_Q_used |   Compton_frequency_ratio_fixed_h |   M_Ch_ratio_fixed_hbar_G |   N_Ch_ratio_fixed_hbar_G |
|-------------:|--------------:|--------------------------------------:|---------------------------:|---------------------------------:|----------------------------------:|--------------------------:|--------------------------:|
|      0       |       1       |                                     1 |                   1        |                          0       |                           1       |                   1       |                   1       |
|      0.1     |       1.1     |                                     1 |                   0.909091 |                          0.21    |                           1.1     |                   1.39596 |                   1.53556 |
|      0.5     |       1.5     |                                     1 |                   0.666667 |                          1.25    |                           1.5     |                   4.13351 |                   6.20027 |
|      1       |       2       |                                     1 |                   0.5      |                          3       |                           2       |                  11.3137  |                  22.6274  |
|      1.5     |       2.5     |                                     1 |                   0.4      |                          5.25    |                           2.5     |                  24.7053  |                  61.7632  |
|      2.26226 |       3.26226 |                                     1 |                   0.306536 |                          9.64234 |                           3.26226 |                  62.7069  |                 204.566   |

The identity `integral(1+z) dz = z(1+z/2)` passes. The failure is not that integral; it is the additional physical identification of redshift, atomic frequency, source luminosity and observed flux.

## Covariance hierarchy sensitivities

| covariance_variant                                   |    N |    chi2 |   delta_BIC_vs_LCDM |
|:-----------------------------------------------------|-----:|--------:|--------------------:|
| released_total_primary                               | 3422 | 3150.15 |             94.345  |
| pantheon_calib_statonly_grouped                      | 3422 | 3232.12 |             98.1247 |
| pantheon_calib_csp_recal_hstcalspec_statonly_grouped | 3422 | 3229.99 |             96.6008 |
| des_allsys_statonly_grouped                          | 3422 | 3192.08 |             86.3631 |
| combined_calib_des_statonly_grouped                  | 3422 | 3274.05 |             89.702  |

The exact-branch sensitivity range is [+86.363125, +98.124671]. The released total covariance remains primary. Grouped reconstructions are diagnostics because they omit substantial covariance structure.

## Release and redshift diagnostics

### Release splits

| scope                     | model               |    N |   parameter_value |      chi2 |   delta_BIC_vs_LCDM |
|:--------------------------|:--------------------|-----:|------------------:|----------:|--------------------:|
| release:PantheonPlusSH0ES | PETER_NOLOSS_LCONST | 1580 |         0         | 1436.66   |             42.2006 |
| release:PantheonPlusSH0ES | COMPOSITE_S_FREE    | 1580 |         0.0488233 | 1425.94   |             38.8405 |
| release:PantheonPlusSH0ES | FLAT_LCDM           | 1580 |         0.33163   | 1387.1    |              0      |
| release:DES_Dovekie_raw   | PETER_NOLOSS_LCONST | 1820 |         0         | 1667.21   |             27.9362 |
| release:DES_Dovekie_raw   | COMPOSITE_S_FREE    | 1820 |         0.0399412 | 1657.28   |             25.5086 |
| release:DES_Dovekie_raw   | FLAT_LCDM           | 1820 |         0.328733  | 1631.77   |              0      |
| release:Union3p1_UNITY1p8 | PETER_NOLOSS_LCONST |   22 |         0         |   46.2706 |             14.4203 |
| release:Union3p1_UNITY1p8 | COMPOSITE_S_FREE    |   22 |         0.0117265 |   45.925  |             17.1658 |
| release:Union3p1_UNITY1p8 | FLAT_LCDM           |   22 |         0.333956  |   28.7592 |              0      |

### Descriptive redshift splits

| scope   | model               |    N |   parameter_value |     chi2 |   delta_BIC_vs_LCDM |
|:--------|:--------------------|-----:|------------------:|---------:|--------------------:|
| z<0.5   | PETER_NOLOSS_LCONST | 2337 |          0        | 2163.12  |            48.9909  |
| z<0.5   | COMPOSITE_S_FREE    | 2337 |          0.105125 | 2101.96  |            -4.4082  |
| z<0.5   | FLAT_LCDM           | 2337 |          0.345514 | 2106.37  |             0       |
| z>=0.5  | PETER_NOLOSS_LCONST | 1085 |          0        |  962.031 |            26.5862  |
| z>=0.5  | COMPOSITE_S_FREE    | 1085 |         -0.138443 |  935.442 |             6.98663 |
| z>=0.5  | FLAT_LCDM           | 1085 |          0.319047 |  928.456 |             0       |

These subset cells profile their own release intercepts and are descriptive, not independent predictive tests.

## Calibration-budget modes

| dataset           | group       |   N |   budget_sigma_mag |   delta_mode_mag |   abs_delta_over_budget | within_budget   |
|:------------------|:------------|----:|-------------------:|-----------------:|------------------------:|:----------------|
| PantheonPlusSH0ES | IDSURVEY_64 |  53 |         0.00321523 |       -0.03681   |                11.4486  | False           |
| PantheonPlusSH0ES | IDSURVEY_57 |  85 |         0.00383506 |       -0.0407742 |                10.632   | False           |
| PantheonPlusSH0ES | IDSURVEY_63 |  28 |         0.00560802 |       -0.038889  |                 6.93453 | False           |
| PantheonPlusSH0ES | IDSURVEY_65 |  35 |         0.00597239 |       -0.0367685 |                 6.15642 | False           |
| PantheonPlusSH0ES | IDSURVEY_10 | 203 |         0.00693693 |        0.0365831 |                 5.27368 | False           |
| PantheonPlusSH0ES | IDSURVEY_15 | 269 |         0.00637122 |        0.0298566 |                 4.68617 | False           |
| DES_Dovekie_raw   | IDSURVEY_64 |  31 |         0.0124574  |       -0.0557632 |                 4.47631 | False           |
| PantheonPlusSH0ES | IDSURVEY_51 |  38 |         0.00960541 |       -0.0395037 |                 4.11266 | False           |

Release-contained calibration modes are not added again to `STAT+SYS`; the table compares model-difference residual modes with the projected release budgets.

## Gates

- `frozen_protocol_hash_verified`: **PASS**
- `analytic_path_identity`: **PASS**
- `known_fixed_likelihood_regression_within_1e-6`: **PASS**
- `exact_path_absolute_delta_BIC_at_most_2`: **FAIL**
- `exact_path_favours_over_LCDM`: **FAIL**
- `exact_path_all_covariance_sensitivities_abs_delta_BIC_at_most_2`: **FAIL**
- `exact_path_every_release_abs_delta_BIC_at_most_2`: **FAIL**
- `asserted_Q_matches_spectroscopic_closure`: **FAIL**
- `fixed_h_Compton_frequency_matches_asserted_clock`: **FAIL**
- `Chandrasekhar_proxy_absolute_delta_BIC_at_most_2`: **FAIL**
- `diagnostic_composite_power_is_interior`: **PASS**
- `diagnostic_composite_within_2sigma_of_registered_fixed_branch`: **FAIL**
- `exact_model_difference_modes_within_release_budgets`: **FAIL**
- `photon_transfer_and_source_luminosity_separately_predicted`: **FAIL**
- `action_level_clock_matter_electromagnetic_model_supplied`: **FAIL**
- `no_release_covariance_prior_double_counting`: **PASS**
- `escalate_to_light_curve_level`: **FAIL**

## Claim boundary

This result rejects promotion of the explicit v1 closure. It does not reject every FR, varying-constant or static-universe theory. A revised theory must derive, before another cosmological fit, the redshift identity, dimensionless atomic-frequency law, photon transfer, Type Ia luminosity and the dynamics of the clock field.

The proposed next stage was a Pantheon+ light-curve-level refit. Because the frozen escalation rule fails, that stage is not run.

## Reproduction

```powershell
python github_export/code/rawmu/audit_peter_clock_law_v1_2026_07_19.py
```

## Primary data references

- [Riess et al. SH0ES distance-ladder equations](https://arxiv.org/abs/2112.04510)
- [Pantheon+ full data and light-curve release](https://arxiv.org/abs/2112.03863)
- [Official Pantheon+ data release](https://github.com/PantheonPlusSH0ES/DataRelease)
