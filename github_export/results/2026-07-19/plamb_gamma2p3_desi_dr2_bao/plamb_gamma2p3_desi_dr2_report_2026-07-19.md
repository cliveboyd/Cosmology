# PLAMB gamma = 2.3 observable adapter and DESI DR2 BAO audit

Generated: 2026-07-19 09:01:57

## Decision

**THE PHYSICAL BAO-ADAPTER GATE FAILS. DO NOT START AN OVERNIGHT SAMPLER.**

A mathematically self-consistent late-time shape adapter can be written, but no tested branch simultaneously preserves the registered alpha=0 supernova distance, standard distance duality, the direct clock-law radial distance and an independently predicted sound horizon. The DESI diagnostic is nevertheless evaluated because its scale can be profiled analytically and it provides a useful falsification readout.

## Registered clock law

```text
c(z) / c0                    = 1 + gamma_c z
|dz/dT| / H0                = (1 + z)^p
gamma_c                     = 2.3
p_SN                        = 0.800469 +/- 0.023641
I(z)                        = integral[0,z] (1+gamma_c u)/(1+u)^p du
K(z)                        = dI/dz = (1+gamma_c z)/(1+z)^p
q                           = c0/(H0 r_d)
```

The stationary-ruler branch uses `D_M/r_d=q I` and `D_H/r_d=q K`. This is the direct clock-path mapping, but with alpha=0 it implies `D_L=D_M`, which is incompatible with standard Etherington duality `D_L=(1+z)D_M`.

The metric-duality control instead uses `D_M/r_d=q I/(1+z)` and differentiates this quantity for `D_H`. It preserves standard duality and radial/transverse derivative consistency, but its radial distance is no longer the direct clock kernel `q K`.

The early-Universe sound horizon is not calculated. A common late-time scale `q` is profiled, so this is a BAO-shape test and does not extrapolate `c/c0=1+2.3z` to the drag epoch.

## Numerical validation

| branch                 |   maximum_abs_dDMdz_minus_DH |   minimum_DM_dimensionless |   minimum_DH_dimensionless | positive_to_z2p5   |   covariance_minimum_eigenvalue |   covariance_condition_number |
|:-----------------------|-----------------------------:|---------------------------:|---------------------------:|:-------------------|--------------------------------:|------------------------------:|
| stationary_path        |                  1.69335e-09 |                1.00001e-05 |                   1.00001  | True               |                      0.00578999 |                       116.824 |
| metric_duality_control |                  5.91842e-10 |                9.99997e-06 |                   0.311654 | True               |                      0.00578999 |                       116.824 |

## DESI DR2 likelihood

| model                                          |          p |   p_shift_from_SN_sigma |   q_c0_over_H0rd |   chi2_data |   chi2_SN_prior |   BIC_data_only |   goodness_upper_tail | at_p_bound   |
|:-----------------------------------------------|-----------:|------------------------:|-----------------:|------------:|----------------:|----------------:|----------------------:|:-------------|
| PLAMB_stationary_path_fixed_sn                 |   0.800469 |                 0       |           7.5861 |  33458.1    |          0      |      33460.6    |          0            | False        |
| PLAMB_stationary_path_sn_gaussian_prior        |   2.28026  |                62.5953  |          25.5457 |    695.648  |       3918.18   |        700.777  |          4.62333e-142 | False        |
| PLAMB_stationary_path_free_bao                 |   2.5      |                71.8906  |          28.763  |     33.1341 |          0      |         38.264  |          0.000500469  | True         |
| PLAMB_metric_duality_control_fixed_sn          |   0.800469 |                 0       |          29.4051 |    176.07   |          0      |        178.635  |          2.72905e-31  | False        |
| PLAMB_metric_duality_control_sn_gaussian_prior |   0.878538 |                 3.30233 |          30.9602 |     26.4131 |         10.9054 |         31.543  |          0.00563096   | False        |
| PLAMB_metric_duality_control_free_bao          |   0.885316 |                 3.58905 |          31.095  |     25.4669 |          0      |         30.5968 |          0.00778433   | False        |
| LCDM_flat_Omega_m_free                         | nan        |               nan       |          29.5246 |     10.271  |          0      |         15.4009 |          0.506184     | False        |

For the external fixed-SN prediction, the stationary branch gives chi2=33458.052 and delta BIC=33445.216 relative to flat Lambda-CDM.
The metric-duality control gives chi2=176.070 and delta BIC=163.234.

Allowing `p` to be selected by DESI is not an independent validation. It is reported only to diagnose whether the failure is a fixed-parameter mismatch; boundary solutions are gating failures.

## SN-prior posterior updates

- `stationary_path`: p = 2.280279 +/- 0.008855, 95% [2.262427, 2.297161].
- `metric_duality_control`: p = 0.878551 +/- 0.006675, 95% [0.864992, 0.891153].

## Largest redshift-block contributions

| model                                 |     z | kinds                 |   N |   chi2_contribution |
|:--------------------------------------|------:|:----------------------|----:|--------------------:|
| PLAMB_stationary_path_fixed_sn        | 2.33  | DH_over_rd+DM_over_rd |   2 |          10600.4    |
| PLAMB_stationary_path_fixed_sn        | 0.934 | DM_over_rd+DH_over_rd |   2 |           6546.78   |
| PLAMB_stationary_path_fixed_sn        | 0.51  | DM_over_rd+DH_over_rd |   2 |           5243.59   |
| PLAMB_stationary_path_fixed_sn        | 0.706 | DM_over_rd+DH_over_rd |   2 |           4991.78   |
| PLAMB_stationary_path_fixed_sn        | 0.295 | DV_over_rd            |   1 |           4577.32   |
| PLAMB_stationary_path_fixed_sn        | 1.321 | DM_over_rd+DH_over_rd |   2 |           1267.1    |
| PLAMB_stationary_path_fixed_sn        | 1.484 | DM_over_rd+DH_over_rd |   2 |            231.095  |
| PLAMB_metric_duality_control_fixed_sn | 2.33  | DH_over_rd+DM_over_rd |   2 |            136.302  |
| PLAMB_metric_duality_control_fixed_sn | 0.51  | DM_over_rd+DH_over_rd |   2 |             15.6868 |
| PLAMB_metric_duality_control_fixed_sn | 0.934 | DM_over_rd+DH_over_rd |   2 |             12.7744 |

## Theoretical gates

- `stationary_path_derivative_identity`: **PASS**
- `stationary_path_positive_to_z2p5`: **PASS**
- `metric_duality_control_derivative_identity`: **PASS**
- `metric_duality_control_positive_to_z2p5`: **PASS**
- `single_common_BAO_scale`: **PASS**
- `no_early_time_gamma_extrapolation`: **PASS**
- `stationary_branch_standard_distance_duality`: **FAIL**
- `duality_branch_direct_clock_radial_mapping`: **FAIL**
- `independent_sound_horizon_prediction`: **FAIL**
- `one_branch_satisfies_all_physical_requirements`: **FAIL**

## Statistical gates

- `fixed_stationary_goodness_p_at_least_0p05`: **FAIL**
- `fixed_duality_goodness_p_at_least_0p05`: **FAIL**
- `fixed_stationary_delta_BIC_vs_LCDM_at_most_2`: **FAIL**
- `fixed_duality_delta_BIC_vs_LCDM_at_most_2`: **FAIL**
- `free_stationary_p_not_at_bound`: **FAIL**
- `covariance_positive_definite`: **PASS**

## Claim boundary

The gamma value is outcome-aware and the p prior comes from the same supernova family being tested. DESI is an independent dataset, but the adapter lacks a complete physical ruler and reciprocity derivation. These results therefore reject the present adapter branches; they do not test every possible PLAMB completion.

## Data source

- [DESI DR2 Gaussian BAO data product](https://github.com/CobayaSampler/bao_data/tree/master/desi_bao_dr2)
