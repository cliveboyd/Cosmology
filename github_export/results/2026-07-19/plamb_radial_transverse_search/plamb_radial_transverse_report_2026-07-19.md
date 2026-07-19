# Targeted PLAMB radial/transverse search

Generated: 2026-07-19 11:26:46

## Decision

**THE RADIAL/TRANSVERSE SPLIT FITS FULL DESI BUT FAILS HIGH-REDSHIFT PREDICTION.**

With the SN-trained `p=0.797428918`, separate powers give `b_perpendicular=1.016666951` and `b_parallel=1.778668238`, with DESI `chi-squared=14.463048` for `10` degrees of freedom and `p=0.152899`.
The `z=2.33` leave-one-bin-out prediction fails with `p=0.001712`; consequently the adapter does not pass its registered predictive gate.

The isotropic restriction is worse by `Delta BIC=-4139.079103`. The exponent difference is `66.268450` local-Hessian standard deviations from zero. This conditional curvature is not a discovery significance because the power-law family fails the high-redshift hold-out.

At `z=2.33`, the implied ruler factors are `3.397440` transverse and `8.496792` radial. The joint model remains `Delta BIC=+14.558420` versus flat Lambda-CDM.

## Targeted model

This is an explicitly post hoc follow-up to the missing-degree suite. It tests

$$
\begin{aligned}
k(z)                         &= \frac{1+2.3z}{(1+z)^p}, \\
\chi(z)                    &= \int_0^z k(u)\,du, \\
D_M(z)                       &= \frac{\chi(z)}{(1+z)^{b_\perp}}, \\
D_H(z)                       &= \frac{k(z)}{(1+z)^{b_\parallel}}.
\end{aligned}
$$

The equality `b_perpendicular=b_parallel` is the isotropic-ruler hypothesis. Independent powers allow the radial/transverse ratio to change and therefore form a symmetry-breaking control, not a completed covariant model.

## Supernova profile

|    p_mle |   p_profile_lower_delta_chi2_1 |   p_profile_upper_delta_chi2_1 |   p_sigma_average |   sn_chi2_minimum | optimisation_success   |
|---------:|-------------------------------:|-------------------------------:|------------------:|------------------:|:-----------------------|
| 0.797429 |                        0.78294 |                       0.811935 |         0.0144976 |           3047.67 | True                   |

## DESI fits at SN-trained p

| model           | label                                 | fit_mode            |        p |   b_perpendicular |   b_parallel |   delta_b_parallel_minus_perpendicular |   q_c0_over_H0rd |      chi2 |   N |   N_parameters |   dof |        BIC |   goodness_upper_tail | parameter_at_boundary   | optimisation_success   | optimisation_message                                 |
|:----------------|:--------------------------------------|:--------------------|---------:|------------------:|-------------:|---------------------------------------:|-----------------:|----------:|----:|---------------:|------:|-----------:|----------------------:|:------------------------|:-----------------------|:-----------------------------------------------------|
| baseline        | No ruler evolution                    | BAO_at_SN_trained_p | 0.797429 |          0        |      0       |                               0        |          7.56134 | 33561.8   |  13 |              1 |    12 | 33564.4    |              0        | False                   | True                   | no active ruler parameter                            |
| isotropic       | Isotropic ruler power                 | BAO_at_SN_trained_p | 0.797429 |          1.66346  |      1.66346 |                               0        |         35.4426  |  4156.11  |  13 |              2 |    11 |  4161.24   |              0        | False                   | True                   | CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL     |
| transverse_only | Transverse-only power                 | BAO_at_SN_trained_p | 0.797429 |         -0.511467 |      0       |                               0.511467 |          6.48746 | 30300.3   |  13 |              2 |    11 | 30305.4    |              0        | False                   | True                   | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| radial_only     | Radial-only power                     | BAO_at_SN_trained_p | 0.797429 |          0        |      1.223   |                               1.223    |         15.2394  |  8173.87  |  13 |              2 |    11 |  8179      |              0        | False                   | True                   | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| separate        | Separate transverse and radial powers | BAO_at_SN_trained_p | 0.797429 |          1.01667  |      1.77867 |                               0.762001 |         30.1698  |    14.463 |  13 |              3 |    10 |    22.1579 |              0.152899 | False                   | True                   | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |

## Joint fits

| model           | label                                 | fit_mode     |          p |   b_perpendicular |   b_parallel |   delta_b_parallel_minus_perpendicular |   q_c0_over_H0rd |   sn_chi2 |   bao_chi2 |     chi2 |    N |   N_parameters |   dof |      BIC |   goodness_upper_tail |   parameter_at_boundary | optimisation_success   | optimisation_message                                 |   omega_m |
|:----------------|:--------------------------------------|:-------------|-----------:|------------------:|-------------:|---------------------------------------:|-----------------:|----------:|-----------:|---------:|-----:|---------------:|------:|---------:|----------------------:|------------------------:|:-----------------------|:-----------------------------------------------------|----------:|
| baseline        | No ruler evolution                    | joint_SN_BAO |   2.06905  |         0         |     0        |                              0         |          22.4962 |   9943.15 |  2324.86   | 12268    | 3435 |              5 |  3430 | 12308.7  |          0            |                       0 | True                   | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH | nan       |
| isotropic       | Isotropic ruler power                 | joint_SN_BAO |   1.17593  |         1.27691   |     1.27691  |                              0         |          33.9381 |   3708.15 |  2781.45   |  6489.6  | 3435 |              6 |  3429 |  6538.45 |          2.72004e-192 |                       0 | True                   | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH | nan       |
| transverse_only | Transverse-only power                 | joint_SN_BAO |   2.05261  |        -0.0868871 |     0        |                              0.0868871 |          21.5243 |   9775.87 |  2413.66   | 12189.5  | 3435 |              6 |  3429 | 12238.4  |          0            |                       0 | True                   | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH | nan       |
| radial_only     | Radial-only power                     | joint_SN_BAO |   1.44837  |         0         |     0.806913 |                              0.806913  |          20.4739 |   4956.23 |  2728.18   |  7684.41 | 3435 |              6 |  3429 |  7733.26 |          0            |                       0 | True                   | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH | nan       |
| separate        | Separate transverse and radial powers | joint_SN_BAO |   0.798319 |         1.01608   |     1.77775  |                              0.761667  |          30.169  |   3047.68 |    14.4555 |  3062.13 | 3435 |              7 |  3428 |  3119.13 |          0.999998     |                       0 | True                   | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH | nan       |
| lcdm_control    | Flat Lambda-CDM control               | joint_SN_BAO | nan        |       nan         |   nan        |                            nan         |          29.8483 |   3050.98 |    12.8792 |  3063.86 | 3435 |              5 |  3430 |  3104.57 |        nan            |                     nan | True                   | nan                                                  |   0.31175 |

## Relation to the Lambda-CDM ratio

|     z |   PLAMB_separate_DH_over_DM |   LCDM_DH_over_DM |   fractional_PLAMB_minus_LCDM |
|------:|----------------------------:|------------------:|------------------------------:|
| 0.51  |                    1.71044  |          1.68645  |                    0.014227   |
| 0.706 |                    1.15249  |          1.14122  |                    0.0098764  |
| 0.934 |                    0.805259 |          0.800219 |                    0.00629857 |
| 1.321 |                    0.503255 |          0.50188  |                    0.00274054 |
| 1.484 |                    0.427013 |          0.426319 |                    0.00162668 |
| 2.33  |                    0.219286 |          0.220377 |                   -0.00494926 |

Across the paired DESI bins, the fitted PLAMB ratio differs from the joint-fit Lambda-CDM ratio by at most `1.422703%`, with RMS fractional difference `0.789806%`. The radial/transverse powers are therefore acting largely as an empirical reparameterisation of conventional distance geometry.

## Local identifiability

| hessian_available   | positive_definite   |   condition_number |   sigma_b_perpendicular |   sigma_b_parallel |   correlation |   difference_b_parallel_minus_perpendicular |   sigma_difference |   difference_significance_sigma |
|:--------------------|:--------------------|-------------------:|------------------------:|-------------------:|--------------:|--------------------------------------------:|-------------------:|--------------------------------:|
| True                | True                |             3.0816 |               0.0121165 |          0.0107414 |      0.499306 |                                    0.762001 |          0.0114987 |                         66.2684 |

## SN-profile sensitivity

| p_case                     | model    | label                                 | fit_mode            |        p |   b_perpendicular |   b_parallel |   delta_b_parallel_minus_perpendicular |   q_c0_over_H0rd |   chi2 |   N |   N_parameters |   dof |     BIC |   goodness_upper_tail | parameter_at_boundary   | optimisation_success   | optimisation_message                                 |
|:---------------------------|:---------|:--------------------------------------|:--------------------|---------:|------------------:|-------------:|---------------------------------------:|-----------------:|-------:|----:|---------------:|------:|--------:|----------------------:|:------------------------|:-----------------------|:-----------------------------------------------------|
| profile_lower_delta_chi2_1 | separate | Separate transverse and radial powers | BAO_at_SN_trained_p | 0.78294  |           1.02626 |      1.79369 |                               0.767432 |          30.1834 | 14.588 |  13 |              3 |    10 | 22.2828 |              0.147821 | False                   | True                   | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| MLE                        | separate | Separate transverse and radial powers | BAO_at_SN_trained_p | 0.797429 |           1.01667 |      1.77867 |                               0.762001 |          30.1698 | 14.463 |  13 |              3 |    10 | 22.1579 |              0.152899 | False                   | True                   | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| profile_upper_delta_chi2_1 | separate | Separate transverse and radial powers | BAO_at_SN_trained_p | 0.811935 |           1.00708 |      1.76364 |                               0.756551 |          30.1563 | 14.342 |  13 |              3 |    10 | 22.0368 |              0.157957 | False                   | True                   | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |

## Leave-one-redshift-bin-out predictions

|   held_z |   N_test |   N_train |   b_perpendicular_train |   b_parallel_train |   delta_b_train |   q_train |   conditional_chi2 |   conditional_p | training_optimisation_success   |
|---------:|---------:|----------:|------------------------:|-------------------:|----------------:|----------:|-------------------:|----------------:|:--------------------------------|
|    0.295 |        1 |        12 |                1.02265  |            1.7835  |        0.760847 |   30.3153 |            1.37272 |      0.241344   | True                            |
|    0.51  |        2 |        11 |                1.01791  |            1.77556 |        0.757653 |   30.1458 |            5.48887 |      0.0642847  | True                            |
|    0.706 |        2 |        11 |                1.02278  |            1.78106 |        0.758278 |   30.3217 |            4.3065  |      0.116106   | True                            |
|    0.934 |        2 |        11 |                1.01228  |            1.77938 |        0.767102 |   30.052  |            3.79793 |      0.149723   | True                            |
|    1.321 |        2 |        11 |                1.01573  |            1.78183 |        0.766099 |   30.1584 |            1.81514 |      0.403503   | True                            |
|    1.484 |        2 |        11 |                1.01879  |            1.77933 |        0.760546 |   30.1878 |            1.00533 |      0.604917   | True                            |
|    2.33  |        2 |        11 |                0.992066 |            1.75242 |        0.760357 |   29.7592 |           12.7404  |      0.00171185 | True                            |

### Held-out observable components

|   held_z | kind       |   observed |   predicted_from_training |   raw_residual_observed_minus_predicted |   conditional_residual |   conditional_sigma_marginal |   conditional_marginal_pull |   cholesky_whitened_component |
|---------:|:-----------|-----------:|--------------------------:|----------------------------------------:|-----------------------:|-----------------------------:|----------------------------:|------------------------------:|
|    0.295 | DV_over_rd |    7.94168 |                   8.03083 |                              -0.0891518 |             -0.0891518 |                     0.076092 |                   -1.17163  |                     -1.17163  |
|    0.51  | DM_over_rd |   13.5876  |                  13.2402  |                               0.347386  |              0.347386  |                     0.168367 |                    2.06327  |                      2.06327  |
|    0.51  | DH_over_rd |   21.8629  |                  22.6872  |                              -0.824266  |             -0.824266  |                     0.428868 |                   -1.92196  |                     -1.10986  |
|    0.706 | DM_over_rd |   17.3507  |                  17.3791  |                              -0.0283707 |             -0.0283707 |                     0.179931 |                   -0.157675 |                     -0.157675 |
|    0.706 | DH_over_rd |   19.4553  |                  20.0691  |                              -0.613785  |             -0.613785  |                     0.33387  |                   -1.83839  |                     -2.06921  |
|    0.934 | DM_over_rd |   21.5756  |                  21.5435  |                               0.0321168 |              0.0321168 |                     0.161782 |                    0.198519 |                      0.198519 |
|    0.934 | DH_over_rd |   17.6415  |                  17.2898  |                               0.351651  |              0.351651  |                     0.201043 |                    1.74913  |                      1.93869  |
|    1.321 | DM_over_rd |   27.6009  |                  27.6788  |                              -0.0779051 |             -0.0779051 |                     0.324556 |                   -0.240036 |                     -0.240036 |
|    1.321 | DH_over_rd |   14.176   |                  13.8815  |                               0.294525  |              0.294525  |                     0.224551 |                    1.31161  |                      1.32572  |
|    1.484 | DM_over_rd |   30.5119  |                  29.8785  |                               0.633439  |              0.633439  |                     0.763558 |                    0.82959  |                      0.82959  |
|    1.484 | DH_over_rd |   12.817   |                  12.7754  |                               0.0416031 |              0.0416031 |                     0.518012 |                    0.080313 |                      0.563125 |
|    2.33  | DH_over_rd |    8.63155 |                   8.80755 |                              -0.176     |             -0.176     |                     0.101062 |                   -1.7415   |                     -1.7415   |
|    2.33  | DM_over_rd |   38.989   |                  40.0853  |                              -1.09635   |             -1.09635   |                     0.531682 |                   -2.06203  |                     -3.11569  |

These are genuine DESI block hold-outs: the two powers and common BAO scale are fitted without the held redshift block, and the test score uses the Gaussian conditional covariance.

## Gates

- `separate_BAO_goodness_p_at_least_0p05`: **True**
- `separate_beats_isotropic_by_10_BIC`: **True**
- `parameters_not_at_boundary`: **True**
- `positive_identifiable_hessian`: **True**
- `stable_across_SN_profile_interval_max_range_at_most_0p15`: **True**
- `all_LOO_conditional_p_at_least_0p01`: **False**
- `joint_within_10_BIC_of_LCDM`: **False**
- `preserves_background_isotropy`: **False**
- `action_level_radial_transverse_source_supplied`: **False**
- `maximum_power_range_across_SN_profile`: **0.030052606707564822**
- `delta_BIC_separate_minus_isotropic_BAO`: **-4139.079103152821**
- `delta_BIC_separate_joint_minus_LCDM`: **14.55841990789213**
- `radial_transverse_adapter_validated`: **False**
- `physical_degree_identified`: **False**
- `further_physical_sampling_authorised`: **False**

## Validation

| check                               |         value |      target |   absolute_tolerance | passed   |
|:------------------------------------|--------------:|------------:|---------------------:|:---------|
| SN p regression                     |    0.797429   |    0.797429 |               5e-05  | True     |
| Rank-two b_perpendicular regression |    1.01667    |    1.01667  |               5e-05  | True     |
| Rank-two b_parallel regression      |    1.77867    |    1.77867  |               5e-05  | True     |
| Rank-two BAO chi2 regression        |   14.463      |   14.463    |               0.0005 | True     |
| LCDM joint BIC regression           | 3104.57       | 3104.57     |               0.02   | True     |
| BAO covariance positive definite    |    0.00578999 |    0        |               0      | True     |

## Interpretation

If interpreted literally as ruler evolution, the unequal powers select the observer's radial direction and violate background isotropy. They may instead be diagnosing an incorrect adapter between clock time, redshift and the two BAO observables. Distances alone cannot distinguish those possibilities. A physical continuation requires an action that predicts the two modes and independent angular or multipole evidence for the associated symmetry breaking.

No SU(2), antimatter or new-field claim follows from this test case.
