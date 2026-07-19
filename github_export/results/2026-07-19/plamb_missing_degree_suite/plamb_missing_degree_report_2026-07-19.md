# PLAMB missing-degree-of-freedom diagnostic suite

Generated: 2026-07-19 10:28:41

## Decision

This suite searches for an **observable signature** of an omitted degree of freedom. It cannot identify a fundamental field without an action and independent coupling law.

The largest numerical one-degree compromise is `redshift_remap` with joint `Delta BIC=-6525.609301` relative to the static PLAMB baseline and DESI `chi-squared=1199.395662` at its joint fit. It is boundary-seeking and fails the DESI goodness gate.

Promoted observational signatures: `0`. Physically identified new degrees of freedom: `0`.
The registered result therefore disfavours the hypothesis that one omitted scalar background variable repairs the existing SN--BAO closure. The post hoc rank test requires two independent observable modes: a common scale evolution and a separate radial/transverse distortion.

## Registered observable basis

The common static path is

$$
\begin{aligned}
C(z)                         &= 1+2.3z, \\
k(z)                         &= \frac{C(z)}{(1+z)^p}, \\
\chi(z)                    &= \int_0^z k(u)\,du, \\
D_M                          &= \chi, \\
D_H                          &= k.
\end{aligned}
$$

Each extension adds at most one shape variable beyond `p`:

| candidate              | parameters      | symmetry                        | interpretation                                                                     |
|:-----------------------|:----------------|:--------------------------------|:-----------------------------------------------------------------------------------|
| plamb_baseline         | p               | homogeneous_isotropic           | Registered gamma=2.3 static clock-path closure.                                    |
| common_conformal_power | p, s_conformal  | homogeneous_isotropic           | Observable proxy for an omitted A(T,X) or common unit-scaling factor.              |
| scalar_ruler_power     | p, b_ruler      | homogeneous_isotropic           | Comoving ruler B(z)=(1+z)^b; invisible to supernovae.                              |
| spatial_curvature      | p, omega_k      | homogeneous_isotropic           | Constant-curvature transverse map S_K(chi).                                        |
| redshift_remap         | p, nu_redshift  | homogeneous_isotropic           | Clock/atomic proxy y=[(1+z)^nu-1]/nu with dy/dz explicit.                          |
| clock_quadratic        | p, gamma_2      | homogeneous_isotropic           | One smooth nonlinear term in C(z), retaining C(0)=1.                               |
| transverse_focusing    | p, f_focus      | homogeneous_isotropic_effective | Jacobi-map proxy that changes transverse and luminosity distance, not radial path. |
| anisotropic_ruler      | p, e_anisotropy | breaks_spatial_isotropy         | B_perp=(1+z)^(-e/3), B_parallel=(1+z)^(2e/3).                                      |
| radial_closure_slip    | p, h_radial     | breaks_distance_closure         | Control allowing D_H to move independently of D_M.                                 |
| lcdm_control           | omega_m         | homogeneous_isotropic_control   | Standard flat Lambda-CDM distance relation.                                        |

The anisotropic-ruler and radial-slip rows are controls. A good fit from either would locate the missing observable freedom in the radial/transverse closure, but would not establish a covariant or isotropic theory.

## Joint ranking

| candidate              |          p | theta_name   |   theta_value |   sn_chi2_at_fit |   bao_chi2_at_fit |      BIC |   delta_BIC_joint_vs_PLAMB |   delta_BIC_joint_vs_LCDM |   bao_p_at_joint | parameter_at_boundary   |
|:-----------------------|-----------:|:-------------|--------------:|-----------------:|------------------:|---------:|---------------------------:|--------------------------:|-----------------:|:------------------------|
| lcdm_control           | nan        | omega_m      |      0.31175  |          3050.98 |           12.8792 |  3104.57 |                   -9204.15 |                      0    |     0.301296     | False                   |
| redshift_remap         |   2.5      | nu_redshift  |      2.5      |          4534.86 |         1199.4    |  5783.1  |                   -6525.61 |                   2678.54 |     1.94519e-251 | True                    |
| transverse_focusing    |   2.23417  | f_focus      |      0.739752 |          4077.83 |         2347.61   |  6474.29 |                   -5834.42 |                   3369.73 |     0            | False                   |
| scalar_ruler_power     |   1.17593  | b_ruler      |      1.27691  |          3708.16 |         2781.45   |  6538.45 |                   -5770.26 |                   3433.88 |     0            | False                   |
| clock_quadratic        |   0.849395 | gamma_2      |     -3        |          5625.1  |         1332.53   |  7006.47 |                   -5302.24 |                   3901.9  |     3.65235e-280 | True                    |
| radial_closure_slip    |   1.44837  | h_radial     |     -0.806913 |          4956.23 |         2728.18   |  7733.26 |                   -4575.46 |                   4628.69 |     0            | False                   |
| common_conformal_power |   2.5      | s_conformal  |      0.498256 |          5725.79 |         3173.79   |  8948.43 |                   -3360.28 |                   5843.86 |     0            | True                    |
| anisotropic_ruler      |   1.74622  | e_anisotropy |      0.544168 |          6998.4  |         3014.6    | 10061.8  |                   -2246.87 |                   6957.27 |     0            | False                   |
| spatial_curvature      |   2.23898  | omega_k      |      0.994978 |          8828.91 |         2356.84   | 11234.6  |                   -1074.11 |                   8130.03 |     0            | False                   |
| plamb_baseline         |   2.06905  |              |    nan        |          9943.14 |         2324.86   | 12308.7  |                       0    |                   9204.15 |     0            | False                   |

## Cross-probe transfer

| candidate              | theta_name   | SN_invisible_BAO_parameters   | conditional_parameters                                         | conditional_optimisation_success   |   bao_chi2_after_SN_training |   bao_p_after_SN_training |   sn_chi2_at_bao_fit |   delta_sn_chi2_at_bao_fit |   sn_bao_theta_tension_sigma |
|:-----------------------|:-------------|:------------------------------|:---------------------------------------------------------------|:-----------------------------------|-----------------------------:|--------------------------:|---------------------:|---------------------------:|-----------------------------:|
| plamb_baseline         |              | []                            | {"p": 0.7974288977143246}                                      | True                               |                   33561.8    |                 0         |             14932.4  |                 11884.8    |                    nan       |
| common_conformal_power | s_conformal  | []                            | {"p": 0.7678604720851119, "s_conformal": -0.01777583846028693} | True                               |                   34052.6    |                 0         |             15946.8  |                 12899.1    |                    nan       |
| scalar_ruler_power     | b_ruler      | ["b_ruler"]                   | {"b_ruler": 1.6634571692936535, "p": 0.7974288977143246}       | True                               |                    4156.11   |                 0         |             14932.4  |                 11884.8    |                    nan       |
| spatial_curvature      | omega_k      | []                            | {"omega_k": 0.0, "p": 0.8004694772058842}                      | True                               |                   33458.1    |                 0         |             15060.4  |                 12012.7    |                    nan       |
| redshift_remap         | nu_redshift  | []                            | {"nu_redshift": 1.5059791152922406, "p": 1.325347316866768}    | True                               |                   30758.8    |                 0         |             16077    |                 13030.8    |                    nan       |
| clock_quadratic        | gamma_2      | []                            | {"gamma_2": 0.37384169285782276, "p": 0.8706596077231136}      | True                               |                   33732.2    |                 0         |             15445.3  |                 12398.3    |                    nan       |
| transverse_focusing    | f_focus      | []                            | {"f_focus": -0.0059616792764757365, "p": 0.7901313452064006}   | True                               |                   33867.3    |                 0         |             14582.8  |                 11535.1    |                    nan       |
| anisotropic_ruler      | e_anisotropy | ["e_anisotropy"]              | {"e_anisotropy": 1.0471550081944505, "p": 0.7974288977143246}  | True                               |                   15986.8    |                 0         |             14932.4  |                 11884.8    |                    nan       |
| radial_closure_slip    | h_radial     | ["h_radial"]                  | {"h_radial": -1.2230021586446178, "p": 0.7974288977143246}     | True                               |                    8173.87   |                 0         |             14932.4  |                 11884.8    |                    nan       |
| lcdm_control           | omega_m      | []                            | {"omega_m": 0.33067964259081406}                               | True                               |                      23.4467 |                 0.0241645 |              3058.11 |                    10.4429 |                      2.43976 |

A candidate that exists only in the BAO likelihood cannot be validated by supernova-to-BAO transfer. For those rows, the suite fixes the supernova-trained `p` and fits only the SN-invisible BAO variable. A scalar candidate visible to supernovae must predict DESI without retuning its extra variable.

## Exploratory rank-two localisation

After every registered one-degree extension failed, a deliberately post hoc rank test combined a common ruler power with an independent radial slip while fixing `p` to the supernova-only result:

$$
\begin{aligned}
D_M(z)                      &= \frac{\chi(z)}{(1+z)^b}, \\
D_H(z)                      &= \frac{k(z)}{(1+z)^b}(1+z)^h, \\
B_\perp(z)                 &= (1+z)^b, \\
B_\parallel(z)             &= (1+z)^{b-h}.
\end{aligned}
$$

| analysis_status                                            |   p_fixed_from_SN |   b_common_ruler |   h_independent_radial |   B_parallel_power_b_minus_h |   B_perpendicular_z2p33 |   B_parallel_z2p33 |   AP_ratio_multiplier_z2p33 |   bao_q_c0_over_H0rd |   bao_chi2 |   bao_dof |   bao_goodness_upper_tail |   bao_BIC |   delta_bao_BIC_vs_LCDM |   sn_chi2 |   joint_chi2 |   joint_N_parameters |   joint_BIC |   delta_joint_BIC_vs_LCDM | optimisation_success   | parameter_at_boundary   | preserves_background_isotropy   | preserves_covariant_distance_closure   | physical_degree_identified   |
|:-----------------------------------------------------------|------------------:|-----------------:|-----------------------:|-----------------------------:|------------------------:|-------------------:|----------------------------:|---------------------:|-----------:|----------:|--------------------------:|----------:|------------------------:|----------:|-------------:|---------------------:|------------:|--------------------------:|:-----------------------|:------------------------|:--------------------------------|:---------------------------------------|:-----------------------------|
| post_hoc_rank2_localisation_after_registered_rank1_failure |          0.797429 |          1.01667 |              -0.762001 |                      1.77867 |                 3.39744 |            8.49679 |                     0.39985 |              30.1698 |     14.463 |        10 |                  0.152899 |   22.1579 |                 6.75695 |   3047.67 |      3062.14 |                    7 |     3119.13 |                   14.5622 | True                   | False                   | False                           | False                                  | False                        |

This reaches DESI `chi-squared=14.463048` with `p=0.797429`, `b=1.016667` and `h=-0.762001`. It shows that two observational modes can span most of the discrepancy. At `z=2.33` they require transverse and radial ruler factors of `3.397440` and `8.496792`. The joint result remains `Delta BIC=+14.562195` versus Lambda-CDM. Because the second mode breaks the covariant radial/transverse closure and has no action-level source, this is localisation evidence, not a new physical model.

## Identifiability

| candidate              | mode   | active_parameters     | positive_definite   |   condition_number |   maximum_absolute_correlation |
|:-----------------------|:-------|:----------------------|:--------------------|-------------------:|-------------------------------:|
| plamb_baseline         | sn     | ["p"]                 | True                |            1       |                    0           |
| plamb_baseline         | bao    | ["p"]                 | False               |          nan       |                  nan           |
| plamb_baseline         | joint  | ["p"]                 | True                |            1       |                    0           |
| common_conformal_power | sn     | ["p", "s_conformal"]  | True                |         1559.55    |                    0.998356    |
| common_conformal_power | bao    | ["p", "s_conformal"]  | False               |          nan       |                  nan           |
| common_conformal_power | joint  | ["p", "s_conformal"]  | False               |          nan       |                  nan           |
| scalar_ruler_power     | sn     | ["p"]                 | True                |            1       |                    0           |
| scalar_ruler_power     | bao    | ["p", "b_ruler"]      | False               |          nan       |                  nan           |
| scalar_ruler_power     | joint  | ["p", "b_ruler"]      | True                |           11.6128  |                    0.835069    |
| spatial_curvature      | sn     | ["p", "omega_k"]      | True                |            6.12818 |                    0.715001    |
| spatial_curvature      | bao    | ["p", "omega_k"]      | False               |          nan       |                  nan           |
| spatial_curvature      | joint  | ["p", "omega_k"]      | True                |           35.7755  |                    0.762224    |
| redshift_remap         | sn     | ["p", "nu_redshift"]  | True                |          978.255   |                    0.99795     |
| redshift_remap         | bao    | ["p", "nu_redshift"]  | False               |          nan       |                  nan           |
| redshift_remap         | joint  | ["p", "nu_redshift"]  | False               |          nan       |                  nan           |
| clock_quadratic        | sn     | ["p", "gamma_2"]      | True                |         1188.94    |                    0.98804     |
| clock_quadratic        | bao    | ["p", "gamma_2"]      | False               |          nan       |                  nan           |
| clock_quadratic        | joint  | ["p", "gamma_2"]      | False               |          nan       |                  nan           |
| transverse_focusing    | sn     | ["p", "f_focus"]      | True                |          125.123   |                    0.983366    |
| transverse_focusing    | bao    | ["p", "f_focus"]      | False               |          nan       |                  nan           |
| transverse_focusing    | joint  | ["p", "f_focus"]      | True                |            2.06515 |                    0.346052    |
| anisotropic_ruler      | sn     | ["p"]                 | True                |            1       |                    0           |
| anisotropic_ruler      | bao    | ["p", "e_anisotropy"] | False               |          nan       |                  nan           |
| anisotropic_ruler      | joint  | ["p", "e_anisotropy"] | True                |            3.59313 |                    0.559851    |
| radial_closure_slip    | sn     | ["p"]                 | True                |            1       |                    0           |
| radial_closure_slip    | bao    | ["p", "h_radial"]     | False               |          nan       |                  nan           |
| radial_closure_slip    | joint  | ["p", "h_radial"]     | True                |            4.67412 |                    0.646274    |
| lcdm_control           | sn     | ["omega_m"]           | True                |            1       |                    2.22045e-16 |
| lcdm_control           | bao    | ["omega_m"]           | True                |            1       |                    2.22045e-16 |
| lcdm_control           | joint  | ["omega_m"]           | True                |            1       |                    0           |

## Largest conditional DESI stress

| candidate              |     z |   N |   conditional_chi2 |   conditional_p |
|:-----------------------|------:|----:|-------------------:|----------------:|
| anisotropic_ruler      | 2.33  |   2 |         1406.84    |    3.22386e-306 |
| clock_quadratic        | 0.295 |   1 |          388.888   |    1.44507e-86  |
| common_conformal_power | 2.33  |   2 |         1625.83    |    0            |
| lcdm_control           | 0.706 |   2 |            4.52507 |    0.104086     |
| plamb_baseline         | 2.33  |   2 |         1326.38    |    9.55918e-289 |
| radial_closure_slip    | 2.33  |   2 |         1206.56    |    9.95724e-263 |
| redshift_remap         | 0.295 |   1 |          453.254   |    1.41238e-100 |
| scalar_ruler_power     | 2.33  |   2 |         1020.59    |    2.4067e-222  |
| spatial_curvature      | 2.33  |   2 |         1327.84    |    4.6161e-289  |
| transverse_focusing    | 2.33  |   2 |          790.964   |    1.7555e-172  |

## Promotion gates

| candidate              | symmetry_class                  |   delta_BIC_joint_vs_PLAMB |   delta_BIC_joint_vs_LCDM |   bao_chi2_at_joint |   bao_p_at_joint | delta_BIC_joint_at_most_minus_10   | joint_BAO_goodness_p_at_least_0p05   | SN_trained_BAO_transfer   | joint_hessian_identifiable   | parameter_not_at_boundary   | preserves_background_symmetry   | within_10_BIC_of_LCDM   | action_level_derivation_supplied   | observational_signature_promoted   | physical_degree_identified   |
|:-----------------------|:--------------------------------|---------------------------:|--------------------------:|--------------------:|-----------------:|:-----------------------------------|:-------------------------------------|:--------------------------|:-----------------------------|:----------------------------|:--------------------------------|:------------------------|:-----------------------------------|:-----------------------------------|:-----------------------------|
| plamb_baseline         | homogeneous_isotropic           |                       0    |                   9204.15 |           2324.86   |     0            | False                              | False                                | False                     | True                         | True                        | True                            | False                   | False                              | False                              | False                        |
| common_conformal_power | homogeneous_isotropic           |                   -3360.28 |                   5843.86 |           3173.79   |     0            | True                               | False                                | False                     | False                        | False                       | True                            | False                   | False                              | False                              | False                        |
| scalar_ruler_power     | homogeneous_isotropic           |                   -5770.26 |                   3433.88 |           2781.45   |     0            | True                               | False                                | False                     | True                         | True                        | True                            | False                   | False                              | False                              | False                        |
| spatial_curvature      | homogeneous_isotropic           |                   -1074.11 |                   8130.03 |           2356.84   |     0            | True                               | False                                | False                     | True                         | True                        | True                            | False                   | False                              | False                              | False                        |
| redshift_remap         | homogeneous_isotropic           |                   -6525.61 |                   2678.54 |           1199.4    |     1.94519e-251 | True                               | False                                | False                     | False                        | False                       | True                            | False                   | False                              | False                              | False                        |
| clock_quadratic        | homogeneous_isotropic           |                   -5302.24 |                   3901.9  |           1332.53   |     3.65235e-280 | True                               | False                                | False                     | False                        | False                       | True                            | False                   | False                              | False                              | False                        |
| transverse_focusing    | homogeneous_isotropic_effective |                   -5834.42 |                   3369.73 |           2347.61   |     0            | True                               | False                                | False                     | True                         | True                        | True                            | False                   | False                              | False                              | False                        |
| anisotropic_ruler      | breaks_spatial_isotropy         |                   -2246.87 |                   6957.27 |           3014.6    |     0            | True                               | False                                | False                     | True                         | True                        | False                           | False                   | False                              | False                              | False                        |
| radial_closure_slip    | breaks_distance_closure         |                   -4575.46 |                   4628.69 |           2728.18   |     0            | True                               | False                                | False                     | True                         | True                        | False                           | False                   | False                              | False                              | False                        |
| lcdm_control           | homogeneous_isotropic_control   |                   -9204.15 |                      0    |             12.8792 |     0.301296     | True                               | True                                 | False                     | True                         | True                        | True                            | True                    | True                               | False                              | False                        |

A physical identification additionally requires an action-level derivation of the candidate and its clock, photon, atomic and ruler couplings. None of the PLAMB proxy extensions supplies that derivation.

## Regression checks

| check                                      |          value |       target |   absolute_tolerance | passed   |
|:-------------------------------------------|---------------:|-------------:|---------------------:|:---------|
| SN sample size                             |  3422          |  3422        |                0     | True     |
| PLAMB SN p regression                      |     0.797429   |     0.797429 |                5e-05 | True     |
| PLAMB SN chi2 regression                   |  3047.67       |  3047.67     |                0.02  | True     |
| LCDM SN Omega_m regression                 |     0.33068    |     0.33068  |                5e-05 | True     |
| LCDM SN chi2 regression                    |  3047.66       |  3047.66     |                0.02  | True     |
| LCDM BAO Omega_m regression                |     0.297462   |     0.297462 |                5e-05 | True     |
| LCDM BAO chi2 regression                   |    10.271      |    10.271    |                0.02  | True     |
| Fixed PLAMB BAO chi2 regression            | 33458.1        | 33458.1      |                0.05  | True     |
| BAO covariance minimum eigenvalue positive |     0.00578999 |     0        |                0     | True     |

## Claim boundary

A lower BIC identifies where the present observable closure is deficient; it does not identify SU(2), antimatter, a new field or broken Lorentz symmetry. Candidates that break isotropy must also predict direction-dependent BAO and sky multipoles before they can be interpreted physically.

## Primary references

- [General disformal scalar-tensor transformations](https://arxiv.org/abs/1412.6210)
- [Homogeneous and isotropic SU(2) cosmology](https://arxiv.org/abs/1012.2861)
- [Etherington distance-duality assumptions](https://arxiv.org/abs/1612.08784)
- [DESI DR2 publications and products](https://data.desi.lbl.gov/doc/papers/dr2/)
