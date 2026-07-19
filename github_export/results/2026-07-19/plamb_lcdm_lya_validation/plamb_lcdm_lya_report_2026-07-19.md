# PLAMB / Lambda-CDM equivalence and Lyman-alpha validation audit

Generated: 2026-07-19 12:42:41

## Decision

**THE POWERS CLOSELY APPROXIMATE LAMBDA-CDM BUT MISS THE STRICT EQUIVALENCE GATE; THE DR2 FAILURE IS CROSS-CORRELATION-LED AND SURVIVES THE MULTIPOLE CONTROL.**

The analytic geometry audit gives `chi-squared=4.534951` for the best constant-power representation of the lower-redshift-trained Lambda-CDM curve, with `b_perpendicular=1.030438429` and `b_parallel=1.804606082`. The maximum absolute observable difference is `1.485301%`.

This is a close approximation with geometry goodness `p=0.920006`, not an exact equivalence: the Lambda-CDM-equivalent powers are `2.416077` fitted-parameter standard deviations from the powers selected by the observed full DESI vector. The strict pre-run equivalence gate therefore fails.

The frozen PLAMB prediction remains rejected by the official compressed DESI DR2 combined posterior at `p=0.00170496`. The auto/cross split localises the stronger stress to `DESI_DR2_cross_plus_0p3pct_theory_budget` with `p=0.00297504`.

The alternative DR2 multipole-method control gives `p=0.00516994`; DESI DR1 gives `p=0.0928341`; and the native eBOSS combined grid gives HPD tail mass `0.079816` (`Wilks p=0.07554`). These are controls, not independent confirmations of a new radial/transverse degree of freedom.

The multipole result shows that the DR2 rejection is stable to a substantially different compression and covariance method. It is not an independent-sky replication because it uses the same DR2 spectra. Conversely, the 1% rejection threshold is not crossed by DESI DR1 or native eBOSS DR16.

Relative to the equally frozen Lambda-CDM control, PLAMB has predictive `Delta chi-squared=+8.755836` for the DR2 combined summary and `+7.399697` for the multipole summary. The external eBOSS native grid instead gives `-2.693436`. The sign change across surveys is further evidence against interpreting the DR2 residual as a universal new degree of freedom.

The native DESI DR2 correlation likelihood is not in the current public release. The DR2 rows below therefore use collaboration-published posterior compressions. The unavailable native rerun is retained as an explicit failed availability gate.

## Frozen protocol

The raw-MU supernova likelihood fixes `p`; all DESI points below `z=2` then fix the two powers and common BAO scale. No Lyman-alpha measurement is used to retune either model before it is scored.

$$
\begin{aligned}
k(z)                         &= \frac{1+2.3z}{(1+z)^p}, \\
\chi(z)                    &= \int_0^z k(u)\,du, \\
D_M(z)                       &= q\,\frac{\chi(z)}{(1+z)^{b_\perp}}, \\
D_H(z)                       &= q\,\frac{k(z)}{(1+z)^{b_\parallel}}.
\end{aligned}
$$

The comparator is flat Lambda-CDM with `Omega_m` fitted jointly to the same raw-MU supernova blocks and lower-redshift DESI BAO data, while its common BAO scale is profiled.

## Training fits

| model                   |          p |   b_perpendicular |   b_parallel |       q |    omega_m |
|:------------------------|-----------:|------------------:|-------------:|--------:|-----------:|
| PLAMB radial/transverse |   0.797429 |          0.992066 |      1.75242 | 29.7592 | nan        |
| flat Lambda-CDM         | nan        |        nan        |    nan       | 29.9316 |   0.318484 |

## Analytic equivalence

|   p_SN_trained |   omega_m_LambdaCDM_lower_z_joint |   q_LambdaCDM_lower_z_joint |   b_perpendicular_equivalent |   b_parallel_equivalent |   delta_b_equivalent |   q_PLAMB_equivalent |   chi2_PLAMB_to_LambdaCDM_geometry_best |   dof_geometry_vector_minus_profiled_parameters |   goodness_p_geometry |   chi2_at_full_data_fitted_powers |   q_at_full_data_fitted_powers |   chi2_at_lower_z_fitted_powers |   q_at_lower_z_fitted_powers |   full_data_b_perpendicular |   full_data_b_parallel |   equivalent_minus_data_parameter_distance_sigma |   maximum_absolute_fractional_geometry_difference |   RMS_fractional_geometry_difference | optimisation_success   | optimisation_message                                 | data_hessian_positive_definite   |
|---------------:|----------------------------------:|----------------------------:|-----------------------------:|------------------------:|---------------------:|---------------------:|----------------------------------------:|------------------------------------------------:|----------------------:|----------------------------------:|-------------------------------:|--------------------------------:|-----------------------------:|----------------------------:|-----------------------:|-------------------------------------------------:|--------------------------------------------------:|-------------------------------------:|:-----------------------|:-----------------------------------------------------|:---------------------------------|
|       0.797429 |                          0.318484 |                     29.9316 |                      1.03044 |                 1.80461 |             0.774168 |              30.5008 |                                 4.53495 |                                              10 |              0.920006 |                           10.3039 |                        30.1017 |                         28.7428 |                      29.5683 |                     1.01667 |                1.77867 |                                          2.41608 |                                          0.014853 |                           0.00685038 | True                   | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH | True                             |

The pointwise powers include the globally profiled scale. The local-derivative powers remove that scale and show whether a constant power is genuinely implied by the two background curves.

|     z |   pointwise_b_perpendicular |   pointwise_b_parallel |   pointwise_delta_b |   local_derivative_b_perpendicular |   local_derivative_b_parallel |
|------:|----------------------------:|-----------------------:|--------------------:|-----------------------------------:|------------------------------:|
| 0.51  |                     1.02791 |                1.82979 |            0.801887 |                           0.994963 |                       1.72588 |
| 0.706 |                     1.02234 |                1.80815 |            0.785818 |                           1.01253  |                       1.74627 |
| 0.934 |                     1.02232 |                1.79893 |            0.776614 |                           1.03194  |                       1.77308 |
| 1.321 |                     1.02733 |                1.79706 |            0.769738 |                           1.05828  |                       1.80529 |
| 1.484 |                     1.02996 |                1.79799 |            0.768033 |                           1.06688  |                       1.81338 |
| 2.33  |                     1.04269 |                1.80383 |            0.761136 |                           1.09526  |                       1.82484 |

## DESI compressed validation

### Frozen PLAMB

| measurement_id                                | classification                                 |     z | model                |   observed_DM_over_rd |   observed_DH_over_rd |   predicted_DM_over_rd |   predicted_DH_over_rd |   residual_DM |   residual_DH |   marginal_pull_DM |   marginal_pull_DH |   principal_mode_pull_1 |   principal_mode_pull_2 |   chi2_2d |       p_2d |   systematic_fraction_added |   covariance_minimum_eigenvalue |
|:----------------------------------------------|:-----------------------------------------------|------:|:---------------------|----------------------:|----------------------:|-----------------------:|-----------------------:|--------------:|--------------:|-------------------:|-------------------:|------------------------:|------------------------:|----------:|-----------:|----------------------------:|--------------------------------:|
| DESI_DR2_combined_stat_sys                    | same_spectra_primary_and_internal_splits       | 2.33  | PLAMB_frozen_lower_z |               38.989  |               8.63155 |                40.0853 |                8.80755 |      1.09635  |      0.176    |          2.06203   |           1.7415   |                -2.9442  |              -2.01992   |  12.7484  | 0.00170496 |                       0     |                      0.00825923 |
| DESI_DR2_auto_plus_0p3pct_theory_budget       | same_spectra_primary_and_internal_splits       | 2.33  | PLAMB_frozen_lower_z |               39.3018 |               8.61888 |                40.0853 |                8.80755 |      0.783542 |      0.188665 |          0.983784  |           1.49842  |                -2.12871 |              -0.963002  |   5.45876 | 0.0652596  |                       0.003 |                      0.012918   |
| DESI_DR2_cross_plus_0p3pct_theory_budget      | same_spectra_primary_and_internal_splits       | 2.33  | PLAMB_frozen_lower_z |               38.5886 |               8.65507 |                40.0853 |                8.80755 |      1.49677  |      0.152473 |          2.31943   |           1.17916  |                -2.53921 |              -2.27759   |  11.635   | 0.00297504 |                       0.003 |                      0.0131775  |
| DESI_DR2_multipoles_stat_only                 | same_DR2_spectra_alternative_method_control    | 2.348 | PLAMB_frozen_lower_z |               39.46   |               8.51    |                40.2655 |                8.74386 |      0.80554  |      0.233865 |          1.23929   |           1.94887  |                -3.16892 |              -1.18917   |  11.4562  | 0.00325324 |                       0     |                      0.00993906 |
| DESI_DR2_multipoles_plus_0p3pct_theory_budget | same_DR2_spectra_alternative_method_control    | 2.348 | PLAMB_frozen_lower_z |               39.46   |               8.51    |                40.2655 |                8.74386 |      0.80554  |      0.233865 |          1.21924   |           1.90621  |                -3.02599 |              -1.17182   |  10.5298  | 0.00516994 |                       0.003 |                      0.0107294  |
| DESI_DR1_combined                             | earlier_DESI_partly_overlapping_survey_control | 2.33  | PLAMB_frozen_lower_z |               39.71   |               8.52    |                40.0853 |                8.80755 |      0.375321 |      0.287546 |          0.395075  |           1.69145  |                -2.14947 |              -0.365624  |   4.75388 | 0.0928341  |                       0     |                      0.0220745  |
| DESI_DR1_auto                                 | earlier_DESI_partly_overlapping_survey_control | 2.33  | PLAMB_frozen_lower_z |               39.9844 |               8.55446 |                40.0853 |                8.80755 |      0.100918 |      0.253087 |          0.0705319 |           0.963221 |                -1.12223 |              -0.0547781 |   1.2624  | 0.531954   |                       0     |                      0.0540333  |
| DESI_DR1_cross                                | earlier_DESI_partly_overlapping_survey_control | 2.33  | PLAMB_frozen_lower_z |               39.3964 |               8.51139 |                40.0853 |                8.80755 |      0.688924 |      0.29616  |          0.585814  |           1.40319  |                -1.96659 |              -0.558     |   4.17885 | 0.123758   |                       0     |                      0.0331369  |

### Frozen Lambda-CDM control

| measurement_id                                | classification                                 |     z | model                    |   observed_DM_over_rd |   observed_DH_over_rd |   predicted_DM_over_rd |   predicted_DH_over_rd |   residual_DM |   residual_DH |   marginal_pull_DM |   marginal_pull_DH |   principal_mode_pull_1 |   principal_mode_pull_2 |   chi2_2d |     p_2d |   systematic_fraction_added |   covariance_minimum_eigenvalue |
|:----------------------------------------------|:-----------------------------------------------|------:|:-------------------------|----------------------:|----------------------:|-----------------------:|-----------------------:|--------------:|--------------:|-------------------:|-------------------:|------------------------:|------------------------:|----------:|---------:|----------------------------:|--------------------------------:|
| DESI_DR2_combined_stat_sys                    | same_spectra_primary_and_internal_splits       | 2.33  | LambdaCDM_frozen_lower_z |               38.989  |               8.63155 |                38.6567 |                8.48569 |    -0.332287  |    -0.145854  |          -0.624973 |          -1.4432   |                1.90668  |                0.597629 |  3.99259  | 0.135838 |                       0     |                      0.00825923 |
| DESI_DR2_auto_plus_0p3pct_theory_budget       | same_spectra_primary_and_internal_splits       | 2.33  | LambdaCDM_frozen_lower_z |               39.3018 |               8.61888 |                38.6567 |                8.48569 |    -0.645092  |    -0.133189  |          -0.809952 |          -1.05782  |                1.55824  |                0.794744 |  3.05972  | 0.216566 |                       0.003 |                      0.012918   |
| DESI_DR2_cross_plus_0p3pct_theory_budget      | same_spectra_primary_and_internal_splits       | 2.33  | LambdaCDM_frozen_lower_z |               38.5886 |               8.65507 |                38.6567 |                8.48569 |     0.0681317 |    -0.169381  |           0.105579 |          -1.30992  |                1.41371  |               -0.129064 |  2.01522  | 0.365091 |                       0.003 |                      0.0131775  |
| DESI_DR2_multipoles_stat_only                 | same_DR2_spectra_alternative_method_control    | 2.348 | LambdaCDM_frozen_lower_z |               39.46   |               8.51    |                38.8088 |                8.42105 |    -0.651153  |    -0.0889483 |          -1.00177  |          -0.741236 |                1.56295  |                0.977104 |  3.39755  | 0.182908 |                       0     |                      0.00993906 |
| DESI_DR2_multipoles_plus_0p3pct_theory_budget | same_DR2_spectra_alternative_method_control    | 2.348 | LambdaCDM_frozen_lower_z |               39.46   |               8.51    |                38.8088 |                8.42105 |    -0.651153  |    -0.0889483 |          -0.985563 |          -0.725009 |                1.48458  |                0.962349 |  3.13009  | 0.209078 |                       0.003 |                      0.0107294  |
| DESI_DR1_combined                             | earlier_DESI_partly_overlapping_survey_control | 2.33  | LambdaCDM_frozen_lower_z |               39.71   |               8.52    |                38.6567 |                8.48569 |    -1.05331   |    -0.0343081 |          -1.10875  |          -0.201813 |                0.851833 |                1.09717  |  1.9294   | 0.381098 |                       0     |                      0.0220745  |
| DESI_DR1_auto                                 | earlier_DESI_partly_overlapping_survey_control | 2.33  | LambdaCDM_frozen_lower_z |               39.9844 |               8.55446 |                38.6567 |                8.48569 |    -1.32772   |    -0.0687672 |          -0.927944 |          -0.261721 |                0.788447 |                0.916962 |  1.46247  | 0.481314 |                       0     |                      0.0540333  |
| DESI_DR1_cross                                | earlier_DESI_partly_overlapping_survey_control | 2.33  | LambdaCDM_frozen_lower_z |               39.3964 |               8.51139 |                38.6567 |                8.48569 |    -0.73971   |    -0.0256934 |          -0.628999 |          -0.121734 |                0.512582 |                0.621797 |  0.649372 | 0.722754 |                       0     |                      0.0331369  |

The 0.3% diagonal theoretical BAO-shift budget is added to the DR2 auto/cross and multipole controls where stated. The published DR2 combined covariance already includes it.

## Native eBOSS DR16 likelihood

| measurement_id             | classification                            |     z | model                    |   predicted_DM_over_rd |   predicted_DH_over_rd |   native_grid_best_DM |   native_grid_best_DH | inside_grid   |   delta_chi2_native |   wilks_p_2d |   HPD_tail_mass |
|:---------------------------|:------------------------------------------|------:|:-------------------------|-----------------------:|-----------------------:|----------------------:|----------------------:|:--------------|--------------------:|-------------:|----------------:|
| eBOSS_DR16_native_auto     | external_survey_native_likelihood_control | 2.334 | PLAMB_frozen_lower_z     |                40.1255 |                8.79333 |               37.7634 |               8.91706 | True          |             2.24671 |    0.325186  |       0.333471  |
| eBOSS_DR16_native_auto     | external_survey_native_likelihood_control | 2.334 | LambdaCDM_frozen_lower_z |                38.6906 |                8.47126 |               37.7634 |               8.91706 | True          |             3.14598 |    0.207424  |       0.216967  |
| eBOSS_DR16_native_cross    | external_survey_native_likelihood_control | 2.334 | PLAMB_frozen_lower_z     |                40.1255 |                8.79333 |               37.4433 |               9.05748 | True          |             3.07426 |    0.214998  |       0.228145  |
| eBOSS_DR16_native_cross    | external_survey_native_likelihood_control | 2.334 | LambdaCDM_frozen_lower_z |                38.6906 |                8.47126 |               37.4433 |               9.05748 | True          |             4.86843 |    0.0876667 |       0.097037  |
| eBOSS_DR16_native_combined | external_survey_native_likelihood_control | 2.334 | PLAMB_frozen_lower_z     |                40.1255 |                8.79333 |               37.4433 |               8.98727 | True          |             5.16619 |    0.07554   |       0.079816  |
| eBOSS_DR16_native_combined | external_survey_native_likelihood_control | 2.334 | LambdaCDM_frozen_lower_z |                38.6906 |                8.47126 |               37.4433 |               8.98727 | True          |             7.85962 |    0.0196474 |       0.0209077 |

The eBOSS release states that its native auto- and cross-correlation likelihood grids may be treated as independent. The combined row uses their product. `HPD_tail_mass` is the integrated grid likelihood at density no greater than the frozen prediction; `Wilks_p_2d` is supplied only as an asymptotic likelihood-ratio readout.

## Product availability

| product                                            | publicly_available   | used   | role                                               |
|:---------------------------------------------------|:---------------------|:-------|:---------------------------------------------------|
| DESI DR2 baseline native 9306-vector likelihood    | False                | False  | pending native rerun                               |
| DESI DR2 combined/auto/cross posterior compression | True                 | True   | primary compressed validation                      |
| DESI DR2 alternative multipole likelihood/data     | False                | False  | published summary method control                   |
| DESI DR1 correlation data and covariance           | True                 | False  | published summary; native Vega refit not performed |
| eBOSS DR16 auto/cross likelihood grids             | True                 | True   | external native non-Gaussian validation            |

## Gates

- `analytic_practical_LambdaCDM_equivalence`: **False**
- `analytic_geometry_goodness_p`: **0.9200056212958674**
- `analytic_maximum_fractional_geometry_difference`: **0.014853006574198968**
- `analytic_equivalent_minus_data_parameter_distance_sigma`: **2.4160766337026542**
- `PLAMB_DESI_DR2_combined_p`: **0.001704959545133816**
- `delta_chi2_PLAMB_minus_LambdaCDM_DESI_DR2_combined`: **8.755836179746254**
- `PLAMB_DESI_DR2_auto_p`: **0.06525959122283458**
- `PLAMB_DESI_DR2_cross_p`: **0.002975038046976611**
- `worst_DR2_split`: **DESI_DR2_cross_plus_0p3pct_theory_budget**
- `worst_DR2_split_p`: **0.002975038046976611**
- `PLAMB_DR2_multipole_p`: **0.00516993850520226**
- `delta_chi2_PLAMB_minus_LambdaCDM_DR2_multipoles`: **7.399696991941371**
- `PLAMB_DESI_DR1_combined_p`: **0.09283413059331302**
- `delta_chi2_PLAMB_minus_LambdaCDM_DESI_DR1`: **2.8244831859505783**
- `PLAMB_eBOSS_combined_HPD_tail`: **0.07981595524891492**
- `PLAMB_eBOSS_combined_Wilks_p`: **0.07554002558991839**
- `delta_native_chi2_PLAMB_minus_LambdaCDM_eBOSS`: **-2.693436176519909**
- `independent_DR1_and_eBOSS_controls_pass_0p01`: **True**
- `DR2_failure_survives_same_spectra_multipole_method`: **True**
- `DR2_cross_fails_while_auto_passes_0p01`: **True**
- `DR2_failure_replicated_by_independent_survey_at_0p01`: **False**
- `all_available_compressed_or_native_controls_pass_0p01`: **False**
- `native_DESI_DR2_likelihood_publicly_available`: **False**
- `native_DESI_DR2_likelihood_rerun_complete`: **False**
- `radial_transverse_adapter_validated`: **False**
- `physical_degree_identified`: **False**
- `further_physical_sampling_authorised`: **False**

## Validation checks

| check                                      |       value |    target |   tolerance | passed   |
|:-------------------------------------------|------------:|----------:|------------:|:---------|
| SN p regression                            |  0.797429   |  0.797429 |      5e-05  | True     |
| lower-z b_perpendicular regression         |  0.992066   |  0.992066 |      5e-05  | True     |
| lower-z b_parallel regression              |  1.75242    |  1.75242  |      5e-05  | True     |
| frozen z2.33 DM regression                 | 40.0853     | 40.0853   |      0.0005 | True     |
| frozen z2.33 DH regression                 |  8.80755    |  8.80755  |      0.0005 | True     |
| all Gaussian covariances positive definite |  0.00825923 |  0        |      0      | True     |
| native eBOSS predictions inside both grids |  1          |  1        |      0      | True     |
| equivalence optimisation succeeded         |  1          |  1        |      0      | True     |

## Interpretation

A successful constant-power reconstruction of Lambda-CDM does not identify a new symmetry. It shows that the two fitted powers can encode an already familiar distance geometry. The DESI DR2 split then tests whether the remaining high-redshift residual is common to the forest auto-correlation and quasar cross-correlation or concentrated in one tracer construction.

The alternative multipole analysis is a valuable method and covariance control, but it uses the same DR2 spectra. DESI DR1 is an earlier, partly overlapping survey control; eBOSS DR16 is the external spectroscopic control. These distinctions are preserved in every output table.

The auto/cross readout is a localisation diagnostic, not a formal test of the difference between the two split estimates: their shared spectra and cross-covariance prevent treating those p-values as independent. It nevertheless identifies the Ly-alpha x quasar construction as the immediate place to inspect redshift errors, proximity modelling and quasar-related nuisance priors.

No SU(2), anisotropic metric or evolving-ruler claim is authorised by this audit. A physical claim would additionally require the native DR2 likelihood, independent high-redshift replication, and an action-level source for the radial/transverse split.

## Primary sources

- DESI DR2 baseline: https://arxiv.org/abs/2503.14739
- DESI DR2 figure data: https://zenodo.org/records/15690869
- DESI DR2 multipoles: https://arxiv.org/abs/2603.04281
- DESI DR1: https://arxiv.org/abs/2404.03001
- eBOSS DR16: https://arxiv.org/abs/2007.08995
- SDSS DR16 likelihood release: https://svn.sdss.org/public/data/eboss/DR16cosmo/tags/v1_0_1/likelihoods/BAO-only/
