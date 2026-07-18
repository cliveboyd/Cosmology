# SPARC FR Environment Result: Post Hoc Robustness Audit

Date: 2026-07-18

## Status

This audit is post hoc. It cannot alter the failed preregistered development gate and cannot authorise access to the reserved replication outcomes.

## Component Readout

| component | partial_spearman_rho | permutation_p_two_sided_unadjusted | holm_p_four_components | cv_rmse_fractional_improvement | coefficient_sign_fraction |
| --- | --- | --- | --- | --- | --- |
| log1p_count_2mpc | -0.298361 | 0.0153492 | 0.0477976 | 0.0255981 | 1 |
| log1p_count_5mpc | 0.0232335 | 0.849258 | 0.849258 | -0.0108551 | 0.8 |
| minus_log10_nearest_mpc | -0.246634 | 0.0451977 | 0.0903955 | 0.00422313 | 1 |
| log10_tidal_strength | -0.313725 | 0.0119494 | 0.0477976 | 0.0246004 | 1 |

The four component p-values are Holm-adjusted as one exploratory family. A component result is descriptive even if its adjusted p-value is small because the decomposition was examined after the composite result was known.

## Leave-One-Out Stability

The full partial Spearman correlation is `-0.307254`. Across all `66` single-galaxy removals, rho ranges from `-0.347159` to `-0.277491`, with the original sign retained in `100.0` per cent of removals.

Largest absolute changes:

| removed_galaxy | leave_one_out_rho | change_from_full_rho | same_sign_as_full |
| --- | --- | --- | --- |
| UGC00634 | -0.347159 | -0.0399052 | True |
| UGC00731 | -0.277491 | 0.0297627 | True |
| UGC11914 | -0.336145 | -0.0288912 | True |
| F583-1 | -0.27972 | 0.0275337 | True |
| ESO079-G014 | -0.28007 | 0.027184 | True |
| NGC5033 | -0.333916 | -0.0266621 | True |
| NGC5371 | -0.280813 | 0.026441 | True |
| UGC12732 | -0.333086 | -0.0258317 | True |
| DDO170 | -0.282561 | 0.0246928 | True |
| D512-2 | -0.28361 | 0.0236438 | True |

## Residual Quartiles

| environment_residual_quartile | N | median_environment_residual | mean_outcome_residual | median_outcome_residual |
| --- | --- | --- | --- | --- |
| Q1_low | 17 | -0.4806 | 0.00932051 | 0.00454562 |
| Q2 | 16 | -0.165437 | -0.00456611 | 0.000132527 |
| Q3 | 16 | 0.107651 | 7.59378e-05 | -0.00104339 |
| Q4_high | 17 | 0.45101 | -0.00809015 | -0.000985992 |

## Host-Covariate Correlations

| host_feature | spearman_rho_with_environment_score |
| --- | --- |
| abs_gal_b_deg | 0.561048 |
| log1p_sideband_count_5mpc | -0.469717 |
| frac_distance_error | -0.45091 |
| Inc_deg | 0.28706 |
| log_n_points | -0.188907 |
| log_L36 | 0.187308 |
| log_SBeff | 0.182406 |
| T | -0.149071 |
| log_Vflat | 0.1466 |
| D_Mpc | 0.125018 |
| bulge_proxy | -0.0884075 |
| log_MHI | 0.0151977 |

## Interpretation

The primary pattern should be carried forward only as a possible conventional environment term. The next defensible analysis would use an external group catalogue or a deeper volume-limited tracer, freeze a single dominant component suggested here, and test genuinely held-out galaxies. No antimatter interpretation follows from these diagnostics.
