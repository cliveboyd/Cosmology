# SPARC Posterior Diagnostic Audit Readout

Date: 2026-07-16

RAR rank-normalised split-Rhat gate: does not pass; max `1.0616`, minimum bulk ESS `76.269`. PLAMB max rank-split Rhat is `1.3082` and its minimum bulk ESS is `8.1147`.

| score_case | model | delta_vs_RAR | galaxies_better_than_RAR | galaxies_worse_than_RAR |
| --- | --- | --- | --- | --- |
| paired_posterior_fit_baseline_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | 13.6018 | 78 | 79 |
| paired_posterior_fit_baseline_all_Q2 | RAR | 0 | 0 | 0 |
| paired_posterior_fit_bulge_removed_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | 25.0818 | 69 | 69 |
| paired_posterior_fit_bulge_removed_all_Q2 | RAR | 0 | 0 | 0 |
| paired_posterior_fit_gas_dominated_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | 4.82809 | 20 | 11 |
| paired_posterior_fit_gas_dominated_all_Q2 | RAR | 0 | 0 | 0 |
| paired_posterior_fit_high_inclination_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | 11.9272 | 43 | 40 |
| paired_posterior_fit_high_inclination_all_Q2 | RAR | 0 | 0 | 0 |
| paired_posterior_fit_low_accel_outer_points_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | 0.88998 | 79 | 76 |
| paired_posterior_fit_low_accel_outer_points_all_Q2 | RAR | 0 | 0 | 0 |
| paired_posterior_fit_low_inclination_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | 1.64359 | 36 | 38 |
| paired_posterior_fit_low_inclination_all_Q2 | RAR | 0 | 0 | 0 |

Claim boundary: subset wins, not a full-sample win. Both convergence and out-of-sample prediction remain open gates; the paired-posterior table is an in-sample fit check only.
