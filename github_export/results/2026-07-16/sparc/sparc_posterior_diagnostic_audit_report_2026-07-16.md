# SPARC Posterior Diagnostic Audit

Date: 2026-07-16

## Bottom Line

The rescued RAR chain does not pass the rank-normalised split-Rhat <= 1.05 gate. PLAMB has maximum rank-split Rhat 1.308 and is not robustly converged. The older stress table was not a true paired posterior-predictive calculation because galaxy nuisances were redrawn from their priors. The replacement table below uses stored paired global and nuisance posterior draws, but is deliberately labelled an in-sample posterior fit-density check rather than out-of-sample predictive evidence.

## Robust Convergence Diagnostics

| model | max_classic_rhat | max_rank_split_rhat | min_bulk_ess |
| --- | --- | --- | --- |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | 1.05424 | 1.30819 | 8.11474 |
| RAR | 1.03071 | 1.06163 | 76.2689 |

| model | parameter | classic_rhat | rank_split_rhat | bulk_ess | mcse_mean_approx |
| --- | --- | --- | --- | --- | --- |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | log_ydisk | 1.05186 | 1.30819 | 8.11474 | 0.00825179 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | log_ybul | 1.05424 | 1.28716 | 8.20626 | 0.00840389 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | log10_kappa | 1.01415 | 1.05227 | 80.559 | 0.00150957 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | bridge_exponent | 1.00035 | 1.03386 | 98.76 | 0.000492067 |
| RAR | log_ydisk | 1.0245 | 1.05688 | 79.7275 | 0.00314553 |
| RAR | log_ybul | 1.03071 | 1.06163 | 76.2689 | 0.0035655 |
| RAR | log10_gdagger | 1.00123 | 1.00398 | 898.87 | 0.000259993 |

## Paired Posterior Fit-Density Checks

| score_case | model | N_galaxies | N_points | delta_vs_RAR | median_galaxy_delta_vs_RAR | galaxies_better_than_RAR | galaxies_worse_than_RAR |
| --- | --- | --- | --- | --- | --- | --- | --- |
| paired_posterior_fit_baseline_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | 157 | 3001 | 13.6018 | -0.00202503 | 78 | 79 |
| paired_posterior_fit_baseline_all_Q2 | RAR | 157 | 3001 | 0 | 0 | 0 | 0 |
| paired_posterior_fit_bulge_removed_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | 138 | 2279 | 25.0818 | 0.000189848 | 69 | 69 |
| paired_posterior_fit_bulge_removed_all_Q2 | RAR | 138 | 2279 | 0 | 0 | 0 | 0 |
| paired_posterior_fit_gas_dominated_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | 31 | 449 | 4.82809 | 0.0114855 | 20 | 11 |
| paired_posterior_fit_gas_dominated_all_Q2 | RAR | 31 | 449 | 0 | 0 | 0 | 0 |
| paired_posterior_fit_high_inclination_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | 83 | 1739 | 11.9272 | 0.0136206 | 43 | 40 |
| paired_posterior_fit_high_inclination_all_Q2 | RAR | 83 | 1739 | 0 | 0 | 0 | 0 |
| paired_posterior_fit_low_accel_outer_points_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | 155 | 1401 | 0.88998 | 0.00348601 | 79 | 76 |
| paired_posterior_fit_low_accel_outer_points_all_Q2 | RAR | 155 | 1401 | 0 | 0 | 0 | 0 |
| paired_posterior_fit_low_inclination_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | 74 | 1262 | 1.64359 | -0.0069264 | 36 | 38 |
| paired_posterior_fit_low_inclination_all_Q2 | RAR | 74 | 1262 | 0 | 0 | 0 | 0 |

## Interpretation Boundary

These scores reuse the fitted observations and therefore measure posterior fit, not generalisation. They can identify severe subset tension, but they cannot promote PLAMB over RAR. Their current positive values are provisional because the PLAMB global chain also fails robust convergence. A publication-grade comparison requires converged chains followed by galaxy-level cross-validation or a genuinely held-out SPARC subset with all selection and nuisance priors frozen in advance.

The standing SPARC conclusion remains: subset wins, not a full-sample win. Convergence rescue removes one gate failure; any negative high-inclination or gas-dominated paired score remains an explicit gating failure.

## Configuration

```json
{
  "date": "2026-07-16",
  "run_dir": "C:\\Users\\clive\\Documents\\Cosmology\\plamb_runs\\diagnostics\\sparc_hierarchical_posterior\\posterior_minus_persistent_negative_rar_rescue_with_existing_plamb_20260714",
  "map_dir": "C:\\Users\\clive\\Documents\\Cosmology\\plamb_runs\\diagnostics\\sparc_hierarchical_map\\optical_depth_hierarchical_20260714",
  "max_paired_draws": 4500,
  "seed": 270716,
  "convergence": "rank-normalised split Rhat plus approximate rank-normalised bulk ESS",
  "fit_score": "log mean likelihood over paired stored global/nuisance posterior draws; in-sample diagnostic only"
}
```
