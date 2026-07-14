# SPARC Posterior Run Analysis

Run directory: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\sparc_hierarchical_posterior\posterior_minus_persistent_negative_full_20260714`

## Chain Diagnostics

| Split | Model | Chains | Saved / chain | Max Rhat | Global accept | Galaxy accept |
|---|---|---:|---:|---:|---:|---:|
| all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | 3 | 1500 | 1.054 | 0.074 | 0.680 |
| all_Q2 | RAR | 3 | 1500 | 1.3 | 0.133 | 0.681 |

## PLAMB Checks

| Split | Model | p median [16,84] | p=0.5 in 68%? | kappa median [16,84] | kappa=1 in 68%? |
|---|---|---:|---:|---:|---:|
| all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | 0.5193 [0.51437, 0.52374] | False | 0.78976 [0.76858, 0.81729] | False |

## Nuisance Pulls

| Split | Model | Distance pull >3 | M/L pull >3 | Max |d pull| | Max |M/L pull| |
|---|---|---:|---:|---:|---:|
| all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | 13 | 14 | 3.961 | 3.867 |
| all_Q2 | RAR | 14 | 14 | 3.961 | 3.849 |

## Chain-Based Predictive Scores

| Case | Model | Delta log score vs RAR | Better / worse galaxies |
|---|---|---:|---:|
| posterior_baseline_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | -368.314 | 80 / 77 |
| posterior_stress_bulge_removed_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | 99.583 | 83 / 55 |
| posterior_stress_gas_dominated_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | 93.318 | 21 / 10 |
| posterior_stress_high_inclination_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | -612.830 | 45 / 38 |
| posterior_stress_high_quality_Q1_all | PLAMB_OPTICAL_DEPTH_KAPPA_P | 274.392 | 51 / 44 |
| posterior_stress_low_accel_outer_points_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | -173.603 | 81 / 74 |
| posterior_stress_low_inclination_all_Q2 | PLAMB_OPTICAL_DEPTH_KAPPA_P | -294.276 | 28 / 46 |

## Output Files

- case diagnostics: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\sparc_hierarchical_posterior\posterior_minus_persistent_negative_full_20260714\analysis_full_nd16_20260714\posterior_case_diagnostics.csv`
- global parameter summary: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\sparc_hierarchical_posterior\posterior_minus_persistent_negative_full_20260714\analysis_full_nd16_20260714\posterior_global_parameter_summary.csv`
- PLAMB checks: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\sparc_hierarchical_posterior\posterior_minus_persistent_negative_full_20260714\analysis_full_nd16_20260714\posterior_plamb_checks.csv`
- nuisance summary: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\sparc_hierarchical_posterior\posterior_minus_persistent_negative_full_20260714\analysis_full_nd16_20260714\posterior_nuisance_pull_summary.csv`
- predictive summary: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\sparc_hierarchical_posterior\posterior_minus_persistent_negative_full_20260714\analysis_full_nd16_20260714\posterior_predictive_summary.csv`
