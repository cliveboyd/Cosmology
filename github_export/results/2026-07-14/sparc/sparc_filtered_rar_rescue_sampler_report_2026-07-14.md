# SPARC Hierarchical Posterior Sampler: all_Q2 / RAR

## Setup

- MAP source objective: `6803.058721507953`
- chains: `5`
- sweeps per chain: `96000`
- burn-in sweeps: `36000`
- thinning: `20`
- galaxy update fraction per sweep: `1.0`

## Global Posterior Summary

| Parameter | median | q16 | q84 | Rhat |
|---|---:|---:|---:|---:|
| log_ydisk | -0.584398 | -0.611985 | -0.557115 | 1.0245 |
| log_ybul | -0.503659 | -0.535388 | -0.473469 | 1.0307 |
| log10_gdagger | -10.0499 | -10.0575 | -10.0422 | 1.0012 |
| ydisk | 0.557441 | 0.542274 | 0.57286 |  |
| ybul | 0.604315 | 0.585442 | 0.622838 |  |

## Chain Diagnostics

| Chain | Global accept | Galaxy accept | Saved draws | Final logpost |
|---|---:|---:|---:|---:|
| 1 | 0.226 | 0.681 | 3000 | -11451.1 |
| 2 | 0.230 | 0.681 | 3000 | -11436.4 |
| 3 | 0.227 | 0.681 | 3000 | -11450.3 |
| 4 | 0.228 | 0.681 | 3000 | -11451 |
| 5 | 0.230 | 0.681 | 3000 | -11457.3 |

## Outputs

- global samples: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\sparc_hierarchical_posterior\posterior_minus_persistent_negative_rar_rescue_20260714\all_Q2\RAR\global_chain_samples.csv`
- posterior summary: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\sparc_hierarchical_posterior\posterior_minus_persistent_negative_rar_rescue_20260714\all_Q2\RAR\posterior_global_summary.csv`
- galaxy nuisance summary: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\sparc_hierarchical_posterior\posterior_minus_persistent_negative_rar_rescue_20260714\all_Q2\RAR\posterior_galaxy_nuisance_summary.csv`
- compressed samples: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\sparc_hierarchical_posterior\posterior_minus_persistent_negative_rar_rescue_20260714\all_Q2\RAR\posterior_samples.npz`
