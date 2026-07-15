# SPARC Non-centred NUTS Pilot Report

**Completed:** `2026-07-15T20:39:28.104617+00:00`  
**Pilot galaxies:** `24`  
**Pilot points:** `398`  
**Chains:** `4`

## Decision

Both models pass every strict preregistered geometry gate. The implementation is ready for a separately registered full filtered all-Q2 run.

This pilot assesses sampler geometry only. It is not held-out prediction and must not be used as a PLAMB-versus-RAR evidence statement.

## Gate Results

| Model | Status | Max R-hat | Min bulk ESS | Min tail ESS | Divergences | Min E-BFMI |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| RAR | PASS | 1.00738 | 816.2 | 1193.4 | 0 | 0.8345 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | PASS | 1.00551 | 730.3 | 1194.7 | 0 | 0.8624 |

## Likelihood Equivalence

| Model | Max velocity difference (km/s) | Chi-squared difference | Pass |
| --- | ---: | ---: | --- |
| RAR | 5.684e-14 | 0.000e+00 | True |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | 8.527e-14 | 1.819e-12 | True |

## Global Posterior Readout

These values are subset-level sampler checks, not full-sample estimates.

| Model | Parameter | Mean | SD | 3% HDI | 97% HDI | R-hat | Bulk ESS | Tail ESS |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| RAR | ydisk | 0.520964 | 0.0391771 | 0.443389 | 0.589227 | 1.00397 | 816.2 | 1193.4 |
| RAR | ybul | 0.0875344 | 0.0171414 | 0.0566636 | 0.118526 | 1.00208 | 2903.0 | 1862.9 |
| RAR | log10_gdagger | -10.0222 | 0.0181796 | -10.0546 | -9.98726 | 1.00113 | 1916.4 | 2354.1 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | ydisk | 0.583691 | 0.0422354 | 0.510581 | 0.66726 | 1.00517 | 730.3 | 1213.3 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | ybul | 0.0914199 | 0.0177225 | 0.0584519 | 0.123989 | 1.00221 | 2081.1 | 1556.9 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | log10_kappa | -0.233418 | 0.0275131 | -0.282116 | -0.180224 | 1.00127 | 1554.2 | 2370.7 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | bridge_exponent | 0.589703 | 0.0133946 | 0.565289 | 0.614878 | 1.00042 | 2019.8 | 2659.4 |

## Reproducibility

- PyMC: `5.28.5`
- nutpie: `0.16.11`
- ArviZ: `0.23.4`
- Python: `3.13.9 | packaged by Anaconda, Inc. | (main, Oct 21 2025, 19:09:58) [MSC v.1929 64 bit (AMD64)]`
- selection seed: `160726`
- sampling seed: `16072601`
