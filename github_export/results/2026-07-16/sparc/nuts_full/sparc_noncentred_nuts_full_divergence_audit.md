# SPARC Full Non-centred NUTS Divergence Audit

This post-run audit preserves the preregistered zero-divergence gate. It locates every divergent transition; it does not redefine convergence after seeing the outcome.

## Counts

| Model | Retained draws | Divergences | Rate | By chain | All draws near boundary |
|---|---:|---:|---:|---|---:|
| RAR | 8,000 | 3 | 0.000375 | 1;0;2;0 | 1.000 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | 8,000 | 6 | 0.000750 | 1;1;2;2 | 1.000 |

A transition is labelled near the retained four-standard-deviation nuisance boundary when either maximum latent absolute value is at least 3.8.

## Transition Locations

| Model | Chain | Draw | Depth | Divergence energy error | Max abs(z_D) galaxy | Max abs(z_eta) galaxy | Latents >= 3.8 |
|---|---:|---:|---:|---:|---|---|---:|
| RAR | 1 | 1701 | 6 | 2133.786 | NGC5371 (-3.970) | NGC4217 (-3.972) | 11 |
| RAR | 3 | 997 | 3 | 1292.145 | NGC5371 (-3.972) | NGC4217 (-3.966) | 5 |
| RAR | 3 | 1892 | 5 | 1051.416 | NGC0801 (-3.927) | IC2574 (-3.993) | 10 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | 1 | 620 | 3 | 2555.044 | NGC2915 (3.991) | UGC00128 (3.991) | 8 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | 2 | 1736 | 6 | 1871.874 | NGC0801 (-3.950) | F571-8 (-3.983) | 6 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | 3 | 416 | 6 | 3915.474 | NGC2403 (3.998) | NGC2915 (3.983) | 8 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | 3 | 637 | 4 | 1050.840 | NGC2403 (3.990) | UGC00128 (3.985) | 10 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | 4 | 487 | 4 | 272385.851 | NGC5371 (-3.996) | NGC4217 (-3.948) | 8 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | 4 | 524 | 1 | 1012.914 | NGC7814 (3.992) | IC2574 (-3.997) | 11 |

## Readout

All 9 divergent transitions lie near the explicit nuisance boundary, but so does every retained draw in both models. The mean numbers of latent variables with absolute value at least 3.8 are RAR: divergent 8.67, non-divergent 8.51; PLAMB_OPTICAL_DEPTH_KAPPA_P: divergent 8.50, non-divergent 8.76. There is therefore no transition-level enrichment of boundary occupancy among the nine divergences.
Repeated nearest-boundary galaxies across divergent draws: PLAMB_OPTICAL_DEPTH_KAPPA_P: NGC2403 (2), PLAMB_OPTICAL_DEPTH_KAPPA_P: NGC2915 (2), PLAMB_OPTICAL_DEPTH_KAPPA_P: UGC00128 (2), RAR: NGC4217 (2), RAR: NGC5371 (2).

The truncation is nevertheless active throughout the posterior. The ten highest marginal boundary occupancies are:

| Model | Latent | Galaxy | Fraction with abs(z) >= 3.8 | Mean z |
|---|---|---|---:|---:|
| RAR | z_logd | NGC5371 | 0.981 | -3.947 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | z_logd | NGC5371 | 0.980 | -3.946 |
| RAR | z_logd | NGC2403 | 0.882 | 3.906 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | z_logd | NGC2403 | 0.879 | 3.906 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | z_logd | NGC0801 | 0.844 | -3.892 |
| RAR | z_logd | NGC0801 | 0.827 | -3.888 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | z_logeta | UGC00128 | 0.647 | 3.819 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | z_logeta | F571-8 | 0.626 | -3.803 |
| RAR | z_logeta | UGC00128 | 0.625 | 3.809 |
| RAR | z_logeta | F571-8 | 0.606 | -3.792 |

The divergence rates are small, and the rank-normalised R-hat, effective sample size, E-BFMI and tree-depth gates pass in the main report. Nevertheless, both models fail the separately preregistered zero-divergence requirement. The appropriate rescue is an unchanged-model rerun at higher target acceptance, registered before execution; these draws cannot be relabelled as a strict convergence pass.
