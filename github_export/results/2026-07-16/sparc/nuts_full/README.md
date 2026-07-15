# Separately Preregistered 157-galaxy SPARC NUTS Run

This directory contains the compact reproducibility record for the full filtered all-Q2 run completed on 16 July 2026 (Australia/Sydney). The sample contains 157 galaxies and 3,001 rotation-curve points. RAR and `PLAMB_OPTICAL_DEPTH_KAPPA_P` each used four chains with 2,000 tuning and 2,000 retained draws per chain.

The exact July likelihood-equivalence checks pass. Rank-normalised R-hat, bulk and tail effective sample size, E-BFMI and tree-depth criteria also pass for both models. The locked zero-divergence gate does not pass: RAR has 3 divergences and PLAMB has 6 among 8,000 retained draws per model.

The post-run audit finds active nuisance-boundary saturation throughout both posteriors. Every retained draw has at least one galaxy-level latent with `abs(z) >= 3.8`; the number is not elevated at divergent transitions. This is a modelling sensitivity to preserve and test, not a reason to waive the registered convergence gate.

The large NetCDF traces remain in the local run directory and are intentionally excluded from Git. Their SHA-256 identifiers are:

| Model | Trace file | SHA-256 |
|---|---|---|
| RAR | `sparc_noncentred_nuts_full_RAR_trace.nc` | `5645d268b6bc51f488767341b6c8e07f9e40b3db495627a5ad1d99b5b7089d9d` |
| PLAMB | `sparc_noncentred_nuts_full_PLAMB_OPTICAL_DEPTH_KAPPA_P_trace.nc` | `b517661935b2b93e37e2b90df9d56fe7574c3d91f8e43e9f5efeb1636c0ae983` |

The standing scientific boundary remains: `subset wins, not a full-sample win`. These files establish strong mixing under the non-centred parameterisation, but they neither supply a strict zero-divergence pass nor replace galaxy-held-out predictive validation.
