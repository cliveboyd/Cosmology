# SU2 / Quaia WISE-Depth and Stellar-Density Gate Readout

Date: July 15, 2026

## Bottom Line

AllWISE W1/W2 tile-depth controls do not absorb the locked angular mode.

Gaia DR3 stellar-density control does not absorb the locked angular mode.

The full external stack including WISE depth and stellar density retains a substantial residual mode.

## Key Numbers

| control_group | amp_ratio_vs_baseline | snr_max | snr_delta_vs_baseline | max_pair_direction_sep_deg | mean_resultant_length | N_min | N_max |
| --- | --- | --- | --- | --- | --- | --- | --- |
| baseline_no_template | 1 | 4.7852 | 0 | 13.4723 | 0.996066 | 118193 | 183527 |
| wise_depth_external | 0.978169 | 4.75555 | -0.0296507 | 15.7795 | 0.994619 | 118193 | 183527 |
| gaia_stellar_density_external | 1.09983 | 5.00263 | 0.217428 | 16.0578 | 0.994655 | 118193 | 183527 |
| wise_plus_stellar_external | 1.06899 | 4.92309 | 0.137892 | 17.965 | 0.993594 | 118193 | 183527 |
| sfd_selection_scanlaw_external | 0.866391 | 3.60328 | -1.18192 | 31.6722 | 0.979236 | 118193 | 183527 |
| all_external_wise_stellar | 0.920458 | 3.55534 | -1.22986 | 35.1206 | 0.975352 | 118193 | 183527 |
| local_catalog_quality_proxy | 0.819012 | 6.22483 | 1.43963 | 58.3541 | 0.917343 | 118193 | 183527 |
| all_external_plus_local_quality | 0.711541 | 3.66543 | -1.11977 | 36.7889 | 0.972961 | 118193 | 183527 |
