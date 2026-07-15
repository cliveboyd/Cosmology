# SU2 / Quaia WISE-Depth and Stellar-Density Gate

Date: July 15, 2026

## Bottom Line

- AllWISE W1/W2 tile-depth controls do not absorb the locked angular mode.
- Gaia DR3 stellar-density control does not absorb the locked angular mode.
- The full external stack including WISE depth and stellar density retains a substantial residual mode.

## Control Summary

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

## Template Correlations In Primary Mask

| template | N | mean | std | corr_z |
| --- | --- | --- | --- | --- |
| wise_w1_medcov_log1p | 183527 | 3.54457 | 0.381855 | 0.00387396 |
| wise_w2_medcov_log1p | 183527 | 3.53523 | 0.393788 | 0.00351598 |
| wise_w1_mincov | 183527 | 0.539046 | 1.54214 | -0.000284963 |
| wise_w2_mincov | 183527 | 2.18897 | 5.48744 | -0.000520654 |
| wise_w1_lowcov_log1p | 183527 | 0.00499437 | 0.0054772 | -0.00328039 |
| wise_w2_lowcov_log1p | 183527 | 0.00519293 | 0.068638 | -0.00153156 |
| wise_w1_w2_medcov_diff | 183527 | 0.168489 | 1.47272 | 0.00174024 |
| gaia_dr3_density_hpx4_log1p | 183527 | 4.53015 | 0.79516 | -0.0015647 |
| ebv_sfd | 183527 | 0.0547752 | 0.0537952 | 0.00119195 |
| selection_T | 183527 | 0.641014 | 0.3726 | -0.00948667 |
| random_density_log1p | 183527 | 5.30697 | 0.228901 | 0.000251233 |
| gaia_scan_count_log1p_dr3 | 183527 | 3.83054 | 0.372824 | 0.00742959 |
| gaia_scan_angle_resultant_dr3 | 183527 | 0.300056 | 0.175469 | 0.000134338 |
| w1_w2 | 183527 | 1.19329 | 0.153055 | 0.348229 |

## External Inputs

- AllWISE Atlas Inventory W1/W2 tile-level depth-of-coverage metadata (`medcov`, `mincov`, `maxcov`, low-coverage fractions), matched to Quaia objects by nearest Atlas tile center.
- Gaia DR3 source-density map from `gaiadr3.gaia_source_lite`, `phot_g_mean_mag < 20.5`, HEALPix level 4, when the TAP aggregate is available; otherwise CDS Gaia EDR3 density HiPS tiles are sampled.
- Existing SFD, Quaia selection/random-density, and Gaia scan-law templates are reused from the prior gate.

## Configuration

```json
{
  "date": "2026-07-15",
  "z_min": 0.95,
  "z_max": 1.45,
  "b_cuts": [
    15.0,
    20.0,
    25.0,
    30.0,
    35.0
  ],
  "wise_inventory": "C:\\Users\\clive\\Documents\\Cosmology\\external_datasets\\allwise_atlas_inventory\\allwise_p3am_cdd_w1w2_coverage.csv",
  "gaia_density": "C:\\Users\\clive\\Documents\\Cosmology\\external_datasets\\gaia_dr3_density\\gaia_dr3_source_density_hpx4_glt20p5.csv",
  "gaia_density_present": true,
  "gaia_density_source": "CDS Gaia EDR3 density HiPS Norder2",
  "gaia_density_hips_base_url": "https://alasky.cds.unistra.fr/ancillary/GaiaEDR3/density-map",
  "gaia_density_hips_order": 2,
  "method": "AllWISE Atlas W1/W2 tile coverage matched by nearest tile center; Gaia DR3 density sampled from TAP HEALPix level 4 when available, otherwise from CDS Gaia EDR3 density HiPS tiles."
}
```
