# SU2 / Quaia External Dust Template Gate

Date: July 15, 2026

## Locked Test

- Redshift window: 0.95 <= z < 1.45
- Latitude cuts: |b|>15, |b|>20, |b|>25, |b|>30, |b|>35
- Response: Quaia redshift fitted with a Galactic dipole plus optional standardized template controls.

## External / Template Inputs

- SFD dust: Schlegel, Finkbeiner & Davis 1998 full-sky E(B-V), Harvard Dataverse DOI 10.7910/DVN/EWCNL5.
- Selection depth proxy: Local Quaia G<20.0 selection_function_NSIDE64 map.
- Random density proxy: Local Quaia G<20.0 10x random catalog binned to the selection HEALPix grid.
- Gaia scan law: Nominal Gaia scanning-law FITS bundle from the gaiascanlaw package, counted to the DR3 end date.
- WISE-depth and stellar-density maps were not present as true external all-sky maps in the local tree; WISE magnitudes/colors and Gaia catalog-quality columns are reported as local proxies only.

## Control Summary

| control_group | amp_ratio_vs_baseline | snr_max | snr_delta_vs_baseline | max_pair_direction_sep_deg | mean_resultant_length | N_min | N_max |
| --- | --- | --- | --- | --- | --- | --- | --- |
| baseline_no_template | 1 | 4.7852 | 0 | 13.4723 | 0.996066 | 118193 | 183527 |
| sfd_dust_external | 1.00644 | 4.66601 | -0.119192 | 12.9543 | 0.996222 | 118193 | 183527 |
| quaia_selection_depth | 1.03478 | 4.79898 | 0.0137761 | 13.5025 | 0.995953 | 118193 | 183527 |
| sfd_plus_selection_depth | 1.02617 | 4.71411 | -0.0710932 | 12.9733 | 0.996255 | 118193 | 183527 |
| gaia_scanlaw_external | 0.783971 | 3.36404 | -1.42116 | 29.8825 | 0.980638 | 118193 | 183527 |
| sfd_selection_gaia_external | 0.866391 | 3.60328 | -1.18192 | 31.6722 | 0.979236 | 118193 | 183527 |
| local_catalog_quality_proxy | 0.819012 | 6.22483 | 1.43963 | 58.3541 | 0.917343 | 118193 | 183527 |
| all_available_templates | 0.499563 | 2.92202 | -1.86318 | 48.2136 | 0.948272 | 118193 | 183527 |

## Template Correlations In Primary Mask

| template | N | mean | std | corr_z |
| --- | --- | --- | --- | --- |
| ebv_sfd | 183527 | 0.0547752 | 0.0537952 | 0.00119195 |
| selection_T | 183527 | 0.641014 | 0.3726 | -0.00948667 |
| random_density_log1p | 183527 | 5.30697 | 0.228901 | 0.000251233 |
| gaia_scan_count_log1p_dr3 | 183527 | 3.83054 | 0.372824 | 0.00742959 |
| gaia_scan_angle_resultant_dr3 | 183527 | 0.300056 | 0.175469 | 0.000134338 |
| zerr | 183527 | 0.139388 | 0.125683 | 0.040653 |
| g | 183527 | 19.3213 | 0.58165 | 0.0216554 |
| w1 | 183527 | 15.5204 | 0.570684 | 0.308456 |
| w2 | 183527 | 14.3271 | 0.61946 | 0.198129 |
| w1_w2 | 183527 | 1.19329 | 0.153055 | 0.348229 |

## Readout

- SFD dust does not absorb the locked angular mode in this regression gate.
- Selection/depth-only amp ratio: 1.035; max direction spread: 13.5 deg.
- Gaia scan-law-only amp ratio: 0.784; max direction spread: 29.88 deg.
- SFD plus local selection/depth proxies still leaves the locked angular mode stable.
- SFD plus selection/depth plus Gaia scan-law controls materially suppresses or rotates the mode.
- All-available proxy amp ratio: 0.4996; max direction spread: 48.21 deg.
- The all-available proxy stack remains a stress test rather than a headline correction, because it mixes true external dust, Quaia selection products, and local catalog-quality columns.

## Next Gate

Acquire true external WISE-depth and stellar-density maps on a documented HEALPix grid, then rerun this exact template gate with those external maps replacing local proxy columns.
