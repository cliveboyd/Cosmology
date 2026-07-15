# SU2 / Quaia External Dust Gate Readout

Date: July 15, 2026

## Bottom Line

SFD dust does not absorb the locked angular mode in this regression gate.

SFD plus the local Quaia selection/depth proxies and bundled Gaia scan-law counts are the best external-template gate available from the current local files. True WISE-depth and stellar-density maps still need acquisition before this becomes the full external-systematics gate.

## Key Numbers

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
