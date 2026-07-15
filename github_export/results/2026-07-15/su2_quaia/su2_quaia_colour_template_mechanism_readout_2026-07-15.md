# SU2 / Quaia Colour-Template Mechanism Gate Readout

Date: July 15, 2026

## Bottom Line

- `gaia_bp_rp_colour` is best carried by `gaia_scanlaw`: score `0.497779`, R2 `0.00640632`, residual amp/raw `0.546831`, SMD absorption `0.835919`.
- `wise_w1_w2_colour` is best carried by `gaia_scanlaw`: score `0.410102`, R2 `0.000266778`, residual amp/raw `0.57196`, SMD absorption `0.62123`.
- `gaia_wise_colour_cross` is best carried by `all_external_templates`: score `0.188208`, R2 `0.00267529`, residual amp/raw `0.582949`, SMD absorption `-0.82724`.

## Top Mechanism Families

| rank_within_response | response_field | template_group | mechanism_score | template_r2_mean | resid_amp_ratio_vs_raw | dipole_absorption_fraction | smd_absorption_fraction | max_resid_direction_sep_from_raw_deg |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | gaia_bp_rp_colour | gaia_scanlaw | 0.497779 | 0.00640632 | 0.546831 | 0.453169 | 0.835919 | 16.776 |
| 1 | wise_w1_w2_colour | gaia_scanlaw | 0.410102 | 0.000266778 | 0.57196 | 0.42804 | 0.62123 | 22.828 |
| 1 | gaia_wise_colour_cross | all_external_templates | 0.188208 | 0.00267529 | 0.582949 | 0.417051 | -0.82724 | 75.4432 |

## Strongest Individual Template Correlations

| response_field | template | corr_response_template | corr_z_template | N |
| --- | --- | --- | --- | --- |
| gaia_bp_rp_colour | sfd_log1p | 0.331477 | 0.00145021 | 183527 |
| gaia_bp_rp_colour | ebv_sfd | 0.323996 | 0.00119195 | 183527 |
| gaia_bp_rp_colour | random_density_log1p | -0.295785 | 0.000251233 | 183527 |
| gaia_bp_rp_colour | gaia_dr3_density_hpx4_log1p | 0.138462 | -0.0015647 | 183527 |
| gaia_bp_rp_colour | sfd_sq | 0.119371 | -0.00068312 | 183527 |
| gaia_bp_rp_colour | wise_w1_lowcov_log1p | 0.0775697 | -0.00328039 | 183527 |
| gaia_wise_colour_cross | ebv_sfd | 0.0454231 | 0.00119195 | 183527 |
| gaia_wise_colour_cross | sfd_log1p | 0.0431382 | 0.00145021 | 183527 |
| gaia_wise_colour_cross | random_density_log1p | -0.0384102 | 0.000251233 | 183527 |
| gaia_wise_colour_cross | sfd_sq | 0.026857 | -0.00068312 | 183527 |
| gaia_wise_colour_cross | gaia_scan_angle_cos2_mean_dr3 | 0.00483436 | -0.00248145 | 183527 |
| gaia_wise_colour_cross | gaia_dr3_density_hpx4_log1p | -0.00469891 | -0.0015647 | 183527 |
| wise_w1_w2_colour | sfd_log1p | 0.0707159 | 0.00145021 | 183527 |
| wise_w1_w2_colour | ebv_sfd | 0.0684412 | 0.00119195 | 183527 |
| wise_w1_w2_colour | random_density_log1p | -0.0650176 | 0.000251233 | 183527 |
| wise_w1_w2_colour | gaia_dr3_density_hpx4_log1p | 0.0303894 | -0.0015647 | 183527 |
| wise_w1_w2_colour | sfd_sq | 0.0216433 | -0.00068312 | 183527 |
| wise_w1_w2_colour | wise_w1_lowcov_log1p | 0.0147701 | -0.00328039 | 183527 |
