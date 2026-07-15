# SU2 / Quaia Colour-Template Mechanism Gate

Date: July 15, 2026

## Purpose

This mechanism gate tests whether the colour fields that suppressed the locked Quaia angular mode are spatially carried by external survey templates. The response fields are standardised Gaia `BP-RP`, standardised WISE `W1-W2`, and their standardised cross-product.

For each latitude cut, the response field is regressed against a template group. The residual colour-field dipole and the residual colour imbalance across the locked redshift-dipole hemispheres are then measured.

## Model

```text
c_i       = alpha + T_i beta + eta_i
eta_i     = c_i - alpha - T_i beta
c_i       = c0 + d_c . n_i + epsilon_i
eta_i     = e0 + d_eta . n_i + epsilon_i
SMD_c     = (mean(c | H=1) - mean(c | H=0)) / sigma_pooled
```

`H` is the hemisphere defined by the redshift dipole axis for that latitude cut.

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

## Full Family Ranking

| rank_within_response | response_field | template_group | mechanism_score | template_r2_mean | resid_amp_ratio_vs_raw | dipole_absorption_fraction | smd_absorption_fraction | max_resid_direction_sep_from_raw_deg |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | gaia_bp_rp_colour | gaia_scanlaw | 0.497779 | 0.00640632 | 0.546831 | 0.453169 | 0.835919 | 16.776 |
| 2 | gaia_bp_rp_colour | all_external_templates | 0.399889 | 0.0791478 | 0.580368 | 0.419632 | 0.557786 | 58.7893 |
| 3 | gaia_bp_rp_colour | dust_selection_scanlaw | 0.215295 | 0.0785742 | 0.754644 | 0.245356 | 0.254773 | 70.0205 |
| 4 | gaia_bp_rp_colour | sfd_dust | 0.0929932 | 0.0784038 | 0.828195 | 0.171805 | -0.0276075 | 65.5799 |
| 5 | gaia_bp_rp_colour | quaia_selection_depth | 0.0512422 | 0.0578869 | 0.911856 | 0.0881441 | -0.337811 | 68.917 |
| 6 | gaia_bp_rp_colour | wise_depth | 0.0154503 | 0.00376069 | 0.967337 | 0.0326625 | -0.0176232 | 15.1807 |
| 7 | gaia_bp_rp_colour | wise_depth_stellar_density | 0.00210676 | 0.0105338 | 1.33532 | -0.335315 | -1.1221 | 62.7433 |
| 8 | gaia_bp_rp_colour | stellar_density | 0.00166199 | 0.00830994 | 1.27132 | -0.271322 | -1.16174 | 61.3639 |
| 1 | gaia_wise_colour_cross | all_external_templates | 0.188208 | 0.00267529 | 0.582949 | 0.417051 | -0.82724 | 75.4432 |
| 2 | gaia_wise_colour_cross | gaia_scanlaw | 0.171019 | 3.7239e-05 | 0.919413 | 0.0805873 | 0.384993 | 18.0424 |
| 3 | gaia_wise_colour_cross | dust_selection_scanlaw | 0.134216 | 0.00255262 | 0.702876 | 0.297124 | -0.0248601 | 47.292 |
| 4 | gaia_wise_colour_cross | wise_depth_stellar_density | 0.10325 | 4.51326e-05 | 0.770575 | 0.229425 | -0.445323 | 18.5381 |
| 5 | gaia_wise_colour_cross | stellar_density | 0.0992833 | 2.38909e-05 | 0.779381 | 0.220619 | -0.42552 | 18.6524 |
| 6 | gaia_wise_colour_cross | quaia_selection_depth | 0.0751217 | 0.00066166 | 1.50195 | -0.501948 | 0.214255 | 19.4437 |
| 7 | gaia_wise_colour_cross | sfd_dust | 0.0466911 | 0.00251408 | 0.897359 | 0.102641 | -0.678798 | 42.2745 |
| 8 | gaia_wise_colour_cross | wise_depth | 0.00978132 | 2.15802e-05 | 0.978273 | 0.0217267 | -0.0117503 | 2.49687 |
| 1 | wise_w1_w2_colour | gaia_scanlaw | 0.410102 | 0.000266778 | 0.57196 | 0.42804 | 0.62123 | 22.828 |
| 2 | wise_w1_w2_colour | all_external_templates | 0.301259 | 0.00393026 | 0.687356 | 0.312644 | 0.456524 | 42.6647 |
| 3 | wise_w1_w2_colour | dust_selection_scanlaw | 0.169632 | 0.00385735 | 0.862545 | 0.137455 | 0.305732 | 48.8697 |
| 4 | wise_w1_w2_colour | wise_depth | 0.00892997 | 0.000142099 | 0.98891 | 0.0110898 | 0.0111747 | 7.87183 |
| 5 | wise_w1_w2_colour | sfd_dust | 0.00208202 | 0.00369397 | 0.997015 | 0.00298494 | -0.0235928 | 36.5183 |
| 6 | wise_w1_w2_colour | quaia_selection_depth | 0.00058827 | 0.00294135 | 1.04087 | -0.0408694 | -0.1085 | 37.9112 |
| 7 | wise_w1_w2_colour | wise_depth_stellar_density | 9.61905e-05 | 0.000480952 | 1.24112 | -0.241121 | -0.384508 | 32.4707 |
| 8 | wise_w1_w2_colour | stellar_density | 8.52448e-05 | 0.000426224 | 1.2532 | -0.253202 | -0.425161 | 32.5351 |

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

## Configuration

```json
{
  "date": "2026-07-15",
  "analysis": "su2_quaia_colour_template_mechanism",
  "z_min": 0.95,
  "z_max": 1.45,
  "b_cuts": [
    15.0,
    20.0,
    25.0,
    30.0,
    35.0
  ],
  "primary_b_cut": 15.0,
  "responses": [
    "gaia_bp_rp_colour",
    "wise_w1_w2_colour",
    "gaia_wise_colour_cross"
  ],
  "template_groups": {
    "sfd_dust": [
      "ebv_sfd",
      "sfd_log1p",
      "sfd_sq"
    ],
    "quaia_selection_depth": [
      "selection_T",
      "random_density_log1p"
    ],
    "gaia_scanlaw": [
      "gaia_scan_count_log1p_dr3",
      "gaia_scan_angle_cos2_mean_dr3",
      "gaia_scan_angle_sin2_mean_dr3",
      "gaia_scan_angle_resultant_dr3"
    ],
    "wise_depth": [
      "wise_w1_medcov_log1p",
      "wise_w2_medcov_log1p",
      "wise_w1_mincov",
      "wise_w2_mincov",
      "wise_w1_lowcov_log1p",
      "wise_w2_lowcov_log1p",
      "wise_w1_w2_medcov_diff",
      "wise_w1_w2_medcov_ratio"
    ],
    "stellar_density": [
      "gaia_dr3_density_hpx4_log1p"
    ],
    "dust_selection_scanlaw": [
      "ebv_sfd",
      "sfd_log1p",
      "sfd_sq",
      "selection_T",
      "random_density_log1p",
      "gaia_scan_count_log1p_dr3",
      "gaia_scan_angle_cos2_mean_dr3",
      "gaia_scan_angle_sin2_mean_dr3",
      "gaia_scan_angle_resultant_dr3"
    ],
    "wise_depth_stellar_density": [
      "wise_w1_medcov_log1p",
      "wise_w2_medcov_log1p",
      "wise_w1_mincov",
      "wise_w2_mincov",
      "wise_w1_lowcov_log1p",
      "wise_w2_lowcov_log1p",
      "wise_w1_w2_medcov_diff",
      "wise_w1_w2_medcov_ratio",
      "gaia_dr3_density_hpx4_log1p"
    ],
    "all_external_templates": [
      "ebv_sfd",
      "sfd_log1p",
      "sfd_sq",
      "selection_T",
      "random_density_log1p",
      "gaia_scan_count_log1p_dr3",
      "gaia_scan_angle_cos2_mean_dr3",
      "gaia_scan_angle_sin2_mean_dr3",
      "gaia_scan_angle_resultant_dr3",
      "wise_w1_medcov_log1p",
      "wise_w2_medcov_log1p",
      "wise_w1_mincov",
      "wise_w2_mincov",
      "wise_w1_lowcov_log1p",
      "wise_w2_lowcov_log1p",
      "wise_w1_w2_medcov_diff",
      "wise_w1_w2_medcov_ratio",
      "gaia_dr3_density_hpx4_log1p"
    ]
  },
  "mechanism_score": "0.45*dipole_absorption + 0.35*hemisphere-SMD absorption + 0.20*template R2, each clipped to [0,1] where applicable.",
  "external_inputs": {
    "wise_inventory": "C:\\Users\\clive\\Documents\\Cosmology\\external_datasets\\allwise_atlas_inventory\\allwise_p3am_cdd_w1w2_coverage.csv",
    "gaia_density": "C:\\Users\\clive\\Documents\\Cosmology\\external_datasets\\gaia_dr3_density\\gaia_dr3_source_density_hpx4_glt20p5.csv",
    "gaia_density_source": "CDS Gaia EDR3 density HiPS Norder2",
    "gaia_density_present": true,
    "gaiascanlaw_data_dir": "C:\\Users\\clive\\anaconda3x\\Lib\\site-packages\\gaiascanlaw\\data"
  }
}
```
