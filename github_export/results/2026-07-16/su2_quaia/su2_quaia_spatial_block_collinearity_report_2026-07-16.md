# SU2 / Quaia Spatial-Block and Template-Collinearity Audit

Date: 2026-07-16

## Bottom Line

The z=1.0-1.5 residual family does not remain below the 1% gate at every spatial-cluster scale. The earlier object-shuffle p-value is therefore not robust enough for promotion.

## Why This Replaces the Simple Target-Family Shuffle Readout

The earlier residual scan used 1,000 object-level shuffles and independently permuted overlapping windows. Its p=0.000999 value was the Monte Carlo floor, not a resolved probability. The present audit addresses the targeted z=1.0-1.5 family with shared spatial blocks and a substantially larger null ensemble.

## Spatial Null

For each latitude cut, redshift was regressed on the seven Gaia scan-law and catalogue-colour controls. Control residuals were multiplied by shared Rademacher signs on an equal-area Galactic longitude x sin(latitude) grid. The same cell signs were used for all nested latitude cuts in each mock, preserving their covariance and within-cell residual structure.

| cell_deg | n_cells | n_mocks | observed_family_max_snr | null_max_snr_p95 | null_max_snr_p99 | family_p | monte_carlo_floor |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 8 | 910 | 50000 | 5.75359 | 4.97344 | 5.9165 | 0.0135197 | 1.99996e-05 |
| 12 | 424 | 50000 | 5.75359 | 6.29349 | 7.51404 | 0.0922982 | 1.99996e-05 |
| 16 | 231 | 50000 | 5.75359 | 7.59853 | 9.14024 | 0.232595 | 1.99996e-05 |

## Observed Controlled Dipoles

| bcut_deg | N | residual_snr | residual_amp | residual_l_deg | residual_b_deg | control_r2_z | canonical_corr | condition_number |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 10 | 194552 | 5.75359 | 0.00335336 | 222.657 | -4.64068 | 0.135974 | 0.610202 | 2.67274 |
| 15 | 185739 | 5.68228 | 0.00346028 | 219.602 | -4.55171 | 0.138497 | 0.637852 | 2.73476 |
| 20 | 172054 | 4.58277 | 0.00297222 | 223.911 | -5.6127 | 0.139857 | 0.664291 | 2.81115 |
| 25 | 155609 | 4.29184 | 0.00305014 | 227.5 | -3.95679 | 0.141136 | 0.684438 | 2.9228 |
| 30 | 137850 | 2.9821 | 0.00236787 | 233.617 | -3.70101 | 0.141786 | 0.697506 | 3.08892 |
| 35 | 119719 | 2.79086 | 0.00252649 | 233.021 | -4.166 | 0.142452 | 0.700607 | 3.21877 |

## Direction-Specific Template Collinearity

| direction_label | bcut_deg | projection_control_r2 | projection_vif | max_corr_control | max_corr_signed | dipole_control_canonical_corr | full_design_condition_number |
| --- | --- | --- | --- | --- | --- | --- | --- |
| locked_b15 | 10 | 0.142733 | 1.1665 | gaia_scan_angle_sin2_mean_dr3 | -0.377628 | 0.610202 | 2.67274 |
| locked_b25 | 10 | 0.112518 | 1.12678 | gaia_scan_angle_sin2_mean_dr3 | -0.3352 | 0.610202 | 2.67274 |
| cmb | 10 | 0.377652 | 1.60682 | gaia_scan_angle_sin2_mean_dr3 | 0.614099 | 0.610202 | 2.67274 |
| anti_cmb | 10 | 0.377652 | 1.60682 | gaia_scan_angle_sin2_mean_dr3 | -0.614099 | 0.610202 | 2.67274 |
| locked_b15 | 15 | 0.163632 | 1.19565 | gaia_scan_angle_sin2_mean_dr3 | -0.404376 | 0.637852 | 2.73476 |
| locked_b25 | 15 | 0.12824 | 1.1471 | gaia_scan_angle_sin2_mean_dr3 | -0.357906 | 0.637852 | 2.73476 |
| cmb | 15 | 0.410972 | 1.69771 | gaia_scan_angle_sin2_mean_dr3 | 0.64069 | 0.637852 | 2.73476 |
| anti_cmb | 15 | 0.410972 | 1.69771 | gaia_scan_angle_sin2_mean_dr3 | -0.64069 | 0.637852 | 2.73476 |
| locked_b15 | 20 | 0.190895 | 1.23593 | gaia_scan_angle_sin2_mean_dr3 | -0.436701 | 0.664291 | 2.81115 |
| locked_b25 | 20 | 0.149313 | 1.17552 | gaia_scan_angle_sin2_mean_dr3 | -0.386152 | 0.664291 | 2.81115 |
| cmb | 20 | 0.440019 | 1.78578 | gaia_scan_angle_sin2_mean_dr3 | 0.663049 | 0.664291 | 2.81115 |
| anti_cmb | 20 | 0.440019 | 1.78578 | gaia_scan_angle_sin2_mean_dr3 | -0.663049 | 0.664291 | 2.81115 |
| locked_b15 | 25 | 0.205892 | 1.25927 | gaia_scan_angle_sin2_mean_dr3 | -0.453431 | 0.684438 | 2.9228 |
| locked_b25 | 25 | 0.160566 | 1.19128 | gaia_scan_angle_sin2_mean_dr3 | -0.400393 | 0.684438 | 2.9228 |
| cmb | 25 | 0.453891 | 1.83114 | gaia_scan_angle_sin2_mean_dr3 | 0.673486 | 0.684438 | 2.9228 |
| anti_cmb | 25 | 0.453891 | 1.83114 | gaia_scan_angle_sin2_mean_dr3 | -0.673486 | 0.684438 | 2.9228 |
| locked_b15 | 30 | 0.21505 | 1.27397 | gaia_scan_angle_sin2_mean_dr3 | -0.463361 | 0.697506 | 3.08892 |
| locked_b25 | 30 | 0.168343 | 1.20242 | gaia_scan_angle_sin2_mean_dr3 | -0.409921 | 0.697506 | 3.08892 |
| cmb | 30 | 0.453453 | 1.82967 | gaia_scan_angle_sin2_mean_dr3 | 0.673195 | 0.697506 | 3.08892 |
| anti_cmb | 30 | 0.453453 | 1.82967 | gaia_scan_angle_sin2_mean_dr3 | -0.673195 | 0.697506 | 3.08892 |
| locked_b15 | 35 | 0.225227 | 1.2907 | gaia_scan_angle_sin2_mean_dr3 | -0.474216 | 0.700607 | 3.21877 |
| locked_b25 | 35 | 0.178452 | 1.21721 | gaia_scan_angle_sin2_mean_dr3 | -0.422076 | 0.700607 | 3.21877 |
| cmb | 35 | 0.445489 | 1.80339 | gaia_scan_angle_sin2_mean_dr3 | 0.667254 | 0.700607 | 3.21877 |
| anti_cmb | 35 | 0.445489 | 1.80339 | gaia_scan_angle_sin2_mean_dr3 | -0.667254 | 0.700607 | 3.21877 |

## Interpretation

For the locked b=15-degree direction at the matching cut, the control-template R2 is 0.1636 (VIF 1.196); its largest single-template correlation is -0.4044 with gaia_scan_angle_sin2_mean_dr3.

The collinearity table is diagnostic rather than a correction. A high canonical correlation or direction-specific VIF explains why unconstrained controls can rotate or inflate recovered injections; it does not prove that the angular mode is physical.

## Claim Boundary

Treat the z=1.0-1.5 mode as a targeted follow-up candidate only. Promotion requires an independently specified catalogue or release, a pre-registered angular statistic, and a survey-selection likelihood whose injection power and false-positive rate are calibrated before examining the new data.

## Configuration

```json
{
  "date": "2026-07-16",
  "analysis": "su2_quaia_spatial_block_collinearity",
  "z_min": 1.0,
  "z_max": 1.5,
  "b_cuts": [
    10.0,
    15.0,
    20.0,
    25.0,
    30.0,
    35.0
  ],
  "cell_degrees": [
    8.0,
    12.0,
    16.0
  ],
  "n_mocks_per_cell_scale": 50000,
  "batch_size": 1000,
  "seed": 200716,
  "controls": [
    "gaia_scan_count_log1p_dr3",
    "gaia_scan_angle_cos2_mean_dr3",
    "gaia_scan_angle_sin2_mean_dr3",
    "gaia_scan_angle_resultant_dr3",
    "bp_rp_z",
    "w1_w2_z",
    "bp_rp_x_w1_w2_z"
  ],
  "null": "Shared Rademacher wild bootstrap on equal-area Galactic longitude x sin(latitude) cells; common signs across nested latitude cuts.",
  "load_meta": {
    "gaiascanlaw_data_dir": "C:\\Users\\clive\\anaconda3x\\Lib\\site-packages\\gaiascanlaw\\data",
    "scanlaw_templates_present": true
  }
}
```
