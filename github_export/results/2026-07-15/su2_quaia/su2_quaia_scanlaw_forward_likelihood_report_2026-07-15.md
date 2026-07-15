# SU2 / Quaia Gaia Scan-Law Forward-Template Likelihood

Date: July 15, 2026

## Purpose

This test treats Gaia scan-law and catalogue-colour fields as an explicit forward-template model, then asks whether a redshift dipole adds explanatory power after those fields are already in the likelihood.

## Model Set

```text
M0: z_i = alpha + epsilon_i
Md: z_i = alpha + d . n_i + epsilon_i
Mt: z_i = alpha + S_i beta_s + C_i beta_c + epsilon_i
Mtd: z_i = alpha + d . n_i + S_i beta_s + C_i beta_c + epsilon_i
```

`S_i` contains Gaia scan count and scan-angle geometry. `C_i` contains Gaia BP-RP, WISE W1-W2, and their colour cross-term.

## Bottom Line

Scan-law plus colour controls are favoured over adding a dipole; the angular mode does not pass this mechanism gate.

## Summary

| comparison | mean_delta_bic_full_minus_restricted | mean_partial_r2 | max_full_dipole_snr | mean_amp_ratio_vs_raw_dipole | max_direction_sep_vs_raw_dipole_deg | bic_readout |
| --- | --- | --- | --- | --- | --- | --- |
| dipole_added_to_all_external_colour | 21.7205 | 9.06266e-05 | 4.26371 | 0.835209 | 110.642 | restricted model favoured |
| dipole_added_to_scanlaw_colour | 16.7231 | 0.00012216 | 5.29092 | 1.12947 | 94.6803 | restricted model favoured |
| dipole_added_to_scanlaw_colour_quality | 22.2655 | 8.59438e-05 | 4.32075 | 0.741865 | 119.676 | restricted model favoured |
| dipole_only_vs_null | 15.4904 | 0.000138258 | 4.7852 | 1 | 0 | restricted model favoured |

## B-Cut Details

| comparison | bcut_deg | N | delta_bic_full_minus_restricted | partial_r2 | lrt_p | full_dipole_snr | amp_ratio_vs_raw_dipole | direction_sep_vs_raw_dipole_deg |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dipole_only_vs_null | 15 | 183527 | 18.423 | 9.77323e-05 | 0.000453123 | 4.23531 | 1 | 0 |
| dipole_added_to_scanlaw_colour | 15 | 183527 | 8.36702 | 0.000152518 | 3.64376e-06 | 5.29092 | 1.35439 | 86.9901 |
| dipole_added_to_scanlaw_colour_quality | 15 | 183527 | 17.6907 | 0.000101722 | 0.000319942 | 4.32075 | 0.839913 | 119.289 |
| dipole_added_to_all_external_colour | 15 | 183527 | 18.1804 | 9.90536e-05 | 0.000403823 | 4.26371 | 0.909876 | 101.625 |
| dipole_only_vs_null | 20 | 169941 | 18.4501 | 0.000104028 | 0.000512123 | 4.20476 | 1 | 0 |
| dipole_added_to_scanlaw_colour | 20 | 169941 | 13.1657 | 0.00013512 | 4.10879e-05 | 4.79208 | 1.28678 | 84.7603 |
| dipole_added_to_scanlaw_colour_quality | 20 | 169941 | 17.9704 | 0.00010685 | 0.000407815 | 4.26126 | 0.903046 | 118.408 |
| dipole_added_to_all_external_colour | 20 | 169941 | 18.6916 | 0.000102607 | 0.000574284 | 4.17579 | 0.93904 | 102.649 |
| dipole_only_vs_null | 25 | 153651 | 13.0928 | 0.000147951 | 4.58691e-05 | 4.76819 | 1 | 0 |
| dipole_added_to_scanlaw_colour | 25 | 153651 | 15.7335 | 0.000130767 | 0.000162313 | 4.48259 | 1.10527 | 83.406 |
| dipole_added_to_scanlaw_colour_quality | 25 | 153651 | 21.7432 | 9.16591e-05 | 0.0027928 | 3.75277 | 0.725024 | 110.96 |
| dipole_added_to_all_external_colour | 25 | 153651 | 20.4327 | 0.000100187 | 0.00150866 | 3.92349 | 0.83299 | 99.6867 |
| dipole_only_vs_null | 30 | 136052 | 15.3425 | 0.000147873 | 0.000160304 | 4.48562 | 1 | 0 |
| dipole_added_to_scanlaw_colour | 30 | 136052 | 24.7009 | 7.90954e-05 | 0.0130879 | 3.28041 | 0.926127 | 94.6803 |
| dipole_added_to_scanlaw_colour_quality | 30 | 136052 | 25.7527 | 7.13652e-05 | 0.0212018 | 3.1159 | 0.693281 | 119.676 |
| dipole_added_to_all_external_colour | 30 | 136052 | 25.3799 | 7.41046e-05 | 0.0178781 | 3.17516 | 0.769895 | 110.642 |
| dipole_only_vs_null | 35 | 118193 | 12.1435 | 0.000193704 | 4.24356e-05 | 4.7852 | 1 | 0 |
| dipole_added_to_scanlaw_colour | 35 | 118193 | 21.6482 | 0.0001133 | 0.00386125 | 3.65944 | 0.974782 | 90.5854 |
| dipole_added_to_scanlaw_colour_quality | 35 | 118193 | 28.1704 | 5.81223e-05 | 0.0761637 | 2.62089 | 0.548063 | 103.397 |
| dipole_added_to_all_external_colour | 35 | 118193 | 25.9176 | 7.71809e-05 | 0.0277045 | 3.02022 | 0.724245 | 95.3741 |

## Configuration

```json
{
  "date": "2026-07-15",
  "analysis": "scanlaw_forward_template_likelihood",
  "z_min": 0.95,
  "z_max": 1.45,
  "b_cuts": [
    15.0,
    20.0,
    25.0,
    30.0,
    35.0
  ],
  "external_inputs": {
    "gaiascanlaw_data_dir": "C:\\Users\\clive\\anaconda3x\\Lib\\site-packages\\gaiascanlaw\\data",
    "scanlaw_templates_present": true
  },
  "models": {
    "dipole_only": "z = alpha + d.n + epsilon",
    "scanlaw_colour": "z = alpha + Gaia_scanlaw beta_s + colour beta_c + epsilon",
    "scanlaw_colour_plus_dipole": "z = alpha + d.n + Gaia_scanlaw beta_s + colour beta_c + epsilon"
  },
  "gate": "The dipole term should add material BIC-supported explanatory power after scan-law and colour controls before the angular mode can be promoted."
}
```
