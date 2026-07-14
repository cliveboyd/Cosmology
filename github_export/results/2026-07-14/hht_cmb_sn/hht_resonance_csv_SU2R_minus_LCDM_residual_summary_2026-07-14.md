# HHT CSV diagnostic

Input: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\hht_resonance\sn_residuals_for_hht.csv`
x column: `log1pz`
y column: `delta_residual_SU2R_baocov_4D_minus_LCDM_baocov`
reference column: ``

## Dominant IMF components

| rank | series | component | energy fraction | weighted frequency | median positive frequency |
|---:|---|---|---:|---:|---:|
| - | no IMF extracted | residue only | 1 | nan | nan |

## Caution

HHT on cosmology residuals is exploratory. Treat any oscillatory component as a candidate feature only after comparing against LCDM/SU2 residuals, redshift shuffles, bootstrap resamples, and survey-window/systematic controls.

## Output files

- hht_resonance_csv_SU2R_minus_LCDM_residual_imf_summary.csv
- hht_resonance_csv_SU2R_minus_LCDM_residual_overview.png
