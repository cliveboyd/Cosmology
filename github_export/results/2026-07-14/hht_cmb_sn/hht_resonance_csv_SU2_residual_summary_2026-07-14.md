# HHT CSV diagnostic

Input: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\hht_resonance\sn_residuals_for_hht.csv`
x column: `log1pz`
y column: `residual_SU2_baocov_4D`
reference column: ``

## Dominant IMF components

| rank | series | component | energy fraction | weighted frequency | median positive frequency |
|---:|---|---|---:|---:|---:|
| 1 | SU2_residual | IMF4 | 0.25027 | 7.908 | 2.7519 |
| 2 | SU2_residual | IMF1 | 0.22115 | 183.51 | 67.301 |
| 3 | SU2_residual | IMF2 | 0.15926 | 68.757 | 12.11 |
| 4 | SU2_residual | IMF3 | 0.064424 | 53.348 | 26.84 |
| 5 | SU2_residual | IMF5 | 0.026585 | 18.179 | 11.132 |

## Caution

HHT on cosmology residuals is exploratory. Treat any oscillatory component as a candidate feature only after comparing against LCDM/SU2 residuals, redshift shuffles, bootstrap resamples, and survey-window/systematic controls.

## Output files

- hht_resonance_csv_SU2_residual_imf_summary.csv
- hht_resonance_csv_SU2_residual_overview.png
