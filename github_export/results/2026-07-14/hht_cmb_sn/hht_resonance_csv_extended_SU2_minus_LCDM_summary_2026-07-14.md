# HHT CSV diagnostic

Input: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\hht_resonance\extended_realdata_pulls_for_hht.csv`
x column: `x_hht`
y column: `joint_delta_pull_SU2_baocov_4D_minus_LCDM_baocov`
reference column: ``

## Dominant IMF components

| rank | series | component | energy fraction | weighted frequency | median positive frequency |
|---:|---|---|---:|---:|---:|
| 1 | extended_SU2_minus_LCDM | IMF1 | 0.6177 | 16.164 | 8.588 |
| 2 | extended_SU2_minus_LCDM | IMF2 | 0.42901 | 8.8001 | 6.1189 |
| 3 | extended_SU2_minus_LCDM | IMF3 | 0.20675 | 5.0914 | 3.6729 |
| 4 | extended_SU2_minus_LCDM | IMF4 | 0.070741 | 1.9735 | 1.6576 |

## Caution

HHT on cosmology residuals is exploratory. Treat any oscillatory component as a candidate feature only after comparing against LCDM/SU2 residuals, redshift shuffles, bootstrap resamples, and survey-window/systematic controls.

## Output files

- hht_resonance_csv_extended_SU2_minus_LCDM_imf_summary.csv
- hht_resonance_csv_extended_SU2_minus_LCDM_overview.png
