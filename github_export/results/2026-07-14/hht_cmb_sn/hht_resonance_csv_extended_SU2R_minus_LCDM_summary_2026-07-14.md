# HHT CSV diagnostic

Input: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\hht_resonance\extended_realdata_pulls_for_hht.csv`
x column: `x_hht`
y column: `joint_delta_pull_SU2R_baocov_4D_minus_LCDM_baocov`
reference column: ``

## Dominant IMF components

| rank | series | component | energy fraction | weighted frequency | median positive frequency |
|---:|---|---|---:|---:|---:|
| 1 | extended_SU2R_minus_LCDM | IMF1 | 0.61539 | 16.169 | 8.6041 |
| 2 | extended_SU2R_minus_LCDM | IMF2 | 0.43131 | 8.7386 | 6.0592 |
| 3 | extended_SU2R_minus_LCDM | IMF3 | 0.2677 | 4.5174 | 2.8201 |
| 4 | extended_SU2R_minus_LCDM | IMF4 | 0.025287 | 2.541 | 1.7448 |

## Caution

HHT on cosmology residuals is exploratory. Treat any oscillatory component as a candidate feature only after comparing against LCDM/SU2 residuals, redshift shuffles, bootstrap resamples, and survey-window/systematic controls.

## Output files

- hht_resonance_csv_extended_SU2R_minus_LCDM_imf_summary.csv
- hht_resonance_csv_extended_SU2R_minus_LCDM_overview.png
