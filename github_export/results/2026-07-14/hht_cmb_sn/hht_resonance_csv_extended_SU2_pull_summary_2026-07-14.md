# HHT CSV diagnostic

Input: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\hht_resonance\extended_realdata_pulls_for_hht.csv`
x column: `x_hht`
y column: `joint_pull_SU2_baocov_4D`
reference column: ``

## Dominant IMF components

| rank | series | component | energy fraction | weighted frequency | median positive frequency |
|---:|---|---|---:|---:|---:|
| 1 | extended_SU2_pull | IMF1 | 0.30946 | 20.184 | 17.682 |
| 2 | extended_SU2_pull | IMF3 | 0.23527 | 3.6816 | 2.2691 |
| 3 | extended_SU2_pull | IMF2 | 0.15498 | 10.502 | 8.5667 |
| 4 | extended_SU2_pull | IMF4 | 0.11342 | 2.6804 | 2.2656 |
| 5 | extended_SU2_pull | IMF5 | 0.0064117 | 2.043 | 1.8543 |

## Caution

HHT on cosmology residuals is exploratory. Treat any oscillatory component as a candidate feature only after comparing against LCDM/SU2 residuals, redshift shuffles, bootstrap resamples, and survey-window/systematic controls.

## Output files

- hht_resonance_csv_extended_SU2_pull_imf_summary.csv
- hht_resonance_csv_extended_SU2_pull_overview.png
