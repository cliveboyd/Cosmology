# HHT CSV diagnostic

Input: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_scan\su2_quaia_hht_series.csv`
x column: `x`
y column: `amp_par_cmb`
reference column: ``

## Dominant IMF components

| rank | series | component | energy fraction | weighted frequency | median positive frequency |
|---:|---|---|---:|---:|---:|
| 1 | SU2_Quaia_amp_par_cmb | IMF1 | 0.44898 | 4.0793 | 4.5314 |
| 2 | SU2_Quaia_amp_par_cmb | IMF2 | 0.12599 | 1.7802 | 1.7983 |

## Caution

HHT on cosmology residuals is exploratory. Treat any oscillatory component as a candidate feature only after comparing against LCDM/SU2 residuals, redshift shuffles, bootstrap resamples, and survey-window/systematic controls.

## Output files

- hht_resonance_csv_SU2_Quaia_amp_par_cmb_imf_summary.csv
- hht_resonance_csv_SU2_Quaia_amp_par_cmb_overview.png
