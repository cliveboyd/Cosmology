# HHT CSV diagnostic

Input: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_scan\su2_quaia_hht_series.csv`
x column: `x`
y column: `frac_par_cmb`
reference column: ``

## Dominant IMF components

| rank | series | component | energy fraction | weighted frequency | median positive frequency |
|---:|---|---|---:|---:|---:|
| 1 | SU2_Quaia_frac_par_cmb | IMF1 | 0.45332 | 3.525 | 4.2808 |
| 2 | SU2_Quaia_frac_par_cmb | IMF2 | 0.095854 | 1.7131 | 1.6565 |

## Caution

HHT on cosmology residuals is exploratory. Treat any oscillatory component as a candidate feature only after comparing against LCDM/SU2 residuals, redshift shuffles, bootstrap resamples, and survey-window/systematic controls.

## Output files

- hht_resonance_csv_SU2_Quaia_frac_par_cmb_imf_summary.csv
- hht_resonance_csv_SU2_Quaia_frac_par_cmb_overview.png
