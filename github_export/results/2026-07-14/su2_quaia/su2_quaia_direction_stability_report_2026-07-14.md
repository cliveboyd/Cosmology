# SU2 / Quaia Direction-Amplitude Stability Grid

Date: July 14, 2026

## Purpose

This gate tests whether the remaining `z ~ 1.0-1.5` SU2/Quaia angular hook behaves like a stable angular mode under nearby redshift-window and Galactic-latitude perturbations.

The reference direction is the strongest original family row: `1.0 <= z < 1.5`, `|b| > 35`.

## Fixed Width 0.50 Window Shifts

| tag | z_lo | z_hi | amp_mean | amp_cv | snr_max | sep_ref_max_deg | max_pair_direction_sep_deg | passes_30deg_ref_gate | passes_45deg_pair_gate |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| shift_m0p10_width0p50 | 0.9 | 1.4 | 0.0023115 | 0.313022 | 3.68867 | 17.3313 | 12.4341 | True | True |
| shift_m0p05_width0p50 | 0.95 | 1.45 | 0.00302017 | 0.185348 | 4.7852 | 11.9959 | 13.4723 | True | True |
| baseline_width0p50 | 1 | 1.5 | 0.00234501 | 0.158962 | 3.70594 | 17.9416 | 17.9416 | True | True |
| shift_p0p05_width0p50 | 1.05 | 1.55 | 0.00139514 | 0.0883359 | 2.75824 | 61.3353 | 34.8113 | False | True |
| shift_p0p10_width0p50 | 1.1 | 1.6 | 0.00104422 | 0.162183 | 2.1157 | 76.0101 | 54.9998 | False | False |

## Most Stable Rows By Reference Separation

| tag | family | z_lo | z_hi | amp_mean | amp_cv | snr_max | sep_ref_max_deg | max_pair_direction_sep_deg |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| shift_m0p05_width0p50 | width0p50 | 0.95 | 1.45 | 0.00302017 | 0.185348 | 4.7852 | 11.9959 | 13.4723 |
| shift_m0p10_width0p50 | width0p50 | 0.9 | 1.4 | 0.0023115 | 0.313022 | 3.68867 | 17.3313 | 12.4341 |
| baseline_width0p50 | width0p50 | 1 | 1.5 | 0.00234501 | 0.158962 | 3.70594 | 17.9416 | 17.9416 |
| slide_c1p20_w0p30 | width0p30 | 1.05 | 1.35 | 0.00158965 | 0.237578 | 3.41209 | 20.8484 | 25.6744 |
| slide_c1p25_w0p20 | width0p20 | 1.15 | 1.35 | 0.00127717 | 0.0359257 | 3.29694 | 30.7813 | 18.5137 |
| slide_c1p15_w0p30 | width0p30 | 1 | 1.3 | 0.00143488 | 0.428232 | 3.3759 | 36.1363 | 12.3231 |
| slide_c1p25_w0p30 | width0p30 | 1.1 | 1.4 | 0.00166316 | 0.0878992 | 3.92753 | 36.8986 | 16.0016 |
| slide_c1p20_w0p20 | width0p20 | 1.1 | 1.3 | 0.000864528 | 0.140452 | 2.0805 | 39.0802 | 18.2819 |
| slide_c1p30_w0p30 | width0p30 | 1.15 | 1.45 | 0.00159905 | 0.0790334 | 3.93256 | 54.748 | 43.4883 |
| shift_p0p05_width0p50 | width0p50 | 1.05 | 1.55 | 0.00139514 | 0.0883359 | 2.75824 | 61.3353 | 34.8113 |

## Stability By Latitude Cut

| bcut_deg | n_windows | baseline_amp | baseline_snr | amp_ratio_min_vs_bcut_baseline | amp_ratio_max_vs_bcut_baseline | sep_from_bcut_baseline_max_deg | passes_30deg_bcut_baseline_gate |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 10 | 27 | 0.00209207 | 3.47955 | 0.0734795 | 1.17709 | 171.013 | False |
| 15 | 27 | 0.00224973 | 3.59778 | 0.0832233 | 1.13289 | 156.612 | False |
| 20 | 27 | 0.00184579 | 3.02083 | 0.152967 | 1.4103 | 171.241 | False |
| 25 | 27 | 0.0023678 | 3.49857 | 0.11505 | 1.35041 | 177.373 | False |
| 30 | 27 | 0.00246501 | 3.38376 | 0.125558 | 1.31616 | 173.31 | False |
| 35 | 27 | 0.00304966 | 3.70594 | 0.128587 | 1.33287 | 160.647 | False |

## Interpretation Gate

- A stable mode should keep direction clustered under small redshift-window shifts and latitude-cut changes.
- Smooth amplitude changes are acceptable; sharp direction flips or large window sensitivity are not.
- This observed-grid test is diagnostic only. Any candidate-stable region still needs mock stability percentiles and external selection templates.

## Outputs

- fit grid: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_direction_stability_20260714\su2_quaia_direction_stability_fit_grid.csv`
- window summary: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_direction_stability_20260714\su2_quaia_direction_stability_window_summary.csv`
- bcut summary: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_direction_stability_20260714\su2_quaia_direction_stability_bcut_summary.csv`
- report: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_direction_stability_20260714\su2_quaia_direction_stability_report.md`
