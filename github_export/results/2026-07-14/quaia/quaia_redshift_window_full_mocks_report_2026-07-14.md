# Quaia Redshift-Window Full Mock Catalogs

Date: July 14, 2026

## Purpose

This run replaces the earlier summary-only Quaia mock control for redshift-window anisotropy. For each selected observed redshift window, it generates object-level mock catalogs by putting the observed window redshift distribution onto random Quaia sky positions passing the same Galactic latitude cut, then refits the same redshift dipole model `z = a0 + D.n`.

The saved mock catalogs are compact `.npz` object catalogs for the first requested mock(s) in each window. All mocks are nevertheless generated and fitted object-by-object; the per-mock fit summaries are retained in CSV.

## Configuration

- observed Quaia table: `C:\Users\clive\Documents\Cosmology\quaia_G20p0_basic_gal.csv`
- random sky catalog: `C:\Users\clive\Documents\Cosmology\quaia_G20p0_randoms_simple.fits`
- scan CSV: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_scan\su2_quaia_scan.csv`
- top window CSV: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_scan\su2_quaia_scan_top.csv`
- window mode: `top`
- top_n: `25`
- n_mocks per window: `200`
- saved full catalogs per window: `1`
- seed: `260714`

## Smallest Amplitude p-values

| window_index | bcut_deg | tag | N | observed_amp | mock_amp_p95 | mock_amp_p99 | p_amp_ge_observed | priority_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | 15 | z1p00_1p50 | 185739 | 0.00224973 | 0.00218972 | 0.00252277 | 0.0497512 | 0.0345398 |
| 4 | 10 | z1p00_1p50 | 194552 | 0.00209207 | 0.00204948 | 0.00257093 | 0.0497512 | 0.0334047 |
| 3 | 25 | z1p00_1p50 | 155609 | 0.0023678 | 0.00284842 | 0.00325062 | 0.134328 | 0.0335872 |
| 1 | 35 | z1p00_1p50 | 119719 | 0.00304966 | 0.00379078 | 0.00462667 | 0.144279 | 0.0355781 |
| 5 | 30 | z1p00_1p50 | 137850 | 0.00246501 | 0.00292915 | 0.00424623 | 0.144279 | 0.0324851 |
| 25 | 10 | z0p60_0p70 | 23898 | 0.000957104 | 0.00116704 | 0.00134554 | 0.144279 | 0.0240407 |
| 9 | 10 | z0p90_1p00 | 36372 | 0.00071159 | 0.000871486 | 0.0010059 | 0.169154 | 0.0266726 |
| 20 | 15 | z1p50_2p50 | 264161 | 0.0029781 | 0.00372779 | 0.00475453 | 0.179104 | 0.0244657 |
| 23 | 20 | z0p80_0p90 | 28027 | 0.00096454 | 0.00122142 | 0.00146808 | 0.179104 | 0.0243314 |
| 7 | 20 | z1p00_1p50 | 172054 | 0.00184579 | 0.00242302 | 0.00273227 | 0.18408 | 0.0290008 |
| 16 | 10 | z1p70_1p80 | 34300 | 0.00075108 | 0.00100197 | 0.00121455 | 0.18408 | 0.0249246 |
| 24 | 20 | z1p70_1p80 | 30496 | 0.000824803 | 0.00115175 | 0.00139295 | 0.21393 | 0.0242936 |

## Smallest CMB-Projected p-values

| window_index | bcut_deg | tag | N | observed_abs_amp_par_cmb | mock_abs_par_p95 | mock_abs_par_p99 | p_abs_par_ge_observed | priority_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 16 | 10 | z1p70_1p80 | 34300 | 0.00068264 | 0.000739142 | 0.001005 | 0.0895522 | 0.0249246 |
| 24 | 20 | z1p70_1p80 | 30496 | 0.000723612 | 0.000906418 | 0.0012732 | 0.124378 | 0.0242936 |
| 18 | 15 | z1p70_1p80 | 32830 | 0.000702641 | 0.000880028 | 0.00115455 | 0.129353 | 0.0245837 |
| 4 | 10 | z1p00_1p50 | 194552 | 0.00113534 | 0.00168389 | 0.00217749 | 0.20398 | 0.0334047 |
| 1 | 35 | z1p00_1p50 | 119719 | 0.00189556 | 0.0032021 | 0.0037465 | 0.228856 | 0.0355781 |
| 2 | 15 | z1p00_1p50 | 185739 | 0.00121995 | 0.00177062 | 0.00234851 | 0.228856 | 0.0345398 |
| 7 | 20 | z1p00_1p50 | 172054 | 0.00122665 | 0.00187744 | 0.00260479 | 0.228856 | 0.0290008 |
| 25 | 10 | z0p60_0p70 | 23898 | 0.000614797 | 0.000884831 | 0.00114845 | 0.233831 | 0.0240407 |
| 23 | 20 | z0p80_0p90 | 28027 | 0.000581435 | 0.000973358 | 0.001281 | 0.263682 | 0.0243314 |
| 5 | 30 | z1p00_1p50 | 137850 | 0.00145793 | 0.00230148 | 0.00363415 | 0.268657 | 0.0324851 |
| 17 | 10 | z1p20_1p30 | 37668 | 0.000439621 | 0.000714654 | 0.000953942 | 0.288557 | 0.0247922 |
| 22 | 15 | z0p90_1p00 | 34529 | 0.000411634 | 0.000758016 | 0.00103864 | 0.313433 | 0.0243769 |

## Interpretation Gate

- Do not interpret a Quaia anisotropy window as physical unless it remains unusual against these object-level, window-matched mocks.
- A p-value near the finite mock floor means the window needs more mocks, not immediate promotion.
- These mocks preserve the observed redshift distribution inside each window and the random-catalog sky footprint, but they are still not a complete Quaia selection-function likelihood.

## Outputs

- per-mock summary: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\quaia_redshift_window_full_mocks_20260714\quaia_redshift_window_mock_summary.csv`
- observed-vs-mock p-values: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\quaia_redshift_window_full_mocks_20260714\quaia_redshift_window_mock_pvalues.csv`
- manifest: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\quaia_redshift_window_full_mocks_20260714\quaia_redshift_window_mock_manifest.csv`
- config: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\quaia_redshift_window_full_mocks_20260714\quaia_redshift_window_mock_config.json`
- saved mock catalogs directory: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\quaia_redshift_window_full_mocks_20260714\mock_catalogs`
