# Quaia Redshift-Window Full Mock Catalogs

Date: July 14, 2026

## Purpose

This run replaces the earlier summary-only Quaia mock control for redshift-window anisotropy. For each selected observed redshift window, it generates object-level mock catalogs by putting the observed window redshift distribution onto random Quaia sky positions passing the same Galactic latitude cut, then refits the same redshift dipole model `z = a0 + D.n`.

The saved mock catalogs are compact `.npz` object catalogs for the first requested mock(s) in each window. All mocks are nevertheless generated and fitted object-by-object; the per-mock fit summaries are retained in CSV.

## Configuration

- observed Quaia table: `C:\Users\clive\Documents\Cosmology\quaia_G20p0_basic_gal.csv`
- random sky catalog: `C:\Users\clive\Documents\Cosmology\quaia_G20p0_randoms_simple.fits`
- scan CSV: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_scan\su2_quaia_scan.csv`
- top window CSV: `C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs\quaia_redshift_window_highstat_focus_windows_2026-07-14.csv`
- window mode: `top`
- top_n: `3`
- n_mocks per window: `5000`
- saved full catalogs per window: `1`
- seed: `270714`

## Smallest Amplitude p-values

| window_index | bcut_deg | tag | N | observed_amp | mock_amp_p95 | mock_amp_p99 | p_amp_ge_observed | priority_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 15 | z1p00_1p50 | 185739 | 0.00224973 | 0.00217232 | 0.00267981 | 0.0415917 | 0.0345398 |
| 2 | 10 | z1p00_1p50 | 194552 | 0.00209207 | 0.00204508 | 0.00252232 | 0.0437912 | 0.0334047 |
| 3 | 10 | z1p70_1p80 | 34300 | 0.00075108 | 0.000974398 | 0.00119261 | 0.180764 | 0.0249246 |

## Smallest CMB-Projected p-values

| window_index | bcut_deg | tag | N | observed_abs_amp_par_cmb | mock_abs_par_p95 | mock_abs_par_p99 | p_abs_par_ge_observed | priority_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3 | 10 | z1p70_1p80 | 34300 | 0.00068264 | 0.000800331 | 0.00106073 | 0.0961808 | 0.0249246 |
| 1 | 15 | z1p00_1p50 | 185739 | 0.00121995 | 0.00181021 | 0.00235965 | 0.189162 | 0.0345398 |
| 2 | 10 | z1p00_1p50 | 194552 | 0.00113534 | 0.00168363 | 0.00220393 | 0.191162 | 0.0334047 |

## Interpretation Gate

- Do not interpret a Quaia anisotropy window as physical unless it remains unusual against these object-level, window-matched mocks.
- A p-value near the finite mock floor means the window needs more mocks, not immediate promotion.
- These mocks preserve the observed redshift distribution inside each window and the random-catalog sky footprint, but they are still not a complete Quaia selection-function likelihood.

## Outputs

- per-mock summary: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\quaia_redshift_window_highstat_mocks_20260714\quaia_redshift_window_mock_summary.csv`
- observed-vs-mock p-values: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\quaia_redshift_window_highstat_mocks_20260714\quaia_redshift_window_mock_pvalues.csv`
- manifest: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\quaia_redshift_window_highstat_mocks_20260714\quaia_redshift_window_mock_manifest.csv`
- config: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\quaia_redshift_window_highstat_mocks_20260714\quaia_redshift_window_mock_config.json`
- saved mock catalogs directory: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\quaia_redshift_window_highstat_mocks_20260714\mock_catalogs`
