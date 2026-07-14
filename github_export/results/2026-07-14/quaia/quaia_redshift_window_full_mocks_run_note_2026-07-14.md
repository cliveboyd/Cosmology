# Quaia Redshift-Window Full Mock Run Note

Date: July 14, 2026

## Purpose

This run closes the previous Quaia mock-control gap: the earlier control used summary-level full-sample mock dipole products and could not validate redshift-window anisotropy targets. The new run generates object-level mocks per redshift window before interpreting the anisotropy scan.

## Method

For each selected observed Quaia redshift window:

1. Select observed Quaia objects in the same redshift window and Galactic latitude cut.
2. Preserve the observed window redshift distribution.
3. Draw the same number of random sky positions from `quaia_G20p0_randoms_simple.fits` after the same Galactic latitude cut.
4. Place the observed redshift distribution on the random sky positions.
5. Fit the same redshift dipole model, `z = a0 + D.n`.
6. Compare the observed window amplitude, CMB-projected component, formal SNR, and CMB alignment against the window-matched mock distribution.

The run keeps per-mock fit summaries for all mocks and saves one compact full mock catalog per selected window as `.npz`.

## Configuration

- selected windows: top 25 SU2-Quaia priority redshift windows from `su2_quaia_scan_top.csv`
- mocks per window: 200
- saved full mock catalogs per window: 1
- random catalog: `C:\Users\clive\Documents\Cosmology\quaia_G20p0_randoms_simple.fits`
- observed catalog: `C:\Users\clive\Documents\Cosmology\quaia_G20p0_basic_gal.csv`
- output directory: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\quaia_redshift_window_full_mocks_20260714`

## Background Job

PID: `4224`

Tail the run:

```powershell
Get-Content -Wait -Tail 80 "C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs\quaia_redshift_window_full_mocks_stdout.log"
```

Check errors:

```powershell
Get-Content -Tail 80 "C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs\quaia_redshift_window_full_mocks_stderr.log"
```

## Expected Final Outputs

- `quaia_redshift_window_mock_summary.csv`: every mock fit for every selected window.
- `quaia_redshift_window_mock_pvalues.csv`: observed-vs-mock p-values per window.
- `quaia_redshift_window_mock_manifest.csv`: window definitions and saved catalog paths.
- `quaia_redshift_window_mock_report.md`: interpretation and smallest p-value tables.
- `mock_catalogs/*.npz`: compact full object-level mock catalogs for the first mock in each window.

## Interpretation Gate

Do not interpret a Quaia anisotropy window as physical unless it remains unusual against these redshift-window-matched object-level mocks. A p-value at the finite mock floor means more mocks are needed, not immediate promotion.
