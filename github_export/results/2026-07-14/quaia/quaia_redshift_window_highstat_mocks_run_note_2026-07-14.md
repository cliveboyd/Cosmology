# Focused Quaia High-Stat Mock Run Note

Date: July 14, 2026

## Purpose

The completed 25-window object-level Quaia mock run found only marginal evidence in a small number of windows. This focused run spends statistics only on those borderline targets rather than rerunning every window.

## Focus Windows

- `z1p00_1p50`, `|b| > 10`: borderline amplitude p-value in the 200-mock run.
- `z1p00_1p50`, `|b| > 15`: borderline amplitude p-value in the 200-mock run.
- `z1p70_1p80`, `|b| > 10`: lowest CMB-projected p-value in the 200-mock run.

## Configuration

- mocks per focused window: 5000
- saved full mock catalogs per window: 1
- observed catalog: `C:\Users\clive\Documents\Cosmology\quaia_G20p0_basic_gal.csv`
- random catalog: `C:\Users\clive\Documents\Cosmology\quaia_G20p0_randoms_simple.fits`
- focus CSV: `C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs\quaia_redshift_window_highstat_focus_windows_2026-07-14.csv`
- output directory: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\quaia_redshift_window_highstat_mocks_20260714`

## Background Job

PID: `16532`

Tail the run:

```powershell
Get-Content -Wait -Tail 80 "C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs\quaia_redshift_window_highstat_mocks_stdout.log"
```

Check errors:

```powershell
Get-Content -Tail 80 "C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs\quaia_redshift_window_highstat_mocks_stderr.log"
```

## Interpretation Gate

If the high-stat p-values move upward, the Quaia anisotropy result remains a non-promoted diagnostic. If a p-value remains small, the next step is not immediate physical interpretation; it is a systematics/selection-function audit of the same window.
