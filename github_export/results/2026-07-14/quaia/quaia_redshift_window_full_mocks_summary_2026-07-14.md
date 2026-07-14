# Quaia Window-Matched Mock Summary

Date: July 14, 2026

## Result

The 25-window object-level Quaia mock run completed cleanly. The earlier summary-only mock limitation is now closed for the selected redshift-window anisotropy targets: each observed window was tested against mocks that preserve the window redshift distribution and use random Quaia sky positions under the same Galactic latitude cut.

The strongest amplitude cases are marginal rather than decisive. The two smallest amplitude p-values are both finite-mock-borderline results at p = 0.04975 from 200 mocks:

- `z1p00_1p50`, `|b| > 10`
- `z1p00_1p50`, `|b| > 15`

The strongest CMB-projected component is not below 0.05; the smallest is `z1p70_1p80`, `|b| > 10`, with p = 0.08955. This keeps Quaia as a diagnostic target, not a promoted physical anisotropy claim.

## Smallest Amplitude p-values

| idx | bcut | window     | N      | obs amp    | mock p95   | mock p99  | p amp    | p CMB proj |
| --- | ---- | ---------- | ------ | ---------- | ---------- | --------- | -------- | ---------- |
| 2   | 15   | z1p00_1p50 | 185739 | 0.0022497  | 0.0021897  | 0.0025228 | 0.049751 | 0.22886    |
| 4   | 10   | z1p00_1p50 | 194552 | 0.0020921  | 0.0020495  | 0.0025709 | 0.049751 | 0.20398    |
| 3   | 25   | z1p00_1p50 | 155609 | 0.0023678  | 0.0028484  | 0.0032506 | 0.13433  | 0.32836    |
| 1   | 35   | z1p00_1p50 | 119719 | 0.0030497  | 0.0037908  | 0.0046267 | 0.14428  | 0.22886    |
| 25  | 10   | z0p60_0p70 | 23898  | 0.0009571  | 0.001167   | 0.0013455 | 0.14428  | 0.23383    |
| 5   | 30   | z1p00_1p50 | 137850 | 0.002465   | 0.0029292  | 0.0042462 | 0.14428  | 0.26866    |
| 9   | 10   | z0p90_1p00 | 36372  | 0.00071159 | 0.00087149 | 0.0010059 | 0.16915  | 0.34328    |
| 23  | 20   | z0p80_0p90 | 28027  | 0.00096454 | 0.0012214  | 0.0014681 | 0.1791   | 0.26368    |
| 20  | 15   | z1p50_2p50 | 264161 | 0.0029781  | 0.0037278  | 0.0047545 | 0.1791   | 0.70149    |
| 7   | 20   | z1p00_1p50 | 172054 | 0.0018458  | 0.002423   | 0.0027323 | 0.18408  | 0.22886    |

## Smallest CMB-Projected p-values

| idx | bcut | window     | N      | obs |CMB|  | mock p95   | mock p99  | p CMB proj | p amp    |
| --- | ---- | ---------- | ------ | ---------- | ---------- | --------- | ---------- | -------- |
| 16  | 10   | z1p70_1p80 | 34300  | 0.00068264 | 0.00073914 | 0.001005  | 0.089552   | 0.18408  |
| 24  | 20   | z1p70_1p80 | 30496  | 0.00072361 | 0.00090642 | 0.0012732 | 0.12438    | 0.21393  |
| 18  | 15   | z1p70_1p80 | 32830  | 0.00070264 | 0.00088003 | 0.0011546 | 0.12935    | 0.23881  |
| 4   | 10   | z1p00_1p50 | 194552 | 0.0011353  | 0.0016839  | 0.0021775 | 0.20398    | 0.049751 |
| 1   | 35   | z1p00_1p50 | 119719 | 0.0018956  | 0.0032021  | 0.0037465 | 0.22886    | 0.14428  |
| 2   | 15   | z1p00_1p50 | 185739 | 0.00122    | 0.0017706  | 0.0023485 | 0.22886    | 0.049751 |
| 7   | 20   | z1p00_1p50 | 172054 | 0.0012267  | 0.0018774  | 0.0026048 | 0.22886    | 0.18408  |
| 25  | 10   | z0p60_0p70 | 23898  | 0.0006148  | 0.00088483 | 0.0011485 | 0.23383    | 0.14428  |
| 23  | 20   | z0p80_0p90 | 28027  | 0.00058144 | 0.00097336 | 0.001281  | 0.26368    | 0.1791   |
| 5   | 30   | z1p00_1p50 | 137850 | 0.0014579  | 0.0023015  | 0.0036342 | 0.26866    | 0.14428  |

## Focused High-Stat Rerun

Prepared focus windows:

| bcut | window     | N      | obs amp    | p amp 200 | obs |CMB|  | p CMB 200 |
| ---- | ---------- | ------ | ---------- | --------- | ---------- | --------- |
| 10   | z1p00_1p50 | 194552 | 0.0020921  | 0.049751  | 0.0011353  | 0.20398   |
| 15   | z1p00_1p50 | 185739 | 0.0022497  | 0.049751  | 0.00122    | 0.22886   |
| 10   | z1p70_1p80 | 34300  | 0.00075108 | 0.18408   | 0.00068264 | 0.089552  |

The follow-up should not rerun all 25 windows. It should spend statistics only on these three windows, using 5000 mocks each. Promotion still requires these high-stat p-values to remain small and a separate physical/systematics argument; otherwise the result remains a window-level diagnostic.

## Copied Outputs

- `quaia_redshift_window_full_mocks_pvalues_2026-07-14.csv`
- `quaia_redshift_window_full_mocks_report_2026-07-14.md`
- `quaia_redshift_window_full_mocks_pvalues_2026-07-14.png`
- `quaia_redshift_window_full_mocks_manifest_2026-07-14.csv`
- `quaia_redshift_window_highstat_focus_windows_2026-07-14.csv`
