# SU2 / Quaia Global Look-Elsewhere Mock Audit

Date: July 14, 2026

## Purpose

This run tests whether the strongest SU2/Quaia redshift-window targets remain unusual after accounting for the fact that the scan searched many redshift windows and Galactic latitude cuts.

For each mock, every non-full redshift-window row from the original SU2/Quaia scan is regenerated using the observed window redshift distribution and random Quaia sky positions under the same latitude cut. The recorded null statistic is the maximum over the scan.

## Configuration

- scan CSV: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_scan\su2_quaia_scan.csv`
- Quaia table: `C:\Users\clive\Documents\Cosmology\quaia_G20p0_basic_gal.csv`
- random sky catalog: `C:\Users\clive\Documents\Cosmology\quaia_G20p0_randoms_simple.fits`
- n_mocks: `300`
- seed: `280714`
- min_count: `50`

## Mock Maxima

- max formal SNR mean / p95 / p99: `3.65876` / `4.22811` / `4.61164`
- max SU2-priority mean / p95 / p99: `0.0326638` / `0.0381545` / `0.0415383`
- max CMB-projected component mean / p95 / p99: `0.00368273` / `0.00535212` / `0.0065612`

## Observed Thresholds

| role | tag | bcut_deg | N | observed_formal_snr | observed_priority | global_p_any_window_snr_ge_observed | global_p_any_window_priority_ge_observed | z1p00_family_p_any_bcut_snr_ge_observed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| observed_max_snr+observed_max_priority | z1p00_1p50 | 35 | 119719 | 3.70594 | 0.0355781 | 0.418605 | 0.166113 | 0.0166113 |
| target | z1p00_1p50 | 10 | 194552 | 3.47955 | 0.0334047 | 0.687708 | 0.365449 | 0.0564784 |
| target | z1p00_1p50 | 15 | 185739 | 3.59778 | 0.0345398 | 0.54485 | 0.239203 | 0.0332226 |
| target | z1p00_1p50 | 25 | 155609 | 3.49857 | 0.0335872 | 0.674419 | 0.33887 | 0.0564784 |

## Interpretation Gate

- If the global p-values are ordinary, the SU2/Quaia result remains a targeting diagnostic, not evidence for anisotropy.
- If a z1p00_1p50 target stays small under the global scan correction, the next gate is a selection-function/template regression.
- CMB-projected components remain controls unless their global p-values become small under the same protocol.

## Outputs

- per-mock maxima: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_global_lookelsewhere_20260714\su2_quaia_global_lookelsewhere_mock_maxima.csv`
- observed thresholds: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_global_lookelsewhere_20260714\su2_quaia_global_lookelsewhere_observed_thresholds.csv`
- config: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_global_lookelsewhere_20260714\su2_quaia_global_lookelsewhere_config.json`
