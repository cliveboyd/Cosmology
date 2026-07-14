# SU2 / Quaia Template Regression Control

Date: July 14, 2026

## Purpose

This is a first-pass selection/systematics stress test for the surviving `z1p00_1p50` Quaia amplitude hook. It refits the redshift dipole after adding internal catalog-quality proxies and smooth Galactic quadrupole terms.

This is not yet a full external dust, scanning-law, or survey-depth likelihood. It is an immediate internal-control regression using object-level columns available in `quaia_G20.0.fits`.

## Controls

- photometry/common controls: redshift error, Gaia G, BP-RP, W1, W2, W1-W2, total proper motion, pmra, pmdec.
- quadrupole/common controls: six smooth Galactic quadrupole basis terms.
- common-sample rows compare against `baseline_common`, so amplitude changes are not just missing-data artifacts.

## Summary

| tag | bcut_deg | model | N | amp | amp_ratio_vs_baseline_common | formal_snr | snr_delta_vs_baseline_common | l_deg | b_deg |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| z1p00_1p50 | 10 | baseline_all | 194552 | 0.00209207 |  | 3.47955 |  | 159.018 | -32.1216 |
| z1p00_1p50 | 10 | baseline_common | 194552 | 0.00209207 | 1 | 3.47955 | 0 | 159.018 | -32.1216 |
| z1p00_1p50 | 10 | photometry_common | 194552 | 0.00236387 | 1.12992 | 5.25302 | 1.77347 | 295.823 | -19.7083 |
| z1p00_1p50 | 10 | quadrupole_common | 194552 | 0.00212861 | 1.01746 | 3.51122 | 0.0316734 | 161.642 | -30.8216 |
| z1p00_1p50 | 10 | photometry_plus_quadrupole_common | 194552 | 0.00251073 | 1.20011 | 5.52826 | 2.04871 | 294.284 | -17.1212 |
| z1p00_1p50 | 15 | baseline_all | 185739 | 0.00224973 |  | 3.59778 |  | 155.979 | -29.0193 |
| z1p00_1p50 | 15 | baseline_common | 185739 | 0.00224973 | 1 | 3.59778 | 0 | 155.979 | -29.0193 |
| z1p00_1p50 | 15 | photometry_common | 185739 | 0.00203881 | 0.906248 | 4.43812 | 0.840332 | 284.598 | -21.649 |
| z1p00_1p50 | 15 | quadrupole_common | 185739 | 0.00227283 | 1.01027 | 3.61603 | 0.018242 | 157.135 | -28.2848 |
| z1p00_1p50 | 15 | photometry_plus_quadrupole_common | 185739 | 0.00216453 | 0.96213 | 4.64687 | 1.04909 | 284.443 | -18.6315 |
| z1p00_1p50 | 20 | baseline_all | 172054 | 0.00184579 |  | 3.02083 |  | 151.066 | -38.0799 |
| z1p00_1p50 | 20 | baseline_common | 172054 | 0.00184579 | 1 | 3.02083 | 0 | 151.066 | -38.0799 |
| z1p00_1p50 | 20 | photometry_common | 172054 | 0.00213023 | 1.1541 | 4.29868 | 1.27786 | 287.911 | -18.2571 |
| z1p00_1p50 | 20 | quadrupole_common | 172054 | 0.00186943 | 1.01281 | 3.03969 | 0.0188674 | 152.221 | -37.2769 |
| z1p00_1p50 | 20 | photometry_plus_quadrupole_common | 172054 | 0.00224066 | 1.21393 | 4.46125 | 1.44042 | 287.942 | -15.8131 |
| z1p00_1p50 | 25 | baseline_all | 155609 | 0.0023678 |  | 3.49857 |  | 165.22 | -31.3107 |
| z1p00_1p50 | 25 | baseline_common | 155609 | 0.0023678 | 1 | 3.49857 | 0 | 165.22 | -31.3107 |
| z1p00_1p50 | 25 | photometry_common | 155609 | 0.00199143 | 0.84105 | 3.73815 | 0.239579 | 276.986 | -18.6875 |
| z1p00_1p50 | 25 | quadrupole_common | 155609 | 0.00243337 | 1.02769 | 3.56746 | 0.0688891 | 166 | -30.4501 |
| z1p00_1p50 | 25 | photometry_plus_quadrupole_common | 155609 | 0.00209439 | 0.884533 | 3.86235 | 0.363788 | 278.812 | -16.1618 |
| z1p00_1p50 | 30 | baseline_all | 137850 | 0.00246501 |  | 3.38376 |  | 152.225 | -30.474 |
| z1p00_1p50 | 30 | baseline_common | 137850 | 0.00246501 | 1 | 3.38376 | 0 | 152.225 | -30.474 |
| z1p00_1p50 | 30 | photometry_common | 137850 | 0.00178301 | 0.72333 | 2.95035 | -0.433411 | 284.846 | -14.7926 |
| z1p00_1p50 | 30 | quadrupole_common | 137850 | 0.00252912 | 1.02601 | 3.44829 | 0.0645255 | 153.869 | -29.8734 |
| z1p00_1p50 | 30 | photometry_plus_quadrupole_common | 137850 | 0.0019138 | 0.776385 | 3.10369 | -0.280078 | 286.05 | -12.2258 |
| z1p00_1p50 | 35 | baseline_all | 119719 | 0.00304966 |  | 3.70594 |  | 145.328 | -26.7809 |
| z1p00_1p50 | 35 | baseline_common | 119719 | 0.00304966 | 1 | 3.70594 | 0 | 145.328 | -26.7809 |
| z1p00_1p50 | 35 | photometry_common | 119719 | 0.00127997 | 0.419709 | 2.10667 | -1.59926 | 280.25 | -24.3177 |
| z1p00_1p50 | 35 | quadrupole_common | 119719 | 0.00308834 | 1.01268 | 3.7324 | 0.0264643 | 146.742 | -26.3787 |
| z1p00_1p50 | 35 | photometry_plus_quadrupole_common | 119719 | 0.0014002 | 0.459132 | 2.18915 | -1.51679 | 281.618 | -19.9528 |

## Interpretation Gate

- A large dipole-amplitude drop under photometry or quadrupole controls points toward selection/depth leakage.
- A stable amplitude under these internal controls still does not promote the result physically; it only justifies moving to external dust/scanning-law/depth templates.
- The global look-elsewhere mock audit remains the first promotion gate.

## Outputs

- summary CSV: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_template_regression_20260714\su2_quaia_template_regression_summary.csv`
- report: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_template_regression_20260714\su2_quaia_template_regression_report.md`
