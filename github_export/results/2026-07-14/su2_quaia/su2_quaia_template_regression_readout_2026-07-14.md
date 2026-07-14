# SU2 / Quaia Template Regression Readout

Date: July 14, 2026

## Result

The first internal template-regression control completed for the surviving `z1p00_1p50` Quaia amplitude hook across `|b| > 10, 15, 20, 25, 30, 35`.

## Main Readout

Smooth Galactic quadrupole controls do not materially reduce the original dipole. Across the latitude cuts, the quadrupole-controlled amplitude stays close to the baseline, with amplitude ratios of about `1.01-1.03`, and the fitted direction remains close to the baseline direction.

Internal photometry/depth/proper-motion controls are much less benign. They rotate the fitted dipole direction from the original `l ~ 145-165 deg, b ~ -27 to -38 deg` region toward roughly `l ~ 276-296 deg, b ~ -12 to -24 deg`, and they suppress the stricter latitude-cut amplitudes:

| bcut_deg | photometry+quadrupole amplitude ratio | SNR change | readout |
| ---: | ---: | ---: | --- |
| 10 | 1.200 | +2.049 | amplitude grows, direction rotates strongly |
| 15 | 0.962 | +1.049 | amplitude stable, direction rotates strongly |
| 20 | 1.214 | +1.440 | amplitude grows, direction rotates strongly |
| 25 | 0.885 | +0.364 | amplitude modestly reduced, direction rotates strongly |
| 30 | 0.776 | -0.280 | amplitude reduced, direction rotates strongly |
| 35 | 0.459 | -1.517 | amplitude strongly reduced, direction rotates strongly |

## Claim Boundary

This does not promote the SU2/Quaia signal. The result says:

- The original `z1p00_1p50` dipole is not simply a smooth Galactic quadrupole artifact.
- The signal is sensitive to internal photometric/depth/proper-motion controls, especially at higher latitude cuts.
- Direction instability under these controls is a gating warning.
- The next selection-function audit should use external dust, depth, scanning-law, and ecliptic/Galactic templates before any physical interpretation.

The global look-elsewhere mock audit remains the first promotion gate.

