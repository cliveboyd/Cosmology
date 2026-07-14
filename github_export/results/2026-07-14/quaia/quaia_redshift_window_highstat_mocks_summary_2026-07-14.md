# Quaia Redshift-Window High-Stat Mock Summary

Date: July 14, 2026

## Run

This follow-up run regenerated object-level Quaia mock catalogs for the three selected redshift-window/bcut cases, using 5000 mocks per case. Each mock preserves the observed redshift distribution inside the window and draws random sky positions from the Quaia random catalog with the same Galactic latitude cut.

## Results

| window | bcut_deg | N | amplitude p | CMB-projected p | readout |
| --- | ---: | ---: | ---: | ---: | --- |
| z1p00_1p50 | 15 | 185739 | 0.04159 | 0.18916 | Marginal amplitude excess; no CMB-direction support |
| z1p00_1p50 | 10 | 194552 | 0.04379 | 0.19116 | Marginal amplitude excess; no CMB-direction support |
| z1p70_1p80 | 10 | 34300 | 0.18076 | 0.09618 | Does not survive high-stat mock control |

## Interpretation Gate

The high-stat mocks leave a narrow, window-level amplitude anomaly in the `z1p00_1p50` Quaia slice, but they do not support a physical CMB-aligned anisotropy interpretation. The correct claim boundary is therefore:

- Treat `z1p00_1p50` as a focused Quaia selection-function/systematics follow-up target.
- Do not promote the result as broad cosmological anisotropy.
- Do not treat the `z1p70_1p80` case as surviving the high-stat mock control.
- Require a documented Quaia selection-function or survey-systematics likelihood before interpreting the residual amplitude excess physically.

## Next Analysis Step

Run a targeted Quaia audit for the surviving `z1p00_1p50` amplitude cases:

- compare bcut 10 and bcut 15 stability;
- add sky-footprint and depth-modulation controls;
- test ecliptic, Galactic, scanning-law, and extinction-aligned templates;
- keep CMB-projected anisotropy as a control unless it becomes significant under the same mock protocol.

