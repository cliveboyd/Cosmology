# SU2 / Quaia Direction-Amplitude Stability Readout

Date: July 14, 2026

## Test

This gate perturbed the surviving `z ~ 1.0-1.5` SU2/Quaia angular hook and refit the object-level redshift dipole across `|b| > 10, 15, 20, 25, 30, 35`.

The reference direction was the strongest original family row: `1.0 <= z < 1.5`, `|b| > 35`.

## Fixed Width 0.50 Result

| window | mean amplitude | max SNR | max separation from reference | readout |
| --- | ---: | ---: | ---: | --- |
| 0.90-1.40 | 0.0023115 | 3.6887 | 17.33 deg | stable direction |
| 0.95-1.45 | 0.0030202 | 4.7852 | 12.00 deg | strongest and stable |
| 1.00-1.50 | 0.0023450 | 3.7059 | 17.94 deg | baseline stable |
| 1.05-1.55 | 0.0013951 | 2.7582 | 61.34 deg | direction breaks |
| 1.10-1.60 | 0.0010442 | 2.1157 | 76.01 deg | direction breaks |

## Interpretation

The remaining SU2/Quaia hook is not a stable broad `1.0 < z < 1.5` angular mode. It is concentrated on the lower side of that window:

- shifting the window down to `0.95-1.45` strengthens the signal and keeps the direction stable;
- the baseline `1.00-1.50` window remains directionally coherent;
- shifting upward to `1.05-1.55` or `1.10-1.60` weakens the amplitude and rotates the direction substantially;
- narrow-window scans show multiple direction changes outside the `~1.15-1.35` subrange.

This fails as a broad symmetry claim. It remains a narrower diagnostic target around approximately `z ~ 1.15-1.35`, with the strongest stable broad control being `0.95-1.45`.

## Claim Boundary

- Do not describe this as a stable SU2 angular symmetry.
- If continuing, focus on `0.95-1.45` and `1.15-1.35` as targeted follow-up windows.
- The next gate should be mock stability percentiles for those targeted windows, followed by external dust/depth/scanning-law templates.
- Because the global look-elsewhere test already failed, this remains exploratory follow-up rather than promotion evidence.

