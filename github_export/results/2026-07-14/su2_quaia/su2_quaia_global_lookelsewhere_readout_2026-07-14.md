# SU2 / Quaia Global Look-Elsewhere Readout

Date: July 14, 2026

## Test

The global look-elsewhere test reran the SU2/Quaia redshift-window scan inside 300 object-level mocks. Each mock preserved the observed redshift distribution in each scan window and drew random Quaia sky positions under the same Galactic latitude cut. The null statistic was the maximum over all 168 non-full redshift-window scan rows.

## Result

The strongest observed SU2/Quaia row is `z1p00_1p50`, `|b| > 35`, with formal SNR `3.70594` and SU2-priority score `0.0355781`.

Against the global scan maxima:

| observed target | global p, max SNR | global p, max SU2-priority | z1p00-family p, max SNR |
| --- | ---: | ---: | ---: |
| z1p00_1p50, |b| > 35 | 0.418605 | 0.166113 | 0.0166113 |
| z1p00_1p50, |b| > 10 | 0.687708 | 0.365449 | 0.0564784 |
| z1p00_1p50, |b| > 15 | 0.544850 | 0.239203 | 0.0332226 |
| z1p00_1p50, |b| > 25 | 0.674419 | 0.338870 | 0.0564784 |

Mock maxima summary:

- max formal SNR mean / p95 / p99: `3.65876` / `4.22811` / `4.61164`
- max SU2-priority mean / p95 / p99: `0.0326638` / `0.0381545` / `0.0415383`
- max CMB-projected component mean / p95 / p99: `0.00368273` / `0.00535212` / `0.0065612`

## Interpretation

Gate 1 does not pass as a full-scan global anisotropy detection. The full look-elsewhere corrected p-values are ordinary: `0.42` for max SNR and `0.17` for max SU2-priority at the strongest observed row.

The one interesting residue is narrower: inside the pre-identified `z1p00_1p50` family, the `|b| > 35` row has a family-level SNR p-value of `0.0166`. That is a targeted follow-up hook, not a promoted SU2 symmetry claim, especially because the internal template regression already showed direction instability under photometry/depth/proper-motion controls.

## Claim Boundary

- Do not claim a stable SU2-associated angular mode from the current global scan.
- Keep `z1p00_1p50` as a focused diagnostic target.
- Move next to direction/amplitude stability and external selection-function templates only as follow-up controls, not as confirmatory promotion work.

