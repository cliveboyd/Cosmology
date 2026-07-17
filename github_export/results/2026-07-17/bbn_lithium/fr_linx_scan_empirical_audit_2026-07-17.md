# FR LINX Lithium Scan: Empirical Audit

Generated: 2026-07-17T04:40:04Z

## Catalogue Integrity

- Successful unique rows read: `30503`; unresolved failures: `0`.
- Reconstructed registered target: `30503` points; missing IDs: `0`; unexpected IDs: `0`.
- Semantic point-ID collision groups: `0`; random rounding collision groups: `0`.
- Exact duplicate rows removed after key-network normalisation: `39`.

## Score Integrity

- Maximum absolute reconstructed total-score discrepancy: `2.22e-12`.
- Maximum absolute reconstructed lower-bound-score discrepancy: `2.16e-12`.
- These are floating-point reconstruction checks, not a validation of the adopted likelihood.

## Proposal Diagnostics

- Nuclear pulls were proposed from `N(0, 1.35^2)` but scored against `N(0, 1.0^2)`.
- Pulls with `|q| <= 0.25` were replaced by exactly zero.
- Random rows with all three available Li/Be pulls exactly zero: `91` / `30000`.
- Such random q=0 rows inside the registered SU2 base cuts: `5`.
- The q=0 expansion-only scenario therefore contains proposal-generated point mass as well as the deliberate fixed-rate grid.

## Current Best Penalised Row

- Total measurement objective: `215.256`.
- Predictions: D/H x1e5 `2.58506`, Yp `0.248691`, Li/H x1e10 `4.69411`.
- Controls: eta factor `0.975073`, S `1.04233`, initial Delta Neff `0.531119`, q-pull penalty `22.8451`.

## Interpretation Boundary

The catalogue arithmetic and ID coverage pass the checks above when the final target is complete. The scan remains a proposal-driven, penalised sensitivity study. It is not posterior sampling, evidence calculation, a goodness-of-fit chi-squared with a fixed number of degrees of freedom, or a self-consistent FR/fundamental-constant BBN model.
