# SPARC Filtered-Family RAR Rescue Readout

Date: July 14, 2026

## Bottom Line

The RAR convergence-rescue run succeeded, but it does not promote the filtered-family SPARC result to a full-sample PLAMB win.

The claim boundary remains: subset wins and an interesting PLAMB posterior shape, not a full-sample victory over RAR.

## Convergence Result

The RAR-only rescue run used the same filtered all-Q2 sample as the completed full posterior.

| Split | Model | Chains | Saved / chain | Max Rhat | Global accept | Galaxy accept |
|---|---|---:|---:|---:|---:|---:|
| all_Q2 | RAR | 5 | 3000 | 1.031 | 0.228 | 0.681 |

This clears the pre-run convergence target of about max Rhat <= 1.05.

## Combined Predictive Stress Readout

The combined analysis used the rescued RAR chains and the existing PLAMB chains.

| Stress case | PLAMB delta log score vs RAR | Better / worse galaxies | Gate readout |
|---|---:|---:|---|
| all-Q2 baseline | 28.941 | 65 / 92 | Aggregate positive, galaxy-count split unfavorable |
| bulge removed | 6.459 | 63 / 75 | Weak aggregate positive |
| gas dominated | -28.175 | 17 / 14 | Failure |
| high inclination | -358.865 | 49 / 34 | Major failure |
| high-quality Q1 | 379.463 | 48 / 47 | Strong subset win |
| low-acceleration outer points | 48.508 | 80 / 75 | Reversed to positive |
| low inclination | 21.586 | 33 / 41 | Aggregate positive, galaxy-count split unfavorable |

## Claim Boundary

What improved:

- RAR convergence is rescued: max Rhat improved from 1.300 in the earlier full posterior run to 1.031.
- The all-Q2 aggregate delta changed from negative to modestly positive.
- The low-acceleration outer-point aggregate delta changed from negative to positive.
- The high-quality Q1 subset remains a strong PLAMB-favorable subset.

What still blocks a full-sample claim:

- High-inclination stress remains strongly negative at delta log score = -358.865.
- Gas-dominated stress is negative at delta log score = -28.175.
- Baseline all-Q2 has more galaxies worse than RAR than better than RAR, 65 / 92, even though the aggregate score is positive.
- Low-inclination has the same asymmetry, 33 / 41, despite positive aggregate score.
- PLAMB posterior checks still do not include p = 0.5 or kappa = 1 inside the 68% interval: p = 0.5193 [0.51437, 0.52374], kappa = 0.78976 [0.76858, 0.81729].

## Operational Readout

The convergence rescue did its job. The next scientific question is no longer "is RAR convergence contaminating the comparison?" It is "why does PLAMB still fail badly in high-inclination and gas-dominated controls?"

Do not headline this as a full-sample PLAMB win. The updated headline should be:

> RAR convergence rescued; PLAMB gains aggregate all-Q2 and low-acceleration support, but high-inclination and gas-dominated failures keep the SPARC result in the subset-win category.

## Source Outputs

RAR-only run:

`C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\sparc_hierarchical_posterior\posterior_minus_persistent_negative_rar_rescue_20260714`

Combined rescued-RAR plus existing-PLAMB analysis:

`C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\sparc_hierarchical_posterior\posterior_minus_persistent_negative_rar_rescue_with_existing_plamb_20260714\analysis_rescued_rar_existing_plamb_nd16_20260714`
