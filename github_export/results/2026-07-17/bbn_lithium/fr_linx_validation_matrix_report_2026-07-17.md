# FR/LINX BBN Targeted Validation Matrix Report

Generated: 2026-07-17T07:31:47Z

## Completion

- Registered evaluations: `92`.
- Successful evaluations: `92`.
- Missing or unresolved evaluations: `0`.

## Strict Baseline

| Network | D/H x 1e5 | Yp | Li-7/H x 1e10 | Runtime (s) |
|---|---:|---:|---:|---:|
| full | 2.434914 | 0.2469325 | 5.637694 | 1.62e+03 |
| key | 2.434293 | 0.2469465 | 5.695731 | 223 |
| small | 2.435004 | 0.2469494 | 5.638004 | 698 |

## Maximum Absolute Differences

| Check | D/H x 1e5 | Yp | Li-7/H x 1e10 |
|---|---:|---:|---:|
| Strict key minus stored loose key | 0.0005507 | 6.532e-05 | 0.00171 |
| Sampling 600 minus 300 | 0.0001157 | 2.017e-05 | 0.0003191 |
| Small minus key | 0.0007384 | 3.743e-06 | 0.05775 |
| Full minus key | 0.0006485 | 1.483e-05 | 0.05804 |

## D-Rate Theory Response

Independent centred unit-pull responses were combined in quadrature. These are local response estimates, not a complete marginal likelihood.

| Network | State | sigma(D/H x 1e5) | sigma(Yp) | sigma(Li x 1e10) |
|---|---|---:|---:|---:|
| full | baseline | 0.02608 | 2.166e-05 | 0.08423 |
| full | best_registered_direct | 0.02704 | 2.14e-05 | 0.06955 |
| key | baseline | 0.02609 | 2.162e-05 | 0.08493 |
| key | best_registered_direct | 0.02706 | 2.138e-05 | 0.07012 |
| small | baseline | 0.02609 | 2.147e-05 | 0.08423 |
| small | best_registered_direct | 0.02706 | 2.129e-05 | 0.06955 |

## Gate and Claim Boundary

- Registered D+He+Li two-sigma passes among successful designed evaluations: `0`.
- Combined-current D+He+Li two-sigma passes: `0`.
- Largest absolute full-network omitted-channel Li response: `0.06079` in Li/H x 1e10 units.

This matrix validates numerical and selected-network stability and quantifies local nuclear-rate sensitivity. It remains a standard-background LINX calculation. It does not validate a non-expanding FR thermal history or identify nuclear pulls with FR or SU2 dynamics.
