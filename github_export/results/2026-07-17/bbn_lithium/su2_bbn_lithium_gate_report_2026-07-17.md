# SU2-Compatible BBN Lithium Gate

Generated: 2026-07-17T14:28:55
Status: **complete**; successful unique LINX points `30503 / 30503`, attempted `30503`, unresolved `0`.

## Scope

This is an effective SU2 gate, not a first-principles SU2 nucleosynthesis calculation. A non-negative SU2-like early-radiation contribution is represented only through an expansion factor and its equivalent initial radiation shift:

\[
\begin{aligned}
S                    &= H_{\rm BBN}/H_{\rm SBBN}, \\
\Delta N_{\rm eff}  &= \frac{43}{7}\left(S^2-1\right), \\
\eta_{\rm fac}      &= \eta_{\rm BBN}/\eta_{\rm CMB}, \\
\tau_{n,{\rm fac}}  &= \tau_n/\tau_{n,0}.
\end{aligned}
\]

The SU2-compatible rows require `S >= 1`, two-sigma CMB consistency in `eta_fac`, and two-sigma neutron-lifetime consistency. D/H and He-4 must each pass their two-sigma abundance gates before lithium is assessed. Nuclear-rate pulls are controls and are not attributed to SU2 without an explicit coupling model.

## Scenario Gates

| Scenario | Rows | D+He pass | D+He+Li pass | Best chi2, Li measurement | Min Li given D+He |
|---|---:|---:|---:|---:|---:|
| su2_expansion_only | 25 | 16 | 0 | 273.8 | 5.555 |
| su2_plus_modest_rate_controls | 1773 | 1224 | 0 | 228 | 5.113 |
| su2_plus_scanned_rate_controls | 2211 | 1505 | 0 | 219 | 4.924 |
| all_scanned_controls | 30503 | 4715 | 0 | 215.3 | 4.741 |

## Readout

- Expansion-only proxy: best measurement chi-squared `273.8`; lithium gate pass `False`.
- With modest joint rate controls (`sum(q_i^2) <= 9`): best chi-squared `228`; lithium gate pass `False`.
- With the full scanned rate-control range: best chi-squared `219`; lithium gate pass `False`.

A lithium pass requires simultaneous two-sigma agreement with D/H, He-4 and the lithium plateau. A lower total chi-squared that still misses this joint gate is not a lithium solution.

## Claim Boundary

The scan can veto an expansion-only SU2 explanation within the tested effective-clock mapping. It cannot validate or exclude a fully specified SU2 model that changes binding energies, neutron-proton mass splitting, weak rates and nuclear Q-values self-consistently. The selected LINX reaction pulls remain survey-style nuisance controls for nuclear physics, not evidence for SU2.
