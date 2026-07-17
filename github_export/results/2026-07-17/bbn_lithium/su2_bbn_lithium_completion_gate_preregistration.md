# SU2 BBN Lithium Completion-Gate Registration

- **Locked:** `2026-07-17T10:30:37.6020806+10:00`
- **Scan state at lock:** `22,171 / 30,503` successful unique LINX points
- **Unresolved points at lock:** `187`
- **Final target:** exactly `30,503` successful unique points and zero unresolved point IDs

## Registration Status

This is a prospective completion gate for the uncomputed final portion of an existing scan. It is informed by the first 22,171 successful evaluations and is therefore not represented as a blind preregistration of the entire analysis. Its purpose is to prevent the final 8,332 points from changing the decision rule after their abundances are known.

## Effective Mapping

The primary SU2 test is an effective expansion-only proxy:

\[
\begin{aligned}
S                    &= H_{\rm BBN}/H_{\rm SBBN}, \\
\Delta N_{\rm eff}  &= \frac{43}{7}\left(S^2-1\right), \\
\eta_{\rm fac}      &= \eta_{\rm BBN}/\eta_{\rm CMB}, \\
\tau_{n,{\rm fac}}  &= \tau_n/\tau_{n,0}.
\end{aligned}
\]

This mapping does not claim that an SU2 Lagrangian has been propagated through the BBN network. It tests only whether a non-negative SU2-like early-radiation contribution can solve lithium through the expansion history while respecting the registered controls.

## Locked Gates

All primary rows must satisfy:

- `S >= 1`, equivalent to non-negative initial `Delta N_eff` in this proxy;
- `abs(eta_fac - 1) <= 2(0.04 / 6.12)`;
- `abs(tau_n_fac - 1) <= 2(0.6 / 879.4)`;
- no selected nuclear-rate pull: `sum(q_i^2) = 0`;
- D/H within two observational standard deviations;
- He-4 mass fraction within two observational standard deviations;
- Li-7/H within two observational standard deviations.

The expansion-only arm receives a candidate gate pass only if at least one row satisfies all seven requirements. A pass would identify a follow-up region, not establish evidence for SU2. Zero joint rows is a veto of the tested expansion-only mapping over the registered scan.

Two control arms are locked:

1. **Modest rate controls:** the same SU2-compatible rows with `sum(q_i^2) <= 9`.
2. **Full scanned rate controls:** the same SU2-compatible rows with the complete selected-rate range.

Rate-control improvements are not attributed to SU2 without a separate physical coupling that predicts the rate shifts.

## Completion Rule

The final report is labelled complete only when the point table contains 30,503 successful unique point IDs and no failed ID lacking a successful retry. Failed attempts remain in the raw table for auditability. The completion analysis must be rerun with the programme hash below without changing thresholds or scenarios.

## Programme Hashes

| File | SHA-256 |
|---|---|
| `analyze_bbn_lithium_linx_fr_network_2026-07-16.py` | `d553fab1dd69c56ff17781cabf1a80b4a0719bc4098fef184ba1ae3f2dcde772` |
| `resume_bbn_lithium_linx_key_chunks_2026-07-17.ps1` | `1e7002109cef683e733302c32cdcdcd21f23a7a3524e0f211b72a9ebaa5c0df8` |
| `analyze_su2_bbn_lithium_gate_2026-07-17.py` | `d82b38fac5c3da154f796e886e1bd54b32ce547ec4015dca9dca656338d88e8b` |

The locked interim audit JSON has SHA-256 `6549b1f44641aeb5f70e5950e3fb9423f6190fd338af1fc2d73cda58fc75cacc`.
