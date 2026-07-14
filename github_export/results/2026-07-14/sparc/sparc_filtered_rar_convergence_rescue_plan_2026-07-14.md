# SPARC Filtered-Family RAR Convergence-Rescue Plan

Date: July 14, 2026

## Purpose

This is a convergence-rescue run for the filtered all-Q2 SPARC RAR posterior. It does not by itself promote PLAMB. Its job is to reduce the RAR max Rhat from the completed full posterior value of 1.300 toward the promotion gate of about 1.05 or better.

## Claim Boundary Locked In

The filtered-family SPARC result should be described as subset wins plus a persistent p near one half, not as a full-sample win over RAR. Current gating failures are:

- all-Q2 baseline: delta log score vs RAR = -368.314
- low-acceleration outer points: delta log score vs RAR = -173.603
- high inclination: delta log score vs RAR = -612.830
- low inclination: delta log score vs RAR = -294.276
- RAR convergence: max Rhat = 1.300

Promotion requires a follow-up run that improves RAR convergence and makes the all-Q2, low-acceleration, and inclination stress scores nonnegative, or else supplies a preregistered physical exclusion rule before looking at the outcome.

## Prepared Background Run

Launch from PowerShell:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File "C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs\start_sparc_filtered_rar_convergence_rescue_background.ps1"
```

Tail the log after launch:

```powershell
Get-Content -Wait -Tail 80 "C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs\sparc_filtered_rar_convergence_rescue_stdout.log"
```

## Sampler Settings

- input sample: filtered persistent-negative-family removal used in the completed full posterior
- split: all_Q2
- model: RAR only
- chains: 5
- sweeps per chain: 96000
- burn-in: 36000
- thin: 20
- saved draws per chain: 3000
- galaxy update fraction: 1.0
- global log M/L step: 0.006
- log10 acceleration-scale step: 0.004
- log distance step: 0.018
- log stellar M/L step: 0.035
- start jitter: 0.015
- seed: 270714

## Outputs

Primary RAR-only run:

```text
C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\sparc_hierarchical_posterior\posterior_minus_persistent_negative_rar_rescue_20260714
```

Primary RAR convergence report after completion:

```text
C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\sparc_hierarchical_posterior\posterior_minus_persistent_negative_rar_rescue_20260714\analysis_rar_convergence_rescue_20260714\posterior_run_analysis_report.md
```

If the prior PLAMB case is present, the pipeline also builds a combined directory using the rescued RAR chains plus the existing PLAMB chains, then runs posterior-predictive stress scoring:

```text
C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\sparc_hierarchical_posterior\posterior_minus_persistent_negative_rar_rescue_with_existing_plamb_20260714
```

## Pass / Fail Readout

Pass for convergence rescue:

- RAR max global Rhat near 1.05 or better.

Pass for any later full-sample promotion:

- RAR convergence passes.
- all-Q2 posterior-predictive delta is nonnegative.
- low-acceleration outer-point delta is nonnegative.
- high- and low-inclination deltas are nonnegative, or an exclusion rule for inclination is preregistered before viewing the outcome.

Fail:

- RAR max Rhat remains high, or any of the negative all-Q2 / low-acceleration / inclination stress deltas remain negative.

## Completed Run Readout

Completed: July 14, 2026 at 21:57 local time.

Convergence rescue passed:

- RAR max Rhat improved from 1.300 to 1.031.
- RAR global acceptance improved to about 0.228.
- RAR galaxy acceptance remained stable at about 0.681.

Promotion gate remains closed:

- all-Q2 aggregate delta is now positive: +28.941, but the galaxy split remains unfavorable at 65 better / 92 worse.
- low-acceleration outer-point aggregate delta is now positive: +48.508.
- high-inclination stress remains strongly negative: -358.865.
- gas-dominated stress is negative: -28.175.
- low-inclination aggregate delta is positive, +21.586, but the galaxy split remains unfavorable at 33 better / 41 worse.

Updated claim boundary: RAR convergence has been rescued, and PLAMB now has aggregate all-Q2 and low-acceleration support in this filtered-family comparison, but the high-inclination and gas-dominated failures keep SPARC in the subset-win category rather than a full-sample PLAMB win.
