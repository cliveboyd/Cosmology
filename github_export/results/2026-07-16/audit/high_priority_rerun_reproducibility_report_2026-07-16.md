# High-priority programme rerun and reproducibility report

**Date:** 16 July 2026  
**Runtime:** `C:\Users\clive\anaconda3x\python.exe`  
**Status:** Complete

## Scope

The recent Model 5.6 Sol audit identified two completed analyses that warranted immediate clean reruns:

1. the SU2/Quaia 50,000-mock spatial-block null and template-collinearity audit; and
2. the corrected SPARC rank-normalised convergence and paired-posterior fit-density audit.

The earlier SPARC posterior sampler was not repeated unchanged. Its principal remaining problem is poor hierarchical mixing, especially for the PLAMB mass-to-light parameters. Additional sweeps under the same parameterisation are not expected to resolve that problem; the appropriate next sampler is a non-centred hierarchy using HMC/NUTS or an equivalently efficient gradient-based method.

## SU2/Quaia rerun

Programme:

`github_export/code/quaia_su2/analyze_su2_quaia_spatial_block_collinearity_2026-07-16.py`

Configuration:

- redshift interval: `1.0 <= z <= 1.5`;
- Galactic latitude cuts: `|b| > 10, 15, 20, 25, 30, 35` degrees;
- equal-area block scales: `8, 12, 16` degrees;
- mocks per block scale: `50,000`;
- seed: `200716`; and
- common Rademacher block signs across the nested latitude-cut family.

The first clean rerun reproduced every substantive archived output byte for byte. It also emitted a NumPy deprecation warning for `np.row_stack`. The programme was updated to the supported alias `np.vstack` and the complete 150,000-mock calculation was repeated. The patched rerun had an empty standard-error log and again reproduced all seven substantive outputs byte for byte.

| Block scale | Observed family maximum SNR | Family p-value |
| ---: | ---: | ---: |
| 8 degrees | 5.7535877126 | 0.0135197296 |
| 12 degrees | 5.7535877126 | 0.0922981540 |
| 16 degrees | 5.7535877126 | 0.2325953481 |

**Interpretation:** the small p-value at the finest block scale is not stable to plausible enlargement of the spatial correlation scale. The result therefore remains a scale-sensitive follow-up signal, not a robust global anisotropy detection.

Fresh result directories:

- `plamb_runs/diagnostics/su2_quaia_spatial_block_collinearity_repro_20260716`
- `plamb_runs/diagnostics/su2_quaia_spatial_block_collinearity_repro_patched_20260716`

## SPARC rerun

Programme:

`github_export/code/sparc/audit_sparc_posterior_diagnostics_2026-07-16.py`

Configuration:

- combined RAR-rescue and existing-PLAMB posterior directory;
- rank-normalised split R-hat and approximate rank-normalised bulk effective sample size;
- at most 4,500 paired stored global and nuisance posterior draws; and
- seed: `270716`.

All six substantive outputs reproduced byte for byte and the standard-error log was empty.

| Model | Maximum rank-split R-hat | Minimum bulk ESS |
| --- | ---: | ---: |
| RAR | 1.0616317502 | 76.2689 |
| PLAMB optical-depth | 1.3081868983 | 8.11474 |

Paired in-sample log-density differences, PLAMB minus RAR:

| Diagnostic subset | Difference |
| --- | ---: |
| Baseline all-Q2 | +13.6018 |
| Bulge removed | +25.0818 |
| Gas dominated | +4.82809 |
| High inclination | +11.9272 |
| Low-acceleration outer points | +0.889980 |
| Low inclination | +1.64359 |

**Interpretation:** the corrected score calculation is reproducible, but it does not remove the convergence gate. The PLAMB posterior remains too poorly mixed for promotion as a full-sample result. The low-acceleration and inclination differences also remain small relative to the stronger subset results. The paper's claim boundary should remain: **subset wins, not a full-sample win**.

Fresh result directory:

`plamb_runs/diagnostics/sparc_posterior_diagnostic_audit_repro_20260716`

## Decision

No scientific conclusion changed. The reruns strengthen provenance by showing exact deterministic reproduction from the Git-tracked programmes. The highest-priority new SPARC computation is a reparameterised non-centred posterior sampler, not a repetition of the existing random-walk sampler. For SU2/Quaia, the stronger next test remains cross-catalogue validation and a spatial null whose correlation scale is estimated or marginalised rather than selected after inspection.
