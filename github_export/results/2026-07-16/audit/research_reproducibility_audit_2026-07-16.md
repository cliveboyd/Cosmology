# Research Reproducibility Audit - 16 July 2026

## Executive readout

Switching the Codex reasoning model to GPT-5.6 Sol changes only subsequent assistant reasoning, code generation, and prose. It does not alter the existing task history, completed Markov chains, random seeds, data files, logs, Word documents, or Git commits. The model change was therefore treated as a provenance event and an opportunity for an independent audit, not as a numerical intervention.

The recent raw-MU supernova policy remains sound. The principal corrections concern the strength of the Quaia angular-mode readout and the labelling and convergence of the SPARC posterior analysis.

## Material findings

### 1. SU(2)/Quaia spatial null

The previous 1,000-draw object-level shuffle reached its Monte Carlo floor and did not preserve coarse angular correlation or covariance between nested latitude cuts. A new 50,000-draw common-sign spatial wild bootstrap was run for the targeted z = 1.00-1.50 family over |b| cuts of 10-35 degrees.

| Equal-area cell scale | Cells | Observed family max SNR | Null p99 | Targeted-family p |
|---:|---:|---:|---:|---:|
| 8 degrees | 910 | 5.7536 | 5.9165 | 0.01352 |
| 12 degrees | 424 | 5.7536 | 7.5140 | 0.09230 |
| 16 degrees | 231 | 5.7536 | 9.1402 | 0.23260 |

The candidate is weakly unusual only at the finest tested block scale. It is not spatial-null-scale robust and cannot be promoted as a discovery or headline anisotropy result. The historical all-window look-elsewhere result, p_global = 0.432713, remains a separate test and remains ordinary.

At the locked |b|>15 degree direction, projection on seven scan-law and colour controls gives R2 = 0.16363, VIF = 1.19565, and maximum single-control correlation r = -0.40438 with the Gaia scan-angle sin(2 theta) term. The full-design condition number is 2.735. This is moderate template alignment rather than catastrophic numerical collinearity. The CMB direction is more strongly aligned with these controls, R2 = 0.41097.

### 2. SU(2)/Quaia injection calibration

The matched injection study recovers approximately half of a one-times injected angular mode on average, with mean recovery 0.5012. The strict controlled physical fit retains 0.28335, while the final Delta BIC < -10 gate passes only 0.0011 of one-times injections.

The BIC rule remains suitable as a stringent promotion gate, but it has poor power near the observed signal scale. Failure of that gate cannot be interpreted as a stand-alone physical rejection test.

### 3. SPARC convergence and nuisance pairing

The previous analyser used classic unsplit Rhat and generated the reported stress scores after drawing galaxy nuisance parameters from their priors. Those values were therefore prior-redrawn in-sample fit checks, not posterior-predictive scores.

The replacement audit computes rank-normalised split Rhat, approximate bulk effective sample size, Monte Carlo standard error, and paired global-plus-galaxy posterior fit densities.

| Model | Max classic Rhat | Max rank-split Rhat | Minimum bulk ESS | Gate |
|---|---:|---:|---:|---|
| RAR | 1.03071 | 1.06163 | 76.27 | Improved; approximately <=1.05 target not passed |
| PLAMB kappa,p | 1.05424 | 1.30819 | 8.11 | Material non-convergence |

The paired posterior fit-density deltas versus RAR are +13.602 for all-Q2, +25.082 with bulges removed, +4.828 for gas-dominated systems, +11.927 at high inclination, +0.890 for low-acceleration outer points, and +1.644 at low inclination. These sign reversals show that nuisance handling mattered, but they are not promotion evidence because PLAMB is not converged and the same galaxies are reused for fitting and scoring.

Corrected claim boundary: **subset wins, not a full-sample win**. Robust convergence and held-out galaxy prediction remain open gates.

### 4. Raw-MU supernova audit

No material defect was found in the completed hierarchical Gaussian random-effect likelihood. The default remains:

- zHD as the preregistered primary frame for Pantheon+SH0ES and DES-Dovekie;
- zHEL and zCMB as controls unless a separate physical reason is specified before analysis;
- explicit dataset and survey calibration offsets;
- a transparent external-budget sensitivity ladder, with 25 mmag survey and 50 mmag dataset scales as the reference readout.

No numerical rerun was justified during this audit.

### 5. HHT and SkyMap archive

The HHT claim remains correctly bounded as a candidate screening workflow: existing surrogate probabilities near 0.97 and 0.96 do not support a confirmed CMB-supernova coupling. Two unterminated plotting-label strings in `plot_imf2_vs_bao_v3.py` were repaired; the change does not alter an HHT result.

The SkyMap archive remains a legacy/exploratory stream. One Python regular-expression string was changed to a raw string so the archive compiles without an invalid-escape warning. No SkyMap scientific result changed.

## Reproducibility repairs

- Added the spatial-block and template-collinearity SU(2)/Quaia analysis programme.
- Added the robust SPARC posterior diagnostic audit programme.
- Added three missing SPARC helper modules required by the archived MAP and sampler programmes.
- Corrected archived SPARC root-path resolution for execution from `github_export/code/sparc`.
- Repaired two HHT plotting syntax errors and one SkyMap regular-expression warning.
- Produced dated Word and PDF papers with corrected claims and retained the earlier versions as historical records.
- Recorded UK English and vertically aligned equals signs as standing document conventions.

## Research-standing priorities

1. **Independent preregistered SU(2)/Quaia validation.** Fix the catalogue, mask, redshift family, angular statistic, spatial-null scale, and promotion threshold before opening the result. A second quasar/LSS catalogue should reproduce the direction and amplitude class.
2. **Reparameterised SPARC sampling.** More sweeps of the same geometry are unlikely to repair a rank-split Rhat of 1.308. Use a non-centred mass-to-light hierarchy and a proven HMC/NUTS implementation where feasible; require rank-split Rhat <= 1.01-1.05 and adequate bulk and tail ESS.
3. **Out-of-sample SPARC comparison.** After convergence, run preregistered galaxy-held-out or grouped cross-validation against RAR, including all-Q2, low-acceleration, inclination, gas-dominated, and bulge controls.
4. **Clone-level reproducibility.** Add a locked environment, dataset acquisition/licence manifest, checksums, one-command runners, and CI smoke tests. Separate small citation artefacts from large local chains and object-level mocks.
5. **Confirmatory discipline.** Keep exploratory and confirmatory directories distinct, archive seeds/configurations/hashes, and preregister claim gates before new catalogue or survey releases are inspected.
6. **HHT replication.** Use probe-separated surrogates, a global multiple-testing correction, and an independent supernova/CMB residual construction before discussing cross-probe coupling.

## Validation completed

- New SU(2)/Quaia and SPARC audit programmes completed and wrote configuration, tables, reports, and manifests.
- The complete `github_export/code` tree compiles under the bundled Python runtime.
- Git whitespace/error checks pass.
- SU(2)/Quaia paper: 22 Letter pages, native Word PDF export, all pages visually reviewed.
- PLAMB paper: 33 Letter pages, native Word PDF export, all pages visually reviewed; blank-page defect removed.
- DOCX structural audits report one portrait Letter section in each paper and no high-severity accessibility findings after adding PLAMB figure descriptions.

## Publication boundary

The audit increases the standing of the project by making the negative and ambiguous results more reproducible and more accurately labelled. It does not strengthen the case for a detected SU(2) angular symmetry or a full-sample PLAMB victory. The strongest current empirical asset remains the calibration-aware, preregistered raw-MU supernova workflow; the most valuable next investment is independent validation rather than another internally tuned claim.
