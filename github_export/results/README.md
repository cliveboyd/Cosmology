# Results

Small result tables and Markdown readouts copied here are intended to be cited by research documents.

Recommended pattern:

- `results/YYYY-MM-DD/topic_name/summary.md`
- `results/YYYY-MM-DD/topic_name/*.csv`
- `results/YYYY-MM-DD/topic_name/manifest.csv`

Large posterior chains, object-level mock catalogs, logs, and generated figure files should stay in local run folders and be listed in a manifest instead.

For SkyMap, a small number of representative figures and extracted DOCX media are committed because the historical reports are plot-heavy; the full Office binaries and large FITS/NPZ products remain outside Git.

## July 16 Audit

- `2026-07-16/su2_quaia/` contains the spatial-block null, template-collinearity, and over-absorption/injection-calibration readouts.
- `2026-07-16/sparc/` contains robust split convergence, paired-posterior nuisance diagnostics, the preregistered 24-galaxy non-centred NUTS pilot under `nuts_pilot/`, and the separately preregistered 157-galaxy run under `nuts_full/`.
- `2026-07-16/rawmu/` contains the FR `c(z)=c0(1+z)` no-loss dimming sequence for raw supernova MU, the alpha=0 calibration-budget ladder and direct LCDM comparison, and the follow-up catalogue-combination, redshift-band, and calibration-pull gates under the July 14 calibration/redshift-frame controls.
- `2026-07-16/bbn_lithium/` contains the first BBN Lithium anomaly scoping gate for FR/no-expansion ideas, including eta, early-clock, Li-production, and stellar-depletion stress tests.
- `2026-07-16/audit/` contains the cross-stream model-upgrade reproducibility audit and research-standing priorities.

## July 17 Audit

- `2026-07-17/bbn_lithium/` contains the completed 30,503-point FR/LINX scan audit, immutable-catalogue score and coverage checks, final registered SU2 lithium gate/frontier readouts, archived run configuration, and a hash manifest for the local catalogue and committed analysis products.
- The detailed audit fixes the claim boundary at an FR-inspired standard-BBN clock/radiation proxy. It also records the priority requirements for an FR-native thermal history, a matched CMB likelihood, full-network and tighter-tolerance validation, deuterium theory-error marginalisation, and hierarchical lithium-depletion modelling.
- The same folder contains the completed, separately registered 62-point literature-defined SU(2)-CMB thermodynamic bridge. Its report, trajectory, eta scan, figure, configuration and audit remain explicitly separate from the repository's late-time SU2/SU2R fluid and from the FR proxy.
- The completed preregistered 92-evaluation FR/LINX validation matrix contains tighter numerical tolerances, key/small/full network comparisons, deuterium-rate theory responses, full-network omitted Li/Be-channel checks, interpolation checks and updated-neutron-lifetime states. All evaluations succeeded; numerical/network corrections are small, while the local D-rate theory uncertainty is comparable to the observational D/H error and must enter the next likelihood.

## July 18 Audit

- `2026-07-18/sparc_am/` contains the preregistered SPARC matter/antimatter-identifiability gate, outcome-blind matched controls, twelve-scenario conventional-systematics matrix, deterministic multi-start profiles, figures and hash manifest.
- The implemented rotation-curve likelihood is matter/antimatter non-identifiable. Three of the six development galaxies retain a model-relative PLAMB tension, but the matched-block direction is not coherent and no galaxy is classified as antimatter.
- `2026-07-18/sparc_am_single_start_invalidated/` preserves the initial single-start bundle and its invalidation notice. It is retained only as a transparent numerical-audit record and must not be used for interpretation.
- `2026-07-18/fr_charge_conjugation/` derives the matter/antimatter transformation from Peter R. Lamb's supplied Full Relativity source and verifies it with synthetic invariance checks. Complete charge conjugation preserves total clout, asymmetry magnitude and galaxy dynamics; it therefore rules out an intrinsic SPARC rotation-curve antimatter classifier.
- The only conditional non-invariance occurs when a host is swapped while a signed asymmetric environment is held fixed. That environment-relative branch remains untestable until an external signed-background estimator and the FR inertia response are specified in advance; the five reserved SPARC galaxies remain unopened.
