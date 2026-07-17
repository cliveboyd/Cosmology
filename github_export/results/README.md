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
