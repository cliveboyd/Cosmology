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
- `2026-07-16/sparc/` contains robust split convergence and paired-posterior nuisance diagnostics.
- `2026-07-16/audit/` contains the cross-stream model-upgrade reproducibility audit and research-standing priorities.
