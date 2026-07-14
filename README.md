# Cosmology Research Archive

This repository is set up as a curated GitHub archive for the PLAMB/SU2/Quaia/SPARC cosmology analyses.

The working directory contains many raw datasets, large FITS files, posterior chains, Office drafts, and generated run folders. Those are intentionally ignored by Git. GitHub should track the curated material under `github_export/`:

- `github_export/code/` - reproducible analysis scripts copied from the working tree.
- `github_export/results/` - small result tables, Markdown readouts, and report-ready summaries.
- `github_export/manifests/` - dataset and result manifests with local/source references.
- `github_export/docs/` - notes for citing repository paths in research documents.
- `github_export/data/` - README/manifests only by default; raw datasets are referenced, not committed.

Large raw data should be referenced by source, checksum, DOI, or local path, and only placed in Git with Git LFS or an external archive workflow.

