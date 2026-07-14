# GitHub Export Layout

This folder is the GitHub-ready subset of the local cosmology workspace.

## Layout

- `code/` contains copied Python programs and launch scripts needed to reproduce the main diagnostics.
- `results/` contains small CSV/Markdown outputs suitable for direct citation in reports.
- `manifests/` records where raw data and larger outputs live.
- `docs/` contains templates for referencing GitHub artifacts in papers/reports.
- `data/` intentionally contains no large datasets unless a separate Git LFS/archive policy is used.

## Rule of Thumb

Commit code, small CSV tables, Markdown reports, and manifests.

Do not normally commit FITS files, MCMC chains, NetCDF traces, huge covariance products, Office drafts, PDFs, generated PNGs, or local logs.

