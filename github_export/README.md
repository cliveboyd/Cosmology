# GitHub Export Layout

This folder is the GitHub-ready subset of the local cosmology workspace.

## Layout

- `code/` contains copied Python programs and launch scripts needed to reproduce the main diagnostics.
- `results/` contains small CSV/Markdown outputs suitable for direct citation in reports.
- `runs/` contains dated launch scripts for longer or background diagnostics.
- `manifests/` records where raw data and larger outputs live.
- `docs/` contains templates and extracted report text/media for referencing GitHub artifacts in papers/reports.
- `data/` contains only small curated samples unless a separate Git LFS/archive policy is used.

## Rule of Thumb

Commit code, small CSV tables, Markdown reports, and manifests.

Do not normally commit FITS files, MCMC chains, NetCDF traces, huge covariance products, Office drafts, PDFs, generated PNGs, or local logs.

Exceptions should be deliberate and documented in manifests, such as small sample catalogs or representative figures needed to preserve a visual result.
