# Data Policy

Raw scientific datasets are not committed here by default.

Use `../manifests/data_manifest.csv` to record:

- dataset name;
- local filename/path;
- source URL or DOI when known;
- approximate size;
- whether it is public, derived, or generated;
- checksum if available;
- notes on how it was used.

If a dataset must be versioned, use Git LFS or an external archive such as Zenodo/OSF and cite that external artifact from the manifest.

Curated exceptions:

- `bao/lya_validation_2026-07-19/` contains two small official eBOSS DR16
  likelihood grids and the DESI DR2 published auto/cross summary readouts. The
  directory README records provenance, licensing and the distinction between
  native likelihoods and published Gaussian compressions.
