# Research Document GitHub Reference Template

Use this pattern in reports:

```text
Code and reproducibility artifacts are archived in the project GitHub repository under:
github_export/code/<topic>/

Small result summaries and tables used in this section are archived under:
github_export/results/<YYYY-MM-DD>/<topic>/

Large raw datasets and generated posterior/mock products are not committed directly.
They are listed with source/local provenance in:
github_export/manifests/data_manifest.csv
github_export/manifests/results_manifest_<YYYY-MM-DD>.csv
```

For a specific claim, cite the exact script and result table, for example:

```text
The SU2/Quaia global look-elsewhere control was run with
github_export/code/quaia_su2/su2_quaia_global_lookelsewhere.py.
The observed-threshold table is archived at
github_export/results/2026-07-14/su2_quaia/su2_quaia_global_lookelsewhere_observed_thresholds_2026-07-14.csv.
```

