# HHT CMB/SN Code Capture

This folder contains the reproducible pieces of the historical HHT investigation.

## Current Scripts

- `diagnose_hht_resonance.py` - EMD/Hilbert diagnostic for synthetic and CSV residual series.
- `diagnose_hht_surrogate_controls.py` - shuffled-order controls for HHT residual features.
- `export_sn_residuals_for_hht.py` - Pantheon+SH0ES residual export for HHT screening.
- `export_extended_realdata_for_hht.py` - combined SN-bin, BAO, and cosmic-chronometer pull export.
- `plot_imf2_vs_bao_v3.py` and `regen_freq_vs_bao_plots.py` - IMF/BAO visualization utilities.

## Legacy Scripts

The `legacy/` folder contains normalized filename copies from `.spyder-py3`, including early
Pantheon+SH0ES/CMB TT-TE-EE model sketches, HHT residual scripts, BAO-lock scripts, and seed-grid
robustness tools.

These legacy scripts preserve provenance. They are not all publication-ready or path-portable.

