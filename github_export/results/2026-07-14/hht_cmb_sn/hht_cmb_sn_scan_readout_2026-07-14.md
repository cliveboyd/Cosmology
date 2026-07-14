# HHT CMB/SN Scan Readout

Date: July 14, 2026

## What Was Found

I found a substantial historical HHT trail spanning:

- modern HHT residual tools in the main Cosmology folder;
- July 10 HHT residual outputs under `plamb_runs/diagnostics/hht_resonance`;
- older `.spyder-py3` HHT/SN/BAO-lock products;
- early Pantheon+SH0ES plus CMB TT/TE/EE model sketches;
- `hht_bridge` posterior-link and seasonal/seed-grid robustness summaries.

## Key Technical Takeaway

The strongest archive-worthy thread is not a direct CMB-SN detection. It is a residual-screening workflow:

- SN residuals were exported against matched LCDM, SU2, and SU2R SN+BAO+Planck branches.
- Extended real-data pulls combined SN bins, BAO rows, and cosmic chronometer rows into HHT screening series.
- HHT found IMF structure in several residual series.
- Shuffled-order surrogate controls then weakened the interpretation: the extended SU2-minus-LCDM residual IMF energy was not unusual against shuffled controls.

## Useful Numbers

- SN residual weighted RMS: LCDM `0.15159`, SU2 `0.15548`, SU2R `0.15542`.
- Extended SU2-minus-LCDM dominant IMF energy fractions: IMF1 `0.6177`, IMF2 `0.4290`, IMF3 `0.2068`.
- Extended SU2-minus-LCDM surrogate control: max-energy upper-tail p `0.9703`, sum-energy upper-tail p `0.9604`.
- Older IMF2-vs-BAO baseline fit: `r_DM=0.1006, p=0.0078`; `r_DH=-0.1037, p=0.0061`.
- HHT bridge link fit: `epsilon_M - kappa*gamma_c` suggested `kappa=-0.3165`, with `sigma_L=0.1846`.

## Claim Boundary

Archive as candidate-feature and method-development evidence. Do not present as a confirmed CMB-SN HHT coupling signal without a new pre-registered run that separates probes, controls survey/window effects, and passes global/surrogate correction.

## GitHub Capture

Added under:

- `github_export/code/hht_cmb_sn/`
- `github_export/results/2026-07-14/hht_cmb_sn/`
- `github_export/manifests/hht_cmb_sn_capture_manifest_2026-07-14.csv`

