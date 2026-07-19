# Lyman-alpha BAO validation inputs

This directory contains the small public inputs used by the 19 July 2026
PLAMB high-redshift validation audit.

- `desi_dr2_lya_alphas.txt` and `desi_dr2_lya_data_splits.txt` preserve the
  numerical contents of the figure-data archive for DESI DR2 Results I,
  Zenodo record `15690869` (CC BY 4.0); one trailing blank line was removed
  from the split-summary text file.
- `sdss_DR16_LYAUTO_BAO_DMDHgrid.txt` and
  `sdss_DR16_LYxQSO_BAO_DMDHgrid.txt` are copied without modification from
  the official SDSS DR16 cosmology likelihood release, tag `v1_0_1`. They are
  the native non-Gaussian auto-correlation and Ly-alpha x quasar BAO
  likelihood grids at `z=2.334`.
- `lya_published_measurements_2026-07-19.json` records the small published
  Gaussian summaries used for DESI DR1 and the alternative DESI DR2
  multipole-method control. It also records the official source locations and
  the evidential classification of each measurement.

The DESI DR2 figure-data archive does not include the native 9,306-element
correlation covariance, likelihood samples or a two-dimensional likelihood
surface. Consequently the DR2 baseline and its auto/cross split are evaluated
from the collaboration's published compressed posterior summaries. They are
not labelled as native-likelihood reruns.

## Primary sources

- DESI DR2 baseline: https://arxiv.org/abs/2503.14739
- DESI DR2 figure data: https://zenodo.org/records/15690869
- DESI DR2 multipoles: https://arxiv.org/abs/2603.04281
- DESI DR1: https://arxiv.org/abs/2404.03001
- SDSS/eBOSS DR16 likelihoods:
  https://svn.sdss.org/public/data/eboss/DR16cosmo/tags/v1_0_1/likelihoods/BAO-only/
