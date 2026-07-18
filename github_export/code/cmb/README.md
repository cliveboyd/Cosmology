# Registered CMB Programme

This directory implements guarded entry points for the five CMB branches registered
on 18 July 2026. It separates data and software validation from physical-model
readiness so that a background-only proxy cannot be reported as an SU(2) perturbation
or photon-thermalisation result.

## Reproduction order

Run commands from the repository root.

```powershell
python github_export/code/cmb/download_planck_pr3_baseline_chains_2026-07-18.py
python github_export/code/cmb/derive_planck_pr3_distance_prior_2026-07-18.py
python github_export/code/cmb/validate_planck_pr3_distance_prior_2026-07-18.py
python github_export/code/cmb/validate_act_dr6_lite_2026-07-18.py
python github_export/code/cmb/download_act_dr6_lensing_data_2026-07-18.py
python github_export/code/cmb/validate_act_dr6_lensing_2026-07-18.py
python github_export/code/cmb/download_firas_monopole_2026-07-18.py
python github_export/code/cmb/validate_firas_ingestion_2026-07-18.py
python github_export/code/cmb/audit_cmb_branch_readiness_2026-07-18.py
python github_export/code/cmb/test_cmb_contracts.py
```

Large public datasets are stored below `external_datasets/current_cosmology_sources`
and are excluded from Git. Their source URLs, hashes and validation metadata are
recorded in local manifests and the dated result reports.

## Theory boundary

`cmb_theory_contract.py` defines the required background, perturbation and
spectral-distortion interfaces. `planck_distance_prior.py` is valid only for the
registered standard-early-physics branch. ACT primary-CMB and lensing interpretation
requires model spectra; parity work requires parity-odd spectra or an explicitly
registered birefringence transform; FIRAS interpretation requires a physical
energy-injection or photon-kinetic prediction.

## Current gates

- Planck three-variable distance prior: implemented and CAMB-validated.
- ACT DR6 primary-CMB likelihood: official fixed-spectrum regression implemented.
- ACT DR6 lensing likelihood: official four-cell regression is the software gate.
- FIRAS: ingestion is implemented; a primary fit is blocked without full covariance.
- Planck low-l: the leave-one-out simulation maximum-statistic engine is implemented;
  map evaluation is blocked pending maps, masks and simulations.
- TB/EB: spectrum conventions and the constant-rotation transform are implemented;
  map evaluation is blocked pending split maps, angle priors and foreground choices.
