# CMB Branch Readiness

Date: 18 July 2026

| branch | data | software | theory | scientific | boundary |
| --- | --- | --- | --- | --- | --- |
| planck_distance_prior | True | True | True | True | standard early physics and recombination only |
| act_dr6_primary_cmb | True | True | True | True | new-model use requires CAMB/CLASS-compatible spectra |
| firas_spectral_distortion | False | True | False | False | blocked by full covariance and registered dQ/dz or photon-kinetic mapping |
| planck_low_l | False | False | True | False | blocked by PR4 maps/masks/simulations and frozen global statistic family |
| tb_eb_parity | False | False | False | False | blocked by split maps, angle priors, foreground model and parity-odd transfer spectra |
| act_dr6_lensing | True | True | False | False | Lambda-CDM regression ready; physical SU2R inference is blocked by perturbation equations |

A false scientific-readiness flag is a hard implementation block. It must not be bypassed with a phenomenological proxy carrying a physical model label.
