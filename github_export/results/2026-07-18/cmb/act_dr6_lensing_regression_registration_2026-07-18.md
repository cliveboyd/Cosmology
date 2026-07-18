# ACT DR6 Lensing Regression Registration

**Locked:** 18 July 2026, before downloading or evaluating the v1.2
likelihood data.

The official `act_dr6_lenslike==1.2.1` package unit tests provide the frozen
fiducial spectra and expected chi-squared values. The registered regression
cells are:

| variant | lens_only | expected chi-squared | tolerance |
| --- | --- | ---: | ---: |
| act_baseline | true | 14.06 | 0.05 |
| act_baseline | false | 14.13 | 0.05 |
| actplanck_baseline | true | 21.07 | 0.05 |
| actplanck_baseline | false | 21.46 | 0.05 |

All four cells must pass. The same packaged fiducial TT, TE, EE, BB and lensing
potential spectra must be used, with the package-prescribed conversion to
`C_L^kappakappa`. No nuisance, multipole or covariance setting may be changed.

This validates software, data placement, units and likelihood corrections. It
is not a cosmological fit. A physical SU2R run remains blocked until the theory
supplies registered scalar perturbations and lensed primary-CMB spectra.
