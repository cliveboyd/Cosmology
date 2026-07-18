# ACT DR6 Lensing Regression

Date: 18 July 2026

| variant | lens only | expected chi-squared | measured chi-squared | difference | pass |
| --- | --- | ---: | ---: | ---: | --- |
| act_baseline | True | 14.06 | 14.057912 | -0.0020882113 | True |
| act_baseline | False | 14.13 | 14.1315 | 0.0015004774 | True |
| actplanck_baseline | True | 21.07 | 21.071631 | 0.0016312736 | True |
| actplanck_baseline | False | 21.46 | 21.461799 | 0.0017993004 | True |

Overall software/data gate: **PASS**.

This is a package regression, not a cosmological fit. SU2R inference remains blocked until registered perturbation equations predict convergence and lensed primary-CMB spectra.
