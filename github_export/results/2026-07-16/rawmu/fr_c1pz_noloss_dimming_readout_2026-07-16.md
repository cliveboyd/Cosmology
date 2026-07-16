# FR c(z)=c0(1+z) No-Loss Dimming Readout

Date: 2026-07-16

## Executive readout

The clarified flat/no-expansion FR path

`r(z) = (c/H0) z (1 + z/2)`

was tested against raw supernova distance moduli with H0 fixed at 67.5 km/s/Mpc and with the July 14 raw-MU calibration/redshift-frame controls. The luminosity-distance sequence was

`DL_alpha = r(z) (1+z)^alpha`

for `alpha = 0, 0.5, 1`.

The result is unambiguous in this implementation: `alpha = 0` is preferred in every audited frame and offset treatment. Adding one or two extra redshift flux factors strongly worsens the fit.

## Primary HD reference

Reference row: HD frame, Pantheon+SH0ES + DES-Dovekie + Union3.1, dataset+survey calibration priors at 50 mmag dataset / 25 mmag survey.

| alpha | data chi2 | prior chi2 | posterior objective | delta vs alpha=0 | profiled RMS mag | max offset pull |
|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 3180.984 | 168.723 | 3349.707 | 0.000 | 0.224997 | 5.450 |
| 0.5 | 5421.150 | 1460.557 | 6881.707 | 3532.001 | 0.272233 | 15.610 |
| 1.0 | 11994.693 | 4473.469 | 16468.162 | 13118.456 | 0.371298 | 31.702 |

## Frame controls

Best alpha by frame and offset treatment:

| frame | offset treatment | best alpha |
|---|---|---:|
| HD | none | 0.0 |
| HD | dataset+survey priors | 0.0 |
| HEL | none | 0.0 |
| HEL | dataset+survey priors | 0.0 |
| CMB_PANTHEON_ONLY | none | 0.0 |
| CMB_PANTHEON_ONLY | dataset+survey priors | 0.0 |

## Interpretation boundary

This supports the no-extra-dimming reading of the clarified FR rule, not a broader FR/PLAMB promotion. The test says that once the path distance is set by `c(z)=c0(1+z)`, the raw-MU data strongly prefer no additional `(1+z)` flux/time-rate factor in the luminosity distance. The result remains a model-consistency diagnostic and should still be interpreted inside the existing calibration and redshift-frame gates.
