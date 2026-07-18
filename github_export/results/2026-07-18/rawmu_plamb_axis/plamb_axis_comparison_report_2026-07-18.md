# PLAMB Three-release Axis Comparison

Date: 18 July 2026

## Main finding

Peter Lamb's corrected-distance display is valid as a visualisation, but it is algebraically the same fixed PLAMB model already used in the raw-MU likelihood. It is not a new likelihood and does not alter the model comparison.

The released-total-covariance primary result has `N=3422` and `Delta BIC(PLAMB-Lambda-CDM)=94.344983`. Positive values favour Lambda-CDM.

## Release-by-release comparison

| release | N | Delta BIC (PLAMB-Lambda-CDM) | LCDM Omega_m |
| --- | ---: | ---: | ---: |
| PantheonPlusSH0ES | 1580 | 42.200623 | 0.331630 |
| DES_Dovekie_raw | 1820 | 27.936226 | 0.328733 |
| Union3p1_UNITY1p8 | 22 | 14.420284 | 0.333956 |

Every release separately favours Lambda-CDM for the fixed PLAMB law under its released total covariance.

## Algebraic identity

- PantheonPlusSH0ES: maximum residual difference `7.105e-15` mag
- DES_Dovekie_raw: maximum residual difference `1.421e-14` mag
- Union3p1_UNITY1p8: maximum residual difference `7.105e-15` mag

All releases pass the registered `1e-10` mag identity gate.

## Joint sensitivities

| likelihood | N | Delta BIC (PLAMB-Lambda-CDM) | LCDM Omega_m | status |
| --- | ---: | ---: | ---: | --- |
| released total covariance in MU | 3422 | 94.344983 | 0.330680 | primary |
| diagonal of total covariance in MU | 3422 | 101.166348 | 0.343915 | sensitivity only |
| diagonal luminosity distance plus available z error | 3422 | 94.340622 | 0.350241 | sensitivity only |
| all positive-z non-calibrators | 3466 | 91.314727 | 0.331250 | low-z stress |

## Row coverage

Primary rows: `{'PantheonPlusSH0ES': 1580, 'DES_Dovekie_raw': 1820, 'Union3p1_UNITY1p8': 22}`.
All-positive-z rows: `{'PantheonPlusSH0ES': 1624, 'DES_Dovekie_raw': 1820, 'Union3p1_UNITY1p8': 22}`.
Distance-space redshift-error coverage: `{'PantheonPlusSH0ES': {'N': 1580, 'source': 'released zHDERR'}, 'DES_Dovekie_raw': {'N': 0, 'source': 'not supplied by compressed release product'}, 'Union3p1_UNITY1p8': {'N': 0, 'source': 'not supplied by compressed release product'}}`.

DES-Dovekie and Union3.1 do not provide complete per-row redshift-error vectors in the compressed products used here. No errors were invented. Union3.1 consists of 22 correlated spline nodes rather than individual raw supernova measurements.
The Pantheon+ distance-space sensitivity may conservatively count some low-redshift velocity uncertainty twice because peculiar-velocity contributions already enter the released distance uncertainty. This is another reason it is not the primary likelihood.

## Statistical interpretation

Type Ia supernovae are standardised candles, not exact identical candles. Their released distance moduli incorporate light-curve shape, colour, host, selection and calibration corrections with associated uncertainties.

The released covariance matrices encode statistical uncertainty, shared calibration, bias-correction, peculiar-velocity and survey systematic modes. They are not energy matrices. Replacing them with independent luminosity-distance errors discards measured correlation structure and is therefore a diagnostic sensitivity, not the preferred likelihood.

The Lambda-CDM comparator here is a Friedmann-equation distance integral with one fitted `Omega_m`, not a polynomial in velocity, deceleration and jerk. The primary sample retains the full released high-redshift range through `z=2.26226`; only Pantheon+ applied the registered `z>0.01` lower cut. Restoring 44 positive-redshift non-calibrator rows leaves the model direction unchanged.

The correct model comparison is a forward prediction for the same observed distance-modulus vector under identical nuisance terms. Calling one plot an x-axis stretch and the other a y-axis division does not change that requirement.
