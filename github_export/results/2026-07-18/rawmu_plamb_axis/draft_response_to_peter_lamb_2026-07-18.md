# Draft response to Peter Lamb

Dear Peter,

Thank you. Your note makes the proposed PLAMB comparison much clearer. I have now implemented the two displays for Pantheon+SH0ES, DES-Dovekie and Union3.1.

One terminology point is useful at the outset: Type Ia supernovae are standardised candles rather than exact identical candles. The published distance moduli include light-curve shape, colour, host, selection and calibration corrections, together with their uncertainties.

The key mathematical point is that dividing the observed luminosity distance by `(1+z/2)` is exactly equivalent to the fixed PLAMB forward model

`D_L(z) = (c/H0) z (1+z/2)`.

In distance-modulus form the correction is

`mu_corr = mu_obs - 5 log10(1+z/2)`.

The direct PLAMB residual and corrected-distance residual agree to better than `1e-10` magnitude in all three releases. The corrected plot is therefore a useful and intuitive representation, but it is not a different fit from the PLAMB curve already tested.

Using the same 3,422 released points, the same redshift choices, one normalisation intercept per release and each release's full covariance, the result is `Delta BIC(PLAMB-Lambda-CDM)=94.344983`. Positive values favour Lambda-CDM for this fixed PLAMB law. The luminosity-distance diagonal sensitivity gives `Delta BIC=94.340622`, but it is not the primary likelihood because it removes shared calibration and survey-systematic correlations.

The three releases also give the same direction separately: Pantheon+ `+42.200623`, DES-Dovekie `+27.936226`, and Union3.1 `+14.420284`. Adding the 44 positive-redshift Pantheon+ non-calibrator rows below the earlier `z>0.01` threshold gives a joint value of `+91.314727`, so the low-redshift extension does not reverse the result. There is no high-redshift cut in this comparison; the released range is retained through `z=2.26226`.

The Lambda-CDM comparator is not fitted as an arbitrary velocity polynomial here. It is the flat Friedmann distance integral with one fitted matter-density parameter. PLAMB has the fixed `z(1+z/2)` shape and no fitted shape parameter.

The covariance matrices are not assumptions about energy. They describe which reported supernova distances share calibration, selection, peculiar-velocity and other systematic uncertainties. Transforming from magnitude to luminosity distance does not make those correlations disappear; the covariance must be transformed with the data. A diagonal distance-error fit is still useful as a transparent sensitivity check, and I have included it as such.

There are two data limitations to the proposed `z` and `LD` error-only comparison. Pantheon+ publishes `zHDERR`, but the compressed DES-Dovekie Hubble diagram and the 22 Union3.1 spline nodes do not supply equivalent complete per-row redshift-error vectors. Union3.1 is also a correlated compression, not a list of raw individual-supernova distances. I have not invented missing redshift errors.

Although spectroscopic redshift measurement errors are usually small, peculiar velocities are not small relative to the Hubble recession velocity for the nearest objects. They are one reason for the low-redshift sensitivity check. In Pantheon+, part of this velocity uncertainty already enters the released distance error, so adding `zHDERR` again in a diagonal distance fit can be conservative or partly double-counted.

The plots now supplied are:

1. calibrated observed distance against the fitted Lambda-CDM stretched-redshift distance;
2. calibrated observed distance divided by `(1+z/2)` against measured redshift; and
3. covariance-weighted binned residuals for both forward models.

To turn the PLAMB factor into a distinct physical derivation rather than a fitted distance law, the next item I need is the explicit chain from clock-rate or light-speed evolution to `1+z/2`: the time variable, the function `c(t)` or wavelength-evolution equation, the integration limits, and the assumptions connecting emitted luminosity to observed flux. That derivation would let us test alternatives to the fixed factor without introducing an arbitrary polynomial.

Regards,
Clive
