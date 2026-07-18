# PLAMB Axis-comparison Preregistration

**Locked:** 18 July 2026, before evaluating the new coordinate plots,
distance-space sensitivity or no-low-redshift-cut cell.

## Motivation

Peter Lamb proposed expressing the fixed PLAMB supernova law by dividing the
observed luminosity-distance scale by `1+z/2` while leaving measured redshift on
the horizontal axis. This follow-up tests that representation in Pantheon+SH0ES,
DES-Dovekie and Union3.1 without changing the established raw-MU likelihood.

The earlier release-grounded result is already known and is not treated as blind:
the 3,422-point total-covariance comparison gave
`Delta BIC(FR-Lambda-CDM)=+94.344983`.

## Algebraic Identity Gate

For the fixed PLAMB law,

\[
\begin{aligned}
D_L^{\rm PLAMB}(z)
  &= \frac{c}{H_0}z\left(1+\frac{z}{2}\right), \\
\mu_{\rm corr}
  &= \mu_{\rm obs}-5\log_{10}\left(1+\frac{z}{2}\right), \\
\mu_{\rm linear}(z)
  &= 5\log_{10}\left(\frac{c}{H_0}z\right)+25.
\end{aligned}
\]

The residuals `mu_obs-mu_PLAMB` and `mu_corr-mu_linear` must agree to an
absolute tolerance of `1e-10` mag for every retained row. Failure invalidates
all derived plots.

## Primary Likelihood

- Use the same 3,422 rows as the completed release-grounded analysis.
- Use `zHD` for Pantheon+ and DES-Dovekie, and released node `z` for Union3.1.
- Use released total covariance or precision for every catalogue.
- Profile one flat magnitude intercept per release for both models.
- Compare fixed PLAMB `z(1+z/2)` with flat Lambda-CDM having one fitted
  `Omega_m` in `[0.05,0.60]`.
- Keep `H0=67.5 km s^-1 Mpc^-1`; release intercepts absorb the absolute scale.
- Report joint and release-by-release chi-squared and BIC values.

## Registered Displays

1. **Distance versus stretched-redshift distance:** use the fitted Lambda-CDM
   coordinate
   `x_LCDM=(1+z) integral_0^z dz'/E(z')` and calibrated
   `y=H0 D_L,obs/c`. The identity line is `y=x_LCDM`.
2. **Corrected distance versus measured redshift:** use
   `y_PLAMB=(H0/c)D_L,obs/(1+z/2)` and the identity line `y_PLAMB=z`.
3. Show all retained points with low-opacity markers and fixed logarithmic
   redshift bins. Do not tune bins after inspecting residuals.

## Sensitivities

- Repeat the comparison with only the diagonal of each released total
  covariance. Label this non-primary because it discards shared calibration and
  systematic modes.
- Repeat in luminosity-distance space using propagated diagonal magnitude
  errors and reported redshift errors where available. Do not invent missing
  redshift errors. Pantheon+ supplies `zHDERR`; the compressed DES-Dovekie and
  Union3.1 products do not supply an equivalent complete redshift-error vector.
- Add all finite, positive-redshift, non-calibrator Pantheon+ rows to a
  no-low-redshift-cut stress cell; retain all DES and Union3.1 rows.
- The primary result remains the released-total-covariance magnitude-space
  likelihood. Sensitivities cannot override it.

## Interpretation Gates

- The corrected-axis display is a re-expression of the fixed PLAMB law, not a
  new model or independent likelihood.
- A covariance matrix is an uncertainty model, not an energy matrix. Removing
  correlations is permissible only as a labelled sensitivity.
- Redshift precision does not remove peculiar-velocity uncertainty at low
  redshift or release-wide calibration correlations in distance modulus.
- Union3.1 contributes 22 correlated spline-compressed distance nodes, not raw
  individual-supernova luminosity distances.
