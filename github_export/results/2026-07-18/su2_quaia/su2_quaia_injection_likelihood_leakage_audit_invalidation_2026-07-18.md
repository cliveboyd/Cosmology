# Outcome-blind invalidation notice: response leakage in the first injection-calibrated angular likelihood

**Date:** 18 July 2026  
**Status:** methodologically invalidated for held-out prediction  
**Scope:** the first `su2_quaia_injection_calibrated_angular_likelihood` branch only

## Reason for invalidation

The first branch defines the object-level SU(2) angular predictor as

$$
x_{\mathrm{SU2},i}=s(z_i)\,n_i,
$$

where the same observed photometric redshift $z_i$ is also the response being predicted. Consequently, a nominally held-out response enters the held-out design matrix through $s(z_i)$. Cross-validation therefore does not reproduce prediction for an unseen response: the test-fold outcome has already been used to construct one of its predictors.

This is target leakage. Its cross-validated log predictive density (CVLPD), injection recovery and any decision rule that depends on those quantities are not valid held-out predictive tests of the SU(2) term. In-sample summaries remain descriptive only and cannot repair the leakage.

## Outcome-blind basis

This notice follows from the predictor-response dependency alone. It does not depend on the sign, size or statistical significance of the first branch's observed result. The branch had already failed its preregistered promotion rule, so this methodological invalidation does not withdraw or weaken any positive claim: no positive claim was made.

## Corrective requirement

A replacement analysis must construct the scalar curve from information independent of the Quaia photometric-redshift response. The corrected design will use a positional, one-to-one match to official SDSS DR16Q v4 quasars and define

$$
y_i=z_{\mathrm{Quaia},i}-z_{\mathrm{spec},i},
$$

while evaluating the scalar curve only at the independent DR16Q spectroscopic redshift $z_{\mathrm{spec},i}$. All predictors must remain unchanged when the Quaia photometric-redshift response is permuted. That invariance will be tested explicitly before the replacement decision is reported.

## Claim boundary

- The first branch is not evidence for held-out SU(2) prediction.
- Its already-negative promotion decision remains negative.
- No physical conclusion about SU(2) follows from the leakage diagnosis alone.
- Only the separately preregistered, cross-matched spectroscopic-redshift replacement may supply a valid predictive readout.
