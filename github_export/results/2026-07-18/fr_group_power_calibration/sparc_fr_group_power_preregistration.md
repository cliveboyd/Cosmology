# SPARC FR Group-Richness Injection-and-Recovery Preregistration

Date: 2026-07-18
Locked: `2026-07-18T08:08:21.338818+00:00`

## Status

This is a post-result power calibration of the failed 65-galaxy published-group test. It is not a new physical-data test and cannot rescue or alter the failed gate.

## Fixed Design

The exact locked primary frame, group-richness predictor, eleven intrinsic host controls, distance-method indicators, ten folds and ridge penalty are reused. The six previous candidates and five reserved galaxies are absent.

A full-sample host-only ridge fit defines the deterministic host signal. Its centred, standardised empirical residuals are permuted without replacement to generate noise. The host-residualised group-richness vector receives a negative injected effect, matching the direction of the earlier 2MRS composite.

Target absolute residual correlations are: `0, 0.1, 0.2, 0.3, 0.307254, 0.4, 0.5, 0.6`. The exact earlier-effect anchor is `r=0.307254`.

For target r, the injected standardised slope is beta=r/sqrt(1-r^2). There are 5,000 replicates at every grid value.

## Recovery Gates

The two-sided p<=0.01 component is represented by an absolute Spearman-rho threshold calibrated from 100,000 null rank permutations of the fixed environment residual. Each simulation also reruns the locked cross-fitted RMSE and coefficient-sign calculations.

A joint primary recovery requires the p component, RMSE improvement>=5 per cent and one coefficient sign in at least 8/10 folds in the same replicate.

## Locked Interpretation

Power is adequate for the earlier-effect anchor only if at least 80 per cent of anchor injections recover the joint primary gate. The detection threshold is the smallest tested absolute residual correlation with at least 80 per cent joint recovery.

A low recovery rate means the group-richness design is too weak or too coarse to reject an effect of that standardised magnitude. A high recovery rate would make the observed null informative against an effect carried by group richness. Neither outcome is evidence for antimatter.

Programme SHA-256: `bc03d77c9a144ba706e53ba50881e556bdc61c1a3706235346036488086aa645`
