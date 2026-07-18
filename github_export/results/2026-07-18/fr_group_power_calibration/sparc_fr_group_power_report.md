# SPARC Group-Richness Injection-and-Recovery Calibration

Date: 2026-07-18
Completed: `2026-07-18T08:10:13.912325+00:00`

## Bottom Line

At the earlier-effect anchor |r|=`0.307254`, the p component recovers in `66.24` per cent of simulations, the 5 per cent CV component in `43.40` per cent, and the joint primary gate in `42.08` per cent. The current group-richness design is **underpowered under the locked 80 per cent joint-recovery rule** for an effect of that standardised magnitude.

The smallest tested effect with at least 80 per cent joint primary recovery is `0.5`.

## Simulation Grid

| target_abs_residual_correlation | median_recovered_rho | p_component_recovery | median_cv_rmse_fractional_improvement | cv_component_recovery | sign_component_recovery | joint_primary_recovery |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | -0.0109266 | 0.0326 | -0.00948678 | 0.0062 | 0.8464 | 0.0054 |
| 0.1 | -0.135686 | 0.1048 | -0.00660041 | 0.0228 | 0.8816 | 0.021 |
| 0.2 | -0.254327 | 0.3228 | 0.0092607 | 0.1422 | 0.957 | 0.1334 |
| 0.3 | -0.362981 | 0.6392 | 0.0402484 | 0.4112 | 0.9898 | 0.3974 |
| 0.307254 | -0.369886 | 0.6624 | 0.0426998 | 0.434 | 0.9914 | 0.4208 |
| 0.4 | -0.462697 | 0.892 | 0.0839854 | 0.7472 | 0.998 | 0.7368 |
| 0.5 | -0.55271 | 0.9902 | 0.143347 | 0.9402 | 1 | 0.9374 |
| 0.6 | -0.636429 | 0.9998 | 0.216388 | 0.9942 | 1 | 0.994 |

## Calibration

The fixed two-sided p<=0.01 critical value is |rho|>=`0.317526`, calibrated with `100000` null rank permutations. Every grid point uses `5000` empirical-residual injections.

The simulation conditions on the observed 65-galaxy design and empirical residual distribution. It does not include group-catalogue misclassification, alternative host populations, cosmic variance or uncertainty in transferring an effect from the 2MRS composite to group richness.

## Claim Boundary

This study calibrates sensitivity only. It cannot turn the observed null into evidence for FR, antimatter or an environmental effect, and the five reserved galaxies remain sealed.
