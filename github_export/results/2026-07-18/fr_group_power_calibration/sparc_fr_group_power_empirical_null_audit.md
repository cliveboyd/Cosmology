# SPARC Group-Power Empirical-Null Audit

Date: 2026-07-18
Completed: `2026-07-18T08:12:42.703236+00:00`

## Status

This is a post hoc calibration diagnostic. It preserves the locked simulation and cannot alter its preregistered result.

## Finding

The locked rank-permutation threshold was |rho|>=`0.317526`. Under the complete empirical-residual null pipeline, its p-component false-positive rate was `3.26` per cent rather than 1 per cent.

The empirical 1 per cent threshold is |rho|>=`0.377185`. At the earlier-effect anchor, p-component recovery changes from `66.24` to `47.58` per cent, and joint recovery changes from `42.08` to `36.50` per cent.

The smallest tested effect with at least 80 per cent empirically recalibrated joint recovery remains `0.5`.

## Recalibrated Grid

| target_abs_residual_correlation | locked_p_recovery | empirical_null_p_recovery | locked_joint_primary_recovery | empirical_null_joint_recovery |
| --- | --- | --- | --- | --- |
| 0 | 0.0326 | 0.01 | 0.0054 | 0.0032 |
| 0.1 | 0.1048 | 0.0438 | 0.021 | 0.0154 |
| 0.2 | 0.3228 | 0.1856 | 0.1334 | 0.11 |
| 0.3 | 0.6392 | 0.4542 | 0.3974 | 0.3528 |
| 0.307254 | 0.6624 | 0.4758 | 0.4208 | 0.365 |
| 0.4 | 0.892 | 0.7694 | 0.7368 | 0.6878 |
| 0.5 | 0.9902 | 0.9546 | 0.9374 | 0.9184 |
| 0.6 | 0.9998 | 0.9982 | 0.994 | 0.9928 |

## Interpretation

The rank-permutation p component is anti-conservative for this empirical-noise and cross-fitting pipeline. The full locked joint gate remains conservative under the null because it also requires a 5 per cent predictive improvement: its locked null pass rate is 0.54 per cent, while the empirically recalibrated joint null pass rate is 0.32 per cent.

The correction strengthens the conclusion that the 65-galaxy group-richness design is underpowered at the earlier-effect anchor. It does not change the observed catalogue result, whose p-value was 0.704, and it creates no evidence for FR or antimatter.
