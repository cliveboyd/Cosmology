# Fixed-gamma=2.3 PLAMB conditional analysis

**Analysis date:** 2026-07-19  
**Status:** outcome-aware conditional target; cannot promote a PLAMB claim  
**Primary split:** fit z<0.5, predict z>=0.5 with released covariance

## Headline result

With gamma_c fixed at 2.3, the full-sample PATH fit gives p=0.797429 and Delta BIC=+0.011958 relative to flat Lambda-CDM. The two curves are statistically indistinguishable by this BIC threshold when gamma_c is treated as externally fixed.

The z<0.5 PATH posterior gives p=0.800469 +/- 0.023641, with 95% interval [0.754062, 0.846731]. Its low-redshift Delta BIC is -0.379408.

After integrating p and each release intercept, the high-redshift predictive Delta log density is -0.240554, where positive values favour PATH. This is a negligible preference for Lambda-CDM. The PATH upper-tail area is 0.999219.

This is a conditional theory target, not independent support: gamma_c=2.3 was chosen after inspecting these same high-redshift outcomes.

## Equal-parameter fits

| Sample | Model | Shape parameter | chi-squared | BIC | Delta BIC vs LCDM |
|---|---|---:|---:|---:|---:|
| full | `PATH_GAMMA2P3` | 0.797429 | 3047.674476 | 3080.226397 | +0.011958 |
| full | `LCDM` | 0.330680 | 3047.662517 | 3080.214439 | +0.000000 |
| z<0.5 | `PATH_GAMMA2P3` | 0.800422 | 2105.992970 | 2137.019464 | -0.379408 |
| z<0.5 | `LCDM` | 0.345514 | 2106.372379 | 2137.398872 | +0.000000 |

Both cells contain one shape parameter and one release intercept per included release. The BIC differences therefore equal the chi-squared differences.

## Low-redshift posteriors

| Model | Mean | Standard deviation | 68% interval | 95% interval | Edge Delta chi-squared |
|---|---:|---:|---|---|---|
| `PATH_GAMMA2P3` | 0.800469 | 0.023641 | [0.776842, 0.823861] | [0.754062, 0.846731] | 64.69, 63.32 |
| `LCDM` | 0.345903 | 0.016442 | [0.329471, 0.362168] | [0.313966, 0.378415] | 73.18, 56.76 |

## Integrated high-redshift prediction

| Model | N high-z | Log predictive density | Mean predictive chi-squared | Upper-tail area | Two-sided tail |
|---|---:|---:|---:|---:|---:|
| `PATH_GAMMA2P3` | 1085 | 211.407273 | 943.382492 | 0.999219 | 0.001562 |
| `LCDM` | 1085 | 211.647826 | 944.604924 | 0.999092 | 0.001816 |

### Release-level predictive readout

| Release | PATH log density | LCDM log density | Delta log density | PATH upper tail | PATH two-sided tail |
|---|---:|---:|---:|---:|---:|
| PantheonPlusSH0ES | 98.506934 | 98.088386 | +0.418548 | 0.999947 | 0.000106 |
| DES_Dovekie_raw | 89.946654 | 90.192923 | -0.246269 | 0.953459 | 0.093082 |
| Union3p1_UNITY1p8 | 22.920311 | 22.967897 | -0.047586 | 0.665595 | 0.668810 |

Joint mixture log density is the registered comparison. Release-level mixture densities share the same joint p posterior and are descriptive; they are not algebraically additive after parameter integration.
The registered gate uses only a minimum upper-tail area. The additional two-sided readout is a post hoc calibration warning: values near zero indicate chi-squared that is unusually high or unusually low. Both cosmologies have unusually low joint quadratic scores, with Pantheon+ providing the strongest effect.

## Release-specific p stability

| Release | N low-z | Mean p | Standard deviation | 95% interval |
|---|---:|---:|---:|---|
| PantheonPlusSH0ES | 1370 | 0.776678 | 0.035928 | [0.706180, 0.847014] |
| DES_Dovekie_raw | 958 | 0.818392 | 0.039673 | [0.740564, 0.896080] |
| Union3p1_UNITY1p8 | 9 | 0.819345 | 0.051373 | [0.718628, 0.920008] |

Maximum pairwise release tension: 0.779355 standard deviations.

## Physical consequences

Using the joint low-redshift posterior mean p=0.800469:

| z | c(z)/c0 | abs(dz/dT)/H0 | PATH integral | Anchored Delta mu vs LCDM (mag) |
|---:|---:|---:|---:|---:|
| 0.500000 | 2.150000 | 1.383425 | 0.652053 | -0.006697 |
| 1.000000 | 3.300000 | 1.741668 | 1.519839 | +0.007602 |
| 1.500000 | 4.450000 | 2.082279 | 2.530728 | +0.036094 |
| 2.262260 | 6.203198 | 2.576653 | 4.268465 | +0.085484 |

The light-speed and clock-rate ratios are implications of the phenomenological parameterisation. They are not independently measured quantities. In particular, the large c(z)/c0 values require a separate covariant derivation and consistency checks against non-supernova observables.

## Conditional gates

| Gate | Pass |
|---|---|
| `full_absolute_delta_BIC_at_most_2` | **True** |
| `low_z_absolute_delta_BIC_at_most_2` | **True** |
| `high_z_delta_log_predictive_density_nonnegative` | **False** |
| `joint_PATH_tail_area_at_least_0p05` | **True** |
| `every_release_PATH_tail_area_at_least_0p05` | **True** |
| `maximum_pairwise_release_p_tension_at_most_2_sigma` | **True** |
| `p_posterior_interior_and_grid_edges_delta_chi2_at_least_30` | **True** |
| `independent_gamma_derivation` | **False** |

Statistical gates passed: 6/7.

**Claim decision: DO NOT PROMOTE.** The independent-gamma gate is false by construction, and the value 2.3 remains an outcome-aware target for theory rather than an externally predicted constant.

## Reproduction

```powershell
python github_export/code/rawmu/analyze_plamb_gamma2p3_2026_07_19.py
```
