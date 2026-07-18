# SPARC FR Published-Group-Catalogue Replication Preregistration

Date: 2026-07-18
Locked: `2026-07-18T07:53:38.567253+00:00`

## Status

Independent-predictor replication on previously generated development outcomes; not a held-out-outcome replication

The earlier 2MRS composite test had rho=-0.307254, p=0.0117994 and 2.398 per cent cross-validated improvement, and failed its locked gate.

The group covariates and sample coverage were constructed without reading outcome values in this programme. The development outcomes themselves were generated previously, so this is not a new held-out-outcome experiment.

## Catalogue

The predictor comes from [Kourkchi and Tully 2017](https://doi.org/10.3847/1538-4357/aa76db); fixed tables are archived by [VizieR](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJ/843/16). The catalogue lists 15,004 galaxies and group assignments within 3500 km/s.

The catalogue uses overlapping nearby-galaxy and Ks-band information, including 2MRS completeness; it is not an independent raw survey.

## Locked Samples

| sample | N_outcome_blind |
| --- | --- |
| primary_10_40_mpc | 65 |
| near_5_40_mpc | 96 |
| core_10_35_mpc | 61 |
| latitude_strict_10_40_mpc | 61 |
| luminous_10_40_mpc | 44 |

The primary sample is Q<=2, 10-40 Mpc and |b|>=10 degrees. Six outcome-selected galaxies and all five reserved galaxies are excluded. Accepted coordinates must satisfy the locked coordinate/velocity match; ambiguous matches are excluded.

## Primary Predictor

The single primary predictor is `group_richness_score = ln(N_group)`, where catalogue singles have N_group=1 and score zero. It avoids selecting a predictor using the development outcome and is less directly tied to Ks luminosity than catalogue halo mass.

Luminosity-derived group mass and grouped/single status are secondary diagnostics and cannot pass the gate. A +90-degree RA shifted group-richness assignment is the negative control.

## Outcome And Controls

The outcome is the frozen signed-log transformed PLAMB-minus-RAR profile-objective contrast per rotation-curve point under the combined-conventional nuisance scenario. The baseline scenario is a sensitivity control. Host morphology, luminosity, surface brightness, gas mass, velocity, inclination, distance uncertainty, bulge proxy, point count, distance, Galactic latitude and distance-method indicators are controlled by cross-fitting.

## Locked Gate

Every condition is required: two-sided permutation p<=0.01, cross-fitted RMSE improvement>=5 per cent, one coefficient sign in at least 8/10 folds, the same partial-correlation sign in all listed distance/latitude/luminosity and baseline-outcome sensitivities, and a smaller absolute shifted-sky rho.

Reserved outcomes remain sealed even after a pass. A separate held-out replication preregistration is mandatory.

## Claim Boundary

A pass would support a conventional group-environment predictor. It would not identify antimatter, estimate a signed matter/antimatter field or authorise matter/antimatter galaxy labels.

Programme SHA-256: `85263afce4ce39bdcba6886b77ac098ed527f1e7e0b32c5eaabb9d6a341641d3`
