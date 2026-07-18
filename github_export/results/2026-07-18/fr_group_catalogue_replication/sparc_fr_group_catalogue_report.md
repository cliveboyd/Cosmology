# SPARC FR Published-Group-Catalogue Replication

Date: 2026-07-18
Completed: `2026-07-18T07:54:56.855133+00:00`

## Bottom Line

The preregistered development gate is **FAIL**. For the primary published group-richness predictor, partial Spearman rho=`0.0476836`, two-sided permutation p=`0.704415`, and cross-fitted RMSE change=`-5.410` per cent.

This is an independent-predictor replication on previously generated development outcomes. It is not an independent held-out-outcome replication. Kourkchi-Tully group construction is independent of this project's neighbour-count algorithm, but its nearby-galaxy and Ks-band inputs overlap 2MRS.

## Locked Test Matrix

| sample | scenario | control | predictor | N | partial_spearman_rho | permutation_p_two_sided | cv_rmse_fractional_improvement | environment_coefficient_sign_fraction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| primary_10_40_mpc | combined_conventional | actual | group_richness_score | 65 | 0.0476836 | 0.704415 | -0.0540986 | 0.8 |
| near_5_40_mpc | combined_conventional | actual | group_richness_score | 96 | 0.000203473 | 0.9986 | -0.027374 | 0.9 |
| core_10_35_mpc | combined_conventional | actual | group_richness_score | 61 | 0.196245 | 0.127844 | 0.0266283 | 1 |
| latitude_strict_10_40_mpc | combined_conventional | actual | group_richness_score | 61 | 0.0845584 | 0.509575 | -0.0611719 | 0.8 |
| luminous_10_40_mpc | combined_conventional | actual | group_richness_score | 44 | 0.0427061 | 0.784211 | -0.0815618 | 0.875 |
| primary_10_40_mpc | baseline | actual | group_richness_score | 65 | 0.0333042 | 0.79236 | -0.040622 | 0.9 |
| primary_10_40_mpc | combined_conventional | shifted_ra90 | shifted_group_richness_score | 65 | 0.126311 | 0.320184 | -0.0123885 | 1 |
| primary_10_40_mpc | combined_conventional | secondary_group_mass | group_log_mass_luminosity | 65 | 0.0245629 | 0.846958 | -0.0346974 | 0.9 |
| primary_10_40_mpc | combined_conventional | secondary_grouped_indicator | grouped_indicator | 65 | 0.0276224 | 0.822909 | -0.0210472 | 0.8 |

## Gate

| condition | passed |
| --- | --- |
| primary_p | False |
| cv_improvement | False |
| coefficient_sign | True |
| sensitivity_signs | True |
| negative_control | False |

The five reserved galaxies remain sealed. No result in this programme identifies antimatter or authorises matter/antimatter labels.

## Sources

[Kourkchi and Tully 2017](https://doi.org/10.3847/1538-4357/aa76db); [VizieR fixed catalogue tables](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJ/843/16).
