# SPARC FR Environmental-Asymmetry Study

Date: 2026-07-18
Completed: `2026-07-18T07:22:27.170305+00:00`

## Bottom Line

The preregistered development gate is **FAIL**. The primary charge-blind 2MRS environment association has partial Spearman rho=`-0.307254`, two-sided permutation p=`0.0117994`, and cross-fitted RMSE change=`2.398` per cent.

This test cannot identify antimatter. It asks whether ordinary observed environment predicts the residual difference between the current PLAMB bridge and RAR after conventional host controls. A positive association would motivate an environment-conditioned FR response; it would not measure a signed matter/antimatter background.

## Primary Test

The primary sample contains `66` previously unselected, non-reserved Q<=2 galaxies at 10-40 Mpc and |b|>=10 degrees. Neighbours are 2MRS galaxies with M_K<=-22 and |Delta v|<=500 km/s.

The RA+90-degree negative control has rho=`0.0308319` and p=`0.80956`.

## Locked Test Matrix

| sample | scenario | control | N | partial_spearman_rho | permutation_p_two_sided | cv_rmse_fractional_improvement | environment_coefficient_sign_fraction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| primary_10_40_mpc_mk22 | combined_conventional | actual | 66 | -0.307254 | 0.0117994 | 0.0239843 | 1 |
| primary_10_40_mpc_mk22 | combined_conventional | shifted_ra90 | 66 | 0.0308319 | 0.80956 | -0.00620592 | 0.6 |
| near_5_40_mpc_mk22 | combined_conventional | actual | 97 | -0.22687 | 0.0245988 | 0.0118278 | 1 |
| wide_10_80_mpc_mk23 | combined_conventional | actual | 80 | -0.134201 | 0.243088 | -0.0412698 | 0.8 |
| latitude_strict_10_40_mpc_mk22 | combined_conventional | actual | 61 | -0.154045 | 0.233388 | -0.0292581 | 0.9 |
| primary_10_40_mpc_mk22 | baseline | actual | 66 | -0.221376 | 0.0751462 | 0.0141707 | 1 |

## Gate Conditions

| condition | passed |
| --- | --- |
| primary_p | False |
| cv_improvement | False |
| coefficient_sign | True |
| sensitivity_signs | True |
| negative_control | True |

Replication was not unsealed. Even a passing development gate requires a separate frozen replication programme.

## Previously Selected Galaxies

The following environment comparison is descriptive only because the six galaxies were selected previously from model outcomes. No candidate outcome enters the primary test.

| candidate | candidate_environment_score | matched_control_mean_environment_score | candidate_minus_control_environment_score |
| --- | --- | --- | --- |
| UGC06787 | 0.0664188 | -0.782241 | 0.84866 |
| UGC03580 | -1.05922 | -0.400735 | -0.658489 |
| UGC02487 | -0.068884 | -0.581377 | 0.512493 |
| NGC5985 | 0.141995 | -0.588562 | 0.730557 |
| UGC09133 | -0.645282 | -0.569472 | -0.0758094 |
| NGC2903 | -0.890017 | -0.089046 | -0.800971 |

## Data and Model Boundary

2MRS is a K_s<=11.75 near-infrared redshift catalogue with 97.6 per cent redshift completeness over 91 per cent of the sky. [2MRS catalogue](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJS/199/26); [Huchra et al. 2012](https://doi.org/10.1088/0067-0049/199/2/26). Central positions and redshifts were resolved independently through [NASA/IPAC NED](https://ned.ipac.caltech.edu/Documents/Guides/Interface/ObjectLookup).

The proxy is charge-blind: it traces luminous galaxy environment, not matter-minus-antimatter sign. Its principal limitations are redshift-space projection, 2MRS luminosity selection, peculiar velocities in the nearby volume and the use of K-band luminosity as a mass tracer.

Deterministic profile optimiser failures: `0` out of `444` model-galaxy-scenario fits.

## Locked Claim Boundary

- Do not label any galaxy as antimatter from this study.
- Do not reinterpret a charge-blind density association as a signed FR background.
- Do not inspect the five reserved replication outcomes under this programme.
- If the gate fails, retain conventional environment only as a possible nuisance/control path.
- If the gate passes, preregister an independent replication before opening any reserved outcome.
