# Priority follow-up protocol: Quaia/SU2 and supernova calibration

**Locked:** 2026-07-18, before inspection of any new SDSS DR16Q, injection-calibrated, or release-budget fit outcome.

## Scope

This protocol joins three separately executable analyses without pooling their evidence. Each branch has its own machine-readable preregistration and result manifest. A failure in one branch cannot be offset by success in another.

## 1. Independent quasar catalogue

- Independent catalogue: SDSS DR16Q v4 spectroscopic quasar catalogue.
- Primary sample: `IS_QSO_FINAL = 1`, finite position and primary spectroscopic redshift, `1.00 <= z < 1.50`, and `|b| >= 35 degrees`.
- The SDSS footprint is retained. Multipole statistics are calibrated on the observed mask; no full-sky interpretation is permitted.
- The original Quaia target direction is fixed at equatorial `(RA, Dec) = (145.3283284733, -26.7809296987) degrees`.
- Redshift-window, latitude-cut and spatial-block-scale families are fixed in the branch preregistration before the catalogue outcome is calculated.
- Promotion is conjunctive: empirical global `p <= 0.01` at every locked spatial-block scale; target-direction separation `<= 30 degrees`; amplitude ratio in `[0.5, 2.0]` across the primary stability family; and mask-calibrated joint `l = 1, 2, 3` coherence `p <= 0.05`.

## 2. Injection-calibrated angular likelihood

- The null model contains the locked dust, depth, Gaia scanning-geometry, WISE-colour and catalogue-quality templates.
- The candidate model adds the locked SU2 scalar-residual/angular component without changing the nuisance basis.
- Training and assessment use fixed spatial sky folds. The primary score is held-out predictive gain, not in-sample BIC.
- The decision threshold is derived from locked null and signal injections. The former fixed `Delta BIC < -10` rule is retained only as a historical diagnostic.
- Promotion is conjunctive: empirical null false-positive probability `<= 0.01`; recovery probability `>= 0.80` at the locked one-times signal; positive pooled held-out predictive gain over template-only; and recovery of the locked amplitude/direction tolerances.
- Failure to attain 80 per cent power is reported as underpowered. It does not authorise moving the threshold after outcomes are known.

## 3. Supernova calibration-budget hierarchy

- `zHD` is the primary readout. `zHEL` and `zCMB` are controls only.
- FR and Lambda-CDM use identical data rows, covariance treatment, offset design, priors, optimisation bounds and hold-out definitions.
- Pantheon+ calibration modes are derived from its statistical-only covariance and published calibration systematic groups. DES modes are derived from its statistical-only and total precision releases. Union3 uses release-provided covariance and documented zero-point budgets only to the extent that separable components are available.
- Released total covariances are reconstruction checks. A calibration mode already present in a total covariance is not added again as an independent prior.
- Dataset, survey and high-redshift hold-outs are fixed before fit outcomes. A robust FR preference must preserve its sign under all primary hold-outs, remain within every external calibration budget, and not require a larger calibration pull than Lambda-CDM under the same nuisance model.

## Primary sources

- SDSS DR16Q catalogue and description: https://www.sdss4.org/dr17/algorithms/qso_catalog/
- SDSS DR16Q data model: https://data.sdss.org/datamodel/files/BOSS_QSO/DR16Q/DR16Q_v4.html
- Pantheon+ official data release: https://github.com/PantheonPlusSH0ES/DataRelease
- Pantheon+ calibration covariance groups: https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/sytematic_groupings
- Union3/UNITY analysis: https://doi.org/10.3847/1538-4357/adc0a5

