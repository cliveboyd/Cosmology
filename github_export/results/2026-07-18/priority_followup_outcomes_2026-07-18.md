# Priority Follow-up Outcomes

Date: 18 July 2026

This readout closes the three branches frozen in
`priority_followup_protocol_2026-07-18.md`. Detailed reports, tables,
configurations and manifests remain in the topic directories.

## 1. Quaia Cross-catalogue Validation

The independent SDSS DR16Q v4 test is **not promoted**.

- The registered primary sample contains 139,305 quasars at
  `1.0 <= z < 1.5` and `|b| >= 35 deg`.
- The primary SDSS dipole amplitude is 0.0481791, but its direction is
  116.205 degrees from the locked Quaia direction; the gate required at
  most 30 degrees.
- Redshift-window perturbations fail direction stability by 165.587 and
  141.879 degrees, and the `1.0 <= z < 1.4` amplitude ratio is 0.106361.
- The global family probability is small at each registered block scale,
  but the individual partial `l = 1`, `2` and `3` probabilities are
  ordinary. The significant joint statistic is therefore compatible with
  footprint-coupled catalogue structure rather than a coherent replicated
  physical multipole pattern.

The fixed promotion rule fails four components. No exploratory choice of
window, mask, block scale or coordinates may be used to recast this as
positive SU(2) evidence.

Detailed report:
`su2_quaia/sdss_dr16q_v4_cross_catalogue_validation_2026-07-18_report.md`.

## 2. Injection-calibrated Angular Likelihood

The first Quaia-only implementation is **methodologically invalidated**
because its candidate predictor was evaluated at the same photometric
redshift used as the response. Its negative outcome is preserved only as an
audit record and is not scientific evidence.

The registered replacement uses reciprocal one-to-one Quaia x DR16Q matches
within 1 arcsec, independent spectroscopic redshift in every predictor, and

`y_i = z_Quaia,i - z_spec,i`.

The candidate angular term is

`y_i = X_i beta + A s(z_spec,i) (n_i dot d) + epsilon_i`,

where `X_i` jointly includes unmodulated sky direction, cubic spectroscopic
redshift, the scalar main effect, dust, depth, Gaia scanning geometry, WISE
colour, Quaia quality and DR16Q quality. Twelve spatial sectors supply
held-out prediction. A response-permutation audit gives bit-identical
predictors.

The full-rank 46,687-object fit has these results:

- frozen 99th-percentile null threshold: 3.38514;
- untouched-null false-positive rate: 0.009;
- predictive detection rate at the locked one-times signal: 0.080;
- joint direction/amplitude recovery at one-times: 0.000;
- observed held-out Delta CVLPD: -4.0084;
- empirical null probability: 0.727273;
- observed direction separation: 57.1066 degrees.

The outcome is **underpowered**, not a physical rejection. The present
crossmatch cannot recover the locked amplitude with the required 80 per cent
power. The observed candidate also fails held-out predictive gain and
direction gates, so it cannot be promoted.

Detailed reports:
`su2_quaia/su2_quaia_injection_likelihood_leakage_audit_invalidation_2026-07-18.md`
and
`su2_quaia/su2_quaia_crossmatched_specz_angular_likelihood_report_2026-07-18.md`.

## 3. Supernova Calibration-budget Hierarchy

The primary raw-MU comparison uses released total `STAT+SYS` covariance,
`zHD` for Pantheon+ and DES, identical nuisance treatment for FR and
Lambda-CDM, and one profiled release intercept. Grouped dataset and survey
offset modes are sensitivity diagnostics only because they do not reproduce
the release covariance.

For 3,422 supernovae, the fixed FR form

`D_L = (c/H_0) z (1 + z/2)`

is compared with flat Lambda-CDM with fitted `Omega_m`. The primary result is
`Delta BIC(FR - Lambda-CDM) = +94.344983` at
`Omega_m = 0.330680`; positive values favour Lambda-CDM. Dataset hold-outs
and all three registered high-redshift hold-outs retain that direction.

The stricter promotion gate nevertheless fails:

- only 10 of 17 eligible survey hold-outs preserve the full-sample sign,
  below the required 80 per cent;
- several FR-minus-Lambda-CDM survey residual modes exceed release-grounded
  calibration budgets by factors from about 2.5 to 11.45;
- dataset/survey common modes leave 94.0 per cent of the Pantheon+ `CALIB`
  covariance norm and 99.9 per cent of the DES all-systematics covariance
  norm unexplained.

Thus the primary likelihood strongly favours Lambda-CDM over this fixed FR
form, but the result is not promoted as a survey-robust model comparison
contained within calibration budgets under the registered conjunctive rule.

Detailed report:
`rawmu/rawmu_release_grounded_analysis_report_2026-07-18.md`.

## Claim Boundary and Next Gate

1. Do not promote the current Quaia angular mode: the independent catalogue
   fails direction/window coherence and the valid spectroscopic-redshift
   likelihood is underpowered at the locked signal.
2. Any next angular study should increase independent spectroscopic overlap
   or combine separately registered spectroscopic catalogues, then calibrate
   power before inspecting the observed score.
3. Keep released total covariance and `zHD` as the default supernova readout.
   A next FR comparison should improve the physical FR luminosity-distance
   prediction before adding nuisance flexibility, and must retain the same
   nuisance model and external calibration budgets for both cosmologies.
