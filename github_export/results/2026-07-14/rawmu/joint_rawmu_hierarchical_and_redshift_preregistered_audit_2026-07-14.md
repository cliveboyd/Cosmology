# Joint Raw-MU Hierarchical Offsets and Pre-Registered Redshift Audit

Date: 2026-07-14

This note packages the July 14 diagnostics for the two requested promotions:

1. Promote the joint DES-Dovekie, Pantheon+SH0ES, and Union3.1 raw-MU likelihood to a hierarchical offset model with explicit calibration-offset priors.
2. Repeat the multiplicative redshift-frame audit with a documented external-flow convention and pre-registered dataset redshift choices.

The underlying run products are in `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\joint_rawmu_hierarchical_offset_20260714` and `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\plamb_redshift_external_flow_preregistered_20260714`. Key files have been copied into this task's `outputs` directory.

## 1. Hierarchical Raw-MU Offset Likelihood

### Model

The promoted likelihood uses the same independent block-diagonal covariance/precision treatment as the July 14 joint raw-MU runner:

`residual = MU_observed - MU_model(H0, beta)`

`residual = X_dataset delta_dataset + X_survey delta_survey + noise`

with explicit zero-mean Gaussian priors:

`delta_dataset ~ N(0, sigma_dataset^2)`

`delta_survey ~ N(0, sigma_survey^2)`

Tested prior scales were:

| prior_config | sigma_dataset_mag | sigma_survey_mag |
|---|---:|---:|
| tight | 0.02 | 0.01 |
| moderate | 0.05 | 0.02 |
| loose | 0.10 | 0.05 |
| very_loose | 0.20 | 0.10 |

The joint model included Pantheon+SH0ES, DES-Dovekie raw, and Union3.1 UNITY1.8 compressed nodes. Union3.1 has no object-level survey labels, so it only receives a dataset-level random offset.

### Headline Rows

Rows below are sorted by the marginal random-effect BIC-like score. Lower is better.

| frame | random_mode | prior | model | N | k_random | H0 | beta | data_chi2 | prior_chi2 | marginal_BIC_like |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| HD | dataset+survey | moderate | FR_BETA05_H0fixed675 | 3422 | 27 | 67.5 | 0.500000 | 3097.754 | 19.907 | 3147.056 |
| HD | dataset+survey | loose | FR_BETA05_H0fixed675 | 3422 | 27 | 67.5 | 0.500000 | 3082.878 | 9.594 | 3150.488 |
| HD | dataset+survey | moderate | FR_BETAfree_H0fixed675 | 3422 | 27 | 67.5 | 0.521314 | 3095.494 | 18.011 | 3151.038 |
| HD | dataset+survey | moderate | FR_BETA05_H0free | 3422 | 27 | 69.226 | 0.500000 | 3097.263 | 16.891 | 3151.687 |
| HD | dataset+survey | moderate | FR_BETAfree_H0free | 3422 | 27 | 69.714 | 0.526515 | 3095.490 | 12.530 | 3153.690 |
| HD | dataset | loose | FR_BETAfree_H0fixed675 | 3422 | 3 | 67.5 | 0.541273 | 3128.862 | 2.261 | 3158.339 |

For comparison, the no-offset HD H0-fixed beta-free row gives beta `0.466532` and BIC-like score `4118.004`. The hierarchical calibration terms therefore remove a large amount of cross-dataset tension, but that improvement is explicitly prior-mediated rather than a free offset fit.

### Interpretation

The promoted model does what we wanted: it turns the raw-MU offset audit into an explicit calibration model with dataset and survey priors. Under the registered HD frame, moderate dataset+survey priors keep beta close to the exact-log value, with beta-free H0-fixed `beta = 0.521314`. The best-scoring row is still the fixed beta=0.5 model, but the evidence is sensitive to allowing survey-level calibration structure.

The largest offset pulls remain calibration-diagnostic rather than cosmology-diagnostic. Under moderate priors, the Pantheon+SH0ES dataset offset is about `-0.135 mag` in the HD beta-free H0-fixed row, a `-2.71 sigma` pull. Tight priors drive larger formal pulls, above `5 sigma` for the Pantheon dataset offset, which argues against treating the raw-MU shape preference as standalone cosmological evidence before the calibration tension is audited.

## 2. External-Flow Redshift Audit

### Pre-Registered Dataset Choices

The external-flow convention is: use each release's published Hubble-diagram redshift as the external peculiar-velocity corrected frame, then treat heliocentric or CMB-frame alternatives as controls.

| dataset | primary redshift | controls | status |
|---|---|---|---|
| PantheonPlusSH0ES | zHD | zHEL, zCMB | audited |
| DES_Dovekie_raw | zHD | zHEL | audited |
| Union3p1_UNITY1p8 | z | none | compressed z only |

The multiplicative blend was pre-specified in `1+z`:

`(1 + z_alpha) = (1 + zHEL) * ((1 + zHD) / (1 + zHEL))^alpha`

where `alpha=0` is zHEL and `alpha=1` is zHD.

### External-Flow Validation

Pantheon+ exposes VPEC, so the audit validates the release relation:

`VPEC_reconstructed_kms = c * ((1 + zCMB) / (1 + zHD) - 1)`

The reconstruction has correlation `0.999951` with the release VPEC column, median absolute reconstruction error `0.788 km/s`, and p90 absolute error `1.772 km/s`.

DES-Dovekie raw exposes zHEL, zHD, and MUERR_VPEC but not explicit object-level VPEC/RA/Dec in the local HD table, so the audit treats release zHD as the documented external-flow output and directly audits the zHEL-to-zHD multiplicative shift.

### Primary z_min = 0.01 Fits

Negative delta AIC favors the exact logarithmic kernel; positive favors LCDM.

| dataset | z_col | N | log_p | lcdm_Om | delta_AIC_log_minus_lcdm | preferred |
|---|---:|---:|---:|---:|---:|---|
| DES_Dovekie_raw | zHD | 1820 | 1.369024 | 0.328734 | 1.483130 | LCDM_Omfree |
| DES_Dovekie_raw | zHEL | 1820 | 1.325037 | 0.346543 | -5.753557 | EXACT_LOG_Pfree |
| PantheonPlusSH0ES | zCMB | 1578 | 1.344532 | 0.346665 | -0.637611 | EXACT_LOG_Pfree |
| PantheonPlusSH0ES | zHD | 1580 | 1.380136 | 0.331631 | 3.202212 | LCDM_Omfree |
| PantheonPlusSH0ES | zHEL | 1576 | 1.317479 | 0.358616 | -3.386572 | EXACT_LOG_Pfree |

### Multiplicative Alpha Crossings

| dataset | alpha where delta_AIC crosses zero |
|---|---:|
| PantheonPlusSH0ES | 0.473871 |
| DES_Dovekie_raw | 0.777855 |

This is the key result: as the redshift basis is moved from heliocentric to release zHD using the multiplicative rule, the exact-log preference smoothly disappears. Under the pre-registered primary frame, both Pantheon+ and DES-Dovekie favor LCDM by AIC. The exact-log preference is therefore redshift-frame sensitive, not robust under the external-flow convention.

## References Used For The External-Flow Convention

- Pantheon+SH0ES public distance table: https://github.com/PantheonPlusSH0ES/DataRelease/blob/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat
- Carrick et al. 2015, 2M++ peculiar-velocity field: https://arxiv.org/abs/1504.04627
- DES-SN5YR / DES-Dovekie public release repository: https://github.com/des-science/DES-SN5YR

## Output Artifacts

- `joint_rawmu_hierarchical_offset_report_2026-07-14.md`
- `joint_rawmu_hierarchical_offset_summary_2026-07-14.csv`
- `joint_rawmu_hierarchical_offset_pulls_2026-07-14.csv`
- `joint_rawmu_hierarchical_best_rows_2026-07-14.csv`
- `external_flow_preregistration_2026-07-14.json`
- `preregistered_dataset_frame_choices_2026-07-14.csv`
- `external_flow_model_validation_2026-07-14.csv`
- `plamb_redshift_external_flow_audit_2026-07-14.md`
- `plamb_redshift_primary_frame_summary_2026-07-14.csv`
- `multiplicative_redshift_alpha_crossings_2026-07-14.csv`

