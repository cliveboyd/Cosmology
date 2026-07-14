# Raw-MU Calibration-Budget Default and zHD Frame Policy

Date: July 14, 2026

## Decision

The default supernova diagnostic is now the completed hierarchical raw-MU offset likelihood in the preregistered `HD` frame. The headline configuration is `dataset+survey` random effects, with the 25 mmag survey / 50 mmag dataset calibration-budget row as the reference case, and the 10, 16, and 50 mmag survey-budget rows retained as a visible sensitivity analysis.

The default redshift-frame readout is the completed preregistered `zHD` external-flow audit. `zHEL` and `zCMB` are retained as controls; they should not be promoted to headline frames unless a separate physical reason is specified before looking at the fit outcome.

## External Calibration-Budget Rationale

The sensitivity ladder is intentionally labelled as calibration budgets, not as a hidden fixed grid. The survey-level widths are 10, 16, 25, and 50 mmag; dataset-level widths are twice those values to cover between-compendium calibration/standardization mismatch. The 25 mmag survey-width reference is tied to the Pantheon+ inter-survey zeropoint sensitivity scale; 10 and 16 mmag are tighter calibration cases, while 50 mmag is a stress case.

Primary references used for this policy package:

- Pantheon+ calibration/zeropoint sensitivity: https://arxiv.org/abs/2110.03486
- DES-SN5YR DES-Dovekie data release: https://github.com/des-science/DES-SN5YR
- Pantheon+SH0ES data release: https://github.com/PantheonPlusSH0ES/DataRelease
- External-flow context for 2M++-style flow fields: https://arxiv.org/abs/1504.04627

## Headline Hierarchical Result

Reference row: `HD`, `dataset+survey`, `budget_025mmag_ds050`, `FR_BETA05_H0fixed675`.

- N = 3422
- calibration random effects = 27
- data chi2 at MAP = 3092.819
- calibration prior chi2 = 19.154
- marginal/BIC-like score = 3146.032

The matching beta-free H0-fixed row gives beta = 0.518263 and BIC-like score = 3151.260. The no-offset HD beta-free control gives beta = 0.466532 and BIC-like score = 4118.004, so the old no-offset readout is no longer an acceptable default diagnostic.

## HD Calibration-Budget Sensitivity

These rows keep H0 fixed at 67.5 and fit beta. The stable feature is the move from the no-offset beta near 0.467 to calibrated-offset beta near 0.515-0.519 for the dataset+survey budget ladder.

| budget                  | offset model   | beta     | beta-0.5  | BIC-like | prior chi2 |
| ----------------------- | -------------- | -------- | --------- | -------- | ---------- |
| no offsets              | none           | 0.466532 | -0.033468 | 4118.004 | 0.000      |
| 10 survey / 20 dataset  | dataset+survey | 0.514866 | 0.014866  | 3187.627 | 49.877     |
| 16 survey / 32 dataset  | dataset+survey | 0.519264 | 0.019264  | 3160.275 | 28.200     |
| 25 survey / 50 dataset  | dataset+survey | 0.518263 | 0.018263  | 3151.260 | 17.855     |
| 50 survey / 100 dataset | dataset+survey | 0.514357 | 0.014357  | 3157.093 | 8.984      |

## zHD Headline Frame, zHEL/zCMB Controls

At the reference 25 mmag survey / 50 mmag dataset budget with fixed beta and H0, `HD` is also the best-scoring frame among the audited choices. The controls remain useful precisely because they show how much the exact-log preference can move when the redshift basis is changed.

| frame             | N    | data chi2 | prior chi2 | BIC-like |
| ----------------- | ---- | --------- | ---------- | -------- |
| HD                | 3422 | 3092.819  | 19.154     | 3146.032 |
| CMB_PANTHEON_ONLY | 3420 | 3122.244  | 16.550     | 3172.834 |
| HEL               | 3418 | 3128.392  | 14.110     | 3176.537 |

Preregistered frame choices:

| dataset           | primary | controls  | status                                      |
| ----------------- | ------- | --------- | ------------------------------------------- |
| PantheonPlusSH0ES | zHD     | zHEL,zCMB | audited                                     |
| DES_Dovekie_raw   | zHD     | zHEL      | audited                                     |
| Union3p1_UNITY1p8 | z       | none      | compressed z only; not object-level audited |

Primary redshift-frame audit rows at z_min = 0.01:

| dataset           | z column | N    | log p  | LCDM Om | Delta AIC log-LCDM | AIC preference  |
| ----------------- | -------- | ---- | ------ | ------- | ------------------ | --------------- |
| DES_Dovekie_raw   | zHD      | 1820 | 1.3690 | 0.3287  | 1.483              | LCDM_Omfree     |
| DES_Dovekie_raw   | zHEL     | 1820 | 1.3250 | 0.3465  | -5.754             | EXACT_LOG_Pfree |
| PantheonPlusSH0ES | zCMB     | 1578 | 1.3445 | 0.3467  | -0.638             | EXACT_LOG_Pfree |
| PantheonPlusSH0ES | zHD      | 1580 | 1.3801 | 0.3316  | 3.202              | LCDM_Omfree     |
| PantheonPlusSH0ES | zHEL     | 1576 | 1.3175 | 0.3586  | -3.387             | EXACT_LOG_Pfree |

## Model-Flexibility Checks

The reference frame and calibration-budget row remain the basis for H0/beta flexibility checks; these are secondary diagnostics, not a frame-selection mechanism.

| model                  | H0     | beta     | data chi2 | prior chi2 | marginal score | BIC-like |
| ---------------------- | ------ | -------- | --------- | ---------- | -------------- | -------- |
| FR_BETA05_H0fixed675   | 67.500 | 0.500000 | 3092.819  | 19.154     | 3146.032       | 3146.032 |
| FR_BETA05_H0free       | 69.278 | 0.500000 | 3092.356  | 15.929     | 3142.345       | 3150.482 |
| FR_BETAfree_H0fixed675 | 67.500 | 0.518263 | 3091.208  | 17.855     | 3143.122       | 3151.260 |
| FR_BETAfree_H0free     | 69.702 | 0.523532 | 3091.173  | 12.507     | 3137.739       | 3154.015 |

## Largest Reference-Row Offset Pulls

These are diagnostics for calibration tension under the reference budget. They are not grounds to drop datasets or surveys unless a new exclusion rule is preregistered.

| offset group                     | type    | sigma mag | offset mag | pull sigma |
| -------------------------------- | ------- | --------- | ---------- | ---------- |
| PantheonPlusSH0ES:DATASET        | dataset | 0.050     | -0.1288    | -2.576     |
| PantheonPlusSH0ES:IDSURVEY_OTHER | survey  | 0.025     | -0.0409    | -1.635     |
| PantheonPlusSH0ES:IDSURVEY_15    | survey  | 0.025     | 0.0348     | 1.394      |
| DES_Dovekie_raw:IDSURVEY_10      | survey  | 0.025     | 0.0339     | 1.358      |
| PantheonPlusSH0ES:IDSURVEY_4     | survey  | 0.025     | 0.0315     | 1.261      |
| PantheonPlusSH0ES:IDSURVEY_150   | survey  | 0.025     | -0.0266    | -1.063     |
| PantheonPlusSH0ES:IDSURVEY_64    | survey  | 0.025     | -0.0241    | -0.965     |
| PantheonPlusSH0ES:IDSURVEY_10    | survey  | 0.025     | 0.0150     | 0.598      |
| DES_Dovekie_raw:IDSURVEY_63      | survey  | 0.025     | -0.0130    | -0.520     |
| PantheonPlusSH0ES:IDSURVEY_100   | survey  | 0.025     | 0.0121     | 0.485      |
| PantheonPlusSH0ES:IDSURVEY_65    | survey  | 0.025     | -0.0117    | -0.469     |
| DES_Dovekie_raw:IDSURVEY_65      | survey  | 0.025     | -0.0105    | -0.421     |

## Output Files

- `rawmu_calibration_budget_default_frame_summary_2026-07-14.csv`: top rows by BIC-like score.
- `rawmu_calibration_budget_hd_sensitivity_2026-07-14.csv`: HD beta-free H0-fixed sensitivity ladder.
- `rawmu_calibration_budget_control_frames_2026-07-14.csv`: HD, HEL, and CMB control-frame sensitivity rows.
- `rawmu_calibration_budget_frame_default_comparison_2026-07-14.csv`: 25 mmag reference frame comparison.
- `rawmu_calibration_budget_default_offsets_2026-07-14.csv`: offset pulls for the reference row.
- `rawmu_calibration_budget_full_summary_2026-07-14.csv`, `rawmu_calibration_budget_full_offsets_2026-07-14.csv`, `rawmu_calibration_budget_blocks_2026-07-14.csv`, and `rawmu_calibration_budget_full_hierarchical_report_2026-07-14.md`: complete underlying run products copied from the diagnostics directory.
