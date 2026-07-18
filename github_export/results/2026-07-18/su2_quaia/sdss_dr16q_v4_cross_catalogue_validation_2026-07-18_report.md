# SDSS DR16Q v4 Independent Quaia Cross-Catalogue Validation

Date: 2026-07-18

## Decision

**NOT PROMOTED.** Promotion is conjunctive; 4 registered gate component(s) failed.

## Explicit Failures

- cross_catalogue_direction: primary_row failed (116.205; required <= 30.0 deg)
- window_direction_stability: z0p95_1p45 failed (165.587; required <= 30.0 deg)
- window_direction_stability: z1p00_1p40 failed (141.879; required <= 30.0 deg)
- amplitude_stability: z1p00_1p40_b35 failed (0.106361; required in [0.5, 2.0])

## Primary Row

| window_tag | bcut_deg | N | dipole_amplitude | dipole_ra_deg | dipole_dec_deg | direction_separation_from_locked_quaia_deg | partial_F_l1 | partial_F_l2 | partial_F_l3 | joint_F_l1_l3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| z1p00_1p50 | 35 | 139305 | 0.0481791 | 59.5866 | 84.2057 | 116.205 | 1.01554 | 0.859078 | 0.795405 | 5.17216 |

## Empirical Global Null

| cell_deg | n_mocks | n_occupied_catalogue_cells | observed_family_max_score | null_max_score_p99 | global_p |
| --- | --- | --- | --- | --- | --- |
| 8 | 20000 | 384 | 30.5982 | 8.15621 | 4.99975e-05 |
| 12 | 20000 | 203 | 30.5982 | 11.0279 | 4.99975e-05 |
| 16 | 20000 | 126 | 30.5982 | 14.6289 | 4.99975e-05 |

## Primary Mask-Calibrated Probabilities

| cell_deg | statistic | observed_F | point_p | exceedances | n_mocks |
| --- | --- | --- | --- | --- | --- |
| 8 | partial_F_l1 | 1.01554 | 0.630318 | 12606 | 20000 |
| 8 | partial_F_l2 | 0.859078 | 0.79236 | 15847 | 20000 |
| 8 | partial_F_l3 | 0.795405 | 0.865407 | 17308 | 20000 |
| 8 | joint_F_l1_l3 | 5.17216 | 4.99975e-05 | 0 | 20000 |
| 12 | partial_F_l1 | 1.01554 | 0.691515 | 13830 | 20000 |
| 12 | partial_F_l2 | 0.859078 | 0.846658 | 16933 | 20000 |
| 12 | partial_F_l3 | 0.795405 | 0.899655 | 17993 | 20000 |
| 12 | joint_F_l1_l3 | 5.17216 | 9.9995e-05 | 1 | 20000 |
| 16 | partial_F_l1 | 1.01554 | 0.741463 | 14829 | 20000 |
| 16 | partial_F_l2 | 0.859078 | 0.906355 | 18127 | 20000 |
| 16 | partial_F_l3 | 0.795405 | 0.946503 | 18930 | 20000 |
| 16 | joint_F_l1_l3 | 5.17216 | 0.00109995 | 21 | 20000 |

## Interpretation

The global family winner is `z0p95_1p45`, `|b| >= 35 deg`, `joint_F_l1_l3`. It is not the locked primary row.

The primary SDSS dipole amplitude is 15.8 times the locked Quaia amplitude. This cross-catalogue ratio is descriptive because the registered amplitude gate compares the fixed SDSS perturbation rows with the SDSS primary row.

Although the registered primary joint l=1-3 coherence statistic passes at every block scale, the individual partial l=1, l=2 and l=3 empirical probabilities span 0.6303 to 0.9465. The joint significance therefore does not show that the three multipoles are separately unusual; footprint-coupled collective structure remains a viable explanation.

## Promotion Gates

| gate | component | value | threshold | pass | failure |
| --- | --- | --- | --- | --- | --- |
| integrity | all_registered_rows | 20 | exactly 20 valid rows | True |  |
| global_rarity | cell_8_deg | 4.99975e-05 | <= 0.01 | True |  |
| global_rarity | cell_12_deg | 4.99975e-05 | <= 0.01 | True |  |
| global_rarity | cell_16_deg | 4.99975e-05 | <= 0.01 | True |  |
| joint_l1_l3_coherence | cell_8_deg | 4.99975e-05 | <= 0.05 | True |  |
| joint_l1_l3_coherence | cell_12_deg | 9.9995e-05 | <= 0.05 | True |  |
| joint_l1_l3_coherence | cell_16_deg | 0.00109995 | <= 0.05 | True |  |
| cross_catalogue_direction | primary_row | 116.205 | <= 30.0 deg | False | cross_catalogue_direction: primary_row failed (116.205; required <= 30.0 deg) |
| window_direction_stability | z0p95_1p45 | 165.587 | <= 30.0 deg | False | window_direction_stability: z0p95_1p45 failed (165.587; required <= 30.0 deg) |
| window_direction_stability | z1p05_1p55 | 6.40657 | <= 30.0 deg | True |  |
| window_direction_stability | z1p00_1p40 | 141.879 | <= 30.0 deg | False | window_direction_stability: z1p00_1p40 failed (141.879; required <= 30.0 deg) |
| window_direction_stability | z1p10_1p50 | 1.09823 | <= 30.0 deg | True |  |
| latitude_direction_stability | bcut_20 | 2.31356 | <= 30.0 deg | True |  |
| latitude_direction_stability | bcut_25 | 2.72728 | <= 30.0 deg | True |  |
| latitude_direction_stability | bcut_30 | 0.428611 | <= 30.0 deg | True |  |
| amplitude_stability | z0p95_1p45_b35 | 1.09677 | in [0.5, 2.0] | True |  |
| amplitude_stability | z1p05_1p55_b35 | 1.28773 | in [0.5, 2.0] | True |  |
| amplitude_stability | z1p00_1p40_b35 | 0.106361 | in [0.5, 2.0] | False | amplitude_stability: z1p00_1p40_b35 failed (0.106361; required in [0.5, 2.0]) |
| amplitude_stability | z1p10_1p50_b35 | 1.0678 | in [0.5, 2.0] | True |  |
| amplitude_stability | z1p00_1p50_b20 | 1.14064 | in [0.5, 2.0] | True |  |
| amplitude_stability | z1p00_1p50_b25 | 1.16044 | in [0.5, 2.0] | True |  |
| amplitude_stability | z1p00_1p50_b30 | 1.07299 | in [0.5, 2.0] | True |  |
| amplitude_stability | z1p00_1p50_b35 | 1 | in [0.5, 2.0] | True |  |

## Quality Controls

These predeclared controls are diagnostic and are not additional promotion gates.

| quality_control | N | dipole_amplitude | amplitude_ratio_to_primary | dipole_ra_deg | dipole_dec_deg | direction_separation_from_primary_deg |
| --- | --- | --- | --- | --- | --- | --- |
| zwarning_zero | 109416 | 0.09426 | 1.95645 | 347.004 | 74.217 | 15.0716 |
| zconf_ge_2 | 26740 | 0.0807934 | 1.67694 | 293.641 | -61.4399 | 154.444 |

## Method

The analysis conditions on the observed SDSS positions, jointly fits the registered real l=1-3 harmonic design and uses restricted-residual Rademacher signs shared by equal-area sky blocks across all overlapping rows. The empirical maximum over the fixed 80-statistic family supplies the global probability independently at 8, 12 and 16 degrees.

## Catalogue Audit

```json
{
  "fits_path": "C:\\Users\\clive\\Documents\\Cosmology\\external_datasets\\sdss_dr16q_v4\\DR16Q_v4.fits",
  "fits_bytes": 1718484480,
  "fits_rows": 750414,
  "primary_quality_rows_before_deduplication": 749749,
  "primary_quality_rows_after_deduplication": 749749,
  "duplicate_rows_removed": 0,
  "finite_z_range": [
    -999.0,
    7.023917198181152
  ]
}
```

## Claim Boundary

This is an independent catalogue validation of a fixed angular target. A non-promotion result remains an explicit catalogue-validation failure and must not be reframed as positive SU(2) evidence through exploratory window, mask, block-scale or coordinate changes.
