# SU2 / Quaia Gaia Scan-Law Follow-Up

Date: July 15, 2026

## Purpose

This follows the external-template gate where the Gaia scan-law controls weakened and rotated the locked `0.95 <= z < 1.45` angular mode.

## Tests

1. Scan-law preserving null mocks: shuffle redshifts inside bins of Gaia scan count, scan-angle resultant, and scan-angle phase.
2. Scan-law strata: fit the locked mode separately in scan-count and scan-angle-resultant quartiles.
3. Residual mode: regress out SFD dust, Quaia selection/depth, and Gaia scan-law templates before refitting the dipole.

## Scan-Law Preserving Null

| null_kind | n_mocks | n_scan_bins | observed_snr_max | p_snr_max_ge_observed | p_pair_sep_le_observed | p_coherence_ge_observed | p_joint_snr_and_pair_sep | mock_snr_max_p99 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| scanlaw_count_resultant_phase_zshuffle | 2000 | 63 | 4.7852 | 0.00949525 | 0.0809595 | 0.0069965 | 0.00349825 | 4.73562 |

- The locked mode remains rare under scan-law-preserving redshift shuffles.

## Scan-Law Strata Summary

| strata_variable | quartile | amp_mean | snr_max | max_pair_direction_sep_deg | N_min | amp_q4_over_q1 | max_fit_direction_sep_deg |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gaia_scan_count_log1p_dr3 | 1 | 0.00374179 | 3.56602 | 29.5412 | 30842 |  |  |
| gaia_scan_count_log1p_dr3 | 2 | 0.00244732 | 2.38604 | 56.5079 | 24031 |  |  |
| gaia_scan_count_log1p_dr3 | 3 | 0.00337874 | 2.54589 | 32.4933 | 29139 |  |  |
| gaia_scan_count_log1p_dr3 | 4 | 0.00421821 | 2.64099 | 16.6812 | 34181 |  |  |
| gaia_scan_count_log1p_dr3 | Q4/Q1 |  |  |  |  | 1.12732 | 106.058 |
| gaia_scan_angle_resultant_dr3 | 1 | 0.00223215 | 1.28731 | 67.8868 | 25810 |  |  |
| gaia_scan_angle_resultant_dr3 | 2 | 0.00438108 | 3.76392 | 21.9315 | 30421 |  |  |
| gaia_scan_angle_resultant_dr3 | 3 | 0.00160545 | 1.76318 | 70.1235 | 29056 |  |  |
| gaia_scan_angle_resultant_dr3 | 4 | 0.0065073 | 4.66115 | 11.2555 | 32906 |  |  |
| gaia_scan_angle_resultant_dr3 | Q4/Q1 |  |  |  |  | 2.91527 | 100.41 |

## Residual Mode Summary

| label | amp_ratio_vs_baseline | snr_max | snr_delta_vs_baseline | max_pair_direction_sep_deg | mean_resultant_length | N_min | N_max |
| --- | --- | --- | --- | --- | --- | --- | --- |
| baseline_no_residualization | 1 | 4.7852 | 0 | 13.4723 | 0.996066 | 118193 | 183527 |
| residual_after_sfd_selection_gaia | 0.805144 | 3.54648 | -1.23872 | 22.1118 | 0.989859 | 118193 | 183527 |
| residual_after_all_available_proxy_stack | 0.47531 | 2.8718 | -1.9134 | 35.4794 | 0.971522 | 118193 | 183527 |

- The residual-mode test still retains most of the baseline amplitude after SFD + selection/depth + Gaia scan-law regression.

## Interpretation

The promotion gate is not simply whether the signal survives dust and selection/depth. It must also survive an explicit Gaia scan-law null and retain a coherent residual mode after scan-law regression.

## Configuration

- mocks: `2000`
- redshift window: `0.95 <= z < 1.45`
- latitude cuts: `15.0, 20.0, 25.0, 30.0, 35.0`
- seed: `150715`
