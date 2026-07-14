# SU2 / Quaia Targeted Stability Mock Percentiles

Date: July 14, 2026

## Purpose

This run checks whether the observed direction/amplitude stability in the two targeted follow-up windows is unusual against nested random-sky mocks.

For each target window, the mock draws random Quaia sky positions at the loosest latitude cut and evaluates the stricter latitude cuts as nested subsets. This preserves the latitude-cut correlation structure needed for a stability test.

## Configuration

- n_mocks: `1500`
- latitude cuts: `10, 15, 20, 25, 30, 35`
- seed: `290714`
- Quaia table: `C:\Users\clive\Documents\Cosmology\quaia_G20p0_basic_gal.csv`
- random sky catalog: `C:\Users\clive\Documents\Cosmology\quaia_G20p0_randoms_simple.fits`

## Observed Stability Metrics

| target | z_lo | z_hi | amp_mean | snr_max | max_pair_direction_sep_deg | mean_resultant_length | coherence_score |
| --- | --- | --- | --- | --- | --- | --- | --- |
| broad_0p95_1p45 | 0.95 | 1.45 | 0.00302017 | 4.7852 | 13.4723 | 0.996101 | 3.66831 |
| narrow_1p15_1p35 | 1.15 | 1.35 | 0.00127717 | 3.29694 | 18.5137 | 0.991724 | 2.31658 |

## Mock Percentiles

| target | n_mocks | observed_snr_max | p_snr_max_ge_observed | observed_max_pair_direction_sep_deg | p_pair_sep_le_observed | observed_coherence_score | p_coherence_ge_observed | p_joint_snr_and_pair_sep |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| broad_0p95_1p45 | 1500 | 4.7852 | 0.000666223 | 13.4723 | 0.00799467 | 3.66831 | 0.000666223 | 0.000666223 |
| narrow_1p15_1p35 | 1500 | 3.29694 | 0.0439707 | 18.5137 | 0.0166556 | 2.31658 | 0.0153231 | 0.00532978 |

## Interpretation Gate

- A targeted window is interesting only if its amplitude/SNR and direction clustering are unusual together.
- Small `p_pair_sep_le_observed` means the observed directions are unusually tightly clustered across latitude cuts.
- Small `p_joint_snr_and_pair_sep` means the mock rarely gets both comparable SNR and comparable direction clustering.
- This remains exploratory because the full global look-elsewhere gate did not pass.

## Outputs

- observed fits: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_targeted_stability_mocks_20260714\su2_quaia_targeted_stability_observed_fits.csv`
- observed summary: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_targeted_stability_mocks_20260714\su2_quaia_targeted_stability_observed_summary.csv`
- mock summary: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_targeted_stability_mocks_20260714\su2_quaia_targeted_stability_mock_summary.csv`
- p-values: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_targeted_stability_mocks_20260714\su2_quaia_targeted_stability_pvalues.csv`
- config: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_targeted_stability_mocks_20260714\su2_quaia_targeted_stability_config.json`
