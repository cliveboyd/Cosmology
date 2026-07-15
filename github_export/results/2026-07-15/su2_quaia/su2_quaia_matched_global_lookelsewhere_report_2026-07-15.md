# SU2 / Quaia Matched-Quality Global Look-Elsewhere Audit

Date: July 15, 2026

## Purpose

This test repeats the all-window SU2/Quaia look-elsewhere scan after all-catalogue-quality overlap weighting. It asks whether the targeted z ~ 1-1.5 mode remains globally unusual once the catalogue-quality mechanism gate is imposed across the whole scan.

## Model

```text
z_i      = c0 + d . n_i + epsilon_i
w_i      = overlap propensity weight for all catalogue-quality variables
SNR_d    = sqrt(d^T Cov(d)^(-1) d)
priority = SNR_d * max_abs_delta_mu_SU2(window)
```

## Configuration

- scan windows prepared: `168`
- redshift-shuffle mocks: `200`
- propensity strata: `10`

## Mock Maxima

- max matched SNR mean / p95 / p99: `10.4253` / `11.6319` / `11.9651`
- max matched priority mean / p95 / p99: `0.0951951` / `0.107132` / `0.11327`

## Observed Matched Thresholds

| role | tag | bcut_deg | N | ess | weighted_snr | weighted_priority | global_p_any_window_snr_ge_observed | global_p_any_window_priority_ge_observed | z1p00_family_p_any_bcut_snr_ge_observed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| matched_observed_max_snr+matched_observed_max_priority | z1p90_2p00 | 35 | 18180 | 18043.2 | 3.5321 | 0.0315634 | 1 | 1 |  |
| target | z1p00_1p50 | 10 | 194552 | 194415 | 0.705141 | 0.00676955 | 1 | 1 | 1 |
| target | z1p00_1p50 | 15 | 185739 | 185585 | 1.46456 | 0.0140602 | 1 | 1 | 1 |
| target | z1p00_1p50 | 25 | 155609 | 155398 | 1.95877 | 0.0188047 | 1 | 1 | 1 |
| target | z1p00_1p50 | 35 | 119719 | 119345 | 1.84756 | 0.0177371 | 1 | 1 | 1 |

## Top Matched Observed Windows

| tag | bcut_deg | N | ess | weighted_snr | weighted_priority | amp_ratio_weighted_vs_baseline | direction_sep_weighted_vs_baseline_deg |
| --- | --- | --- | --- | --- | --- | --- | --- |
| z1p90_2p00 | 35 | 18180 | 18043.2 | 3.5321 | 0.0315634 | 1.03833 | 1.68929 |
| z1p90_2p00 | 30 | 20791 | 20553.4 | 3.06786 | 0.0274148 | 1.03045 | 0.945282 |
| z1p90_2p00 | 20 | 25800 | 25592.4 | 3.05483 | 0.0272984 | 1.03301 | 1.65294 |
| z1p90_2p00 | 10 | 28768 | 28493.7 | 2.96097 | 0.0264597 | 1.08391 | 4.09115 |
| z1p90_2p00 | 15 | 27642 | 27398.5 | 2.95715 | 0.0264256 | 1.06267 | 4.54585 |
| z0p60_0p70 | 10 | 23898 | 23724.9 | 3.48602 | 0.0259756 | 1.09078 | 1.39785 |
| z1p90_2p00 | 25 | 23377 | 23123.1 | 2.89735 | 0.0258912 | 1.06983 | 2.36692 |
| z0p80_0p90 | 20 | 28027 | 27982.6 | 2.94153 | 0.0257905 | 1.06152 | 0.423149 |
| z1p20_1p30 | 20 | 33196 | 32897.4 | 2.634 | 0.0252872 | 0.956533 | 2.28122 |
| z1p70_1p80 | 10 | 34300 | 34101.6 | 2.7263 | 0.0251661 | 1.01833 | 1.505 |
| z1p20_1p30 | 10 | 37668 | 37319.9 | 2.6061 | 0.0250193 | 1.01367 | 2.74299 |
| z1p20_1p30 | 30 | 26478 | 26226.3 | 2.56149 | 0.024591 | 0.936543 | 8.84713 |

## Interpretation Gate

- Passing condition: the locked z ~ 1-1.5 family remains small after global look-elsewhere correction under the matched-quality null.
- Failing condition: global p-values are ordinary or the matched weighted mode is no longer the scan maximum.
- This is a promotion gate, not a discovery test by itself.
