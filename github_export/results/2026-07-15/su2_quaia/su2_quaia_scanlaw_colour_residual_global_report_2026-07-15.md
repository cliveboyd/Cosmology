# SU2 / Quaia Scan-Law + Colour Residual Global Look-Elsewhere Search

Date: July 15, 2026

## Purpose

This promotion gate asks whether any scanned redshift-window / Galactic-latitude-cut angular mode remains globally unusual after regressing Quaia redshift against Gaia scan-law geometry and catalogue-colour controls.

## Residual Model

```text
z_i      = alpha + S_i beta_s + C_i beta_c + epsilon_i
r_i      = z_i - fitted(z_i | S_i, C_i)
r_i      = a0 + d . n_i + eta_i
priority = residual_SNR_d * max_abs_delta_mu_SU2(window)
```

`S_i` contains Gaia scan count log1p, scan-angle cos2 mean, scan-angle sin2 mean, and scan-angle resultant. `C_i` contains Gaia BP-RP z-score, WISE W1-W2 z-score, and their cross term.

## Null

For each window/cut, residuals are shuffled within quantiles of the scan-law/colour fitted redshift. Each mock records the maximum statistic over all prepared windows.

## Bottom Line

A residual angular mode is globally unusual at the 1% scan gate and should be treated as a follow-up candidate only.

## Configuration

- prepared windows: `168`
- residual-shuffle mocks: `1000`
- strata per window: `10`
- seed: `190715`

## Mock Maxima

- max residual SNR mean / p95 / p99: `3.76301` / `4.47673` / `4.89758`
- max residual priority mean / p95 / p99: `0.0330881` / `0.0391785` / `0.0426926`

## Observed Thresholds

| role | tag | bcut_deg | N | residual_formal_snr | residual_priority | global_p_any_window_snr_ge_observed | global_p_any_window_priority_ge_observed | z1p00_family_p_any_bcut_snr_ge_observed | z1p00_family_p_any_bcut_priority_ge_observed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| observed_max_residual_snr | z1p00_1p50 | 10 | 194552 | 5.75359 | 0.0552361 | 0.000999001 | 0.000999001 | 0.000999001 | 0.000999001 |
| observed_max_residual_priority | z1p00_1p50 | 10 | 194552 | 5.75359 | 0.0552361 | 0.000999001 | 0.000999001 | 0.000999001 | 0.000999001 |
| z1p00_1p50_family_max_residual_snr | z1p00_1p50 | 10 | 194552 | 5.75359 | 0.0552361 | 0.000999001 | 0.000999001 | 0.000999001 | 0.000999001 |
| z1p00_1p50_family_max_residual_priority | z1p00_1p50 | 10 | 194552 | 5.75359 | 0.0552361 | 0.000999001 | 0.000999001 | 0.000999001 | 0.000999001 |
| z1p00_1p50_bcut_10 | z1p00_1p50 | 10 | 194552 | 5.75359 | 0.0552361 | 0.000999001 | 0.000999001 | 0.000999001 | 0.000999001 |
| z1p00_1p50_bcut_15 | z1p00_1p50 | 15 | 185739 | 5.68228 | 0.0545516 | 0.000999001 | 0.000999001 | 0.000999001 | 0.000999001 |
| z1p00_1p50_bcut_20 | z1p00_1p50 | 20 | 172054 | 4.58277 | 0.0439959 | 0.031968 | 0.002997 | 0.000999001 | 0.000999001 |
| z1p00_1p50_bcut_25 | z1p00_1p50 | 25 | 155609 | 4.29184 | 0.041203 | 0.0859141 | 0.021978 | 0.002997 | 0.002997 |
| z1p00_1p50_bcut_30 | z1p00_1p50 | 30 | 137850 | 2.9821 | 0.028629 | 0.999001 | 0.944056 | 0.18981 | 0.18981 |
| z1p00_1p50_bcut_35 | z1p00_1p50 | 35 | 119719 | 2.79086 | 0.0267931 | 1 | 0.997003 | 0.293706 | 0.293706 |

## Top Observed Residual Windows

| tag | bcut_deg | N | control_r2 | raw_formal_snr | residual_formal_snr | residual_priority | residual_l_deg | residual_b_deg | residual_sep_cmb_deg |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| z1p00_1p50 | 10 | 194552 | 0.135974 | 3.47955 | 5.75359 | 0.0552361 | 222.657 | -4.64068 | 64.0406 |
| z1p00_1p50 | 15 | 185739 | 0.138497 | 3.59778 | 5.68228 | 0.0545516 | 219.602 | -4.55171 | 65.4894 |
| z1p00_1p50 | 20 | 172054 | 0.139857 | 3.02083 | 4.58277 | 0.0439959 | 223.911 | -5.6127 | 64.2893 |
| z1p00_1p50 | 25 | 155609 | 0.141136 | 3.49857 | 4.29184 | 0.041203 | 227.5 | -3.95679 | 61.1618 |
| z1p50_2p50 | 10 | 275246 | 0.224318 | 2.66667 | 3.68628 | 0.032652 | 133.013 | -24.2538 | 134.815 |
| z1p90_2p00 | 35 | 18180 | 0.00214225 | 3.47876 | 3.39697 | 0.0303558 | 305.373 | -35.7866 | 91.767 |
| z1p00_1p50 | 30 | 137850 | 0.141786 | 3.38376 | 2.9821 | 0.028629 | 233.617 | -3.70101 | 58.3373 |
| z1p00_1p50 | 35 | 119719 | 0.142452 | 3.70594 | 2.79086 | 0.0267931 | 233.021 | -4.166 | 59.0006 |
| z1p90_2p00 | 20 | 25800 | 0.00146369 | 2.96668 | 2.89085 | 0.025833 | 257.163 | -50.476 | 98.9048 |
| z1p50_2p50 | 15 | 264161 | 0.22416 | 2.76208 | 2.8945 | 0.0256386 | 132 | -25.806 | 136.557 |
| z1p20_1p30 | 30 | 26478 | 0.0131769 | 2.83648 | 2.67002 | 0.025633 | 282.26 | -27.897 | 77.8885 |
| z1p90_2p00 | 30 | 20791 | 0.00172738 | 2.98274 | 2.8671 | 0.0256208 | 305.301 | -51.7604 | 106.039 |
| z1p90_2p00 | 15 | 27642 | 0.00164275 | 2.79567 | 2.64784 | 0.0236615 | 262.018 | -53.1988 | 101.466 |
| z1p90_2p00 | 10 | 28768 | 0.00161984 | 2.74888 | 2.59209 | 0.0231633 | 266.936 | -57.7709 | 106.051 |
| z1p90_2p00 | 25 | 23377 | 0.00171338 | 2.73637 | 2.56078 | 0.0228835 | 302.282 | -76.2159 | 126.872 |

## Claim Boundary

This is not a discovery run. Ordinary residual global p-values favour the scan-law/colour explanation. A surviving residual mode is a follow-up candidate only.

## Files

- observed windows: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_scanlaw_colour_residual_global_20260715\su2_quaia_scanlaw_colour_residual_global_observed_windows.csv`
- mock maxima: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_scanlaw_colour_residual_global_20260715\su2_quaia_scanlaw_colour_residual_global_mock_maxima.csv`
- thresholds: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\su2_quaia_scanlaw_colour_residual_global_20260715\su2_quaia_scanlaw_colour_residual_global_thresholds.csv`
