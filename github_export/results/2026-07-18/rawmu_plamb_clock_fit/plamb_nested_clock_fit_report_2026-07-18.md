# PLAMB nested clock-law fit

**Analysis date:** 2026-07-18  
**Protocol:** outcome-blind flexible-law protocol written before execution  
**Primary sample:** 3,422 supernovae, released total covariance, one profiled intercept per release

## Headline result

The preregistered best PLAMB-family cell is `PATH_FREE`. Its full-sample Delta BIC relative to flat Lambda-CDM is 13.345052. This does not clear the registered in-sample BIC gate.

The fit is a phenomenological curve test. Supernovae constrain the combined flux exponent alpha, but do not separately identify photon-energy evolution, time dilation and intrinsic standardised-luminosity evolution.

The fixed Peter-law Delta BIC reproduces the locked prior result to better than 1e-12, providing an exact regression check on the likelihood implementation.

## Primary model ladder

| Model | chi-squared | Shape k | BIC | Delta BIC vs LCDM | Parameters |
|---|---:|---:|---:|---:|---|
| `LCDM` | 3047.662517 | 1 | 3080.214439 | +0.000000 | `{"Omega_m": 0.33067962104536625}` |
| `PATH_FREE` | 3052.869589 | 2 | 3093.559491 | +13.345052 | `{"alpha": 0.0, "gamma_c": 2.0, "p": 0.6281917761590368}` |
| `GENERAL_FREE` | 3047.174295 | 3 | 3096.002177 | +15.787738 | `{"alpha": 0.7015657756158429, "gamma_c": 1.8319888929719106, "p": 1.753281643154555}` |
| `P_ALPHA_FREE` | 3063.499775 | 2 | 3104.189677 | +23.975238 | `{"alpha": 1.1409362553052282, "gamma_c": 1.0, "p": 2.000000049458507}` |
| `GAMMA_FREE` | 3128.846794 | 1 | 3161.398716 | +81.184277 | `{"alpha": 0.0, "gamma_c": 1.0845990114931208, "p": 0.0}` |
| `ALPHA_FREE` | 3131.423438 | 1 | 3163.975360 | +83.760921 | `{"alpha": 0.037604901114341276, "gamma_c": 1.0, "p": 0.0}` |
| `P_FREE` | 3136.721888 | 1 | 3169.273810 | +89.059371 | `{"alpha": 0.0, "gamma_c": 1.0, "p": -0.05256061428060205}` |
| `PETER_FIXED` | 3150.145481 | 0 | 3174.559422 | +94.344983 | `{"alpha": 0.0, "gamma_c": 1.0, "p": 0.0}` |
| `FRACTIONAL_FIXED` | 8081.254481 | 0 | 8105.668422 | +5025.453983 | `{"alpha": 0.0, "gamma_c": 1.0, "p": 1.0}` |

BIC counts the three release intercepts in every cell and adds the listed shape parameters. The common intercept count therefore cancels in model differences.

## Component attribution

| Release from fixed Peter law | Delta chi-squared | Delta BIC |
|---|---:|---:|
| `P_FREE` | -13.423593 | -5.285612 |
| `ALPHA_FREE` | -18.722043 | -10.584062 |
| `GAMMA_FREE` | -21.298687 | -13.160707 |
| `PATH_FREE` | -97.275892 | -80.999931 |
| `P_ALPHA_FREE` | -86.645706 | -70.369745 |
| `GENERAL_FREE` | -102.971186 | -78.557245 |

Among one-component releases, freeing gamma_c gives the largest improvement, followed by alpha and then p. None is sufficient by itself. Joint gamma_c-p path freedom removes most of the fixed-law discrepancy, but selects the gamma_c upper bound. The fully general curve lowers chi-squared by only 0.488223 relative to Lambda-CDM while using two additional shape parameters, so its BIC remains 15.787738 higher.

## Identifiability audit

- Selected-candidate central Hessian available: False.
- Selected-candidate Hessian rank: None of 2.
- Hessian condition number: None.
- Minimum fractional distance from a registered bound: 0.0.
- Boundary fractions: `{"gamma_c": 0.0, "p": 0.3760639253863456}`.
- Approximate standard errors: `{}`.

| Model | Hessian rank | Condition number | Minimum bound fraction |
|---|---:|---:|---:|
| `P_FREE` | 1 / 1 | 1.0 | 0.149146461906466 |
| `ALPHA_FREE` | 1 / 1 | 1.0 | 0.2688024505571706 |
| `GAMMA_FREE` | 1 / 1 | 1.0 | 0.4577004942534396 |
| `PATH_FREE` | None / 2 | None | 0.0 |
| `P_ALPHA_FREE` | 2 / 2 | 1863.609575559605 | 0.16666665018049764 |
| `GENERAL_FREE` | 3 / 3 | 44610.39022753036 | 0.0840055535140447 |
| `LCDM` | 1 / 1 | 1.0 | 0.4896734162811523 |

The general curve is formally full rank, but p and alpha have correlation 0.997363; their local-Hessian standard errors are large. This is a strong clock-rate/flux degeneracy, not a clean measurement of either mechanism.

## Hold-out readout

Complete-release shape wins over Lambda-CDM: 0/3.
Survey plus high-redshift wins: 10/20 (0.500).
Largest held-out loss: 321.926337 chi-squared.

| Type | Held-out set | N test | Delta chi-squared candidate - LCDM | Candidate wins |
|---|---|---:|---:|---|
| release_shape | `PantheonPlusSH0ES` | 1580 | +3.367234 | False |
| release_shape | `DES_Dovekie_raw` | 1820 | +1.658285 | False |
| release_shape | `Union3p1_UNITY1p8` | 22 | +0.423853 | False |
| survey_conditional | `PantheonPlusSH0ES:survey=1` | 321 | -0.837224 | True |
| survey_conditional | `PantheonPlusSH0ES:survey=4` | 160 | +1.625264 | False |
| survey_conditional | `PantheonPlusSH0ES:survey=5` | 74 | +0.241678 | False |
| survey_conditional | `PantheonPlusSH0ES:survey=10` | 203 | -1.724638 | True |
| survey_conditional | `PantheonPlusSH0ES:survey=15` | 269 | +3.956694 | False |
| survey_conditional | `PantheonPlusSH0ES:survey=50` | 27 | -0.523616 | True |
| survey_conditional | `PantheonPlusSH0ES:survey=51` | 38 | -0.951767 | True |
| survey_conditional | `PantheonPlusSH0ES:survey=56` | 36 | +0.103116 | False |
| survey_conditional | `PantheonPlusSH0ES:survey=57` | 85 | -0.356875 | True |
| survey_conditional | `PantheonPlusSH0ES:survey=63` | 28 | -0.345612 | True |
| survey_conditional | `PantheonPlusSH0ES:survey=64` | 53 | +0.433904 | False |
| survey_conditional | `PantheonPlusSH0ES:survey=65` | 35 | +0.162652 | False |
| survey_conditional | `PantheonPlusSH0ES:survey=150` | 173 | +2.101961 | False |
| survey_conditional | `DES_Dovekie_raw:survey=10` | 1623 | -16.066281 | True |
| survey_conditional | `DES_Dovekie_raw:survey=64` | 31 | -1.352565 | True |
| survey_conditional | `DES_Dovekie_raw:survey=65` | 22 | -0.405605 | True |
| survey_conditional | `DES_Dovekie_raw:survey=150` | 117 | -2.183650 | True |
| high_z_conditional | `z>=0.5` | 1085 | +321.926337 | False |
| high_z_conditional | `z>=0.8` | 167 | +20.055476 | False |
| high_z_conditional | `z>=1.0` | 46 | +25.498202 | False |

Survey and redshift tests use the exact Gaussian conditional covariance. Complete-release tests profile a test-release intercept and therefore test shape transfer, not absolute cross-release calibration. These are point-predictive scores: fitted-parameter uncertainty is not integrated into the conditional likelihood.

## Calibration readout

Maximum candidate-versus-Lambda-CDM release-intercept difference: 0.018000 mag.

## Flux-exponent interpretation

The fully general diagnostic gives alpha = 0.701566 with local-Hessian standard error 0.803913. With the external DES time-dilation summary b = 1.003 +/- 0.011180, this implies:

\[
\begin{aligned}
e\; (s=0) &= 2\alpha-b = 0.400132 \pm 1.607866,\\
s\; (e=1) &= 1+b-2\alpha = 0.599868 \pm 1.607866.
\end{aligned}
\]

These are assumption-conditional translations, not separate measurements of e or s. The time-dilation input is the DES supernova result reported in [arXiv:2406.05050](https://arxiv.org/abs/2406.05050).

## Registered promotion gates

| Gate | Pass |
|---|---|
| `primary_delta_BIC_below_minus_10` | **False** |
| `no_shape_parameter_within_one_percent_of_bound` | **False** |
| `full_rank_hessian_condition_below_1e8` | **False** |
| `wins_all_three_complete_release_shape_holdouts` | **False** |
| `wins_at_least_75_percent_survey_and_redshift_holdouts` | **False** |
| `no_single_holdout_loss_above_delta_chi2_10` | **False** |
| `release_intercept_differences_within_0p050_mag` | **True** |
| `external_e_b_s_physical_identification` | **False** |

**Overall promotion decision: DO NOT PROMOTE.**

The external physical-identification gate is deliberately false unless an independent model or measurement fixes enough of e, b and s to identify the claimed mechanism.

## Reproduction

```powershell
python github_export/code/rawmu/run_plamb_nested_clock_fit_2026_07_18.py
```

The CSV tables, JSON summary, plot and SHA-256 manifest in this directory are generated by that command.
