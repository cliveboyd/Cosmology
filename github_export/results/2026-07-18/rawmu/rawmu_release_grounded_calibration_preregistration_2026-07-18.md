# Release-Grounded Raw-MU Calibration Preregistration

Date: 18 July 2026

Status at registration: methods, models, frames, nuisance treatment, hold-outs and decision gates fixed before inspecting any new fit outcomes from this analysis.

## Scope and Claim Boundary

This analysis replaces the heuristic 10/16/25/50 mmag ladder used on 14 July with release-grounded covariance projections and documented zero-point budgets. It is a model-comparison diagnostic, not a measurement of `H0` and not a discovery claim.

The primary redshift frame is `zHD` for Pantheon+SH0ES and DES-Dovekie. Union3.1 exposes only its compressed release redshift `z`. `zHEL` and `zCMB` are controls and are outside the primary run.

## Locked Data Products

| Dataset | Distance product | Primary covariance | Calibration information |
| --- | --- | --- | --- |
| Pantheon+SH0ES | `Pantheon+SH0ES.dat`, `MU_SH0ES`, non-calibrators, `zHD > 0.01` | `Pantheon+SH0ES_STAT+SYS.cov` | release `STATONLY` and calibration-only `Pantheon+SH0ES_122221_CALIB.cov` products |
| DES-Dovekie | `DES-Dovekie_HD.csv`, `MU`, `zHD > 0.01` | `STAT+SYS.npz` | `STATONLY.npz`; their covariance difference is all-systematics, not calibration-only |
| Union3.1 | UNITY1.8 compressed 22-node FITS product | precision block embedded in FITS | calibration derivatives and priors are already propagated by UNITY; no local calibration-only block |

Source metadata and URLs are fixed in the companion JSON. File hashes and realised dimensions will be recorded by the prior-derivation program.

The official `2_CALIBRATION/FRAGILISTIC_COVARIANCE.npz` calibration-parameter covariance is archived and hashed. It is not mapped directly to SN-distance offsets without the corresponding release Jacobian; distance-covariance grouping products define this analysis.

## Calibration Hierarchy

Release normalisation and calibration uncertainty are distinct nuisance layers. DES explicitly states that its released distance moduli assume a conventional normalisation and that the SN likelihood marginalises the absolute magnitude. Consequently, the primary comparison profiles one unpenalised intercept per release. These intercepts are not interpreted as calibration shifts and are not tested against zero-point budgets.

For a group indicator `x_g`, define `w_g = x_g / sum(x_g)`. The common-mode covariance projection and residual external budget are

```math
\begin{aligned}
\sigma_{\mathrm{proj},g}(C) &= \sqrt{w_g^{\mathsf T} C w_g}, \\
C_{\mathrm{sys}}             &= C_{\mathrm{tot}}-C_{\mathrm{stat}}, \\
\sigma_{\mathrm{add},g}      &= \sqrt{\max\!\left(\sigma_{\mathrm{doc},g}^{2}-\sigma_{\mathrm{included},g}^{2},0\right)}.
\end{aligned}
```

The locked interpretation is:

- Pantheon+: project the release calibration-only covariance for dataset and `IDSURVEY` groups. These are release-internal calibration priors already represented in `STAT+SYS`; they are reported but not added again.
- DES-Dovekie: project `C_sys` only as an upper bound because it contains all systematics. It is not promoted to a calibration-only Gaussian prior.
- Union3.1: use the documented 0.01 mag scale for important dataset zero points as metadata, not as an extra prior, because the released compressed covariance already propagates calibration derivatives and uncertainties.
- An additive Gaussian residual prior is permitted only when a documented external component is demonstrably absent from the covariance. Otherwise `sigma_add = 0` and `prior_action = covariance_contained`.

The explicit hierarchy uses a dataset common mode plus survey residual common modes. The dataset variance is the mean cross-survey covariance; each survey variance is the non-negative remainder of its projected common-mode variance after subtracting the dataset mode. Small surveys are pooled into `OTHER` for hierarchy construction.

Locked reconstruction sensitivities are:

- Pantheon+ stat-only grouping reconstruction: each official grouping file includes the statistical covariance, so define `Delta C_group = C_group - C_STATONLY`; use `C_base = C_STATONLY`, then add grouped common modes derived from `Delta C_group`. `CALIB` is primary; repeat separately for `CSP_RECAL`, `HSTCALSPEC`, and their summed increments. Compare each reconstruction with both its official grouping covariance and released `C_total`.
- DES stat-only hierarchy: use `C_base = inv(P_STATONLY)` and grouped caps derived from `C_sys = inv(P_STAT+SYS) - inv(P_STATONLY)`. Because `C_sys` contains all systematics, this is an all-systematics common-mode sensitivity, not a calibration-only reconstruction.
- Union3.1: no explicit reconstruction unless separable stat-only and calibration components become available; retain the released total precision.

These stat-only reconstruction sensitivities never supplement `STAT+SYS`. A separate residual-offset sensitivity may supplement `STAT+SYS` only for a demonstrably external, non-overlapping `sigma_add` component. Every sensitivity must use identical nuisance parameterisation and priors for FR and Lambda-CDM, report covariance reconstruction error, and cannot supersede the primary `STAT+SYS` result or be counted as independent information.

## Same-Nuisance Models

Both models use the same rows, covariance, offset design matrix, offset parameter count, offset prior widths, profiling or marginalisation rule and hold-out scoring in every comparison cell. No model-specific offset is permitted. A fixed `H0 = 67.5 km s^-1 Mpc^-1` is retained for continuity, but the profiled release intercepts remove SN-only sensitivity to this normalisation.

The primary cell uses covariance-contained calibration and one flat profiled normalisation intercept per release. A residual-offset sensitivity cell may add only the non-overlapping `sigma_add` component defined above, with the identical Gaussian residual-offset design and priors for FR and Lambda-CDM. This sensitivity is correlated with, and subordinate to, the primary release likelihood; it is not new independent information.

```math
\begin{aligned}
D_{L,\mathrm{FR}}(z)        &= \frac{c}{H_0}z\left(1+\frac{z}{2}\right), \\
E_{\Lambda\mathrm{CDM}}(z) &= \sqrt{\Omega_m(1+z)^3+1-\Omega_m}, \\
D_{L,\Lambda\mathrm{CDM}}  &= \frac{c}{H_0}(1+z)\int_0^z\frac{\mathrm{d}z'}{E_{\Lambda\mathrm{CDM}}(z')}, \\
\mu(z)                      &= 5\log_{10}\!\left(D_L/\mathrm{Mpc}\right)+25.
\end{aligned}
```

The FR shape has no fitted cosmological parameter. Flat Lambda-CDM fits `Omega_m` on `[0.05, 0.60]`. For residual vector `r`, covariance precision `P` and release-intercept design `X`, the profiled intercepts and score are

```math
\begin{aligned}
\widehat a &= (X^{\mathsf T}PX)^{-1}X^{\mathsf T}Pr, \\
\chi^2     &= (r-X\widehat a)^{\mathsf T}P(r-X\widehat a), \\
\mathrm{BIC} &= \chi^2+k\ln N.
\end{aligned}
```

Define `Delta_BIC = BIC_FR - BIC_LCDM`; negative values favour FR.

## Locked Hold-Out Matrix

1. Full primary fit: all three releases, full release `STAT+SYS` covariance, one profiled intercept per release and no re-added release-internal calibration prior.
2. Dataset hold-outs: omit each of Pantheon+, DES-Dovekie and Union3.1 in turn. Fit `Omega_m` on retained releases. Score the omitted release after profiling one test-only intercept, identically for both models.
3. Survey hold-outs: for every Pantheon+ or DES `IDSURVEY` group with at least 20 selected objects, fit on the complement and score the omitted group after profiling one test-only intercept. Covariance subsets are marginal covariance subsets, never raw precision submatrices.
4. High-redshift hold-outs: primary threshold `z >= 0.50`; controls `z >= 0.80` and `z >= 1.00`. Fit on lower-redshift rows and score held-out rows with one test-only intercept per represented release.

No hold-out is deleted because of its outcome. Union3.1 has compressed nodes and no survey labels, so it participates in dataset and high-redshift hold-outs only.

## External-Budget Diagnostic and Robustness Gate

For each survey, the program reports the covariance-weighted residual common mode after removing its release intercept. It also reports the change in that mode between FR and Lambda-CDM, `Delta a_g`. The calibration-budget gate applies to `Delta a_g`, not to the arbitrary release intercept:

```math
|\Delta a_g| \leq \max(\sigma_{\mathrm{doc},g},\sigma_{\mathrm{proj},g}).
```

A model preference may be called robust only if all of the following hold:

- the full-data `Delta_BIC` has the same sign as the primary `z >= 0.50` held-out score;
- no dataset hold-out gives a strong opposite test preference (`|Delta_chi2_test| >= 4`);
- at least 80% of eligible survey hold-outs preserve the full-data sign and none gives a strong opposite test preference;
- every reported model-difference survey mode remains inside its release-grounded external budget;
- conclusions do not depend on adding a calibration prior already present in `STAT+SYS`.

Failure of any gate is reported explicitly and blocks promotion of a full-sample model preference.

## Primary Sources

- Pantheon+SH0ES DataRelease, distance and covariance products: https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR
- Pantheon+ systematic groupings, including calibration-only covariance: https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/sytematic_groupings
- Pantheon+ calibration analysis: https://arxiv.org/abs/2110.03486
- DES-SN5YR release: https://github.com/des-science/DES-SN5YR
- Union3 analysis: https://doi.org/10.3847/1538-4357/adc0a5
