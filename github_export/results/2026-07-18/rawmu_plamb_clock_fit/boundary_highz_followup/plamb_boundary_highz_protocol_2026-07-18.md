# Post hoc PLAMB boundary and high-redshift audit protocol

**Protocol date:** 18 July 2026  
**Status at registration:** no result with \(\gamma_c>2\) has been run or inspected  
**Claim status:** diagnostic only; widening an outcome-motivated bound cannot promote the model

## Known motivating results

The separately preregistered nested clock-law fit selected `PATH_FREE` as the lowest-BIC PLAMB-family cell, but it remained worse than flat Lambda-CDM by \(\Delta\mathrm{BIC}=+13.345052\). It reached the registered upper bound \(\gamma_c=2\), lost all three complete-release shape hold-outs and gave \(\Delta\chi^2_{\rm PATH-LCDM}=+321.926337\) for the conditional \(z\geq0.5\) prediction. The fully general model had a \(p\)-\(\alpha\) local-Hessian correlation of 0.997363.

## Questions

1. Does the PATH likelihood find a stable interior minimum when the diagnostic \(\gamma_c\) range is extended?
2. Does boundary extension improve genuine \(z\geq0.5\) transfer, or only the fitted low-redshift region?
3. Which releases, redshift bands, surveys and whitened observations account for the high-redshift excess?

## Locked likelihood

- Same 3,422-supernova sample and row order as the preceding nested fit.
- Pantheon+ and DES use `zHD`; Union3.1 uses released `z`.
- Released total covariance matrices.
- One analytically profiled magnitude intercept per release.
- \(H_0=67.5\ {\rm km\,s^{-1}\,Mpc^{-1}}\), absorbed by those intercepts.
- Primary boundary model: `PATH_FREE`, with \(\alpha=0\).
- Controls: `GENERAL_FREE`, with \(\alpha\) profiled, and spatially flat Lambda-CDM.
- Existing bounds remain \(-0.5\leq p\leq2.5\) and \(-0.5\leq\alpha\leq1.5\).

The clock family remains

\[
\begin{aligned}
\frac{c(z)}{c_0}                    &= 1+\gamma_c z,\\
\left|\frac{dz}{dT}\right|         &= H_0(1+z)^p,\\
D_{\rm path}(z)                    &= \frac{c_0}{H_0}\int_0^z
\frac{1+\gamma_c u}{(1+u)^p}\,du,\\
D_L(z)                             &= D_{\rm path}(z)(1+z)^\alpha.
\end{aligned}
\]

## Registered boundary map

Evaluate fixed-\(\gamma_c\) profiles for both `PATH_FREE` and `GENERAL_FREE` on:

- the full sample;
- the training sample \(z<0.5\).

At every grid point, profile \(p\), release intercepts and, for `GENERAL_FREE`, \(\alpha\). The fixed grid is:

- \(0\leq\gamma_c\leq2\) in steps of 0.05;
- \(2.1\leq\gamma_c\leq4\) in steps of 0.1;
- \(4.25\leq\gamma_c\leq8\) in steps of 0.25.

Separately polish the minima under \(\gamma_{c,\max}=2,4,8\) using deterministic bounded multi-start L-BFGS-B. The best two finite solutions must agree within \(10^{-4}\) in \(\chi^2\).

## Exact high-redshift prediction

For each release, partition residuals and covariance into low- and high-redshift blocks. The primary conditional score is

\[
\begin{aligned}
\widetilde r_H
    &= r_H-C_{HL}C_{LL}^{-1}r_L,\\
S_H
    &= C_{HH}-C_{HL}C_{LL}^{-1}C_{LH},\\
\chi^2_{H|L}
    &= \widetilde r_H^{\mathsf T}S_H^{-1}\widetilde r_H.
\end{aligned}
\]

Models are fitted only to \(z<0.5\), then held fixed for this score. The high-redshift conditional covariance is Cholesky-whitened after sorting by increasing redshift and source index. Squared innovations are additive and are grouped into:

- \(0.50\leq z<0.65\);
- \(0.65\leq z<0.80\);
- \(0.80\leq z<1.00\);
- \(z\geq1.00\);
- release and survey strata;
- individual source-index leverage rows.

The sum of all grouped innovations must reproduce the joint conditional score to \(10^{-8}\) relative tolerance.

## Covariance sensitivities

For the \(\gamma_{c,\max}=8\) low-redshift solution, compare:

1. primary released-covariance conditional prediction;
2. released-covariance marginal high-redshift score, without the conditional mean adjustment;
3. diagonal-covariance low-redshift refit and diagonal high-redshift score.

These are sensitivity directions, not additional searches.

## Diagnostic gates

The boundary is considered resolved only if:

1. the full-sample cap-8 PATH optimum is at least 5% of the \([0,8]\) interval from either boundary;
2. the fixed-grid profile at \(\gamma_c=8\) is at least \(\Delta\chi^2=3.84\) above the cap-8 minimum;
3. the two best deterministic cap-8 solutions agree within \(10^{-4}\) in \(\chi^2\).

High-redshift transfer is considered acceptable only if:

4. cap-8 PATH has \(\Delta\chi^2_{H|L}\leq0\) relative to Lambda-CDM;
5. no release/redshift-band contribution loses by more than 10;
6. no single ordered innovation contributes more than 25% of a positive net excess;
7. no survey contributes more than 50% of a positive net excess;
8. conditional, marginal and diagonal sensitivities retain the same sign.

The overall diagnostic passes only if all eight gates pass. Regardless of the result, this post hoc audit cannot itself promote a PLAMB law. If the cap-8 solution remains unresolved or high-redshift transfer fails, the phenomenological PATH law is retired as a supernova explanation unless a new physical derivation fixes the parameterisation before another data test.
