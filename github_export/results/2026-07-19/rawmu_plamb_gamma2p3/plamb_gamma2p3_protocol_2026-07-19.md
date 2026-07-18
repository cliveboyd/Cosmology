# Outcome-aware fixed-gamma PLAMB protocol

**Protocol date:** 19 July 2026  
**Claim status:** post hoc conditional target; cannot promote a PLAMB claim  
**Selection disclosure:** \(\gamma_c=2.3\) was selected after inspecting the 18 July full-sample and high-redshift profile

## Known motivation

The diagnostic boundary scan found an interior full-sample PATH minimum at \(\gamma_c=2.276517\), \(p=0.784634\). On the registered fixed grid, \(\gamma_c=2.3\) gave the best outcome-visible \(z\geq0.5\) conditional score, \(\Delta\chi^2_{\rm PATH-\Lambda CDM}=-5.490541\), while costing \(\Delta\chi^2=4.597673\) relative to the independently selected low-redshift PATH optimum. This analysis asks what the fixed value would imply if a future physical derivation specified it independently. It is not a new held-out discovery test.

## Locked data and models

- Same 3,422-supernova sample, row order, redshift frames and released total covariance as the 18 July nested analysis.
- Training set: \(z<0.5\); predictive set: \(z\geq0.5\).
- Pantheon+ and DES use `zHD`; Union3.1 uses released `z`.
- One magnitude intercept per release with a flat prior, updated only by the low-redshift data.
- \(H_0=67.5\ {\rm km\,s^{-1}\,Mpc^{-1}}\), absorbed by the intercepts.
- PATH target: \(\gamma_c=2.3\), \(\alpha=0\), \(p\) free.
- Control: spatially flat Lambda-CDM with \(\Omega_m\) free.

The PATH model is

\[
\begin{aligned}
\frac{c(z)}{c_0}                    &= 1+2.3z,\\
\left|\frac{dz}{dT}\right|         &= H_0(1+z)^p,\\
D_L(z)                             &= \frac{c_0}{H_0}\int_0^z
\frac{1+2.3u}{(1+u)^p}\,du.
\end{aligned}
\]

## Parameter priors and integration

Use uniform priors

\[
\begin{aligned}
-0.5 &\leq p        \leq 2.5,\\
0.05 &\leq \Omega_m \leq 0.60.
\end{aligned}
\]

For each model:

1. find the bounded low-redshift maximum-likelihood point;
2. estimate local curvature and construct a deterministic 1,601-point grid over at least \(\pm8\) local standard deviations, clipped to the registered bounds;
3. use the full bound with 3,001 points if either grid edge is less than \(\Delta\chi^2=30\) above the minimum;
4. normalise \(\exp[-\Delta\chi^2/2]\) with the uniform prior;
5. report posterior mean, median, standard deviation and equal-tail 68% and 95% intervals.

The same procedure is repeated separately for each release's \(z<0.5\) subset to audit parameter stability.

## Integrated conditional prediction

For a release, let \(A=C_{HL}C_{LL}^{-1}\), and let \(\widehat\delta\) and \(V_\delta\) be the low-redshift posterior mean and variance of its magnitude intercept. Then

\[
\begin{aligned}
v                              &= \mathbf 1_H-A\mathbf 1_L,\\
\widetilde r_H(\theta)         &= r_H(\theta)-Ar_L(\theta)-\widehat\delta(\theta)v,\\
S_H                            &= C_{HH}-C_{HL}C_{LL}^{-1}C_{LH},\\
V_H                            &= S_H+V_\delta vv^{\mathsf T},\\
\log p(y_H\mid y_L,\theta)    &= -\frac12\left[
\widetilde r_H^{\mathsf T}V_H^{-1}\widetilde r_H+
\log|V_H|+N_H\log(2\pi)\right].
\end{aligned}
\]

The final predictive density integrates this expression over the low-redshift posterior grid. The reported tail area is the posterior average of the corresponding \(\chi^2_{N_H}\) survival probability. Release-level predictive densities use the same joint low-redshift parameter posterior but are descriptive and do not add to the joint mixture log density.

## Physical readout

At \(z=0.5,1.0,1.5\) and the Union3.1 endpoint \(z=2.262260\), report:

- \(c(z)/c_0=1+2.3z\);
- \(|dz/dT|/H_0=(1+z)^p\) using the joint posterior mean of \(p\);
- the dimensionless PATH integral;
- the distance-modulus shape difference from best-fit Lambda-CDM after anchoring both curves at \(z=0.1\).

These are consequences of the fitted phenomenology, not established measurements of light-speed evolution.

## Conditional gates

1. Full-sample \(|\Delta\mathrm{BIC}_{\rm PATH-\Lambda CDM}|\leq2\).
2. Low-redshift \(|\Delta\mathrm{BIC}_{\rm PATH-\Lambda CDM}|\leq2\).
3. Integrated high-redshift \(\Delta\log p_{\rm PATH-\Lambda CDM}\geq0\).
4. Joint PATH posterior-predictive tail area is at least 0.05.
5. Every release-level PATH tail area is at least 0.05.
6. Maximum pairwise release-specific \(p\)-posterior tension is at most 2 standard deviations.
7. The joint \(p\) posterior is at least 1% of its prior range from either bound and both integration edges satisfy \(\Delta\chi^2\geq30\).
8. An independent pre-data physical derivation specifies \(\gamma_c=2.3\) to the required precision.

The first seven gates describe conditional statistical performance. Gate 8 is deliberately false for the present analysis. Therefore this run cannot promote the model even if every statistical gate passes.
