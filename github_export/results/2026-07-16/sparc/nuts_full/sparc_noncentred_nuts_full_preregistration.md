# SPARC Full Filtered All-Q2 Non-centred NUTS Preregistration

- **Written before full sampling:** `2026-07-15T20:50:34.117599+00:00`
- **Galaxies:** `157`
- **Rotation points:** `3001`
- **Models:** `RAR, PLAMB_OPTICAL_DEPTH_KAPPA_P`

## Fixed Scope

Every galaxy remaining after the previously documented six-galaxy persistent-negative-family exclusion is included. There is no further selection, outcome ranking, stress-subset filtering or residual-based removal.

The July 14 likelihood, fixed H0, distance treatment, velocity-error floor, stellar mass-to-light priors and physical parameter bounds are unchanged. The verified non-centred coordinates are

\[
\begin{aligned}
z_{D,g}        &\sim \mathcal{N}(0,1), & \log D_g        &= \mu_{D,g} + \sigma_{D,g}z_{D,g}, \\
z_{\eta,g}     &\sim \mathcal{N}(0,1), & \log\eta_g     &= \sigma_{\log\eta}z_{\eta,g}, \\
z_{\Upsilon d} &\sim \mathcal{N}(0,1), & \log\Upsilon_d &= \log(0.5)+\sigma_{\Upsilon}z_{\Upsilon d}, \\
z_{\Upsilon b} &\sim \mathcal{N}(0,1), & \log\Upsilon_b &= \log(0.7)+\sigma_{\Upsilon}z_{\Upsilon b}.
\end{aligned}
\]

The same four-standard-deviation nuisance truncation is retained.

## Sampling Budget

- four chains, `2000` warm-up and `2000` retained draws per chain;
- target acceptance `0.9`;
- maximum tree depth `12`; and
- fixed seed `16072611`.

## Locked Gates

- likelihood equivalence within `1e-08` km/s;
- zero post-warm-up divergences;
- maximum rank-normalised split R-hat `<= 1.01`;
- absolute R-hat failure ceiling `1.05`;
- minimum bulk and tail ESS `>= 400.0`; and
- minimum chain E-BFMI `>= 0.3`.

A strict pass clears the convergence gate only. It is not a full-sample PLAMB win, a Bayes factor, or held-out predictive validation. Model promotion remains blocked until galaxy-held-out comparison is completed.
