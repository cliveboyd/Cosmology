# Preregistered PLAMB nested clock-law fit

**Protocol date:** 18 July 2026  
**Status at registration:** no flexible clock-law fit in this protocol has been run or inspected  
**Known prior result:** the fixed Peter-law cell has \(\Delta\mathrm{BIC}_{\rm PLAMB-\Lambda CDM}=+94.344983\) on the locked 3,422-supernova sample.

## Question

Which term is responsible for the failure of the fixed Peter-law raw-distance curve: the light-speed coefficient \(\gamma_c\), the redshift-rate power \(p\), the observed-flux exponent \(\alpha\), or a combination of them?

## Locked data and likelihood

- Pantheon+SH0ES: non-calibrators with \(z_{\rm HD}>0.01\), `MU_SH0ES`, released total covariance.
- DES-Dovekie raw MU: finite rows with \(z_{\rm HD}>0.01\), `MU`, released total covariance.
- Union3.1 UNITY1.8: released compressed \(z\), `MU` and precision matrix.
- Expected total: 3,422 supernovae.
- Redshift frames: `zHD` for Pantheon+ and DES; released `z` for Union3.1.
- One unpenalised, analytically profiled magnitude intercept per release in every model.
- \(H_0=67.5\ {\rm km\,s^{-1}\,Mpc^{-1}}\); its scale is absorbed by the release intercepts.
- No row, covariance or redshift-frame change is permitted after inspecting the flexible-fit results.

The clock family is

\[
\begin{aligned}
\frac{c(z)}{c_0}                    &= 1+\gamma_c z,\\
\left|\frac{dz}{dT}\right|         &= H_0(1+z)^p,\\
D_{\rm path}(z)                    &= \frac{c_0}{H_0}\int_0^z
\frac{1+\gamma_c u}{(1+u)^p}\,du,\\
D_L(z)                             &= D_{\rm path}(z)(1+z)^\alpha,\\
\mu_{\rm model}(z)                &= 5\log_{10}\!\left(\frac{D_L}{\rm Mpc}\right)+25.
\end{aligned}
\]

The integral is evaluated analytically, with logarithmic limits at \(p=1\) and \(p=2\). The physical decomposition of the fitted flux exponent is

\[
\begin{aligned}
\alpha &= \frac{e+b-s}{2},\\
e      &= \text{photon-energy exponent},\\
b      &= \text{arrival-rate or time-dilation exponent},\\
s      &= \text{standardised intrinsic-luminosity evolution exponent}.
\end{aligned}
\]

Supernova distances identify \(\alpha\), not \(e\), \(b\) and \(s\) separately.

## Registered model ladder

| Model key | Free shape parameters | Fixed values |
|---|---:|---|
| `PETER_FIXED` | 0 | \(\gamma_c=1,p=0,\alpha=0\) |
| `FRACTIONAL_FIXED` | 0 | \(\gamma_c=1,p=1,\alpha=0\) |
| `P_FREE` | 1 | \(\gamma_c=1,\alpha=0\) |
| `ALPHA_FREE` | 1 | \(\gamma_c=1,p=0\) |
| `GAMMA_FREE` | 1 | \(p=0,\alpha=0\) |
| `PATH_FREE` | 2 | \(\alpha=0\) |
| `P_ALPHA_FREE` | 2 | \(\gamma_c=1\) |
| `GENERAL_FREE` | 3 | none |
| `LCDM` | 1 | spatially flat, \(\Omega_m\) free |

Registered bounds are

\[
\begin{aligned}
0.0  &\leq \gamma_c \leq 2.0,\\
-0.5 &\leq p        \leq 2.5,\\
-0.5 &\leq \alpha   \leq 1.5,\\
0.05 &\leq \Omega_m \leq 0.60.
\end{aligned}
\]

## Optimisation and checks

- Profile \(\alpha\) analytically whenever it is free, subject to its registered bounds.
- Use bounded scalar minimisation for one nonlinear parameter.
- Use two deterministic differential-evolution searches, seeds 20260718 and 20260719, followed by L-BFGS-B polishing for two nonlinear parameters.
- Require the two searches to agree to \(10^{-4}\) in \(\chi^2\); otherwise flag optimisation instability.
- Count all release intercepts and free shape parameters in BIC.
- Report boundary distances, numerical Hessian rank and condition number for the best PLAMB-family model.

## Outcome-independent hold-out rule

Select the PLAMB-family cell with the lowest full-sample BIC. Compare that cell, `PETER_FIXED`, `FRACTIONAL_FIXED` and `LCDM` in:

1. each complete-release shape hold-out, profiling one test-release intercept;
2. every non-compressed survey stratum with at least 20 objects and at least 20 training objects left in its release;
3. high-redshift hold-outs at \(z\geq0.5,0.8,1.0\).

For survey and redshift partitions, train on the complement and score the held-out residual with the exact Gaussian conditional covariance. Free shape parameters are re-fitted on each training set. A full-release hold-out is explicitly a shape-transfer test because its otherwise unknown magnitude intercept is profiled on the test release.

## Promotion gates

A flexible PLAMB clock law is not promoted as a physical result unless all of the following pass:

1. \(\Delta\mathrm{BIC}_{\rm best\ PLAMB-\Lambda CDM}<-10\) on the primary sample.
2. No fitted shape parameter lies within 1% of its registered bound.
3. The profiled Hessian has full rank and condition number below \(10^8\).
4. The candidate beats ΛCDM in all three complete-release shape hold-outs.
5. It beats ΛCDM in at least 75% of eligible survey and redshift hold-outs, with no single held-out loss exceeding \(\Delta\chi^2=10\).
6. Its release intercepts differ from the corresponding ΛCDM intercepts by no more than 0.050 mag.
7. A separately specified external constraint breaks the \(e,b,s\) degeneracy; the supernova fit alone cannot pass this physical-identification gate.

Failure of a promotion gate is reported directly. A statistically preferred flexible curve remains phenomenological unless the physical-identification gate also passes.
