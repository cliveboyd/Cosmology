# Peter sparse-coupling closure search: frozen protocol

**Registration date:** 19 July 2026  
**Status:** prospectively frozen before sparse enumeration; local and not third-party timestamped

## Purpose and prior knowledge

The existing released-covariance audit already found `Delta BIC=+94.344983`
for the exact Peter path and `Delta BIC=+83.760921` for its best free
composite photon/source power, relative to flat Lambda-CDM (`N=3422`). Those
results are disclosed prior knowledge and are not counted as new evidence.

This screen asks a narrower question: how many additional power-law degrees of
freedom are required to make the redshift, atomic-frequency, photon-transfer,
supernova-luminosity and time-dilation equations mutually consistent?

For every quantity `X`, define an exponent `x_X` by

$$
\begin{aligned}
Z                              &= 1+z, \\
X(z)/X(0)                      &= Z^{x_X}.
\end{aligned}
$$

## Level 1: stated phenomenology

The frozen baseline is `x_a=0`, `x_c=1`, Peter's asserted `x_Q=-1`, no extra
spectral coupling, constant standardised luminosity and no photon-distance
transfer. The allowed perturbations obey

$$
\begin{aligned}
x_Q                            &= -1+\delta_Q, \\
x_z                            &= x_c+\sigma-x_a-x_Q=1, \\
s                              &= \tau-\frac{\ell}{2}=\hat{s}, \\
\hat{s}                        &= 0.0376049011\pm0.0086909617.
\end{aligned}
$$

All active subsets of `(delta_Q, sigma, tau, ell)` up to size three are
enumerated. Exact central-value closure is required; candidates are ranked by
active count, full-column-rank identification, largest endpoint change at
`z_max=2.26226`, exponent norm and name. Rank-deficient exact families are
retained as diagnostics but are not presented as identified point candidates.

## Level 2: microscopic exponent closure

The registered illustrative closure uses

$$
\begin{aligned}
Q                              &\propto \frac{\alpha^2m_ec^2}{h}, \\
M_{\rm Ch}                     &\propto
                                 \left(\frac{\hbar c}{G}\right)^{3/2}m_p^{-2}, \\
L_{\rm Ia}                     &\propto M_{\rm Ch}Z^{\eta_{\rm Ni}}, \\
s                              &= \tau-\frac{x_{L_{\rm Ia}}}{2}.
\end{aligned}
$$

The baseline exponents are `x_a=0`, `x_c=1`, `x_me=x_mp=-1`, with fixed
`h`, `G` and `alpha`. Perturbations to `h`, `G`, `m_p`, `m_e`, `alpha`, the
spectral law, photon transfer and nickel/luminosity conversion are enumerated
up to size four. The two frozen linear constraints are

$$
\begin{aligned}
\delta_h-\delta_{m_e}-2\delta_\alpha+\sigma
                               &= 1, \\
\tau-\frac{3}{4}\delta_h+\frac{3}{4}\delta_G
 +\delta_{m_p}-\frac{1}{2}\eta_{\rm Ni}
                               &= 1.7876049011.
\end{aligned}
$$

Changes in dimensional constants are parametrisation diagnostics, not
separately observable claims. Interpretation is restricted to `alpha`,
`m_e/m_p`, `Gm_p^2/(hbar c)` and observable redshift, flux and duration laws.

## Level 3: light-curve diffusion control

This deliberately simplified control adds

$$
\begin{aligned}
\kappa                         &\propto
                                 \frac{\alpha^2h^2}{m_e^2c^2m_p}, \\
t_{\rm diff}                   &\propto
                                 \left(\frac{\kappa M_{\rm Ch}}{vc}\right)^{1/2}, \\
v                              &\propto c.
\end{aligned}
$$

Together with a unit arrival-stretch exponent, the baseline predicts a
duration exponent `b=2.25`. The DES measurement `b=1.003 +/- 0.005 (stat) +/-
0.010 (sys)` therefore adds

$$
\begin{aligned}
\frac{7}{4}\delta_h-\frac{3}{4}\delta_G
-\frac{3}{2}\delta_{m_p}-\delta_{m_e}+\delta_\alpha
+\delta_{\rm duration}         &= -1.247.
\end{aligned}
$$

This is a sensitivity calculation, not a complete Type Ia explosion model.

## Gates and claim boundary

No candidate is promoted unless it also supplies an action-level scalar,
matter and electromagnetic model, passes the same-nuisance released-covariance
supernova gate `abs(Delta BIC)<=2`, passes time dilation, and survives external
dimensionless-constant and photon-transfer constraints. Atomic-clock limits are
not converted into an invented exponent prior until a scalar history maps
redshift to a local time derivative.

The search identifies algebraic rank and observational consequences. It does
not establish that any enumerated coupling is physically viable, and failure
constrains this power-law closure rather than every possible FR or
varying-constant theory.
