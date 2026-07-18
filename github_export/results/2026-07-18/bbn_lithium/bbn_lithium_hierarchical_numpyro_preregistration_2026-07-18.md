# Hierarchical LINX BBN Lithium Posterior: Preregistration

Registered: {payload['registered_at']}

## Registered Question

Can homogeneous standard BBN, after marginalising all 12 `key_PRIMAT_2023`
nuclear-rate pulls and the laboratory neutron-lifetime uncertainty, predict the
observed lithium abundance? If an explicit positive stellar-depletion nuisance
is admitted, what depletion is required?

## Likelihood And Priors

\[
egin{aligned}
y_{m D}        &= 10^5({m D/H}), \
Y_{m p}        &= 4Y_{{}^4{m He}}, \
y_{{m Li},0}   &= 10^{10}[({}^7{m Li}+{}^7{m Be})/{m H}], \
y_{{m Li},\star} &= y_{{m Li},0}\,10^{-\delta_{m Li}}.
\end{aligned}
\]

The abundance likelihoods are Normal summaries with D/H = 2.533 +/- 0.024,
Yp = 0.2458 +/- 0.0013 and Li/H = 1.45 +/- 0.25 in the scaled units above.
The CMB information is a five-sigma truncated Gaussian
`eta_fac = 1.0 +/- 0.00653595`
summary, not a matched map-level CMB likelihood. Neutron lifetime is
`878.4 +/- 0.5 s`, truncated at five sigma; each LINX nuclear rate pull is a
standard Normal truncated at +/-4. The depletion
arm uses `delta_Li ~ HalfNormal(0.30 dex)`.

## Arms

1. `predictive`: condition on D, He and the CMB eta summary; do not condition on Li.
2. `depletion`: add the Li likelihood and positive stellar depletion.
3. Reweight the depletion posterior to 0.15 and 0.60 dex HalfNormal controls,
   reporting importance effective sample size before interpreting either control.

## Numerical Gates

The headline posterior requires maximum split Rhat <= 1.01, minimum bulk ESS >=
400 and zero divergent transitions. Failed gates make the output exploratory.

Maximum tree depth is 4, matching LINX's supplied NumPyro example. Four chains
run sequentially to control peak memory. An initial depth-8 smoke attempt was
stopped during XLA compilation at 0/16 transitions and supplied no posterior
draws. A second smoke attempt failed at warm-up transition 2 when automatic
initial step-size behaviour drove the implicit solver to a non-finite result.
Before retaining any posterior draw, the upstream LINX initial step size 0.01
and negligible-tail truncations described above were adopted. No likelihood
anchor was changed after observing results. The bounded reverse-mode smoke then
failed on its first trajectory, also before retaining a draw. The final gradient
trial uses forward-mode differentiation, Kvaerno3, `max_steps = 8092 x 12` and
`rtol = 1e-4`, following the upstream numerical envelope. A persistent run-local
JAX cache is enabled before the final smoke and reused by the overnight process.

## Matter-Antimatter Boundary

This run assumes homogeneous standard BBN with one net baryon-density field.
It contains no matter-antimatter domains, boundary annihilation, entropy
injection, spectral-distortion or antimatter-galaxy population parameter.
Consequently it can test standard-BBN and nuisance explanations of lithium, but
cannot validate or exclude an FR matter-antimatter formation mechanism.
