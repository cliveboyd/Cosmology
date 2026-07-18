# Hierarchical LINX BBN Posterior: Detailed Interpretation

Date: 18 July 2026

## Scope

This is the registered gradient-free replacement for the direct LINX NUTS
programme that failed its preposterior gradient smoke tests. It uses 4,800 exact
`key_PRIMAT_2023` LINX evaluations spanning the CMB baryon-density summary,
neutron lifetime and all 12 key-network rate pulls. The fixed split contains
4,000 training and 800 holdout evaluations.

Two 64-walker ensemble posteriors were sampled:

1. `predictive`: D/H, helium and the CMB baryon-density summary determine a
   posterior-predictive primordial lithium abundance without using lithium data.
2. `depletion`: the lithium likelihood is included with a positive stellar
   depletion parameter.

The abundance definitions are

\[
\begin{aligned}
y_{\rm D}          &= 10^5({\rm D/H}), \\
Y_{\rm p}          &= 4Y_{{}^4{\rm He}}, \\
y_{{\rm Li},0}     &= 10^{10}[({}^7{\rm Li}+{}^7{\rm Be})/{\rm H}], \\
y_{{\rm Li},\star} &= y_{{\rm Li},0}10^{-\delta_{\rm Li}}, \\
f_{\rm dep}        &= 10^{\delta_{\rm Li}}.
\end{aligned}
\]

## Numerical Validation

All registered numerical gates passed.

| diagnostic | predictive | depletion |
|:--|--:|--:|
| Walkers x production steps | 64 x 12,000 | 64 x 12,000 |
| Mean acceptance fraction | 0.3561 | 0.3435 |
| Maximum split Rhat | 1.0200 | 1.0203 |
| Minimum autocorrelation ESS | 4,582.5 | 3,641.7 |
| Exact posterior validations | 160/160 | 160/160 |
| Exact-correction importance ESS | 159.98 | 159.98 |

The selected quadratic surrogate passed the independent 800-row holdout gate.
Its holdout root-mean-square errors were `0.000148` for D/H, `0.0000124` for
helium and `0.000570` for lithium in the scaled abundance units. Across the
selected exact posterior validations, the 95th-percentile absolute surrogate
errors were below `0.020` observational standard deviations for every
abundance. Exact LINX therefore confirms that surrogate error is negligible for
the reported inference.

## Posterior Abundances

The intervals below are 95% intervals from the retained full surrogate chains.
The exact-validation importance correction changes the lithium medians by less
than `0.01` in the scaled units.

| quantity | predictive median | predictive 95% interval | depletion median | depletion 95% interval |
|:--|--:|--:|--:|--:|
| D/H x 1e5 | 2.5041 | [2.4640, 2.5440] | 2.5047 | [2.4645, 2.5447] |
| Yp | 0.246682 | [0.246462, 0.246903] | 0.246677 | [0.246454, 0.246898] |
| Primordial Li/H x 1e10 | 5.4897 | [5.1627, 5.8462] | 5.4797 | [5.1451, 5.8210] |
| Surface Li/H x 1e10 | not conditioned | not conditioned | 1.5188 | [1.0595, 1.9870] |
| Lithium depletion, dex | 0 | fixed | 0.5577 | [0.4366, 0.7145] |
| Lithium depletion factor | 1 | fixed | 3.611 | [2.733, 5.182] |

The exact-corrected medians are `5.4949` for predictive primordial lithium,
`5.4859` for depletion-arm primordial lithium and `0.5519 dex` for depletion.
The latter corresponds to an exact-corrected median depletion factor of about
`3.56`.

## Lithium Gate

Using the predictive posterior standard deviation `0.1745` and the lithium
plateau summary `1.45 +/- 0.25`, the approximate posterior-predictive discrepancy
is

\[
\begin{aligned}
z_{\rm Li}
  &= \frac{5.4897-1.45}{\sqrt{0.1745^2+0.25^2}} \\
  &\simeq 13.3.
\end{aligned}
\]

Marginalising all 12 key-network rate pulls therefore does not produce a nuclear
solution to the lithium problem. Agreement with the plateau occurs only after
introducing substantial positive stellar depletion.

## Parameter Readout

The predictive arm gives `eta_fac = 0.99109` with a 95% interval
`[0.98064, 1.00140]`, corresponding to approximately
`Omega_b h^2 = 0.02222 [0.02199, 0.02245]`. Neutron lifetime remains
prior-dominated at `878.40 s [877.42, 879.37]`.

The largest median nuclear pulls are `dpHe3g = -0.90`, `ddHe3n = -0.74` and
`ddtp = -0.63`. Their 95% intervals all include zero. The posterior uses a
combination of a lower baryon-density ratio and moderate deuterium-rate shifts;
it does not identify a single anomalous nuclear channel, and the Li/Be-specific
rate pulls remain close to their priors.

## Claim Boundary

`Promotable = true` in the gate JSON means that the numerical posterior is fit
for scientific interpretation. It does **not** mean that homogeneous standard
BBN solves the lithium discrepancy.

The inference remains conditional on homogeneous standard BBN, a Gaussian CMB
baryon-density summary and the registered abundance summaries. It contains no
matter-antimatter domains, annihilation transport, entropy injection or
antimatter-galaxy population parameter. An FR matter-antimatter cosmology would
require a separate early-universe likelihood before these conclusions could be
transferred to that sector.
