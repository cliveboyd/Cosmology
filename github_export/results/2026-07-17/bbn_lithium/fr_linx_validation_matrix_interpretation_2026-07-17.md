# FR/LINX BBN Validation Matrix: Interpretation

Date: 17 July 2026

## Status

The preregistered targeted validation matrix completed all **92/92** LINX evaluations with zero missing or failed rows. It used

\[
\begin{aligned}
\mathrm{rtol}             &= 10^{-7},\\
\mathrm{atol}             &= 10^{-11},\\
N_{\rm interp,primary}    &= 300,\\
N_{\rm interp,check}      &= 600,
\end{aligned}
\]

across `key_PRIMAT_2023`, `small_PRIMAT_2023` and `full_PRIMAT_2023`. The matrix included ten interpretation-bearing states, centred unit-pull responses for the three principal D-sensitive rates, full-network responses for five Li/Be channels absent from the key network, and sampling-600 checks at the baseline and best registered direct-lithium state.

## Main Findings

### 1. Numerical Convergence Passes

The maximum change between the completed scan's loose key-network settings and the strict key reruns is

\[
\begin{aligned}
\max|\Delta[10^5({\rm D/H})]|           &= 5.51\times10^{-4},\\
\max|\Delta Y_p|                         &= 6.53\times10^{-5},\\
\max|\Delta[10^{10}({}^{7}{\rm Li/H})]| &= 1.71\times10^{-3}.
\end{aligned}
\]

The maximum sampling-600 minus sampling-300 differences over key, small and full networks are respectively `0.000116`, `2.02e-5` and `0.000319` in the same D, helium and lithium units. Numerical tolerance and interpolation error are therefore much smaller than the adopted observational uncertainties and do not drive the lithium conclusion.

### 2. Network Completeness Is a Precision Correction

At the strict baseline:

| Network | D/H x 1e5 | Yp | Li-7/H x 1e10 |
|---|---:|---:|---:|
| key | 2.434293 | 0.2469465 | 5.695731 |
| small | 2.435004 | 0.2469494 | 5.638004 |
| full | 2.434914 | 0.2469325 | 5.637694 |

Across all ten core states, full minus key changes D/H by at most `0.000649`, helium by `1.48e-5`, and lithium by `0.0580`. The full-network correction is important for publication-grade numbers but is far too small to bridge the direct-lithium discrepancy.

For the old best registered direct-lithium controls, the strict full result is

\[
\begin{aligned}
10^5({\rm D/H})                    &= 2.585689,\\
Y_p                                &= 0.248639,\\
10^{10}({}^{7}{\rm Li/H})         &= 4.647982,\\
\chi^2_{\rm registered,direct}    &= 210.57.
\end{aligned}
\]

The corresponding strict key lithium prediction is `4.693676`. The full network improves the penalised score but does not approach a joint abundance pass.

### 3. Deuterium Theory Uncertainty Must Enter the Likelihood

Centred unit-pull responses for `dpHe3g`, `ddHe3n` and `ddtp`, combined in quadrature under independent standard-normal priors, give

| Network | State | Local sigma(D/H x 1e5) | Local sigma(Yp) | Local sigma(Li x 1e10) |
|---|---|---:|---:|---:|
| key | baseline | 0.02609 | 2.16e-5 | 0.08493 |
| small | baseline | 0.02609 | 2.15e-5 | 0.08423 |
| full | baseline | 0.02608 | 2.17e-5 | 0.08423 |
| full | best direct | 0.02704 | 2.14e-5 | 0.06955 |

The baseline D theory scale, `0.02608`, is approximately 90% of the registered observational error `0.029` and exceeds the newer D observational error `0.024`. A provisional independent-error combination would give

\[
\begin{aligned}
\sigma_{\rm D,registered,combined}
  &= \sqrt{0.029^2+0.02608^2}
   = 0.03900,\\
\sigma_{\rm D,current,combined}
  &= \sqrt{0.024^2+0.02608^2}
   = 0.03544.
\end{aligned}
\]

Simple quadrature is only an initial diagnostic. The preferred analysis should marginalise the rate pulls jointly so that their non-linear and cross-abundance responses are retained. Nevertheless, the matrix establishes that scoring D/H with observational uncertainty alone materially understates the likelihood width.

### 4. Omitted Li/Be Channels Do Not Rescue Lithium

All five previously unavailable key-network channels were evaluated in the full network with centred unit pulls. The largest lithium response is from `Be7daap`:

\[
\begin{aligned}
\left|\frac{\partial[10^{10}({}^{7}{\rm Li/H})]}{\partial q_{\rm Be7daap}}\right|_{\rm baseline}
  &= 0.06079,\\
\left|\frac{\partial[10^{10}({}^{7}{\rm Li/H})]}{\partial q_{\rm Be7daap}}\right|_{\rm best}
  &= 0.04640.
\end{aligned}
\]

The next-largest omitted-channel response is only `0.00353`. These are useful completeness corrections but are much smaller than the residual lithium excess of roughly 3-4 in the same units.

## Gate Readout

- Registered D+He+Li two-standard-deviation passes: **0**.
- Combined-current D+He+Li two-standard-deviation passes: **0**.
- The result remains negative under strict tolerances, larger networks, updated neutron-lifetime states, D-rate perturbations and omitted Li/Be channels.

## Consequences

1. Retain the completed key scan for broad exploration, but use small or full networks for publication-grade interpretation-bearing rows.
2. Replace the observational-only D term with a joint nuclear-rate likelihood or validated abundance covariance. This is the main statistical correction exposed by the matrix.
3. Do not launch another broad blind nuclear scan. The full-network and omitted-channel responses are too small to provide a conventional nuclear solution.
4. The next conventional analysis should be a hierarchical posterior over baryon density, neutron lifetime, all key-network rate pulls and stellar-lithium depletion, with matched CMB and abundance likelihoods.
5. The FR claim boundary is unchanged: this remains a standard expanding-background LINX calculation. An FR-native test still requires an explicit non-expanding thermal history and coherent weak, mass, binding-energy and Q-value evolution.
