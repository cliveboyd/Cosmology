# Peter clock-law sparse-coupling closure search

Generated: 2026-07-19 18:35:18

## Decision

**THE STATED EQUATIONS REQUIRE AT LEAST TWO INDEPENDENT REPAIRS; THE ILLUSTRATIVE ATOMIC/CHANDRASEKHAR/DIFFUSION CLOSURE REQUIRES AT LEAST THREE. NO CANDIDATE PASSES THE EXISTING SUPERNOVA OR ACTION-LEVEL PROMOTION GATES.**

The smallest stated-phenomenology closure uses `2` active exponents. The microscopic atomic plus luminosity closure also needs `2`, while adding the simplified diffusion-time equation raises the minimum to `3`. This is evidence for missing structure in the explicit closure, not evidence that any enumerated power law is physical.

Every candidate inherits the best possible Peter-path composite supernova result `Delta BIC=+83.760921` relative to flat Lambda-CDM. The frozen promotion threshold is `abs(Delta BIC)<=2`, so no scalar-field integration or new supernova sampling is warranted from this screen.

## Frozen status

Protocol SHA-256: `2a3287de46ea6cda9e4905bb2575ff32f70f20e74124d50767272db0c18b3530`.

The earlier supernova results were disclosed before enumeration. The new output is the sparse rank, endpoint consequences and explicit separation of identified solutions from rank-deficient families.

## Exponent convention

$$
\begin{aligned}
Z                              &= 1+z, \\
X(z)/X(0)                      &= Z^{x_X}, \\
Z_{\rm max}                  &= 3.262260.
\end{aligned}
$$

## Stated phenomenology

The baseline assertion `x_Q=-1` is allowed either to change directly or to be accompanied by an extra spectral coupling. Photon transfer and source luminosity are kept separate even though supernova magnitudes identify only their combination.

$$
\begin{aligned}
x_Q                            &= -1+\delta_Q, \\
-\delta_Q+\sigma              &= -1, \\
\tau-\frac{\ell}{2}           &= 0.0376049011.
\end{aligned}
$$

There are `4` full-rank, minimum-count stated solutions. The best-ranked endpoint repair is `delta_Q_atomic;tau_photon_distance` with values `{"delta_Q_atomic": 1.0, "tau_photon_distance": 0.037604901114341276}`.

| active_couplings | active_values_json | maximum_active_endpoint_change_factor | atomic_Q_final_exponent | supernova_composite_s | duration_b_predicted |
| --- | --- | --- | --- | --- | --- |
| delta_Q_atomic;tau_photon_distance | {"delta_Q_atomic": 1.0, "tau_photon_distance": 0.037604901114341276} | 3.262260 | 0 | 0.037605 | 1.000000 |
| sigma_spectral;tau_photon_distance | {"sigma_spectral": -1.0, "tau_photon_distance": 0.037604901114341276} | 3.262260 | -1.000000 | 0.037605 | 1.000000 |
| delta_Q_atomic;ell_source | {"delta_Q_atomic": 1.0, "ell_source": -0.07520980222868255} | 3.262260 | 0 | 0.037605 | 1.000000 |
| sigma_spectral;ell_source | {"ell_source": -0.07520980222868255, "sigma_spectral": -1.0} | 3.262260 | -1.000000 | 0.037605 | 1.000000 |

The economical branches either replace the asserted atomic clock law with `Q=constant`, or retain it by inserting a new spectral factor proportional to `(1+z)^-1`. A second independent photon/source law is then still required to reproduce the measured supernova composite.

## Microscopic closure

$$
\begin{aligned}
Q                              &\propto \frac{\alpha^2m_ec^2}{h}, \\
M_{\rm Ch}                     &\propto \left(\frac{\hbar c}{G}\right)^{3/2}m_p^{-2}, \\
L_{\rm Ia}                     &\propto M_{\rm Ch}Z^{\eta_{\rm Ni}}, \\
\delta_h-\delta_{m_e}-2\delta_\alpha+\sigma &= 1, \\
\tau-\frac{3}{4}\delta_h+\frac{3}{4}\delta_G+\delta_{m_p}-\frac{1}{2}\eta_{\rm Ni} &= 1.7876049011.
\end{aligned}
$$

There are `19` full-rank minimum-count microscopic solutions. The best-ranked active-set endpoint factor is `8.27882` for `delta_alpha;tau_photon_distance`. This ranking is descriptive: a small dimensional exponent norm is not a physical prior.

| active_couplings | active_values_json | maximum_active_endpoint_change_factor | maximum_dimensionless_endpoint_change_factor | atomic_Q_exponent | M_Ch_exponent | supernova_composite_s |
| --- | --- | --- | --- | --- | --- | --- |
| delta_alpha;tau_photon_distance | {"delta_alpha": -0.5, "tau_photon_distance": 1.7876049011143413} | 8.278818 | 1.806173 | 0 | 3.500000 | 0.037605 |
| delta_m_p;delta_alpha | {"delta_alpha": -0.5, "delta_m_p": 1.7876049011143413} | 8.278818 | 68.538821 | 0 | -0.075210 | 0.037605 |
| delta_m_e;tau_photon_distance | {"delta_m_e": -1.0, "tau_photon_distance": 1.7876049011143413} | 8.278818 | 3.262260 | 0 | 3.500000 | 0.037605 |
| delta_m_p;delta_m_e | {"delta_m_e": -1.0, "delta_m_p": 1.7876049011143413} | 8.278818 | 68.538821 | 0 | -0.075210 | 0.037605 |
| delta_m_p;sigma_spectral | {"delta_m_p": 1.7876049011143413, "sigma_spectral": 1.0} | 8.278818 | 68.538821 | 1.000000 | -0.075210 | 0.037605 |
| sigma_spectral;tau_photon_distance | {"sigma_spectral": 1.0, "tau_photon_distance": 1.7876049011143413} | 8.278818 | 1.000000 | 1.000000 | 3.500000 | 0.037605 |
| delta_G;delta_alpha | {"delta_G": 2.3834732014857885, "delta_alpha": -0.5} | 16.747799 | 16.747799 | 0 | -0.075210 | 0.037605 |
| delta_G;delta_m_e | {"delta_G": 2.3834732014857885, "delta_m_e": -1.0} | 16.747799 | 16.747799 | 0 | -0.075210 | 0.037605 |
| delta_G;sigma_spectral | {"delta_G": 2.3834732014857885, "sigma_spectral": 1.0} | 16.747799 | 16.747799 | 1.000000 | -0.075210 | 0.037605 |
| delta_h;delta_alpha | {"delta_alpha": -1.691736600742895, "delta_h": -2.3834732014857893} | 16.747799 | 16.747799 | -8.8818e-16 | -0.075210 | 0.037605 |
| delta_h;delta_m_p | {"delta_h": 0.9999999999999993, "delta_m_p": 2.5376049011143404} | 20.095890 | 123.792951 | 6.6613e-16 | -0.075210 | 0.037605 |
| delta_h;tau_photon_distance | {"delta_h": 0.9999999999999993, "tau_photon_distance": 2.5376049011143404} | 20.095890 | 3.262260 | 6.6613e-16 | 5.000000 | 0.037605 |

Order-unity changes in `alpha`, `m_e/m_p` or `G m_p^2/(hbar c)` are warnings, not accepted solutions. A local atomic-clock bound cannot be translated into a cosmological exponent without the scalar field's time history and environmental response.

## Diffusion-time control

$$
\begin{aligned}
\kappa                         &\propto \frac{\alpha^2h^2}{m_e^2c^2m_p}, \\
t_{\rm diff}                   &\propto \left(\frac{\kappa M_{\rm Ch}}{vc}\right)^{1/2}, \\
v                              &\propto c, \\
\frac{7}{4}\delta_h-\frac{3}{4}\delta_G-\frac{3}{2}\delta_{m_p}-\delta_{m_e}+\delta_\alpha+\delta_{\rm duration} &= -1.247.
\end{aligned}
$$

The simplified light-curve system has `65` full-rank solutions at the minimum active count of `3`. Its best-ranked member is `delta_m_p;sigma_spectral;tau_photon_distance` with values `{"delta_m_p": 0.8313333333333331, "sigma_spectral": 0.9999999999999994, "tau_photon_distance": 0.956271567781008}`.

| active_couplings | active_values_json | maximum_active_endpoint_change_factor | maximum_dimensionless_endpoint_change_factor | duration_b_predicted | duration_tension_sigma | time_dilation_gate_pass |
| --- | --- | --- | --- | --- | --- | --- |
| delta_m_p;sigma_spectral;tau_photon_distance | {"delta_m_p": 0.8313333333333331, "sigma_spectral": 0.9999999999999994, "tau_photon_distance": 0.956271567781008} | 3.262260 | 7.141842 | 1.003000 | 3.9721e-14 | PASS |
| delta_G;delta_alpha;tau_photon_distance | {"delta_G": 0.9959999999999993, "delta_alpha": -0.49999999999999956, "tau_photon_distance": 1.0406049011143415} | 3.422709 | 3.246867 | 1.003000 | 9.9301e-14 | PASS |
| delta_G;delta_m_p;delta_m_e | {"delta_G": 1.158419604457365, "delta_m_e": -0.9999999999999987, "delta_m_p": 0.9187901977713165} | 3.934323 | 34.554263 | 1.003000 | 1.9860e-14 | PASS |
| delta_m_p;delta_alpha;tau_photon_distance | {"delta_alpha": -0.49999999999999956, "delta_m_p": 0.4979999999999997, "tau_photon_distance": 1.2896049011143411} | 4.594477 | 3.246867 | 1.003000 | 7.9441e-14 | PASS |
| delta_m_p;delta_m_e;eta_Ni_luminosity | {"delta_m_e": -0.9999999999999993, "delta_m_p": 1.4979999999999987, "eta_Ni_luminosity": -0.5792098022286835} | 5.878287 | 34.554263 | 1.003000 | 1.3902e-13 | PASS |
| delta_m_p;delta_m_e;tau_photon_distance | {"delta_m_e": -0.9999999999999992, "delta_m_p": 1.4979999999999993, "tau_photon_distance": 0.28960490111434134} | 5.878287 | 34.554263 | 1.003000 | 1.9860e-14 | PASS |
| delta_h;delta_alpha;tau_photon_distance | {"delta_alpha": -0.666, "delta_h": -0.3319999999999998, "tau_photon_distance": 1.5386049011143408} | 6.167401 | 2.197877 | 1.003000 | 3.9721e-14 | PASS |
| delta_G;sigma_spectral;eta_Ni_luminosity | {"delta_G": 1.662666666666666, "eta_Ni_luminosity": -1.0812098022286825, "sigma_spectral": 0.9999999999999993} | 7.141842 | 7.141842 | 1.003000 | 5.9581e-14 | PASS |
| delta_G;sigma_spectral;tau_photon_distance | {"delta_G": 1.6626666666666667, "sigma_spectral": 0.9999999999999993, "tau_photon_distance": 0.540604901114341} | 7.141842 | 7.141842 | 1.003000 | 0 | PASS |
| delta_h;sigma_spectral;tau_photon_distance | {"delta_h": -0.7125714285714282, "sigma_spectral": 1.7125714285714286, "tau_photon_distance": 1.2531763296857705} | 7.575952 | 2.322302 | 1.003000 | 5.9581e-14 | PASS |
| delta_m_p;delta_m_e;delta_duration_extra | {"delta_duration_extra": 0.4344073516715123, "delta_m_e": -0.9999999999999987, "delta_m_p": 1.7876049011143404} | 8.278818 | 68.538821 | 1.003000 | 3.9721e-14 | PASS |
| delta_m_p;sigma_spectral;delta_duration_extra | {"delta_duration_extra": 1.4344073516715117, "delta_m_p": 1.7876049011143411, "sigma_spectral": 0.9999999999999993} | 8.278818 | 68.538821 | 1.003000 | 1.9860e-14 | PASS |

This branch is only an opacity/diffusion sensitivity. A credible Type Ia calculation must replace it with radiative transfer, nickel yield, ejecta structure, colour standardisation and selection effects.

## Gate matrix

| Gate | Result | Reason |
| --- | --- | --- |
| Algebraic stated closure | PASS | Exact two-coupling solutions exist. |
| Algebraic microscopic closure | PASS | Exact two-coupling solutions exist. |
| DES time dilation in diffusion control | PASS by construction | The third equation targets the measured central value. |
| Released-covariance supernova gate | FAIL | Best Peter-path composite has Delta BIC `+83.760921`. |
| Dimensionless external constraints | NOT YET TESTED | Requires an action and scalar history; dimensional exponents alone are not observables. |
| Action-level closure | FAIL | No scalar, matter and electromagnetic action was supplied. |
| Promotion | FAIL | Statistical and theory gates remain open. |

## What the search identifies

1. A single missing scalar power cannot repair both the redshift/atomic identity and the supernova photon/source relation.
2. At least two independent combinations are required before light-curve physics; the simplified diffusion equation requires a third.
3. The least elaborate phenomenological branch is `Q=constant` plus a small composite flux correction. Retaining `Q=(1+z)^-1` instead needs an additional inverse spectral factor, which is observationally degenerate at the redshift-identity level.
4. The supernova radial-shape failure remains. Algebraic closure does not make the existing distance law competitive.

## Required next derivation

Before further cosmological sampling, specify one covariant action with a dimensionless scalar and derive from it: (i) the observable spectroscopic redshift; (ii) atomic transition ratios; (iii) photon geodesics, arrival rates and opacity; (iv) dimensionless gravitational couplings; and (v) Type Ia source and diffusion laws. The resulting functions must be fixed before fitting and must reproduce local clocks, supernova time dilation and distance-duality controls.

## Reproduction

```powershell
python github_export/code/rawmu/search_peter_sparse_couplings_2026_07_19.py
```

## Primary research links

- [Bimetric varying-speed-of-light cosmology](https://arxiv.org/abs/gr-qc/0202012)
- [Varying-alpha action and dimensionless observables](https://arxiv.org/abs/gr-qc/0208081)
- [Dimensionless-constant interpretation caveat](https://arxiv.org/abs/hep-th/0208093)
- [DES supernova time-dilation measurement](https://arxiv.org/abs/2406.05050)
- [Type Ia radiative-transfer modelling](https://arxiv.org/abs/astro-ph/0609562)
