# Covariant PLAMB clock-distance-ruler closure

Generated: 2026-07-19 09:28:43

## Decision

**KINEMATIC DERIVATION COMPLETE; PREDICTIVE PHYSICAL CLOSURE FAILS.**

A single covariant optical and ruler relation can be derived, but the current PLAMB inputs do not close it without additional clock, luminosity and ruler functions. Under the intended static background, the required functions are fixed algebraically and are large. DESI's radial/transverse ratio cannot be repaired by any isotropic ruler evolution.

This is a covariant kinematic completion, not an action-level FR field theory. A full theory must still derive the clock scalar dynamics, stress-energy and electromagnetic coupling.

## Covariant setup

Let `T(x)` be a scalar clock field with unit timelike flow `u_mu`. Matter measures proper time on the FLRW metric `g_mu_nu`, while geometric-optics photons follow the disformal metric

```text
g_gamma(mu,nu) = g(mu,nu) + [1-C(T)^2] u_mu u_nu
ds_gamma^2    = -C(T)^2 c0^2 dT^2 + a(T)^2 dSigma_K^2.
```

For a geometric-optics wave covector `k_mu`, the frequency measured by a comoving matter observer is `omega=-u^mu k_mu`, and spatial homogeneity gives `omega proportional to C/a`. Comparison with the local atomic transition therefore introduces the additional ratio `Q(T)=nu_atom(T)/nu_atom(T0)`.

Here `C=c_gamma/c0`. Let `Q(T)=nu_atom(T)/nu_atom(T0)` describe the evolution of the atomic transition used to define spectroscopic redshift. With present values `a0=C0=Q0=1`,

$$
\begin{aligned}
1+z                         &= \frac{C(T)}{a(T)Q(T)}, \\
R(z)\equiv-\frac{dz}{dT} &= (1+z)\left[H+\frac{\dot Q}{Q}-\frac{\dot C}{C}\right], \\
T_0-T(z)                    &= \int_0^z\frac{du}{R(u)}.
\end{aligned}
$$

The registered phenomenology supplies only `C(z)=1+gamma_c z` and `R(z)=H0(1+z)^p`. It does not separately determine `a(T)` and `Q(T)`.

## Optical and BAO distances

Let `B(z)` be the evolution of the comoving BAO ruler, `r_BAO(z)=r_d B(z)`, and let `Xi(z)` collect photon survival and standardised intrinsic luminosity evolution. The flat-background equations are

$$
\begin{aligned}
\chi(z)                   &= \int_0^z\frac{c_0 C(u)}{a(u)R(u)}\,du, \\
D_A(z)                      &= a(z)S_K[\chi(z)], \\
\frac{D_M^{\rm BAO}}{r_d} &= \frac{S_K(\chi)}{r_d B(z)}, \\
\frac{D_H^{\rm BAO}}{r_d} &= \frac{c_0 C(z)}{a(z)R(z)r_d B(z)}, \\
D_L^{\rm SN}(z)           &= (1+z)^2D_A(z)\,\Xi(z) \\
                              &= (1+z)^2a(z)B(z)D_M^{\rm BAO}(z)\,\Xi(z).
\end{aligned}
$$

This is the requested single relation connecting `T`, `D_L`, `D_M`, `D_H` and ruler evolution. For curvature, `S_K` is the usual transverse-curvature map. Differentiation gives

$$
\begin{aligned}
\frac{dD_M^{\rm BAO}}{dz}+\frac{d\ln B}{dz}D_M^{\rm BAO}
&= C_K(\chi)D_H^{\rm BAO}.
\end{aligned}
$$

Thus `D_H=dD_M/dz` is valid only for flat geometry and a non-evolving comoving ruler.

## Relation to the fitted PLAMB path

The supernova programme fitted

$$
\begin{aligned}
D_{\rm path}(z) &= \int_0^z\frac{c_0C(u)}{R(u)}\,du, \\
D_L^{\rm fit}   &= D_{\rm path}.
\end{aligned}
$$

`D_path` equals the comoving optical distance `chi` only when `a=1`. Under that static closure with a fixed BAO ruler, consistency forces

$$
\begin{aligned}
Q(z)                            &= \frac{1+2.3z}{1+z}, \\
\Xi(z)                         &= (1+z)^{-2}, \\
P_\gamma(z)\frac{L_*(z)}{L_{*,0}} &= \Xi(z)^{-2}=(1+z)^4.
\end{aligned}
$$

Consequently, alpha=0 is not photon-conserving Etherington propagation. It requires a specific composite luminosity or photon-transfer law. At fixed atomic standards (`Q=1`), the same gamma law instead requires `a=C/(1+z)` and gives `H(0)/H0=1-gamma_c=-1.3`, a contracting present background.

## Endpoint requirements at z=2.33

| closure               |       a |   Q_atomic_frequency_ratio |   H_over_H0 |   Xi_luminosity_transfer |   J_photon_survival_times_luminosity_ratio |
|:----------------------|--------:|---------------------------:|------------:|-------------------------:|-------------------------------------------:|
| static_fr             | 1       |                    1.90961 |    0        |                0.0901803 |                                   122.964  |
| fixed_atomic_clock    | 1.90961 |                    1       |   -0.160809 |                0.0776504 |                                   165.849  |
| standard_scale_factor | 0.3003  |                    6.359   |    0.786605 |                0.130603  |                                    58.6267 |

## Numerical identity checks

| closure               |   maximum_redshift_identity_error |   maximum_rate_identity_error |   maximum_dchi_dz_minus_DH_error |   maximum_luminosity_identity_error |   minimum_a |   maximum_a |   H_over_H0_at_z0 |
|:----------------------|----------------------------------:|------------------------------:|---------------------------------:|------------------------------------:|------------:|------------:|------------------:|
| static_fr             |                       2.22045e-16 |                   8.88178e-16 |                      1.8656e-07  |                         4.44089e-16 |    1        |    1        |               0   |
| fixed_atomic_clock    |                       2.22045e-16 |                   8.88178e-16 |                      1.33008e-08 |                         8.88178e-16 |    1        |    1.92857  |              -1.3 |
| standard_scale_factor |                       4.44089e-16 |                   8.88178e-16 |                      6.31513e-08 |                         8.88178e-16 |    0.285714 |    0.999999 |               1   |

## DESI ruler tests

| model                                 |   ruler_power_b |   B_z2p33 |   q_c0_over_H0rd |      chi2 |        BIC |   delta_BIC_vs_LCDM |   goodness_upper_tail |
|:--------------------------------------|----------------:|----------:|-----------------:|----------:|-----------:|--------------------:|----------------------:|
| static_fr_constant_ruler              |         0       |   1       |          7.5861  | 33458.1   | 33460.6    |            33445.2  |              0        |
| static_fr_power_law_ruler             |         1.66032 |   7.36929 |         35.4304  |  4144.57  |  4149.7    |             4134.3  |              0        |
| fixed_atomic_clock_constant_ruler     |         0       |   1       |         13.6112  | 20452.1   | 20454.7    |            20439.3  |              0        |
| fixed_atomic_clock_power_law_ruler    |         1.21842 |   4.33068 |         39.4881  |  2684.67  |  2689.8    |             2674.4  |              0        |
| standard_scale_factor_constant_ruler  |         0       |   1       |          2.30395 | 65508.8   | 65511.4    |            65496    |              0        |
| standard_scale_factor_power_law_ruler |         2.70658 |  25.944   |         39.4077  |  8124.03  |  8129.16   |             8113.75 |              0        |
| LCDM_flat                             |         0       |   1       |         29.5246  |    10.271 |    15.4009 |                0    |              0.506184 |

A power-law ruler `B(z)=(1+z)^b` cannot rescue any closure. More generally, allowing an independent isotropic ruler amplitude in every DESI redshift bin leaves the following irreducible Alcock-Paczynski radial/transverse scores:

| closure               |   irreducible_AP_chi2 |   maximum_B_DM_over_B_DH |   minimum_B_DM_over_B_DH |   B_relative_at_z2p33 |
|:----------------------|----------------------:|-------------------------:|-------------------------:|----------------------:|
| fixed_atomic_clock    |               2217.32 |                 0.789996 |                 0.469407 |               3.0584  |
| standard_scale_factor |               7708.33 |                 0.579676 |                 0.279088 |              12.6731  |
| static_fr             |               3881.95 |                 0.687562 |                 0.404186 |               4.73756 |

The ratio `B_DM/B_DH` would equal one if a single isotropic ruler could fit both directions. Its departure from one is independent of the absolute sound horizon and cannot be corrected by choosing another `r_d`.

## Gates

- `static_fr_redshift_identity`: **PASS**
- `static_fr_rate_identity`: **PASS**
- `static_fr_radial_derivative_identity`: **PASS**
- `static_fr_luminosity_identity`: **PASS**
- `fixed_atomic_clock_redshift_identity`: **PASS**
- `fixed_atomic_clock_rate_identity`: **PASS**
- `fixed_atomic_clock_radial_derivative_identity`: **PASS**
- `fixed_atomic_clock_luminosity_identity`: **PASS**
- `standard_scale_factor_redshift_identity`: **PASS**
- `standard_scale_factor_rate_identity`: **PASS**
- `standard_scale_factor_radial_derivative_identity`: **PASS**
- `standard_scale_factor_luminosity_identity`: **PASS**
- `static_background_matches_FR_intent`: **PASS**
- `static_atomic_clock_evolution_is_derived_not_supplied`: **FAIL**
- `static_luminosity_transfer_is_derived_not_supplied`: **FAIL**
- `static_Xi_close_to_unity_at_z2p33`: **FAIL**
- `fixed_atomic_clock_closure_is_noncontracting_today`: **FAIL**
- `power_law_ruler_any_closure_goodness_p_at_least_0p05`: **FAIL**
- `arbitrary_isotropic_ruler_any_closure_AP_p_at_least_0p05`: **FAIL**
- `independent_drag_epoch_ruler_derivation`: **FAIL**
- `action_level_clock_and_photon_dynamics`: **FAIL**

## Research implication

Further sampling is not justified for the present `gamma_c=2.3`, `p=0.800469`, alpha-zero branch. The next theory submission must specify, before seeing another distance dataset: (1) the action or field equations for `T`; (2) whether matter and photons share a metric; (3) the atomic-frequency law `Q(T)`; (4) photon conservation and intrinsic SN luminosity, fixing `Xi(T)`; and (5) the drag-epoch calculation and late-time transport of `B(T)`. These are physical predictions, not nuisance parameters that may all be fitted to the same SN and BAO data.

## Primary references

- [Moffat bimetric VSL construction](https://arxiv.org/abs/gr-qc/0202012)
- [Disformal invariance of cosmological observables](https://arxiv.org/abs/2003.10633)
- [Etherington distance-duality assumptions](https://arxiv.org/abs/1612.08784)
- [Joint-constant requirements in a covariant VSL model](https://arxiv.org/abs/2011.09274)
- [DESI DR2 publications and products](https://data.desi.lbl.gov/doc/papers/dr2/)
- [Peter R. Lamb, Making Sense of Gravity, v8.10](https://www.fullyrelative.com/wp-content/uploads/2023/01/Making-Sense-of-Gravity-book-Vs8-10.pdf)
