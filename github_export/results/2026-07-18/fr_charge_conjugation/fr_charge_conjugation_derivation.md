# Full Relativity Charge-Conjugation Derivation

Date: 2026-07-18
Generated: `2026-07-18T06:57:22.561008+00:00`

## Bottom Line

The supplied Full Relativity (FR) theory fixes a **charge-even galaxy-dynamics rule**, not an antigravity sign. Matter and antimatter make same-sign, non-negative contributions to total clout. Charge conjugation swaps the two chiral components, reverses the signed asymmetry, and leaves both total clout and the asymmetry magnitude unchanged. Because FR assigns stored-energy response to total clout and inertia to the magnitude of the asymmetry, a completely charge-conjugated antimatter galaxy and environment have the same predicted rotation curve as their matter counterpart.

Therefore, the current SPARC rotation-curve likelihood cannot identify an antimatter galaxy. The earlier excluded or persistent-tension galaxies remain dynamical anomalies only; an antimatter interpretation is neither selected nor made more probable by this derivation.

## Source Audit

The governing source used here is Peter R. Lamb's *Making Sense of Gravity*, version 8.10 (12 January 2023). Its component and fractional-asymmetry discussion is on PDF pages 82-85; the rotation-curve interpretation is on pages 96-99; total-clout dynamics are developed on pages 101-104; and the selected attractive-antimatter branch is on page 107. It states that clout arises from stored energy in both matter and antimatter, that light speed responds to their total contribution, that inertia responds to asymmetry, and that antihydrogen should fall rather than exhibit repulsive antigravity. The later 2024 introduction retains the same two-component construction. [Full FR source](https://www.fullyrelative.com/wp-content/uploads/2023/01/Making-Sense-of-Gravity-book-Vs8-10.pdf); [2024 FR introduction](https://www.fullyrelative.com/wp-content/uploads/2024/07/Making-Sense-of-Gravity-intro-to-Full-Relativity-Astro4.pdf).

This is an authorial theory document, not a peer-reviewed derivation of a completed field theory. It explicitly describes parts of the component-combination law as hypotheses requiring further modelling. The present result is consequently a faithful symmetry audit of the supplied FR proposal, not external validation of FR.

No FR action, Lagrangian, covariant field equation or stress-energy source term assigning charge-conjugation parity was found in the supplied research corpus or governing source. The rule below is therefore derived from FR's declared two-component source phenomenology. It must not be represented as a Noether-, CPT- or action-level result.

As an independent empirical boundary condition, ALPHA-g observed antihydrogen behaviour consistent with downward attraction to Earth and ruled out repulsive antigravity for that experiment. [ALPHA-g result](https://www.nature.com/articles/s41586-023-06527-1).

Local source copy used for the audit:

- SHA-256: `21c20bc61c060db1c0dca7478ab1b8cdedc218360c4bbd4f10fe647de3e91061`
- Bytes: `2261469`
- The source PDF is not committed to Git; the URL and checksum provide provenance.

## Component Definitions

Let the non-negative clout contributions at location $x$ be $B_M(x)$ from matter and $B_{\bar M}(x)$ from antimatter. The minimal source construction is

$$
\begin{aligned}
B_M(x)         &= B_{M,0}(x) + \int K(x,x')\,\rho_M(x')\,\mathrm{d}^3x', \\
B_{\bar M}(x) &= B_{\bar M,0}(x) + \int K(x,x')\,\rho_{\bar M}(x')\,\mathrm{d}^3x', \\
B_+(x)         &= B_M(x) + B_{\bar M}(x), \\
B_-(x)         &= B_M(x) - B_{\bar M}(x), \\
\chi(x)        &= \frac{B_-(x)}{B_+(x)}.
\end{aligned}
$$

FR proposes a point-source clout kernel proportional to $M/r$; its spatial derivative supplies the familiar $M/r^2$ acceleration scale. The derivation below does not require the kernel normalisation.

## Charge-Conjugation Operation

For a complete charge conjugation $\mathcal{C}$,

$$
\begin{aligned}
\mathcal{C}:\quad \rho_M              &\longleftrightarrow \rho_{\bar M}, \\
B_M                                             &\longleftrightarrow B_{\bar M}, \\
B_+                                             &\longrightarrow B_+, \\
B_-                                             &\longrightarrow -B_-, \\
\chi                                            &\longrightarrow -\chi, \\
|\chi|                                          &\longrightarrow |\chi|.
\end{aligned}
$$

The FR source assigns the speed of light and stored-energy response to $B_+$, while oscillation frequency and inertia depend on the magnitude of the component asymmetry. The most general response consistent with those statements is

$$
\begin{aligned}
c                  &= c(B_+), \\
m_g                &= m_g(B_+), \\
\eta                &\equiv \frac{m_i}{m_g} = f(|\chi|,B_+,v), \\
g_{\mathrm{dyn}} &= \frac{g_{\mathrm{bar}}}{\eta}.
\end{aligned}
$$

It follows directly that

$$
\begin{aligned}
\mathcal{C}[c]                  &= c, \\
\mathcal{C}[m_g]                &= m_g, \\
\mathcal{C}[\eta]               &= \eta, \\
\mathcal{C}[g_{\mathrm{dyn}}] &= g_{\mathrm{dyn}}.
\end{aligned}
$$

This is the definite FR charge-conjugation rule for galaxy dynamics: **rotation curves are invariant under complete matter-antimatter conjugation**.

## Fixed-Environment Exception

A host-only swap is not a complete charge conjugation. Write the environmental components as $B_M^{\mathrm{env}}$ and $B_{\bar M}^{\mathrm{env}}$, and let the host add non-negative clout $S(r)$. For host label $s=+1$ (matter) or $s=-1$ (antimatter),

$$
\begin{aligned}
B_+(r;s) &= B_M^{\mathrm{env}} + B_{\bar M}^{\mathrm{env}} + S(r), \\
B_-(r;s) &= B_M^{\mathrm{env}} - B_{\bar M}^{\mathrm{env}} + sS(r), \\
\chi(r;s) &= \frac{B_M^{\mathrm{env}}-B_{\bar M}^{\mathrm{env}}+sS(r)}
                  {B_M^{\mathrm{env}}+B_{\bar M}^{\mathrm{env}}+S(r)}.
\end{aligned}
$$

If the external background is exactly balanced, $|\chi(r;+1)|=|\chi(r;-1)|$. If a signed external imbalance is held fixed while only the host is swapped, the two magnitudes can differ. That is an **environment-relative interaction**, not an intrinsic antimatter rotation-curve signature. Testing it requires an external, preregistered estimator of the signed surrounding matter/antimatter field and a specified response $f$; SPARC baryonic and kinematic data alone provide neither.

## Relation to the Current PLAMB Likelihood

The implemented optical-depth bridge is

$$
\begin{aligned}
g_0                    &= \kappa\,\frac{cH_0}{2\pi}, \\
\tau                   &= \left(\frac{g_{\mathrm{bar}}}{g_0}\right)^p, \\
\eta_{\mathrm{PLAMB}} &= 1-\exp(-\tau), \\
g_{\mathrm{pred}}       &= \frac{g_{\mathrm{bar}}}{\eta_{\mathrm{PLAMB}}}.
\end{aligned}
$$

No signed component or host label occurs in these equations. The mapping is therefore charge-even and is compatible with the symmetry above as a phenomenological $\eta(g_{\mathrm{bar}})$ bridge. It is not, however, a derivation of $f(|\chi|,B_+,v)$ from the FR component field. Programme SHA-256: `0d976ff33b12a3fb80752c8719f8a7b6a488e457ef118c22ffd9d6ae270b4fbe`.

All `4` synthetic PLAMB evaluations were exactly charge-invariant.

## Deterministic Checks

| check | expected | observed | result |
| --- | --- | --- | --- |
| global_charge_conjugation | invariant | invariant | PASS |
| host_swap_in_balanced_background | invariant | invariant | PASS |
| host_swap_in_fixed_asymmetric_background | non-invariant | non-invariant | PASS |
| host_and_environment_charge_conjugation | invariant | invariant | PASS |

Overall expectation check: **PASS**.

The deliberately non-invariant row is the fixed asymmetric-background exception. It passes because non-invariance is the expected outcome when the environment is not conjugated.

## What Is Fixed and What Is Missing

Fixed by the supplied FR source and this derivation:

- matter and antimatter contribute with the same sign to total clout;
- complete charge conjugation reverses signed asymmetry but preserves its magnitude;
- attractive rather than repulsive matter-antimatter gravity is the selected branch;
- isolated rotation-curve dynamics are charge-even;
- an intrinsic SPARC antimatter classifier is prohibited.

Not fixed by the supplied FR source:

- the normalised cosmological kernel for $B_M$ and $B_{\bar M}$;
- the exact response $f(|\chi|,B_+,v)$;
- a map of the signed external antimatter background;
- a contact-annihilation likelihood at matter-antimatter domain boundaries.

These missing elements cannot be replaced by a freely fitted galaxy sign without making the model non-identifiable.

## Locked Decision Rule

1. Do not classify any of the six development galaxies or five reserved galaxies as antimatter from rotation curves.
2. Do not open the reserved replication set for an intrinsic matter-versus-antimatter SPARC test; the derived null prediction is exact under complete charge conjugation.
3. Permit an environment-conditioned test only after independently defining $B_M^{\mathrm{env}}-B_{\bar M}^{\mathrm{env}}$, the response $f$, nuisance priors, and a signed directional prediction before examining galaxy residuals.
4. Keep gamma-ray/contact-annihilation searches statistically separate. They test interfaces or mixing, not the charge-even rotation-curve mapping.

## Recommended Next Test

The defensible next step is an **FR environmental-asymmetry identifiability study**, not an antimatter label fit. Build a charge-blind large-scale environment proxy from group catalogues, neighbour mass and tidal fields; determine whether it predicts PLAMB residuals out of sample; then ask whether any additional signed matter/antimatter component is both externally observable and necessary. Until such a signed proxy exists, the antimatter hypothesis has no rotation-curve likelihood ratio.

No SPARC catalogue or reserved-galaxy record was opened by the derivation programme.
