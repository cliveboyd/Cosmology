# FR Matter-to-Antimatter Identifiability Preregistration

Date: 2026-07-18

## Question

Can the currently implemented SPARC RAR/PLAMB rotation-curve likelihood identify whether a galaxy is constructed from matter or antimatter?

## Implemented Equations

\[
\begin{aligned}
\tau_g       &= \left(\frac{g_{{\rm bar},g}}{g_0}\right)^p, \\
g_{{\rm pred},g} &= \frac{g_{{\rm bar},g}}{1-\exp(-\tau_g)}.
\end{aligned}
\]

Every input is positive and the code contains no matter/antimatter state. Therefore, under the implemented likelihood,

\[
\begin{aligned}
\mathcal{L}_{\rm RC}(D_g\mid M,\theta)       &= \mathcal{L}_{\rm RC}(D_g\mid \bar M,\theta), \\
B^{\rm RC}_{\bar M,M}                         &= 1, \\
P(\bar M\mid D_{\rm RC})                     &= P(\bar M).
\end{aligned}
\]

The final equality holds only for the present rotation-curve likelihood. It does not rule out an FR matter-antimatter theory; it shows that an additional, independently derived transformation law is required before the likelihood can test one.

## Locked Transformation Matrix

| observable_or_term | matter_to_antimatter_transformation | implemented_in_sparc_likelihood | identifies_antimatter |
| --- | --- | --- | --- |
| Hydrogen/antihydrogen atomic spectrum | Invariant under CPT to current experimental precision | No particle-sign variable | No |
| Attractive gravitational rotation curve | No sign reversal in standard gravity; ALPHA-g is consistent with downward attraction | Positive gbar and positive predicted speed | No |
| PLAMB optical depth tau=(gbar/g0)^p | Unspecified by the supplied FR material | Same positive tau for every galaxy | No |
| PLAMB predicted acceleration gbar/(1-exp(-tau)) | Unspecified by the supplied FR material | Single branch with no charge-conjugation state | No |
| Matter-antimatter boundary annihilation | Produces an independent gamma-ray and secondary-particle signature at contact | Absent | Potentially, with external data |

## Identifiability Decision

**Gate 1 result: NOT IDENTIFIABLE in the current SPARC likelihood.**

No arbitrary sign flip will be fitted. A future FR antimatter branch must specify before data fitting:

1. which field or source variable changes under charge conjugation;
2. whether the change affects gbar, g0, kappa, p or a separate field term;
3. the predicted radial shape and sign;
4. its matter-matter, antimatter-antimatter and matter-antimatter limits; and
5. a compatible annihilation/contact model for external gamma-ray testing.

## External Physics Boundary

Antihydrogen spectroscopy is consistent with hydrogen, and ALPHA-g observes behaviour consistent with attractive terrestrial gravity. The independent astronomical discriminator is therefore expected at matter-antimatter contact, not in an otherwise isolated optical or 21-cm rotation curve.

References:

- https://www.nature.com/articles/nature21040
- https://www.nature.com/articles/s41586-023-06527-1
- https://arxiv.org/abs/2103.10073
- https://arxiv.org/abs/0808.1122
