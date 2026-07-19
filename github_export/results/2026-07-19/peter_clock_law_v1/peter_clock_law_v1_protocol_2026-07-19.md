# Peter clock-law v1: frozen audit protocol

**Registration date:** 19 July 2026  
**Status:** prospectively frozen before execution of the new physical-consistency matrix; local and not third-party timestamped

## Prior knowledge

The earlier fixed-path analysis was already known to give
`Delta BIC(Peter path - flat Lambda-CDM)=+94.344983` for `N=3422`.
This protocol therefore does not treat repetition of that result as new evidence.
The newly frozen tests concern the atomic-clock closure, photon/source split,
Chandrasekhar consequence, covariance sensitivities and escalation decision.

## Exact clock law

$$
\begin{aligned}
a(z)                         &= 1, \\
C(z)\equiv c_\gamma(z)/c_0 &= 1+z, \\
R(z)\equiv-dz/dT             &= H_0, \\
\chi(z)                      &= \frac{c_0}{H_0}
                                  z\left(1+\frac{z}{2}\right), \\
1+z                           &= \frac{C(z)}{a(z)Q(z)}, \\
Q_{\rm required}(z)           &= 1.
\end{aligned}
$$

Peter's asserted clock control is `Q_claim=(1+z)^-1`. It is tested, not
silently substituted into the spectroscopic-redshift identity.

## Flux and source branches

All deterministic branches use

$$
\begin{aligned}
D_{\rm eff}(z) &= \chi(z)(1+z)^s, \\
s               &= \tau-\frac{\ell}{2}, \\
L_{\rm Ia}(z)/L_{\rm Ia,0} &= (1+z)^\ell.
\end{aligned}
$$

Here `tau=0` is the stated no-post-emission-loss branch and `tau=2` is the
static Etherington-style distance-duality control. The registered luminosity
powers are `ell=0`, Peter's illustrative particle-count `ell=1`, and
`ell=7/2` from the Chandrasekhar scaling with fixed `hbar`, `G` and particle
mass ratio `m(z)/m0=(1+z)^-1`.

A free `s` is fitted only to diagnose the unidentified composite photon/source
law. Because supernova magnitudes identify `s`, not `tau` and `ell` separately,
that fit cannot validate either physical interpretation.

## Likelihood

The primary comparison uses Pantheon+ and DES `zHD`, the released Union3.1
redshift, released total covariance or precision, `z>0.01`, and one identical
profiled magnitude intercept per release for every model. Flat Lambda-CDM has
one fitted `Omega_m`; the exact Peter branches have no fitted shape parameter.

Grouped offset covariance reconstructions, release splits, redshift splits and
release-derived calibration-budget modes are sensitivities. Calibration modes
already contained in `STAT+SYS` are never added to that covariance again.

## Escalation gate

Pantheon+ light-curve-level refitting proceeds only if

1. `abs(Delta BIC)` for the exact no-loss, constant-luminosity branch is at most 2;
2. the asserted and derived atomic-clock laws agree;
3. the Compton-frequency and Chandrasekhar consequences are internally consistent;
4. photon transfer and source luminosity are separately predicted; and
5. an action-level clock, matter and electromagnetic model is supplied.

Failure constrains this explicit closure only. It does not exclude every
possible FR, varying-constant or static-universe construction.
