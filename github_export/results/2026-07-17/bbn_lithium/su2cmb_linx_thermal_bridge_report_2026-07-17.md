# SU(2)-CMB to LINX Thermodynamic Bridge

Generated: 2026-07-17T14:35:19

## Identity and Scope

This calculation is a literature-defined high-temperature SU(2)-CMB benchmark. It is not the repository's late-time SU2/SU2R effective dark-sector fluid, and no equivalence between those models is assumed.

The implemented asymptotic relations are

\[
\begin{aligned}
g_\gamma                         &= 2, \\
g_{\rm YM}                       &= 8, \\
\xi                              &= g_{\rm YM}/g_\gamma = 4, \\
\nu_{\rm CMB}                    &= \xi^{-1/3}=(1/4)^{1/3}, \\
T_\gamma(a)                      &= \nu_{\rm CMB}T_0/a, \\
\rho_{\rm YM}(T)                 &= \xi\rho_\gamma(T), \\
H^2(T)                           &= \frac{8\pi G}{3}\left(\rho_{\rm EM}+\rho_\nu+\rho_{\rm YM,extra}\right).
\end{aligned}
\]

The six additional gauge polarisations are ideal ultra-relativistic modes at BBN temperature. Their energy density and heat capacity enter the LINX background integration. Standard neutrino-electron collision terms then predict a changed T_nu/T_gamma history, which is passed to the standard LINX neutron-proton weak rates.

All PRIMAT nuclear-rate pulls are fixed to zero. Nuclear masses, binding energies and Q-values are standard because the cited SU(2)-CMB thermodynamic model supplies no predicted modifications to them.

## Background Readout

- High-temperature coefficient: `nu_CMB = 0.629960525`.
- Expansion ratio at T_gamma = 1 MeV: `H_SU2CMB/H_SBBN = 1.25326`.
- Final LINX diagnostic N_eff: SBBN `3.0453`; SU2-CMB `20.373`.

## Abundance Readout

| Case | eta_fac | post-e-annihilation baryon factor | D/H x 1e5 | Yp | Li-7/H x 1e10 | chi2 D+He | chi2 all | depletion |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| SBBN reference | 1 | 1 | 2.4337 | 0.24699 | 5.6973 | 6.9948 | 295.62 | 3.9292 |
| SU2-CMB, eta_fac=1 | 1 | 4 | 0.64441 | 0.3165 | 26.476 | 4697.6 | 14719 | 18.26 |
| SU2-CMB, eta_fac near 1/4 | 0.25 | 1 | 7.3229 | 0.30765 | 1.9191 | 28003 | 28006 | 1.3235 |
| SU2-CMB, best D+He | 0.48 | 1.92 | 2.5288 | 0.31214 | 7.6474 | 501.41 | 1115.9 | 5.274 |
| SU2-CMB, best D+He+Li | 0.46 | 1.84 | 2.7138 | 0.31187 | 6.9664 | 547.23 | 1034.1 | 4.8044 |
| SU2-CMB, best in published omega_b prior | 0.62444 | 2.4978 | 1.6109 | 0.31376 | 12.962 | 1482.4 | 3602.7 | 8.9391 |

## Gate Interpretation

- Simultaneous D, He and Li two-sigma passes in the SU2-CMB eta scan: `0`.
- At fixed present reference eta, SU2-CMB predicts D/H x 1e5 `0.64441`, Yp `0.3165`, and Li/H x 1e10 `26.476`.
- Cancelling the fourfold entropy-induced baryon-density increase requires eta_fac near one quarter; that control gives D/H x 1e5 `7.3229` and Li/H x 1e10 `1.9191`.
- The unrestricted best D+He row occurs at eta_fac `0.48`; the unrestricted direct-lithium best row occurs at eta_fac `0.46`.
- Inside the published SU2-CMB omega_b prior, the best abundance chi-squared is `3602.7` at eta_fac `0.62444`.

## Claim Boundary

This bridge implements the published asymptotic SU(2)-CMB equation-of-state limit and propagates it through LINX background and standard weak-rate machinery. It is not a complete first-principles SU(2) gauge-matter BBN calculation: the cited model does not provide SU(2)-specific electron radiative corrections, nuclear binding shifts, Q-value shifts, or reaction-rate changes. Those inputs remain standard and are not fitted.

The repository's SU2/SU2R late-time fluid must remain scientifically separate unless an action-level derivation explicitly connects it to SU(2)-CMB thermodynamics.

## Sources

- High-temperature SU(2)-CMB relation: https://arxiv.org/abs/1712.08561
- Eight-mode cosmological implementation and omega_b prior: https://arxiv.org/abs/1810.01253
- LINX BBN framework: https://arxiv.org/abs/2408.14538
