# SU(2)-CMB LINX Thermal-Bridge Registration

- **Locked:** `2026-07-17T11:45:18.3984414+10:00`
- **Status at lock:** no SU(2)-CMB LINX abundance result had been computed
- **Model identity:** literature-defined high-temperature SU(2)-CMB thermodynamic benchmark
- **Identity exclusion:** not the repository's late-time SU2/SU2R effective dark-sector fluid

## Fixed Thermodynamic Mapping

The benchmark implements the published asymptotic deconfining SU(2)-CMB relations

\[
\begin{aligned}
g_\gamma                         &= 2, \\
g_{\rm YM}                       &= 8, \\
\xi                              &= g_{\rm YM}/g_\gamma = 4, \\
\nu_{\rm CMB}                    &= \xi^{-1/3}=(1/4)^{1/3}, \\
T_\gamma(a)                      &= \nu_{\rm CMB}T_0/a, \\
\rho_{\rm YM}(T)                 &= \xi\rho_\gamma(T), \\
\rho_{\rm YM,extra}(T)           &= (\xi-1)\rho_\gamma(T).
\end{aligned}
\]

At Big Bang Nucleosynthesis temperatures, the six additional gauge polarisations are treated as ideal ultra-relativistic modes sharing the electromagnetic-sector temperature. Their energy density and heat capacity enter the LINX background differential equations. Standard neutrino-electron collision terms then determine the modified neutrino-to-photon temperature ratio used by the LINX neutron-proton weak rates.

## Fixed Microphysical Inputs

- LINX network: `key_PRIMAT_2023`.
- Weak rates: standard LINX weak rates evaluated on the computed SU(2)-CMB temperature trajectory.
- Neutron lifetime factor: fixed to one.
- Nuclear reaction pulls: all fixed to zero.
- Nuclear masses, binding energies and Q-values: standard LINX/PRIMAT values.
- Numerical controls: `sampling_nTOp=150`, `rtol=1e-5`, `atol=1e-9`.

The cited SU(2)-CMB thermodynamic model does not supply predicted shifts to electron radiative corrections, nuclear masses, binding energies, Q-values or reaction rates. Those quantities will not be fitted or reinterpreted as SU2 effects.

## Baryon-Normalisation Scan

The present baryon normalisation is scanned from `eta_fac=0.08` to `1.22` using 58 equally spaced values, supplemented by exact values at `0.25`, `1.0`, and the boundaries implied by the published SU(2)-CMB prior `0.014 <= omega_b h^2 <= 0.027` relative to the LINX reference `omega_b h^2=0.02242`. The factor `eta_fac/nu_CMB^3 = 4 eta_fac` describes the post-electron-positron-annihilation baryon density at fixed photon temperature. During weak freeze-out, the programme integrates the electron entropy explicitly rather than imposing this simple factor.

The primary physically anchored arm is the published `omega_b h^2` prior band. The unrestricted eta scan and the entropy-compensated `eta_fac=0.25` row are controls. They cannot establish support for SU(2)-CMB if the published-prior arm fails.

## Locked Abundance Gate

The observational anchors are

\[
\begin{aligned}
10^5({\rm D/H})                  &= 2.508 \pm 0.029, \\
Y_p                              &= 0.245 \pm 0.003, \\
10^{10}({}^7{\rm Li/H})         &= 1.45 \pm 0.25.
\end{aligned}
\]

A candidate row must place D/H, He-4 and Li-7/H each within two observational standard deviations. A lower total chi-squared without this joint pass is a near miss, not a lithium solution. Zero joint rows inside the published baryon prior vetoes this SU(2)-CMB asymptotic benchmark as a solution to the lithium problem.

## Programme Hashes

| File | SHA-256 |
|---|---|
| `run_su2cmb_linx_thermal_bridge_2026-07-17.py` | `2c0af9bf2adf750921bf594c64fd8b273088ee1ba04c6ec5674e6219cb1e1674` |
| `queue_su2cmb_linx_bridge_after_fr_2026-07-17.ps1` | `75a342cd83c52d8056b08c6f8fa19d234a4b577aee2792a5619d9cca7b8b00b8` |

## Sources

- Hahn and Hofmann, asymptotic temperature-redshift relation: https://arxiv.org/abs/1712.08561
- Hahn, Hofmann and Kramer, eight-mode SU(2)-CMB cosmological implementation: https://arxiv.org/abs/1810.01253
- Giovanetti et al., LINX: https://arxiv.org/abs/2408.14538

## Claim Boundary

This is a literature-defined thermodynamic benchmark propagated through a modern BBN network. It is not a complete first-principles SU(2) gauge-matter nucleosynthesis calculation. The result must remain separate from the late-time SU2/SU2R fluid unless an action-level derivation connects the two models and predicts the missing charged-matter and nuclear corrections.
