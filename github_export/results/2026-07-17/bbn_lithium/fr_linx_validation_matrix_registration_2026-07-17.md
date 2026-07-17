# FR/LINX BBN Targeted Validation Matrix Registration

Date: 17 July 2026
Status: fixed before validation abundances were evaluated

## Objective

Validate the interpretation-bearing rows of the completed 30,503-point key-network scan against tighter numerical settings, larger networks, released deuterium-rate uncertainties, omitted lithium/beryllium channels and the updated neutron-lifetime centre.

## Fixed Numerical Settings

\[
\begin{aligned}
\mathrm{rtol}             &= 1e-07,\\
\mathrm{atol}             &= 1e-11,\\
N_{\rm interp,primary}    &= 300,\\
N_{\rm interp,check}      &= 600.
\end{aligned}
\]

The networks are `key_PRIMAT_2023`, `small_PRIMAT_2023` and `full_PRIMAT_2023`. The full network retains the LINX small-to-full switch at its upstream default temperature.

## Matrix

Total evaluations: **92**.

| Network | Core | D-rate response | Omitted Li/Be response | Sampling check | Total |
|---|---:|---:|---:|---:|---:|
| key | 10 | 12 | 0 | 2 | 24 |
| small | 10 | 12 | 0 | 2 | 24 |
| full | 10 | 12 | 20 | 2 | 44 |

The ten core states are the baseline; registered direct-lithium and lower-bound optima; expansion-only, modest-rate and unrestricted SU2-compatible optima; minimum lithium after the registered D+He gate; combined-current-anchor optimum; and baseline/direct-optimum states at the 2025 PDG neutron-lifetime centre.

The D-rate response uses independent unit pulls for `dpHe3g`, `ddHe3n` and `ddtp` around both the baseline and registered direct-lithium optimum in every network. The full-only completeness response uses unit pulls for `Li7paag`, `Be7naa`, `Be7daap`, `Be7pB8g` and `Li7daan`.

## Readout Rules

1. Numerical convergence compares strict key-network core rows with their stored `rtol=1e-5`, `atol=1e-9`, sampling-150 values.
2. Interpolation convergence compares sampling 300 and 600 at the baseline and registered direct-lithium optimum in each network.
3. Network sensitivity compares key, small and full predictions at identical strict settings.
4. Nuclear theory sensitivity is estimated from centred unit-pull responses. It is descriptive and does not assume that the response is globally linear.
5. Registered and combined-current abundance gates are reported, but this designed matrix is not treated as a new parameter search.
6. No FR or SU2 mechanism is inferred from a nuisance-rate response unless an explicit model predicts that coupled pull pattern.
