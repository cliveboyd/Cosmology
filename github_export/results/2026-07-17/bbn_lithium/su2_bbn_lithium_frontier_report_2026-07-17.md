# SU2 BBN Lithium Near-Miss Frontier

Generated: 2026-07-17T14:29:03
Status: **complete**; successful unique points `30503 / 30503`, attempted `30503`, unresolved `0`.

## Scope

This is an exploratory diagnostic downstream of the separately registered SU2 completion gate. It does not alter that gate. It asks how close the sampled effective expansion model can move lithium while retaining deuterium, helium, Cosmic Microwave Background baryon-density and neutron-lifetime consistency.

The standardised residuals and frontier quantities are

\[
\begin{aligned}
z_i                         &= (X_i-\mu_i)/\sigma_i, \\
\chi^2_{\rm D+He}           &= z_{\rm D}^{2}+z_{\rm He}^{2}, \\
f_{\rm dep}                 &= ({}^7{\rm Li/H})_{\rm BBN}/\mu_{\rm Li}, \\
f_{2\sigma}                 &= ({}^7{\rm Li/H})_{\rm BBN}/(\mu_{\rm Li}+2\sigma_{\rm Li}).
\end{aligned}
\]

The abundance tolerances below require D/H and He-4 each, rather than only their combined chi-squared, to lie inside the stated interval. The uncertainties are the fixed observational errors used by the registered LINX scan; this table does not add theory-error floors.

## Two-Sigma Frontier

| Scenario | Rows passing D+He | Minimum Li/H x 1e10 | Li tension | Depletion to centre | Suppression to 2-sigma upper |
|---|---:|---:|---:|---:|---:|
| Expansion only | 16 | 5.555 | 16.42 sigma | 3.831 | 2.849 |
| Expansion + modest rate controls | 1224 | 5.113 | 14.65 sigma | 3.526 | 2.622 |
| Expansion + all scanned rate controls | 1505 | 4.924 | 13.9 sigma | 3.396 | 2.525 |
| All scanned controls | 4715 | 4.741 | 13.16 sigma | 3.27 | 2.431 |

## Interpretation

- Expansion only bottoms out at Li/H x 1e10 = `5.555` after D and He pass, still `16.42` observational sigma above the plateau centre.
- Modest selected-rate controls lower that floor to `5.113`, requiring a residual depletion factor `3.526` to reach the plateau centre.
- The unrestricted selected-rate scan reaches `4.924`, but still requires suppression by a factor `2.525` merely to enter the upper two-sigma lithium interval.
- A point on this frontier is a near miss, not a model preference. Nuclear-rate pull improvements remain controls unless a specified SU2 microphysical model predicts their signs and magnitudes.

## Claim Boundary

The result tests an effective non-negative expansion shift only. A first-principles SU2 claim would need a self-consistent temperature history and predicted changes to weak rates, binding energies, mass splittings and nuclear Q-values, followed by a new network calculation. The present frontier identifies the size and direction of that burden; it does not substitute for it.
