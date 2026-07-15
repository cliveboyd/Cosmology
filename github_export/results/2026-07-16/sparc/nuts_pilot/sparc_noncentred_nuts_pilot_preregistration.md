# SPARC Non-centred NUTS Pilot Preregistration

**Written before sampling:** `2026-07-15T20:33:53.860371+00:00`  
**Selected galaxies:** `24`  
**Models:** `RAR, PLAMB_OPTICAL_DEPTH_KAPPA_P`

## Fixed Scope

This is a sampler-geometry pilot, not a PLAMB-versus-RAR evidence calculation. The filtered all-Q2 inputs, physical likelihood, fixed H0, and July prior scales are unchanged.

The nuisance coordinates are transformed as

\[
\begin{aligned}
z_{D,g}       &\sim \mathcal{N}(0,1), & \log D_g       &= \mu_{D,g} + \sigma_{D,g} z_{D,g}, \\
z_{\eta,g}    &\sim \mathcal{N}(0,1), & \log \eta_g   &= \sigma_{\log\eta} z_{\eta,g}, \\
z_{\Upsilon d}&\sim \mathcal{N}(0,1), & \log\Upsilon_d&= \log(0.5)+\sigma_{\Upsilon}z_{\Upsilon d}, \\
z_{\Upsilon b}&\sim \mathcal{N}(0,1), & \log\Upsilon_b&= \log(0.7)+\sigma_{\Upsilon}z_{\Upsilon b}.
\end{aligned}
\]

The same four-standard-deviation nuisance truncation and original global physical bounds are retained.

## Outcome-blind Subset

The pilot subset is selected by SHA-256 ranking with seed `160726`, stratified only by published SPARC distance method and inclination below/above 60 degrees. No PLAMB-RAR residual or fit score enters selection.

Selected galaxies:

`DDO154`, `DDO161`, `ESO444-G084`, `ESO563-G021`, `F568-1`, `F568-V1`, `F571-8`, `KK98-251`, `NGC0300`, `NGC2998`, `NGC3949`, `NGC4013`, `NGC4157`, `NGC4183`, `NGC4217`, `NGC5055`, `NGC7331`, `UGC00731`, `UGC02023`, `UGC05750`, `UGC05764`, `UGC05999`, `UGC07608`, `UGC08837`

## Locked Gates

- likelihood-equivalence check must pass at absolute velocity tolerance `1e-08` km/s;
- zero post-warmup divergences;
- rank-normalised split R-hat target `<= 1.01` and absolute ceiling `<= 1.05`;
- minimum bulk and tail ESS `>= 400.0`; and
- minimum chain E-BFMI `>= 0.3`.

A model between the R-hat target and ceiling is labelled `CONDITIONAL_EXTEND_PILOT`, not ready for the full run. Any failure of equivalence, divergences, ESS, energy, or the R-hat ceiling is a pilot failure.
