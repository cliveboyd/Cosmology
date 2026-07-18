# CMB Programme Implementation Report

**Date:** 18 July 2026
**Status:** staged implementation with hard physical-readiness gates

## Purpose

This programme turns five proposed CMB investigations into separate,
reproducible likelihood branches. The branches share typed theory interfaces but
do not share promotion criteria. A successful software regression is not treated
as evidence for SU(2), SU2R, magnetic energy injection, parity violation or FR.

## Common theory interface

Every physical model must declare which predictions it supplies:

1. background expansion and distances;
2. recombination and sound-horizon quantities;
3. scalar perturbation spectra;
4. tensor and parity-odd spectra; and
5. photon thermalisation or a distortion spectrum.

The interface refuses to infer a missing perturbation or distortion calculation
from a background curve. This prevents a phenomenological proxy from silently
acquiring a physical-model label.

## 1. Spectral distortions

The official 43-channel COBE/FIRAS monopole table is downloaded, hashed and
parsed. The baseline linear nuisance and distortion model is

\[
\begin{aligned}
\Delta I_\nu
  &= \Delta T\,G_\nu + \mu\,M_\nu + y\,Y_\nu
     + A_{\rm gal}F^{\rm gal}_\nu + A_Q Q_\nu, \\
-2\ln\mathcal{L}_{\rm FIRAS}
  &= (\mathbf d-\mathbf m)^{\mathsf T}
     \mathbf C^{-1}(\mathbf d-\mathbf m) + \text{constant}.
\end{aligned}
\]

Here `G_nu`, `M_nu` and `Y_nu` are the temperature-shift, chemical-potential and
Compton-y templates; `F_gal` is the released Galactic-pole template; and `Q_nu`
must be a registered physical prediction. A generic early-energy-injection run
must first map `dQ/dz` through a published thermalisation Green's function.
Magnetic or SU(2) labels require their own mapping to `dQ/dz` or to the photon
occupation equation.

The released channel-error column is sufficient for ingestion tests but is not a
full covariance. The likelihood therefore refuses a primary fit unless a
positive-definite 43 by 43 covariance is supplied. The diagonal approximation is
available only through an explicitly named diagnostic switch.

## 2. Low-l anomalies

The registered primary range is `2 <= l <= 30`, using Planck PR4/NPIPE Commander
and SEVEM products. The implementation sequence is:

1. freeze maps, confidence masks, split maps and simulation identifiers;
2. reproduce each map's low-l power spectrum and mask-transfer behaviour;
3. evaluate preferred-axis, dipole-modulation, hemispherical-power and even/odd
   parity statistics on every simulation;
4. define one maximum statistic over all maps, masks, multipole ranges and
   statistics; and
5. quote the global simulation rank, with component-separation and mask changes
   as stability gates.

No map-level statistic has yet been evaluated. The branch is blocked until the
PR4 assets, simulations and HEALPix analysis dependency are present. HHT remains
exploratory and cannot replace the registered map-level global null.

## 3. TB/EB polarisation parity

The map-level likelihood must use cross-spectra of independent splits. For a
constant rotation angle `alpha`, the software transform implements

\[
\begin{aligned}
C_\ell^{TB,\rm obs}
  &= C_\ell^{TE}\sin(2\alpha), \\
C_\ell^{TE,\rm obs}
  &= C_\ell^{TE}\cos(2\alpha), \\
C_\ell^{EB,\rm obs}
  &= \tfrac12(C_\ell^{EE}-C_\ell^{BB})\sin(4\alpha).
\end{aligned}
\]

The full fit must distinguish a cosmic angle from per-frequency instrument-angle
offsets and foreground TB/EB. A chiral-tensor hypothesis is separate and requires
parity-odd transfer spectra, conventionally parameterised by

\[
\chi(k)=\frac{P_R(k)-P_L(k)}{P_R(k)+P_L(k)}.
\]

This branch remains blocked pending split maps, polarisation-angle priors,
foreground choices, pseudo-`C_l` coupling matrices and physical parity-odd spectra.

## 4. CMB lensing

The official `act_dr6_lenslike==1.2.1` likelihood is the primary software gate.
The preregistered regression contains ACT-only and ACT+Planck baseline variants,
both with CMB-marginalised covariance and with likelihood corrections. All four
official fiducial chi-squared targets must be reproduced within 0.05.

The completed local regression gives

\[
\begin{aligned}
\chi^2_{\rm ACT,lens\ only}
  &= 14.057912, \\
\chi^2_{\rm ACT,corrected}
  &= 14.131500, \\
\chi^2_{\rm ACT+Planck,lens\ only}
  &= 21.071631, \\
\chi^2_{\rm ACT+Planck,corrected}
  &= 21.461799.
\end{aligned}
\]

All four cells pass the registered software/data gate.

The convergence convention is

\[
\begin{aligned}
\kappa_{LM}
  &= \tfrac12 L(L+1)\phi_{LM}, \\
C_L^{\kappa\kappa}
  &= \tfrac14[L(L+1)]^2 C_L^{\phi\phi}.
\end{aligned}
\]

An `A_kappa` rescaling is diagnostic only. Physical SU2R interpretation requires
registered scalar perturbation equations that predict `C_L^kappakappa` and lensed
TT, TE, EE and BB spectra. Subsequent ACT-only, Planck-only, combined and low/high-L
hold-outs must avoid double-counting Planck reconstruction information.

### SU2R physical-equation gate

The repository now contains an executable SU2R adapter contract and a dated
equation registry. The adapter requires, as one internally consistent system:

1. the covariant action and field content;
2. the background equations and parameter mapping;
3. scalar evolution and constraint equations;
4. perturbed energy density, pressure, momentum and anisotropic stress;
5. a gauge dictionary and regular initial conditions;
6. the metric/Weyl-potential relation used by lensing; and
7. a hash-registered Boltzmann-backend implementation.

The current registry intentionally leaves these entries empty because they were
not found in the supplied SU2R material. Its audit therefore returns `BLOCK
PHYSICAL SU2R CMB/LENSING INTERPRETATION`, and the adapter raises a hard exception
before spectra can be requested.

This gate also clarifies the status of the earlier code. The 66,150-row SU2R
growth scan is a phenomenological quasi-static closure in which clustering,
slip, sound-speed, gauge-range and helicity filters are assumed. Its best raw
improvement was only `Delta chi2=-0.0971103`, while its AIC-like score retained
the smooth branch. The historical CAMB PPF run is an effective-fluid proxy rather
than a non-Abelian perturbation implementation and was not matched to the source
background. Neither result is a physical SU2R CMB or lensing prediction.

## 5. Planck and ACT constraints

The Planck PR3 baseline chain has been compressed to

\[
\begin{aligned}
R
  &= \sqrt{\Omega_m}\,H_0D_M(z_*)/c, \\
\ell_A
  &= \pi D_M(z_*)/r_s(z_*), \\
\mathbf v
  &= (R,\ell_A,\omega_b h^2).
\end{aligned}
\]

The 24,497 retained weighted samples give

\[
\begin{aligned}
R
  &= 1.750638459 \pm 0.004572828, \\
\ell_A
  &= 301.7607006 \pm 0.0883883, \\
\omega_bh^2
  &= 0.0223597502 \pm 0.000148796.
\end{aligned}
\]

The maximum-posterior chain sample was independently recomputed with CAMB. The
largest discrepancy was 0.047574 prior standard deviations, passing the registered
validation gate. This compression is valid only when early radiation physics,
photon thermodynamics and recombination remain standard.

ACT DR6-lite has separately reproduced the collaboration's fixed-spectrum
log-likelihood: measured `-395.48307451` against `-395.48 +/- 0.01`. ACT is retained
as a full TT/TE/EE likelihood, not relabelled as a Planck-style distance prior. An
ACT distance compression would require a separate chain-level derivation and
out-of-sample validation.

## Current interpretation

The implementation makes background tests immediately usable and gives the
perturbation branches reproducible entry gates. It does not yet create an SU(2) or
FR CMB detection. The SU2R adapter and refusal gate are complete; its physical
backend is not. The highest-value theoretical task is therefore to supply the
registered field equations and implement their scalar hierarchy in a Boltzmann
solver using the same parameters as the late-time model. The ACT lensing software
gate is complete. The highest-value data task is acquisition of the registered
PR4 low-l map and simulation assets.

## Sources

- Planck PR3 ancillary products: https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/
- ACT DR6 products: https://act.princeton.edu/act-dr6-data-products
- ACT DR6-lite: https://github.com/ACTCollaboration/DR6-ACT-lite
- ACT DR6 lensing: https://github.com/ACTCollaboration/act_dr6_lenslike
- COBE/FIRAS monopole: https://lambda.gsfc.nasa.gov/product/cobe/firas_monopole_spect.html
- CLASS: https://class-code.net/
