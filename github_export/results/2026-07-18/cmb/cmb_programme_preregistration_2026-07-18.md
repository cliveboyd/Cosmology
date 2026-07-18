# CMB Programme Preregistration

**Locked:** 18 July 2026, before downloading the Planck chain bundle or
evaluating the ACT DR6-lite regression target.

## Scope

The five branches remain statistically separate. Passing one branch cannot
compensate for failure in another. A phenomenological template is not labelled
an SU(2), SU2R or FR prediction unless the relevant theory adapter supplies the
required background, perturbation, parity or photon-thermalisation equations.

## Immediate Branch A: Planck Distance Compression

- Source: official Planck PR3 baseline chain bundle
  `COM_CosmoParams_base-plikHM-TTTEEE-lowl-lowE_R3.00.zip`.
- Registered vector order: `(R, l_A, omega_b h^2)`.
- Use every finite, positive-weight sample in every chain text file belonging
  to the bundle. Do not tune burn-in or thin after inspecting the vector.
- Use the exact derived `thetastar` chain column, not the approximate sampled
  `theta_MC` parameter.
- For each sample,

  ```text
  theta_star = thetastar / 100
  l_A        = pi / theta_star
  D_M(z*)    = rstar / theta_star
  R          = sqrt(Omega_m) H0 D_M(z*) / c
  omega_bh2  = omegabh2
  ```

- Weighted means and the weighted posterior second central moment define the
  Gaussian compression. The output must be finite, symmetric, positive
  definite and accompanied by source/member hashes.
- This compression is valid only for branches retaining standard early
  physics and recombination. It is diagnostic-only for any model changing the
  photon temperature law, radiation content, sound speed or recombination.

## Immediate Branch B: ACT DR6-lite Regression

- Use the unmodified collaboration likelihood and local official v1.0 SACC
  data, with CAMB and the collaboration's test parameters:
  `ombh2=0.022`, `omch2=0.117`, `ns=0.96`, `As=2e-9`, `tau=0.065`,
  `cosmomc_theta=0.010409`, `A_act=1`, and `P_act=1`.
- Required implementation gate: total ACT TT/TE/EE log-likelihood
  `-395.48 +/- 0.01`.
- This is a software/data regression, not a cosmological fit or evidence claim.

## Deferred Branch C: FIRAS Spectral Distortions

- Freeze the FIRAS monopole residual vector, covariance treatment and Galactic
  spectral nuisance basis before fitting.
- Baseline parameters are temperature shift, `mu`, `y` and registered
  foreground/calibration modes.
- Generic injection uses fixed logarithmic redshift bins and a published
  thermalisation response. Magnetic and SU(2) labels require explicit
  mappings to `dQ/dz` or the photon occupation equation.
- No scientific evaluation occurs in this implementation phase because the
  FIRAS likelihood assets and model mappings are not yet registered locally.

## Deferred Branch D: Low-l Isotropy

- Primary range: `2 <= l <= 30`; Planck PR4/NPIPE Commander and SEVEM maps.
- Register masks, component maps, split maps, simulation ensemble and the
  complete statistic family before evaluating data.
- Preferred-axis, dipole-modulation, hemispherical-power and even/odd-parity
  statistics share a maximum-statistic global null.
- HHT remains exploratory and cannot replace the map-level null.

## Deferred Branch E: TB/EB Parity

- Use cross-spectra of independent map splits and a registered pseudo-`C_l`
  estimator.
- Jointly fit cosmic rotation, per-frequency instrument-angle offsets,
  polarized dust/synchrotron TB/EB, beams and transfer functions.
- Constant birefringence and chiral-tensor templates are distinct hypotheses.
  A chiral claim requires registered parity-odd tensor transfer spectra.

## Deferred Branch F: CMB Lensing

- First reproduce the official ACT DR6 lensing likelihood under Lambda-CDM.
- A phenomenological `A_kappa` branch is diagnostic only.
- A physical SU2R comparison requires registered scalar perturbation equations
  sufficient to predict `C_L^kappakappa` and lensed TT/TE/EE/BB spectra.
- ACT-only, Planck-only and ACT+Planck variants and low/high-`L` hold-outs must
  be fixed without double-counting Planck reconstruction information.

## Primary Sources

- Planck PR3 ancillary products:
  https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/
- ACT DR6 products: https://act.princeton.edu/act-dr6-data-products
- ACT DR6-lite: https://github.com/ACTCollaboration/DR6-ACT-lite
- ACT DR6 lensing: https://github.com/ACTCollaboration/act_dr6_lenslike
- COBE/FIRAS monopole:
  https://lambda.gsfc.nasa.gov/product/cobe/firas_monopole_spect.html
- CLASS: https://class-code.net/
