# Results

Small result tables and Markdown readouts copied here are intended to be cited by research documents.

Recommended pattern:

- `results/YYYY-MM-DD/topic_name/summary.md`
- `results/YYYY-MM-DD/topic_name/*.csv`
- `results/YYYY-MM-DD/topic_name/manifest.csv`

Large posterior chains, object-level mock catalogs, logs, and generated figure files should stay in local run folders and be listed in a manifest instead.

For SkyMap, a small number of representative figures and extracted DOCX media are committed because the historical reports are plot-heavy; the full Office binaries and large FITS/NPZ products remain outside Git.

## July 16 Audit

- `2026-07-16/su2_quaia/` contains the spatial-block null, template-collinearity, and over-absorption/injection-calibration readouts.
- `2026-07-16/sparc/` contains robust split convergence, paired-posterior nuisance diagnostics, the preregistered 24-galaxy non-centred NUTS pilot under `nuts_pilot/`, and the separately preregistered 157-galaxy run under `nuts_full/`.
- `2026-07-16/rawmu/` contains the FR `c(z)=c0(1+z)` no-loss dimming sequence for raw supernova MU, the alpha=0 calibration-budget ladder and direct LCDM comparison, and the follow-up catalogue-combination, redshift-band, and calibration-pull gates under the July 14 calibration/redshift-frame controls.
- `2026-07-16/bbn_lithium/` contains the first BBN Lithium anomaly scoping gate for FR/no-expansion ideas, including eta, early-clock, Li-production, and stellar-depletion stress tests.
- `2026-07-16/audit/` contains the cross-stream model-upgrade reproducibility audit and research-standing priorities.

## July 17 Audit

- `2026-07-17/bbn_lithium/` contains the completed 30,503-point FR/LINX scan audit, immutable-catalogue score and coverage checks, final registered SU2 lithium gate/frontier readouts, archived run configuration, and a hash manifest for the local catalogue and committed analysis products.
- The detailed audit fixes the claim boundary at an FR-inspired standard-BBN clock/radiation proxy. It also records the priority requirements for an FR-native thermal history, a matched CMB likelihood, full-network and tighter-tolerance validation, deuterium theory-error marginalisation, and hierarchical lithium-depletion modelling.
- The same folder contains the completed, separately registered 62-point literature-defined SU(2)-CMB thermodynamic bridge. Its report, trajectory, eta scan, figure, configuration and audit remain explicitly separate from the repository's late-time SU2/SU2R fluid and from the FR proxy.
- The completed preregistered 92-evaluation FR/LINX validation matrix contains tighter numerical tolerances, key/small/full network comparisons, deuterium-rate theory responses, full-network omitted Li/Be-channel checks, interpolation checks and updated-neutron-lifetime states. All evaluations succeeded; numerical/network corrections are small, while the local D-rate theory uncertainty is comparable to the observational D/H error and must enter the next likelihood.

## July 18 Audit

- `2026-07-18/cmb/` contains the registered five-branch CMB programme, a typed
  theory boundary, the official Planck PR3 three-variable distance compression,
  independent CAMB validation, ACT DR6-lite fixed-spectrum regression, ACT DR6
  lensing regression registration, FIRAS ingestion gate and branch-readiness audit.
- The Planck compression uses 24,497 weighted samples and gives
  `R=1.750638459`, `l_A=301.7607006` and `omega_b h^2=0.0223597502`; independent
  CAMB reconstruction agrees within 0.047574 prior standard deviations. ACT
  DR6-lite reproduces the collaboration target with log-likelihood
  `-395.48307451`. Both are implementation results, not SU(2) or FR evidence.
- The official ACT DR6 lensing v1.2 four-cell regression also passes. ACT-only
  lens-only/corrected chi-squared values are `14.057912` and `14.131500`; the
  ACT+Planck values are `21.071631` and `21.461799`. Physical SU2R lensing remains
  blocked until registered scalar perturbations predict convergence and lensed
  primary-CMB spectra.
- The CMB programme explicitly blocks physical lensing, low-l, TB/EB and FIRAS
  claims until the corresponding perturbation, map/simulation, polarisation-angle,
  covariance and photon-thermalisation requirements are supplied.
- A hash-registered SU2R physical-perturbation adapter now makes that boundary
  executable. The current equation registry correctly blocks primary-CMB and
  lensing spectra because the covariant action, scalar field equations, perturbed
  stress-energy, gauge/initial-condition prescription and Weyl relation have not
  been supplied. Earlier quasi-static growth and CAMB PPF effective-fluid scans
  remain labelled phenomenological proxies, not physical SU2R exclusions.
- `2026-07-18/sparc_am/` contains the preregistered SPARC matter/antimatter-identifiability gate, outcome-blind matched controls, twelve-scenario conventional-systematics matrix, deterministic multi-start profiles, figures and hash manifest.
- The implemented rotation-curve likelihood is matter/antimatter non-identifiable. Three of the six development galaxies retain a model-relative PLAMB tension, but the matched-block direction is not coherent and no galaxy is classified as antimatter.
- `2026-07-18/sparc_am_single_start_invalidated/` preserves the initial single-start bundle and its invalidation notice. It is retained only as a transparent numerical-audit record and must not be used for interpretation.
- `2026-07-18/fr_charge_conjugation/` derives the matter/antimatter transformation from Peter R. Lamb's supplied Full Relativity source and verifies it with synthetic invariance checks. Complete charge conjugation preserves total clout, asymmetry magnitude and galaxy dynamics; it therefore rules out an intrinsic SPARC rotation-curve antimatter classifier.
- The only conditional non-invariance occurs when a host is swapped while a signed asymmetric environment is held fixed. That environment-relative branch remains untestable until an external signed-background estimator and the FR inertia response are specified in advance; the five reserved SPARC galaxies remain unopened.
- `2026-07-18/fr_environment_asymmetry/` contains the independently committed preregistration, charge-blind NED/2MRS environmental covariates, deterministic profile fits, cross-fitted tests, gate decision, post hoc robustness audit, figure and hash manifest.
- The preregistered 66-galaxy primary development test found a suggestive negative partial association between local environment and the PLAMB-minus-RAR profile outcome (rho=-0.3073, two-sided permutation p=0.01180), but failed the p<=0.01 and >=5 per cent predictive-improvement gates (observed improvement 2.40 per cent). Sensitivity weakening prevents promotion, and the reserved galaxies remain sealed.
- Post hoc decomposition indicates that the pattern is local rather than broad: the 2 Mpc count and softened tidal proxy each have exploratory Holm-adjusted p=0.0478, while the 5 Mpc count is null. Leave-one-out rho remains negative for all 66 removals. These results motivate a preregistered conventional-environment nuisance study with an independent group catalogue; they do not identify antimatter or revise the failed gate.
- `2026-07-18/fr_group_catalogue_replication/` contains the separately committed independent-predictor preregistration, Kourkchi-Tully/VizieR group matching audit, 65-galaxy primary test, sensitivity matrix, figure and hash manifest. It is explicitly not a new held-out-outcome replication because it uses the previously frozen development outcomes.
- Published group richness is null (rho=0.0477, p=0.7044) and worsens cross-fitted RMSE by 5.41 per cent. Catalogue group mass and grouped/single status are also null, while the shifted-sky control is larger than the actual association. The gate fails, broad group membership does not replicate the earlier local-density hint, and the five reserved galaxies remain sealed.
- `2026-07-18/fr_group_power_calibration/` contains the separately committed 40,000-draw injection-and-recovery calibration, 100,000-permutation rho threshold, complete compressed simulation draws, power curves, empirical-null audit and hash manifest.
- At the earlier |rho|=0.307254 anchor, the locked joint primary recovery probability is 42.08 per cent and the empirically recalibrated value is 36.50 per cent, well below the preregistered 80 per cent standard. Reliable joint recovery begins at about |rho|=0.5. The nominal rank-permutation p component is anti-conservative under the complete empirical-noise pipeline (3.26 per cent rather than 1 per cent null firing), although the full joint null rate remains conservative at 0.54 per cent. The group null therefore weakens but does not exclude an earlier-scale effect carried by this coarse proxy.
- `2026-07-18/su2_quaia/` contains the preregistered SDSS DR16Q v4 independent-catalogue validation, 20,000-draw spatial-block global nulls at three scales, multipole and stability gates, and the injection-calibrated held-out angular likelihood.
- The DR16Q angular family is globally unusual but fails the locked cross-catalogue direction by 116.205 degrees, two redshift-window direction gates and one amplitude-stability gate. It is not promoted; ordinary individual `l=1-3` partial probabilities leave survey-footprint coupling as the leading interpretation.
- The first injection-likelihood branch is retained with a formal invalidation notice because it evaluated the candidate predictor at the response photometric redshift. Its leakage-free 46,687-object Quaia x DR16Q replacement passes the exact predictor-invariance and false-positive audits, but is underpowered: one-times predictive detection is 8.0 per cent, joint recovery is zero, and the observed held-out `Delta CVLPD=-4.0084` has empirical `p=0.7273`.
- `2026-07-18/rawmu/` contains the release-grounded Pantheon+, DES-Dovekie and Union3.1 raw-MU hierarchy, covariance reconstruction audit, external calibration-budget modes and dataset, survey and high-redshift hold-outs with `zHD` primary.
- The 3,422-supernova primary total-covariance comparison gives `Delta BIC(FR-Lambda-CDM)=+94.344983`, favouring Lambda-CDM for the fixed FR distance law. The strict robustness gate still fails because only 58.82 per cent of eligible survey hold-outs preserve the sign and several model-difference modes exceed external calibration budgets; this blocks a survey-robust promotion without reversing the primary likelihood direction.
- `2026-07-18/rawmu_plamb_axis/` contains Peter Lamb's requested three-release
  axis comparison, its outcome-aware preregistration, coordinate and residual
  data, plots, fit table, hash manifest and a draft technical reply. Dividing
  luminosity distance by `(1+z/2)` is verified to be exactly the same fixed PLAMB
  forward law, to `1.421e-14` mag or better. The released-total-covariance result
  remains `Delta BIC(PLAMB-Lambda-CDM)=+94.344983`; every release separately and
  all registered diagonal-distance and all-positive-redshift sensitivities retain
  the same direction.
- The clock-law consistency audit in the same folder replaces the previously
  separate distance integrands with one analytic kernel and explicitly registers
  the fixed raw-MU law as `p=0` and the historical no-expansion pilot as `p=1`.
  All seven analytic, quadrature and wiring gates pass. The full likelihood
  regression changes `Delta BIC` by only `9.095e-13`, so the refactor changes no
  scientific result and does not select either clock law physically.
- `2026-07-18/rawmu_plamb_clock_fit/` contains the separately preregistered nested
  clock-law attribution ladder, deterministic full-covariance fits, complete-
  release shape tests, exact conditional survey/high-redshift tests,
  identifiability audit, figure and hash manifest. The best PLAMB-family BIC cell
  (`PATH_FREE`) remains worse than flat Lambda-CDM by `Delta BIC=+13.345052`,
  reaches the registered `gamma_c=2` boundary, loses all three release-shape
  hold-outs and has a `+321.926337` high-redshift conditional loss at `z>=0.5`.
  The fully general curve has an extreme `p-alpha` correlation of `0.997363` and
  is not a physical identification of clock-rate or flux evolution. Seven of the
  eight promotion gates fail; the result is explicitly **do not promote**.
- The post hoc `boundary_highz_followup/` subfolder resolves the original
  `gamma_c=2` numerical boundary: the full-sample PATH minimum is interior at
  `gamma_c=2.276517`, `p=0.784634`. Its raw chi-squared is nearly identical to
  Lambda-CDM, but it still loses by `Delta BIC=+8.111201`. Independently fitting
  `z<0.5` selects `gamma_c=1.524995`, `p=0.246775` and then loses the exact
  conditional `z>=0.5` prediction by `Delta chi-squared=+321.928234`.
- Union3.1's final compressed node at `z=2.262260` contributes `+246.280504` of
  that excess, but an explicitly post hoc exclusion/refit still leaves
  Pantheon+ plus DES at `Delta chi-squared=+38.606432`. All covariance
  sensitivities retain the negative transfer direction. The boundary checks
  pass, while the primary transfer, redshift-band and concentrated-innovation
  gates fail; the diagnostic remains **fail and cannot promote**.
- `2026-07-18/priority_followup_outcomes_2026-07-18.md` is the integrated claim-boundary readout for these three branches.

## July 19 Audit

- `2026-07-19/rawmu_plamb_gamma2p3/` contains the explicitly outcome-aware
  fixed-`gamma_c=2.3` conditional analysis, including equal-parameter PATH versus
  Lambda-CDM fits, deterministic low-redshift posterior grids, release-intercept-
  integrated high-redshift prediction, release-specific `p` stability, physical
  ratio readout, figure and hash manifest.
- With `gamma_c=2.3` treated as externally fixed, PATH and Lambda-CDM are
  statistically tied: full-sample `Delta BIC=+0.011958`, low-redshift
  `Delta BIC=-0.379408`, and integrated high-redshift
  `Delta log predictive density=-0.240554`. The PATH posterior is
  `p=0.800469+/-0.023641`; release-specific values agree within `0.779355`
  standard deviations.
- Six of seven registered statistical gates pass, but the high-redshift log-score
  gate narrowly fails. A post hoc two-sided calibration readout also finds
  unusually low joint quadratic scores for both cosmologies, driven chiefly by
  Pantheon+. The phenomenology implies `c(z)/c0=6.203198` at `z=2.262260`.
  Because `gamma_c=2.3` was selected after viewing the same outcomes and has no
  independent covariant derivation, the claim decision remains **do not promote**.
- `2026-07-19/rawmu_plamb_gamma2p3_hierarchy/` applies the same outcome-aware
  fixed-gamma cell to the release-grounded calibration hierarchy. The primary
  released-total result remains `Delta BIC=+0.011958`; grouped calibration
  sensitivities span `-0.816882` to `+0.565488`, and every model-difference
  survey mode lies within its external budget. This establishes calibration
  stability of the near-tie but supplies no independent gamma derivation.
- `2026-07-19/plamb_gamma2p3_desi_dr2_bao/` derives and numerically validates
  stationary-path and standard-distance-duality BAO adapters, profiles only the
  common scale `q=c0/(H0 r_d)`, and tests the official 13-element DESI DR2 BAO
  vector with its covariance. No branch simultaneously preserves the alpha-zero
  supernova distance, standard distance duality, the direct clock radial mapping
  and an independent sound horizon, so the physical adapter gate fails and no
  overnight sampler is launched. The fixed supernova-posterior stationary branch
  gives `chi-squared=33458.052`; the duality control gives `chi-squared=176.070`
  and `Delta BIC=+163.234` versus flat Lambda-CDM. The present adapters are
  rejected, without excluding an as-yet unspecified covariant PLAMB completion.
- `2026-07-19/plamb_covariant_distance_closure/` derives a unified disformal
  clock-distance-ruler relation for `T`, `D_L`, `D_M`, `D_H` and an evolving BAO
  ruler, then validates its numerical identities. The kinematic closure passes,
  but physical prediction fails: the static FR branch requires large unprovided
  atomic-clock and luminosity-transfer laws, while even an arbitrary isotropic
  ruler in every DESI bin cannot repair the radial/transverse ratio. The best
  power-law-ruler closure has `chi-squared=2684.670025`, versus `10.271041` for
  flat Lambda-CDM. No further distance sampling is authorised until an
  action-level model independently fixes the clock, photon and ruler evolution.
- `2026-07-19/plamb_missing_degree_suite/` tests nine registered one-degree
  observable extensions of the fixed-`gamma_c=2.3` closure against the 3,422-point
  raw-MU supernova likelihood and the 13-element DESI DR2 BAO likelihood. No
  conformal, clock, curvature, focusing, scalar-ruler or anisotropic rank-one
  proxy passes the cross-probe, DESI-goodness, identifiability and symmetry gates.
  A clearly labelled post hoc rank-two control does fit DESI at the SN-trained
  `p=0.797429`, with `chi-squared=14.463048`, using common ruler power
  `b=1.016667` plus radial slip `h=-0.762001`. It requires transverse and radial
  ruler factors `3.397440` and `8.496792` by `z=2.33`, breaks the covariant
  radial/transverse closure, and remains `Delta BIC=+14.562195` behind flat
  Lambda-CDM jointly. The suite therefore localises a two-mode observable
  deficiency but identifies no new physical degree of freedom and authorises no
  further distance sampling.
- `2026-07-19/plamb_radial_transverse_search/` performs the targeted post hoc
  follow-up with `D_M=chi/(1+z)^b_perpendicular` and
  `D_H=k/(1+z)^b_parallel`, propagating the exact released-total SN profile and
  full DESI covariance. At the SN-trained `p=0.797429`, the full DESI fit is
  acceptable (`chi-squared=14.463048`, 10 degrees of freedom, `p=0.152899`) at
  `b_perpendicular=1.016667` and `b_parallel=1.778668`; the isotropic restriction
  is worse by `Delta BIC=4139.079103`. The powers are stable across the SN
  profile interval, but a genuine leave-`z=2.33`-out prediction fails at
  `p=0.001712`, led by a covariance-decorrelated transverse component of
  `-3.115694` standard deviations. The fitted `D_H/D_M` curve lies within
  `1.422703%` of the joint flat-Lambda-CDM curve across paired bins and remains
  `Delta BIC=+14.558420` behind it jointly. The adapter is therefore best read as
  a post hoc reparameterisation of conventional distance geometry with a
  high-redshift break; it is not validated and identifies no symmetry or field.
- `2026-07-19/plamb_lcdm_lya_validation/` freezes the raw-MU-trained `p` and
  fits the two PLAMB powers and BAO scale using only DESI rows below `z=2`, then
  audits Lambda-CDM equivalence and high-redshift transfer without retuning.
  Constant powers approximate the lower-redshift-trained Lambda-CDM geometry
  with `chi-squared=4.534951` for 10 residual degrees of freedom
  (`p=0.920006`) and a maximum `1.485301%` observable difference, but miss the
  strict equivalence gate because the equivalent powers are `2.416077` fitted
  standard deviations from the observed-data powers. The official compressed
  DESI DR2 hold-out still fails at `p=0.001705`; the auto split passes at
  `p=0.065260`, while the Ly-alpha x quasar split fails at `p=0.002975`.
  The same-spectra alternative multipole analysis also fails at `p=0.005170`,
  establishing method robustness rather than independent replication. DESI
  DR1 passes at `p=0.092834`, and the native eBOSS DR16 auto x cross grid gives
  HPD tail mass `0.079816`. The DR2 native 9,306-element likelihood is not yet
  public, remains an explicit pending gate, and no physical radial/transverse
  or symmetry claim is authorised.
- `2026-07-19/peter_sparse_coupling_search/` prospectively freezes and
  enumerates sparse power-law repairs to the explicit Peter clock, atomic,
  photon-transfer and Type Ia source closure. The stated phenomenology and the
  illustrative microscopic atomic/Chandrasekhar system each require at least
  two independent active couplings; adding the simplified opacity/diffusion
  time-dilation equation raises the minimum to three. The best-ranked
  microscopic two-coupling branch requires an active endpoint change factor of
  `8.278818`, while the best-ranked diffusion branch changes a dimensionless
  gravitational combination by `7.141842` at `z=2.262260`. All algebraic
  candidates remain behind flat Lambda-CDM by the inherited best-path
  `Delta BIC=+83.760921`, and no action-level model is supplied. The result is
  therefore a localisation of missing structure, not a viable revised redshift
  law, and authorises no additional supernova sampling.
