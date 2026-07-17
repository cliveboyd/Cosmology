# FR/LINX Lithium Analysis: Detailed Programme and Assumption Audit

Date: 17 July 2026
Status: final independent audit of the completed 30,503-point key-network catalogue
Language: UK English

## Executive Assessment

The implementation is internally coherent as a **selected-network, penalised sensitivity scan**. All 30,503 registered points completed successfully with zero unresolved evaluations. The raw append-only CSV has 30,700 records because 197 earlier failed attempts were retained before successful retry; reducing it to the latest successful status per point produces exactly the registered catalogue. Catalogue reconstruction finds no semantic point-ID collisions and no random rounding collisions. Recalculation of the stored score components agrees with the programme to approximately `2.3e-12` in absolute chi-squared.

The principal limitation is not a coding failure. It is a mismatch between the physical model named in the research question and the model actually solved. The present programme uses the standard LINX hot-Big-Bang thermal and expansion background, with an effective constant radiation perturbation represented by `Delta Neff`. Your FR framework instead describes a non-expanding universe with evolving speed of light and mass/energy scales. The scan can therefore test an **FR-inspired early-clock proxy**, but it is not a calculation of FR nucleosynthesis and cannot validate or exclude the full FR theory.

The strongest current scientific statement remains negative and useful:

> Within the tested standard-BBN background, selected key network and effective early-clock mapping, none of the 30,503 scanned points reconciles the adopted deuterium, helium and direct stellar-lithium anchors inside the registered two-standard-deviation gates. Treating stellar lithium as a lower bound removes the lithium penalty and favours an astrophysical depletion interpretation, but does not supply positive evidence for FR.

The final D+He-constrained lithium floor is $10^{10}({}^{7}{\rm Li/H})=5.555$ for the nominal expansion-only selection, $5.113$ with modest selected-rate controls and $4.741$ over all scanned controls. These correspond to residual depletion factors of approximately $3.83$, $3.53$ and $3.27$ relative to the adopted plateau centre. The nominal expansion-only selection is itself proposal-contaminated, as documented below, so its row count must not be interpreted as a sample from a continuous prior.

## Audited Components

- Primary scan: [`analyze_bbn_lithium_linx_fr_network_2026-07-16.py`](../../../code/bbn_lithium/analyze_bbn_lithium_linx_fr_network_2026-07-16.py)
- Chunked resume wrapper: [`resume_bbn_lithium_linx_key_chunks_2026-07-17.ps1`](../../../code/bbn_lithium/resume_bbn_lithium_linx_key_chunks_2026-07-17.ps1)
- Registered gate: [`analyze_su2_bbn_lithium_gate_2026-07-17.py`](../../../code/bbn_lithium/analyze_su2_bbn_lithium_gate_2026-07-17.py)
- Exploratory frontier: [`analyze_su2_bbn_lithium_frontier_2026-07-17.py`](../../../code/bbn_lithium/analyze_su2_bbn_lithium_frontier_2026-07-17.py)
- Independent empirical audit: [`audit_bbn_lithium_linx_fr_scan_2026-07-17.py`](../../../code/bbn_lithium/audit_bbn_lithium_linx_fr_scan_2026-07-17.py)
- Final registered gate: [`su2_bbn_lithium_gate_report_2026-07-17.md`](su2_bbn_lithium_gate_report_2026-07-17.md)
- Final exploratory frontier: [`su2_bbn_lithium_frontier_report_2026-07-17.md`](su2_bbn_lithium_frontier_report_2026-07-17.md)
- Updated-anchor reweight: [`fr_linx_anchor_sensitivity_report_2026-07-17.md`](fr_linx_anchor_sensitivity_report_2026-07-17.md)
- LINX source revision: `ec2e9d2ca455e8204137e884da29f5dd13a638fa`
- LINX network: `key_PRIMAT_2023`, 12 reactions
- Registered numerical settings: `rtol=1e-5`, `atol=1e-9`, `sampling_nTOp=150`
- Numerical runtime: Python `3.11.15`, NumPy `2.4.6`, JAX/JAXLIB `0.4.28`, Diffrax `0.6.0`, Equinox `0.11.10`, Lineax `0.0.7`, Optimistix `0.0.10`, Interpax `0.3.1`

## Model Actually Evaluated

The scan defines an effective clock factor and maps it to an initial excess-radiation parameter:

\[
\begin{aligned}
S                         &\equiv H_{\rm BBN}/H_{\rm SBBN},\\
\Delta N_{\rm eff}^{\rm init}
                          &= \frac{43}{7}\left(S^2-1\right),\\
\Omega_bh^2               &= 0.02242\,f_\eta,\\
\eta_{10}                 &\simeq 6.12\,f_\eta,\\
\tau_n                    &= 879.4\,f_\tau\ {\rm s}.
\end{aligned}
\]

For each available reaction $i$, LINX uses a log-normal rate nuisance parameter:

\[
\begin{aligned}
q_i                       &\sim \mathcal N(0,1),\\
R_i(T,q_i)                &= R_{i,0}(T)\exp\!\left[q_i\ln \sigma_{i,\exp}(T)\right].
\end{aligned}
\]

The reported direct-lithium objective is

\[
\begin{aligned}
\chi^2_{\rm scan}
 &= \chi^2_{\rm D}+\chi^2_{Y_p}+\chi^2_{\rm Li}
  +\chi^2_{\eta,\rm CMB}+\chi^2_{\tau_n}+\sum_i q_i^2,\\
\chi^2_X
 &= \left(\frac{X_{\rm pred}-X_{\rm obs}}{\sigma_X}\right)^2.
\end{aligned}
\]

This is a penalised profile objective evaluated over a designed catalogue. It is not a posterior sample, marginal likelihood, Bayesian evidence or chi-squared goodness-of-fit statistic with a single well-defined number of degrees of freedom.

## Findings by Priority

### P0. The Current Background Is Not an FR Nucleosynthesis Model

The programme calls standard LINX `BackgroundModel(delta_neff)` and retains the standard radiation-era relation among temperature, scale factor, density and time. The added component redshifts as radiation. This tests a conventional constant dark-radiation/clock perturbation.

For a genuinely non-expanding FR cosmology, reaction freeze-out cannot be defined by simply setting or perturbing standard (H(T)). FR must provide its own thermal clock and conservation laws. At minimum, a calculable model must specify

\[
\begin{aligned}
T(t)                      &= T_{\rm FR}(t),\\
\rho_i(T)                 &= \rho_{i,\rm FR}(T),\\
c(T)                      &= c_0 F_c(T),\\
m_j(T)                    &= m_{j,0}F_{m_j}(T),\\
\eta(T)                   &= n_b(T)/n_\gamma(T),\\
Q_{np}(T)                 &= m_n(T)-m_p(T),\\
B_D(T)                    &= \text{deuteron binding energy},\\
\Gamma_{n\leftrightarrow p}(T)
                          &= \Gamma_{\rm FR}(T),\\
R_k(T)                    &= R_{k,\rm FR}(T).
\end{aligned}
\]

Without these relations, the result must be labelled **FR-inspired proxy**, not **FR lithium test**. This is the single most important claim-boundary correction.

### P0. The CMB Baryon Prior Is Not Matched to the Scanned Cosmology

The code applies a one-dimensional Gaussian prior $f_\eta=1\pm0.04/6.12$, derived from a Planck-like standard-cosmology result, while independently varying `Delta Neff`. Baryon density, $N_{\rm eff}$, helium, acoustic scale and other CMB parameters are correlated. Under a non-standard FR background, the standard Planck posterior may not be applicable at all.

The immediate conventional correction is a joint $p(\Omega_bh^2,N_{\rm eff},Y_p\mid{\rm CMB})$ prior or a matched CMB likelihood. The eventual FR correction is a forward CMB calculation under the same FR thermal and perturbation model used for BBN. Until then, the existing eta term is a **standard-CMB compatibility stress**, not an FR-native prior.

### P1. The Key Network Omits Most Requested Lithium/Be Destruction Channels

Eight target reactions were requested, but only three are present in `key_PRIMAT_2023`:

- Available: `He3aBe7g`, `Be7nLi7p`, `Li7paa`.
- Unavailable: `Li7paag`, `Be7naa`, `Be7daap`, `Be7pB8g`, `Li7daan`.

Forty individual pull definitions for the five unavailable reactions collapse to one fixed-rate row after network normalisation; 39 exact duplicates are removed. This behaviour is deterministic and the final target count of 30,503 is correct, but the run does not test the omitted channels.

Direct experiments already indicate that enhanced ${}^{7}{\rm Be}(n,p)$, ${}^{7}{\rm Be}(n,\alpha)$ and $d+{}^{7}{\rm Be}$ reactions give only modest lithium reductions. A **targeted full-network rerun** remains necessary for completeness, especially at the baseline, best penalised rows and D+He-constrained lithium frontier. It need not repeat all 30,503 points.

### P1. Deuterium Theory Uncertainty Is Not Marginalised

The key network contains 12 nuclear reactions, but the random and ladder searches vary only the three available targeted lithium/beryllium pulls. The remaining nine pulls are fixed at zero. These include the deuterium-sensitive `dpHe3g`, `ddHe3n` and `ddtp` rates.

The objective nevertheless scores D/H using only the adopted observational uncertainty, $\sigma_{10^5({\rm D/H})}=0.029$. It therefore treats the LINX prediction as exact conditional on the three lithium pulls and can make both the D veto and total $\chi^2$ too restrictive. This matters because deuterium is the principal constraint preventing more extreme clock or baryon shifts.

The next likelihood must either marginalise every available key-network rate pull with its calibrated LINX prior or include a validated theory covariance for D, He and Li. A targeted check should first quantify the change in the best score and gate membership when the three main D-rate pulls are released.

### P1. Proposal Sampling and Scoring Use Different Nuclear-Pull Distributions

Random pulls are proposed as

\[
\begin{aligned}
q_i^{\rm proposal}       &\sim \mathcal N(0,1.35^2),\\
q_i^{\rm stored}         &=0\quad\text{when }|q_i|\leq0.25,\\
-4                       &\leq q_i^{\rm stored}\leq4,
\end{aligned}
\]

but the objective penalises them as unit-Gaussian pulls, $\sum_i q_i^2$. Oversampling the tails can be sensible for discovering a boundary, but it must be described as an exploration proposal. It is not sampling from the stated prior unless importance weights or a separate target distribution are used.

The hard `|q|<=0.25` zeroing also creates an artificial point mass at $q=0$. The final audit found 91 random rows with all three available Li/Be pulls exactly zero, including five inside the registered SU2 base cuts. Consequently, the "expansion-only" subset mixes deliberate fixed-rate grid rows with chance proposal rows. The gate remains a valid searched-support diagnostic, but its row count and implied sampling measure are not those of a continuous prior.

Recommended correction: keep a separate deterministic `q=0` family; draw all stochastic $q_i\sim\mathcal N(0,1)$ without thresholding; and use either posterior sampling, weighted importance sampling or an explicit optimiser for the intended inferential task.

### P1. The Scan Finds Minima but Does Not Establish Posterior Support

The 30,000 random rows plus fixed ladders improve exploration, but they do not provide:

- convergence diagnostics;
- posterior credible intervals;
- marginal likelihood or Bayes factors;
- calibrated profile-likelihood intervals;
- proof that a reported boundary minimum is global;
- a prior-volume penalty beyond the local quadratic pull terms.

The best direct-lithium row has approximately

\[
\begin{aligned}
\chi^2_{\rm total}       &=215.26,\\
10^5({\rm D/H})          &=2.5851,\\
Y_p                       &=0.24869,\\
10^{10}({}^{7}{\rm Li/H})&=4.6941,\\
\text{depletion factor}  &=3.237,\\
\sum_i q_i^2             &=22.85.
\end{aligned}
\]

Its lowest production pull, `He3aBe7g=-3.98`, is effectively on the imposed `-4` boundary. This is evidence that even an extreme selected-rate adjustment does not approach the direct lithium measurement, not evidence for a fitted FR parameter value.

### P1. The Observational Likelihood Needs an Explicit Sensitivity Ladder

The adopted anchors are transparent but not unique:

\[
\begin{aligned}
10^5({\rm D/H})           &=2.508\pm0.029,\\
Y_p                       &=0.245\pm0.003,\\
10^{10}({}^{7}{\rm Li/H})&=1.45\pm0.25.
\end{aligned}
\]

Relevant newer work includes a 2024 combined deuterium value $2.533\pm0.024$, and the 2026 LBT helium result $Y_p=0.2458\pm0.0013$. The latter roughly halves the present helium uncertainty. A May 2026 NLTE study also reports evolutionary-state and metallicity dependence in lithium, reinforcing that one Gaussian plateau anchor is not a complete observation model.

A likelihood-only sensitivity reweight using those D and He anchors together with the 2025 PDG neutron lifetime leaves the conclusion unchanged. The updated expansion-only selection has no D+He two-standard-deviation pass. The modest, complete SU2-compatible and all-scanned selections have 99, 124 and 2,872 D+He passes respectively, but zero D+He+Li passes. Their D+He-constrained lithium floors are respectively $5.122$, $5.009$ and $4.741$ in units of $10^{10}({}^{7}{\rm Li/H})$. Because this reweight does not release the fixed D-sensitive rates or rerun LINX, it is a sensitivity direction rather than a new default likelihood.

The direct-lithium and lower-bound cases are valuable bracketing controls, but a publishable analysis should add a predeclared ladder:

1. PDG/conservative anchors, as now.
2. Current D/H and LBT $Y_p$ anchors.
3. Stellar plateau with an explicit depletion distribution.
4. Lower-red-giant lithium as a separate likelihood.
5. Interstellar SMC lithium as an independent, chemically evolved constraint.

The likelihood should also distinguish observational errors from LINX nuclear-theory uncertainty instead of silently treating the profiled pull penalty as a complete theory covariance.

The neutron-lifetime anchor also requires an update. The bundled LINX constant and scan prior use $879.4\pm0.6$ s, inherited from an older PDG input. The 2025 PDG review gives $878.4\pm0.5$ s. The newer centre corresponds to `tau_n_fac=0.998863` in the existing LINX normalisation and lies inside the scanned range, so catalogue reweighting is possible; however, abundance predictions at the new centre and the gate's prior cuts should be rerun rather than treating the old prior as current.

### P1. Nuclear Pulls Are Nuisance Controls, Not FR Predictions

LINX correctly interprets each $q_i$ as a log-normal rate uncertainty. The present scan assumes independent unit Gaussians and varies only selected rates. A changing $c$, mass scale, fine-structure constant or quark-mass scale would instead move several quantities coherently: neutron-proton mass splitting, weak rates, Coulomb barriers, reaction Q-values, binding energies and reverse rates.

Research on varying constants demonstrates why isolated reaction multipliers are inadequate. Quark-mass variation is constrained at roughly the percent level by BBN, while unified varying-constant analyses propagate coupled changes through nuclear and cosmological parameters. These form a useful technical template for an FR-native response model, even if their particle-physics assumptions differ from FR.

### P1. Numerical and Network Convergence Has Not Yet Been Demonstrated

The chosen tolerances are reasonable for a broad scan, and process recycling successfully removed prior JAX memory failures. Nevertheless, the interpretation-bearing rows should be rerun with:

- `rtol=1e-7` and `atol=1e-11`;
- increased weak-rate interpolation sampling, for example 300 and 600;
- `key_PRIMAT_2023`, `small_PRIMAT_2023` and `full_PRIMAT_2023`;
- direct abundance differences and score differences tabulated;
- finite-value and species-index checks.

This should cover baseline, every registered best row, all gate-nearest rows and a small stratified sample. It is a much higher-value calculation than another untargeted 30,000-point key-network scan.

### P2. Statistical Gates Are Transparent but Not Calibrated Joint Tests

The registered gate requires each of D and He, then lithium, to lie within two marginal standard deviations. This rectangular rule is easy to audit, but it does not account for covariance, joint coverage, selection over 30,503 points or posterior volume. Similarly, $\sum_i q_i^2\leq9$ allows several very different physical pull combinations and is not tied to an external reaction-calibration budget.

Retain the gate as a veto summary, but supplement it with:

- a joint multivariate abundance likelihood;
- posterior-predictive checks;
- a profile-likelihood ratio with simulated calibration where regularity is doubtful;
- prior-predictive probability for the required nuclear pull vector;
- a reaction-by-reaction external uncertainty budget.

### P2. The Frontier and Influence Diagnostics Are Descriptive

The current Pareto frontier optimises D+He disagreement against lithium abundance, using the nuclear penalty only as a tie-break. It is not a three-objective frontier. The linear influence fit standardises a broad random proposal and is therefore proposal-dependent; its coefficients are not local physical derivatives or global Sobol indices.

Recommended additions are local JAX Jacobians at the baseline and selected frontier points, a prior-weighted Sobol analysis, and a three-objective frontier in D+He error, lithium error and nuclear pull cost.

### P2. Provenance Can Be Strengthened

The running config records paths, command arguments and Python version. The original manifest records file size and modification time but not hashes. The independent audit has added SHA-256 hashes for the programme, catalogue snapshot and config, plus the LINX Git revision.

Future runs should record dependency lock files, NumPy/JAX/Diffrax/LINX versions, Git revisions, input hashes, random seed, hardware backend and final output hashes automatically.

## What the Current Result Does Establish

Subject to the targeted numerical and full-network checks still listed below:

1. The key-network score arithmetic is correct.
2. The registered 30,503-point catalogue construction is deterministic and has no detected semantic ID collision.
3. Standard-BBN expansion/clock changes cannot independently suppress lithium while preserving D/H and helium in the tested range.
4. Even extreme selected-rate pulls leave the best lithium prediction far above the adopted stellar plateau.
5. Treating stellar lithium as a lower bound makes a conventional BBN-like point acceptable, pointing towards stellar/observational processing rather than an early-clock solution.
6. None of these statements constitutes a test of a complete non-expanding FR nucleosynthesis model.

## Recommended Next Actions

### Immediate: Validation Before More Breadth

1. **Completed:** finish the registered catalogue and regenerate the gate/frontier once.
2. **Completed:** rerun this independent audit on the immutable final CSV.
3. Run a 20-100 point numerical/network validation matrix at tighter tolerances across key, small and full networks.
4. Regenerate the expansion-only gate from a deterministic `q=0` family, excluding proposal-generated zeros.
5. **Reweight completed:** repeat the numerical validation and formal gate under the current D/H and 2026 LBT helium anchors.
6. **Reweight completed:** rerun the selected numerical validation points at the 2025 PDG neutron-lifetime centre.
7. Release the non-targeted key-network pulls, especially `dpHe3g`, `ddHe3n` and `ddtp`, and quantify the deuterium theory covariance before interpreting its veto as exact.

### Next Statistical Analysis

Construct a hierarchical posterior with parameters

\[
\begin{aligned}
\theta
 &= \{\Omega_bh^2,\Delta N_{\rm eff},\tau_n,\boldsymbol q,
       \delta_{\rm Li,depl},\boldsymbol\psi_{\rm obs}\},\\
p(\theta\mid\mathcal D)
 &\propto p(\mathcal D_{\rm D},\mathcal D_{Y_p},\mathcal D_{\rm Li}\mid\theta)
           p(\Omega_bh^2,\Delta N_{\rm eff},Y_p\mid\mathcal D_{\rm CMB})
           p(\tau_n)p(\boldsymbol q)p(\delta_{\rm Li,depl})p(\boldsymbol\psi_{\rm obs}).
\end{aligned}
\]

LINX is differentiable and was designed to support gradient-based inference, so NUTS or another validated sampler is preferable to interpreting random-catalogue minima.

### Highest-Value FR Investigation

Before a large FR-labelled BBN run, derive and preregister the FR thermal model. The first falsification gate should be whether that model simultaneously produces:

- a thermal photon distribution and an explicit $T(t)$;
- weak freeze-out with a viable neutron fraction;
- a finite deuterium bottleneck epoch;
- D/H and $Y_p$ near observation without post hoc reaction pulls;
- a baryon-to-photon mapping compatible with the FR CMB calculation;
- lithium as a prediction rather than an independently suppressed output.

If FR does not predict these ingredients yet, the honest and scientifically strong position is that BBN is an open model-construction requirement, not a solved anomaly test.

## Alternate Investigation Paths

### A. Hierarchical Stellar-Depletion Model

Model individual stellar lithium abundances using effective temperature, metallicity, evolutionary state, NLTE correction and population-level depletion. This directly tests the path favoured by the present lower-bound result and can incorporate plateau, lower-red-giant and globular-cluster samples.

### B. Model-Independent BBN Response Basis

Use LINX/JAX derivatives to build a response matrix of D, He and Li to $H(T)$, $\eta$, $\tau_n$ and every nuclear pull. Singular-value decomposition can identify which combinations are constrained by D+He and whether any lithium-reducing direction remains. An FR model can later be projected onto this basis without rerunning a blind high-dimensional search.

### C. Full-Network Targeted Nuclear Audit

Use measured external rate uncertainties and full-network reactions, especially ${}^{7}{\rm Be}$ production/destruction. This tests completeness and should reproduce the established conclusion that plausible nuclear changes alone do not remove the discrepancy.

### D. Coupled Fundamental-Constant Response

Translate a specified FR $c(T)$ and mass law into dimensionless variations such as $\alpha(T)$, $m_q/\Lambda_{\rm QCD}$, $Q_{np}$ and $B_D$, then propagate them coherently through weak and nuclear physics. Existing varying-constant BBN calculations provide a methodological starting point, not an FR result.

### E. Independent-Code Replication

Repeat a small benchmark set in PRIMAT or PArthENoPE/AlterBBN. Agreement within declared numerical and nuclear-input differences would substantially strengthen the research standing.

### F. Late-Time Non-Thermal Mechanisms

Only pursue neutron injection, decays or photon cooling if FR supplies a specific mechanism and prior. Such mechanisms are strongly constrained because lithium reduction often overproduces deuterium. An unconstrained scan would add flexibility without increasing explanatory standing.

## Relevant Research Material

- [LINX: fast, differentiable BBN and gradient-based inference](https://arxiv.org/abs/2408.14538). Directly supports moving from a random scan to a joint posterior.
- [PRIMAT precision BBN](https://arxiv.org/abs/1801.08023). Independent benchmark for thermal history, weak rates and full nuclear networks.
- [Planck 2018 cosmological parameters](https://arxiv.org/abs/1807.06209). Shows that the quoted baryon constraint is conditional on a cosmological model and that $N_{\rm eff}$ is jointly constrained.
- [2024 primordial deuterium determination](https://arxiv.org/abs/2401.12797). Reports combined $10^5({\rm D/H})=2.533\pm0.024$ and discusses a moderate standard-model tension.
- [2026 LBT primordial helium](https://arxiv.org/abs/2601.22238). Reports $Y_p=0.2458\pm0.0013$, a high-priority sensitivity update.
- [2026 NLTE metal-poor-star lithium study](https://arxiv.org/abs/2605.19334). Supports evolutionary and metallicity-dependent observation modelling.
- [Interstellar lithium in the SMC](https://arxiv.org/abs/1207.3081). Provides a non-stellar-atmosphere constraint with chemical-evolution caveats.
- [Lower-red-giant lithium diagnostic](https://arxiv.org/abs/1205.2389). Offers an observation channel less sensitive to atomic diffusion than dwarf-star surfaces.
- [Quark-mass variation constraints from BBN](https://arxiv.org/abs/1012.3840). Demonstrates coherent propagation from fundamental parameters to binding energies and abundances.
- [Varying-constant BBN with cosmological degeneracies](https://arxiv.org/abs/2012.10505). Relevant template for joint constant/cosmology priors, while remaining model-dependent.
- [UCN-tau neutron lifetime](https://arxiv.org/abs/2106.10375). Reports $877.75$ s with sub-second uncertainty and motivates a beam/bottle sensitivity case rather than a single uncomplicated Gaussian.
- [2025 PDG neutron properties](https://pdg.lbl.gov/2025/AtomicNuclearProperties/neutron.html). Gives the current scaled world average $878.4\pm0.5$ s used for the recommended default-prior update.
- [LUNA deuterium-burning rate](https://www.nature.com/articles/s41586-020-2878-4). Important for the D/H theory budget.
- [n_TOF ${}^{7}{\rm Be}(n,p){}^{7}{\rm Li}$](https://arxiv.org/abs/1803.05701) and [n_TOF ${}^{7}{\rm Be}(n,\alpha)$](https://arxiv.org/abs/1606.09420). Both find only minor lithium improvement in the BBN energy range.
- [Measured $d+{}^{7}{\rm Be}$ cross-sections](https://pubmed.ncbi.nlm.nih.gov/31144906/). The measured resonance reduces lithium but does not solve the problem.

## Recommended Claim Boundary

> The completed calculation is a reproducible selected-network LINX sensitivity scan around the standard radiation-dominated BBN background. It finds no lithium resolution from the tested effective clock, baryon-density, neutron-lifetime and selected nuclear-rate controls while preserving deuterium and helium. Because the FR framework requires a distinct non-expanding thermal history and coherent evolution of dimensionless particle and nuclear quantities, this scan is not yet a full FR nucleosynthesis calculation. FR remains subject to a separate model-construction and validation gate.
