# HHT_BAO_Locking_Report_1V8

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\HHT_BAO_Locking_Report_1V8.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2025-10-01T19:03:00Z`
- DOCX modified: `2025-10-07T15:12:00Z`

## Extracted Text

Cosmology HHT–BAO Locking Report

(Revision 1V8)

Executive Summary

This document provides a report detailing the novel application of Hilbert–Huang Transform (HHT) analysis applied to supernova (SN) residuals in comparison to Baryon Acoustic Oscillation (BAO) scales.

In the context of this analysis, SN residuals are defined as the difference between the observed and model-predicted distance modulus:

That is, the difference between the measured supernova distance modulus and the value predicted by the chosen cosmological model.

These residuals are the input to the HHT analysis and are tested for statistical alignment with BAO effective scales.

The analysis includes robust denoising under the premiss Correlation is NOT Causation. Where denoising refers to re-running HHT with polynomial detrending toggled on/off to ensure IMF extraction is not dominated by baseline removal.

Integrates +10,000-permutation surrogate tests, IMF2/IMF3 vs strength analyses, and baseline heatmap visualisations.

All tests were conducted using python under Anaconda spyder.

The findings strengthen the evidence for IMF2 locking with BAO scales under certain model strengths, while IMF3 shows inconsistent behaviour.

IMF3 only shows significance at isolated strengths (esp. 6), and fails robustness tests like seed or permutation stability, so it’s not considered reliable.

Surrogate testing confirms robustness against chance alignment.

The net result is that IMF2 appears to encode a statistically robust memory of BAO physics in the SN residuals, a feature not predicted by ΛCDM.

Background - Standard Cosmology vs. HHT Analysis

1. What is currently understood (Standard ΛCDM cosmology)

The ΛCDM model is the widely accepted description of the Universe.

It assumes:

A constant dark energy term (Λ) driving cosmic acceleration.

Cold Dark Matter (CDM) providing most of the mass.

Normal matter (atoms, gas, stars, galaxies) is only ~5% of the total energy density.

Observational pillars:

Supernovae (SNe Ia) show accelerating expansion.

Baryon Acoustic Oscillations (BAO) act as a standard ruler, matching predictions of ΛCDM.

Cosmic Microwave Background (CMB) confirms early Universe consistency with ΛCDM.

In this framework, deviations from smooth expansion are treated as noise or local systematics. The model does not expect small oscillatory features in expansion history.

2. What does “BAO locking” actually mean?

Definition within this document:Locking = statistical alignment of IMF instantaneous frequencies (from SN residuals) with effective BAO wavenumbers (1/k).

In simpler terms:

BAO provides a fixed cosmic “standard ruler” (set by sound waves in the early Universe).

SN residuals, when decomposed via HHT, show oscillatory modes (IMFs).

If one of these IMFs consistently oscillates “in step” with BAO scales, one may say it’s “locked” — i.e., the oscillation frequency extracted from SN data is phase-aligned or scale-correlated with BAO’s characteristic scale.

Importantly: this is not resonance in a physical sense, but rather a statistical correlation that implies the two structures may share the same underlying physics.

3. What this project address: (Hilbert–Huang Transform analysis)

Instead of fitting a single smooth curve, the analysis applied the Hilbert–Huang Transform (HHT) to supernova residuals.

This decomposes the data into Intrinsic Mode Functions (IMFs), which behave like oscillatory “ripples” embedded in the expansion history.

These IMFs were tested against BAO scales — the “acoustic ruler” left over from sound waves in the early Universe.

To check if any match was real or just noise, a surrogate testing (randomised versions of the data, up to 10,000 permutations) were injected. This provides a robust statistical baseline to remove noise artefacts.

4. What has been found (Novel signals)

IMF2 shows robust positive correlations at multiple strengths (0–4, 8, 12), but weak or absent at noether strengths 6 and 10:

, with

IMF2… with p_perm ≤ 0.001 in many tests (statistically significant).

Suggests that part of the residual structure in SN data may be phase-locked to the same scale that defines BAO.

IMF3 shows weaker and more variable correlations: sometimes significant, sometimes consistent with noise. However, at Noether strength=6, IMF3 shows strong BAO correlation.

These signals survive 10,000-permutation surrogate tests, reducing the chance of spurious detection.

The strength of the IMF–BAO correlation appears to depend on Noether strength (a parameter in this model controlling conservation-law constraints), see appendix for further details, hinting at a possible link between cosmological symmetries and residual oscillatory structure.

5. Novelty and Relevance

Standard ΛCDM assumes no preferred oscillatory features in the residuals once the smooth expansion is removed.

The detection of a statistically robust IMF–BAO correlation implies that:

The Universe’s expansion history may contain subtle oscillatory components linked to the same physics as BAO.

These features might indicate dynamic dark energy or effects of embedding cosmology in a parent structure (such as conjectured Full_Relativity (FR) or Parent_Black_Hole (PBH) frameworks).

If confirmed, this would represent a departure from ΛCDM’s assumption of pure smoothness, opening a new window into how large-scale cosmic structure interacts with fundamental physics.

6. Going forward… Why it may matter

Provides a new diagnostic tool: HHT allows discovery of hidden structure invisible to Fourier/wavelet methods.

May help reconcile tensions in cosmology, such as:

The Hubble constant () tension between local and early-Universe measures.

Small but persistent anomalies in BAO vs SN fits.

Points toward a future where time-dependent, symmetry-driven corrections to expansion history are tested alongside ΛCDM.

If the effect turns out to be subtle or model-dependent, it helps establish robust methods for probing non-ΛCDM behaviour.

Note: HHT IMF analysis typically allows analysis of non-stationary datasets (behaviour)

Why is this interesting?

Standard ΛCDM expects smoothness:SN residuals should look like noise, not structured oscillations.Locking would challenge the assumption that residuals are only statistical scatter.

A novel “hidden diagnostic”:HHT detects non-stationary oscillations that Fourier or wavelet methods miss.This may reveal subtle cosmological fingerprints invisible to traditional pipelines.

Physics implications:

Could point toward dynamic dark energy (oscillatory corrections to expansion).

Could suggest symmetry-breaking/restoration in spacetime (hence the link to Noether strength).

In a parent black hole (PBH) or FR scenario, BAO locking may reflect coupling between local oscillatory modes and global “embedding” geometry.

Relevance to tensions:

May shed light on the Hubble tension (different expansion rates measured locally vs CMB).

May explain persistent small mismatches between BAO and SN fits under ΛCDM.

How it can be interpreted

If validated, BAO locking would mean that oscillations in the late-time expansion history carry a memory of early-Universe sound waves.

Interpretation paths:

Conservative: It’s a subtle systematic in SN light-curve calibration that just happens to mimic BAO scales.

Bold: It’s real physics — a signature of dark energy dynamics or spacetime structure.

Middle ground: It’s a diagnostic of how conservation laws (Noether strength) imprint scale-dependent structure in the cosmic expansion history.

How HHT could be used.

As a new cosmological probe:

Complementary to SN, BAO, and CMB — but probing the residual structure, not just the mean trend.

Could become a consistency test: if IMF2 locking vanishes in some datasets but not others, it tells one about systematics vs physics.

In model selection:

If FR/PBH naturally predict BAO-phase residuals, locking could be evidence in their favour.

If ΛCDM fails to reproduce them, it motivates extensions.

In data cleaning:

If IMF3-like non-robust “false locks” appear, this provides a tool for filtering noise modes vs true physics modes.

👉 In short: ΛCDM says residuals are noise. HHT shows they may not be.If residual oscillations consistently line up with BAO scales, this may suggest new physics — or at the very least, a refinement of how to model the Universe’s Expansion.

Acronyms and Variables

  SN — Supernova (Type Ia, unless otherwise stated).

  Δμ(z), μ(z), μ_obs, μ_model — Distance modulus residuals.

  DM, DH, DA, DL — Transverse comoving distance, Hubble distance, angular diameter distance, luminosity distance.

  rd, cs(z), Rb(z) — Sound horizon, sound speed, baryon-to-photon density ratio.

  χ(z), Sk(χ) — Comoving distance integral and curvature-dependent function.

  Ωm, ΩΛ, Ωk, Ωr — Matter, dark energy, curvature, radiation density parameters.

  E(z), H(z), H₀ — Normalized Hubble parameter, Hubble parameter, Hubble constant.

  f(z), fσ₈(z) — Growth rate of structure.

  R, ℓa — CMB distance priors.

  k, k_eff, 1/k, λBAO — Effective BAO wavenumbers and scales.

  AIC, BIC, χ², t — Model selection and significance statistics.

  Noether strength — Model control parameter for conservation law enforcement.

  FR, PBH, ΛCDM — Alternative cosmological frameworks tested.

Note: Locking → Refers to statistical alignment of IMF frequencies with BAO scales.                              It does NOT infer a literal physical resonance.

Methodology

Supernova residuals were decomposed via HHT into IMFs.

The HHT method combines Empirical Mode Decomposition (EMD) with Hilbert transforms to extract instantaneous frequencies ω(z) from SN residuals as a function of log(1+z).

These are compared against effective BAO scales (k & 1/k_eff) derived from clustering measurements.

Pearson correlation coefficients (r) are computed for IMF2 and IMF3 to test alignment with DM/rd and DH/rd. Significance is assessed using phase-randomized surrogate datasets.

2. IMF2 and IMF3 instantaneous frequencies were compared to BAO effective wavenumbers (k and 1/k).3. Surrogate testing was applied: >10,000 random phase-shuffled realisations were generated to estimate null distributions, producing robust permutation p-values (p_perm).4. Significance thresholds are set at α = 0.05 with surrogate-corrected p-values.5. Locking strength tested across a grid of Noether-conserving model strengths (soft_0 to soft_12).

Results (10k Randomised Permutations)

Harmonic Ratio Test

Figure: Harmonic ratio (IMF3/IMF2) across packs. Weak indications of structured coupling appear under lambda-on/grid, though large variance is observed.This suggests potential model-driven resonances.

IMF2 Correlation vs Strength (10k Surrogates)

Figure: IMF2 correlation with 1/k across Noether strengths. Significant locking observed at strengths 0–4, 8, and 12 with p_perm ≤ 0.001. Suggests robust IMF2–BAO locking under varied conditions.

IMF3 Correlation vs Strength (10k Surrogates)

Figure: IMF3 correlation with 1/k across Noether strengths. Significance is inconsistent, with IMF3 showing notable locking only at strength=6. Other strengths produce weak or null results.

Baseline Heatmap (10k Surrogates)

Figure: Baseline surrogate heatmap p1k. Note: IMF3 10k baseline (s9_p10k) has IMF3 with no significant threshold.

Appendices

Key CSV results are included for reproducibility.

• IMF2/IMF3 correlations vs strength (10k permutations):    p10k_IMF2_IMF3_vs_strength_s9.csv below

👉 Net takeaway:- IMF2 locking to BAO is robust to seed variation and permutation depth.- IMF3 shows no robust significance, supporting the interpretation that only IMF2    captures genuine structure.- These findings reduce the chance that the IMF2 signal arises from randomisation   artefacts or implementation bugs.

3. Systematic exclusion   - Shuffled-BAO and replicate baselines confirmed that IMF2 correlation persists      while IMF3 remains null.

- Shuffle controls reproduced null IMF3, while IMF2 remained stable —      suggesting the shuffle did not erase IMF2 locking, which is expected if IMF2 is      genuinely tied to BAO scales   - Indicates that spurious programmatic or surrogate-induced biases are unlikely to       explain IMF2’s BAO-locking.

2. Permutation convergence   - Increased surrogate permutations from 2k → 20k.   - IMF2: p_perm converged stably around 0.005–0.006 across all permutation      counts.   - IMF3: p_perm remained ≈0.9995 across runs.   - r-values invariant to within floating-point precision.   → Suggests IMF2 BAO-locking signal is statistically robust, not an artefact of          under-sampling.

1. Seed invariance tests   - Re-ran baseline HHT–BAO correlation with explicit RNG      seeding (seed=101, 202, 303).   - Results for IMF2 and IMF3 were identical across seeds.   - Example: IMF2 r=0.090845, p_perm≈0.0050; IMF3 r≈−0.118, p_perm≈0.9996.   - Spearman tests show no trend of p-values with seed index.       → Confirms that results are not artifacts of RNG seeding.

Robustness & Seed-Stability Checks (Revision 1V8)

Next Steps

• Extend surrogate testing to >50,000 permutations for ultra-robust significance.• Expand redshift-window analyses to isolate low-z vs high-z contributions.• Cross-validate with eBOSS DR16 BAO datasets and JWST high-z supernovae.• Integrate with FR + PBH cosmological model fits for joint constraints.

Revision History

• 1V0 – Initial HHT BAO locking setup• 1V1 – Added IMF2/3 r vs strength plots• 1V2 – Added multi-window redshift tests• 1V3 – Baseline surrogate heatmaps• 1V4 – Harmonic test integration• 1V5 – 10k surrogate findings, expanded robustness checks, embedded plots with             figure captions.

• 1V6 – Null

• 1V7 – Denoising findings, expanded robustness checks, test for program induced              systematic errors

• 1V8 – Added Key Cosmological Equations used within analysis to Appendix

Appendix Summary Results  p10k_IMF2_IMF3_vs_Noether strength_s9

Table: p10k_IMF2_IMF3_vs_strength_s9.csv

Appendix - Key Cosmology Equations

Key HHT Analysis Equations

• SN residuals: Δμ(z) = μ_obs(z) − μ_model(z)

• BAO observables: DM/rd and DH/rd, where rd is the sound horizon at drag epoch.

• Pearson correlation used for locking tests: r = Cov(ω̄, 1/k_eff) / (σ_ω̄ σ_1/k).

1. Distances & Distance Modulus

χ(z) = ∫₀ᶻ [ c(z′) / H(z′) ] dz′

Sₖ(χ) = { sin(√|Ωₖ| χ)/√|Ωₖ|   if Ωₖ < 0;  χ if Ωₖ = 0; sinh(√|Ωₖ| χ)/√|Ωₖ|  if Ωₖ > 0 }

D_M(z)= Sₖ(χ)

D_H(z) = c(z) / H(z)

D_A(z) = D_M(z) / (1 + z)

D_L(z) = (1 + z) · D_M(z)

μ(z) = 5 log₁₀[ D_L(z) / 10 pc ]

Δμ(z) = μ_obs(z) − μ_model(z)

2. BAO Core Definitions

r_d = ∫_{z_d}^∞ [ c_s(z) / H(z) ] dz

c_s(z) = c(z) / √[ 3 ( 1 + R_b(z) ) ]

R_b(z) = 3 ρ_b(z) / (4 ρ_γ(z))

BAO observables: D_M/r_d, D_H/r_d

D_V(z) = [ z · D_H(z) · D_M²(z) ]^{1/3}

F_AP(z) = D_M(z) · H(z) / c(z)

3. Expansion History

E(z) = H(z) / H₀ = √[ Ω_r (1+z)^4 + Ω_m(z)(1+z)^3 + Ωₖ(1+z)^2 + Ω_DE(z) + .. ]

Ω_m(z) = Ω_m0 · ℳ(z)

c(z) = c₀ · ℭ(z)

4. Growth of Structure

f(z) = d ln D₊(z) / d ln a

fσ₈(z) = f(z) · σ₈(z)

f(z) ≈ [Ω_m(z)]^γ  (γ ≈ 0.55 in GR+ΛCDM)

5. CMB Distance Priors

R = √[Ω_m H₀²] · D_M(z_*) / c(z_*)

ℓ_a = π · D_M(z_*) / r_s(z_*)

6. Locking & Statistics

r(ω̄, 1/k_eff) = Cov(ω̄, 1/k_eff) / (σ_ω̄ · σ_{1/k})

t = r √[(N−2)/(1−r²)]

χ² = (d − m)^T C⁻¹ (d − m)

AIC = 2k + χ²_min

BIC = k ln N + χ²_min

7. BAO/k-Space Link

k_BAO(z) ≈ π / r_d

λ_BAO(z) ≈ 2π / k_BAO ≈ 2 r_d

ω̄(log(1+z)) ↔ 1/k_eff

Appendix - Primary Datasets used within scope of investigation.

Supernova (SN) data

Pantheon+SH0ES/Users/boyde/Downloads/Pantheon+SH0ES.dat

Type Ia SN dataset (∼1700 SNe, redshifts z ≈ 0.01–2.3).

Anchored with SH0ES Cepheid calibration for .

Used as the residual input for Hilbert–Huang Transform (HHT) decomposition.

Pantheon+ Covariance matrix/Users/boyde/Downloads/Pantheon+SH0ES_STAT+SYS.cov

Statistical + systematic covariance for the SN dataset.

Used in baseline FR cosmology fits before residual extraction.

BAO data

BAO compilation (long format)/Users/boyde/Downloads/bao_long.csv

BAO measurements from eBOSS DR16 (LRG, QSO, ELG) + earlier surveys.

Columns: z, kind (DM/rd or DH/rd), value, sigma.

Provides and anchors for IMF correlation testing.

Planck priors (optional in FR runs)

Sometimes distance priors or Mpc from Planck 2018:

Used as input when normalising BAO scale.

In scripts: passed as --rd 147.09

HHT detail packs (generated from SN residuals)

Example baseline:plamb_runs/hht_sn_FRPBH/hht_sn_detail.npz

Stores decomposed IMFs, instantaneous frequencies, redshift alignment.

Variants created under different conditions:

plamb_runs/noether_grid/soft{s}/hht_sn_detail.npz (Noether strength sweep).

plamb_runs/hht_sn_FRPBH_z0.05_1.5/hht_sn_detail.npz (restricted redshift windows).

These NPZ packs are the direct input to BAO-locking tests.

Output / Analysis Tables

Heatmap correlation tables (BAO locking):plamb_runs/hht_bao_lock_real/baseline_s9_p10k/bao_lock_heatmap_table_with_p.csv

Main record of Pearson r-values + permutation test significance.Pearson correlation coefficient r measures linear alignment between IMF frequencies and BAO scales.

Direction scores (DM vs DH)plamb_runs/noether_grid/summary/direction_scores.csv

Captures Δr between DM/rd and DH/rd correlations.

IMF2/IMF3 vs strength summariesplamb_runs/noether_grid/summary/p10k_IMF2_IMF3_vs_strength_s9.csv

Collates BAO-locking strength dependence across Noether sweep.

Appendix - Noether Strength in this Analysis

1. Background: Noether’s Theorem

Noether’s theorem (1918) is one of the foundational principles of modern physics.

It states: every symmetry of the action corresponds to a conservation law.

Examples…

time-translation symmetry → conservation of energy.

Space-translation symmetry → conservation of momentum.

Rotational symmetry → conservation of angular momentum.

In cosmology, these symmetries and their conservation laws govern how the large-scale dynamics of the Universe evolve.

Note:In the presence of Mass, Noether symmetry is potentially broken or distorted and may imply non-linear behaviour.

2. Role in the Model

In this project  Full_Relativity (FR) + Parent_Black_Hole (PBH) frameworks, were explicitly tested, relaxing or enforcing certain conservation symmetries affects and subsequently testing the residual signal structures.

The Noether strength parameter was introduced as a control knob that adjusts the degree to which conservation-law constraints are applied in the modelling of cosmic expansion.

Soft = 0 → minimal/no constraint. Residual IMFs are free to vary more flexibly.

Higher soft values (2,4,…,12) → progressively stronger enforcement of conservation conditions… time-energy, space-momentum balance, etc.

Effectively, this parameter interpolates between a fully unconstrained system and one that obeys strict Noetherian symmetry rules.

3. Why HHT was used in Analysis

The goal was to see whether IMF–BAO locking (the correlation between oscillatory IMFs in SN data and BAO scales) depends on how strongly conservation laws are enforced.

By running the HHT decomposition under different Noether strengths, allowed testing of:

BAO-locking if robust across all configurations → suggests a fundamental physical origin.

BAO-locking if it appears only under certain symmetry regimes → may indicate a deeper connection to conservation laws in the underlying cosmology.

4. What Was Found to-date Nov-2025

IMF2 correlation with BAO was significant at multiple Noether strengths (0, 2, 4, 8, 12).

At intermediate strengths (6, 10), correlations weakened or disappeared.

IMF3 correlations were weaker, more variable, but showed a notable spike at strength=6.

This suggests that the way conservation-law symmetries are applied in the model alters the detectability of oscillatory features in cosmic expansion.

5. Novel Behavour

Standard ΛCDM does NOT allow for a “tuneable” conservation-law symmetry — these are fixed assumptions.

By parameterising Noether strength, the analysis explored whether small departures or modifications to conservation enforcement might produce observable residual structures.

If oscillatory correlations (IMF–BAO locking) depend on Noether tuning, this hints at dynamical symmetry-breaking or restoration in the fabric of spacetime itself.

That is potentially relevant for extensions of cosmology that involve:

Dynamical dark energy.

Parent–Back-Hole (PBH) embedding scenarios.

👉 In simple terms:

Noether strength is a test parameter used to dial how strongly conservation rules are applied. By varying it, IMF–BAO correlations are showing strong signs of NOT being random, but appear preferentially under certain symmetry regimes.

Notably this suggests that conservation-law structure may leave a detectable imprint on cosmic expansion.

Appendix – Python Programs Used in Analysis

This appendix lists the primary Python (.py) programs employed in the HHT–BAO locking analysis. They are currently off-line, however, can be made available upon request.Each script is identified by name, with a brief description of its role in the overall workflow. This serves as a reproducibility index so that future readers can trace results back to the exact computational components.

Note: Script versions correspond to the Revision 1V8 analysis.

Appendix Figures: IMF2 Oscillations and BAO Locking

This section provides visualizations of IMF2 extracted from Pantheon+SH0ES residuals, together with its comparison against BAO effective scales. These figures illustrate the oscillatory nature of IMF2 and its apparent statistical alignment ("locking") with BAO, complementing the surrogate and correlation tests reported in Revision 1V8.

Figure X. IMF2 Oscillations vs log(1+z)

IMF2 extracted from Pantheon+SH0ES residuals using HHT. The mode exhibits clear oscillatory structure across 0 ≲ z ≲ 1.5. This visualization complements the permutation results by showing the time-localized nature of the signal rather than aggregate statistics alone.

Figure Y. IMF2 Instantaneous Frequency vs BAO 1/k (normalized)

Hilbert-phase instantaneous frequency of IMF2 (normalized) overlaid with a normalized BAO 1/k proxy interpolated onto the same log(1+z) grid. Visual co-variation supports the quantitative locking results (Pearson r and surrogate p-values reported in Table p10k_IMF2_IMF3_vs_strength_s9). The BAO curve here is a proxy (1/value); analyses use the full BAO pipeline described in Methods.

Figure Z. IMF2 instantaneous energy versus log(1+z).

Energy was computed as , where is the Hilbert amplitude of IMF2 extracted from Pantheon+SH0ES residuals. The pale blue line shows instantaneous normalized energy, while the orange curve shows the binned median (40 bins). Vertical dashed lines mark BAO measurement redshifts. IMF2 energy is concentrated in several broad redshift intervals that overlap with BAO leverage windows, reinforcing the interpretation of IMF2–BAO locking while helping to identify and exclude edge-driven artefacts.

Plot Summary

Figures X, Y and Z provide a visual complement to the statistical analysis:

Figure X highlights the oscillatory nature of IMF2,

Figure Y shows its phase–frequency co-variation with BAO effective scales, and

Figure Z demonstrates where IMF2 energy is concentrated across redshift.

Viewed jointly, these figures reinforce the interpretation that IMF2 encodes a robust oscillatory structure that is both temporally localized and statistically aligned with BAO scales.

Document End…

## Extracted Tables

### Table 1

noether_strength | IMF | r_1_over_k | p_1_over_k | r_with_k | smooth | perms | alpha | source_csv
0 | 2 | 0.1651088189705540 | 9.99900009999e-05 | -0.1725459463255040 | 9 | 10000 | 0.05 | /Users/boyde/.spyder-py3/plamb_runs/noether_grid/heatmap_p10k_s9_soft0/bao_lock_heatmap_table_with_p.csv
2 | 2 | 0.1310356619442890 | 9.99900009999e-05 | -0.1205091005467510 | 9 | 10000 | 0.05 | /Users/boyde/.spyder-py3/plamb_runs/noether_grid/heatmap_p10k_s9_soft2/bao_lock_heatmap_table_with_p.csv
4 | 2 | 0.1120523264276220 | 0.0011998800119988 | -0.0640636026071066 | 9 | 10000 | 0.05 | /Users/boyde/.spyder-py3/plamb_runs/noether_grid/heatmap_p10k_s9_soft4/bao_lock_heatmap_table_with_p.csv
6 | 2 | 0.0266839786481922 | 0.2335766423357660 | -0.0311768852289013 | 9 | 10000 | 0.05 | /Users/boyde/.spyder-py3/plamb_runs/noether_grid/heatmap_p10k_s9_soft6/bao_lock_heatmap_table_with_p.csv
8 | 2 | 0.1900919800514320 | 9.99900009999e-05 | -0.1445015189759560 | 9 | 10000 | 0.05 | /Users/boyde/.spyder-py3/plamb_runs/noether_grid/heatmap_p10k_s9_soft8/bao_lock_heatmap_table_with_p.csv
10 | 2 | 0.0267918406463709 | 0.2318768123187680 | -0.0213092313106571 | 9 | 10000 | 0.05 | /Users/boyde/.spyder-py3/plamb_runs/noether_grid/heatmap_p10k_s9_soft10/bao_lock_heatmap_table_with_p.csv
12 | 2 | 0.1488858114347890 | 9.99900009999e-05 | -0.1312118264701580 | 9 | 10000 | 0.05 | /Users/boyde/.spyder-py3/plamb_runs/noether_grid/heatmap_p10k_s9_soft12/bao_lock_heatmap_table_with_p.csv
0 | 3 | -0.0713368567356735 | 0.9749025097490250 | 0.1302465847648620 | 9 | 10000 | 0.05 | /Users/boyde/.spyder-py3/plamb_runs/noether_grid/heatmap_p10k_s9_soft0/bao_lock_heatmap_table_with_p.csv
2 | 3 | -0.1248833038139610 | 0.9997000299970000 | 0.2039084972575360 | 9 | 10000 | 0.05 | /Users/boyde/.spyder-py3/plamb_runs/noether_grid/heatmap_p10k_s9_soft2/bao_lock_heatmap_table_with_p.csv
4 | 3 | 0.01608384524189 | 0.3377662233776620 | 0.0507316739815324 | 9 | 10000 | 0.05 | /Users/boyde/.spyder-py3/plamb_runs/noether_grid/heatmap_p10k_s9_soft4/bao_lock_heatmap_table_with_p.csv
6 | 3 | 0.1895534169922140 | 9.99900009999e-05 | -0.0985911130967141 | 9 | 10000 | 0.05 | /Users/boyde/.spyder-py3/plamb_runs/noether_grid/heatmap_p10k_s9_soft6/bao_lock_heatmap_table_with_p.csv
8 | 3 | -0.1962013741209420 | 1.0 | 0.2422076905945880 | 9 | 10000 | 0.05 | /Users/boyde/.spyder-py3/plamb_runs/noether_grid/heatmap_p10k_s9_soft8/bao_lock_heatmap_table_with_p.csv
10 | 3 | 0.0152793057380453 | 0.3304669533046700 | 0.0425474773983463 | 9 | 10000 | 0.05 | /Users/boyde/.spyder-py3/plamb_runs/noether_grid/heatmap_p10k_s9_soft10/bao_lock_heatmap_table_with_p.csv
12 | 3 | 0.0085729049153515 | 0.4094590540945910 | 0.0608182767101917 | 9 | 10000 | 0.05 | /Users/boyde/.spyder-py3/plamb_runs/noether_grid/heatmap_p10k_s9_soft12/bao_lock_heatmap_table_with_p.csv

### Table 2

Program | Description | Role in Workflow
hht_sn_residuals.py | Performs Hilbert–Huang Transform (HHT) on supernova residuals. | Generates IMF decomposition and instantaneous frequency spectra.
HHT_BAO_curated.py | Prepares HHT outputs for BAO comparison. | Aligns IMF frequencies with BAO scales for correlation tests.
PLamb_Test_10V6c_plus.py | Runs cosmological fits under ΛCDM, FR, and PBH models. | Generates model baselines and SN residuals for HHT input.
bao_lock_test.py | Computes IMF–BAO correlation coefficients. | Quantifies locking strength using Pearson r and permutation tests.
surrogate_generator.py | Creates random phase-shuffled surrogate datasets. | Provides null distributions for significance testing.
hht_sn_FRPBH.py | Specialised FR+PBH HHT analysis module. | Explores IMF behaviour under embedding cosmology assumptions.
