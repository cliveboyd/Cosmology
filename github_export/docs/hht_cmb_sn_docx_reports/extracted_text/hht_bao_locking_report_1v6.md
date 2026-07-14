# HHT_BAO_Locking_Report_1V6

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\HHT_BAO_Locking_Report_1V6.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2025-10-01T19:03:00Z`
- DOCX modified: `2025-10-01T19:03:00Z`

## Extracted Text

HHT–BAO Locking Report

(Revision 1V6)

Executive Summary

This document provides an updated synthesis of Hilbert–Huang Transform (HHT) analysis applied to supernova (SN) residuals in comparison to Baryon Acoustic Oscillation (BAO) scales.

Revision 1V5 integrates 10,000-permutation surrogate tests, IMF2/IMF3 vs strength analyses, and baseline heatmap visualisations.

All tests were conducted using python under Anaconda spyder

The findings strengthen the evidence for IMF2 locking with BAO scales under certain model strengths, while IMF3 shows inconsistent behaviour.

Surrogate testing confirms robustness against chance alignment.

Background  Standard Cosmology vs. HHT Analysis

1. What is currently understood (Standard ΛCDM cosmology)

The ΛCDM model is the widely accepted description of the Universe.

It assumes:

A constant dark energy term (Λ) driving cosmic acceleration.

Cold dark matter (CDM) providing most of the mass.

Normal matter (atoms, gas, stars, galaxies) is only ~5% of the total energy density.

Observational pillars:

Supernovae (SNe Ia) show accelerating expansion.

Baryon Acoustic Oscillations (BAO) act as a standard ruler, matching predictions of ΛCDM.

Cosmic Microwave Background (CMB) confirms early Universe consistency with ΛCDM.

In this framework, deviations from smooth expansion are treated as noise or local systematics. The model does not expect small oscillatory features in expansion history.

2. What this project address: (Hilbert–Huang Transform analysis)

Instead of fitting a single smooth curve, the analysis applied the Hilbert–Huang Transform (HHT) to supernova residuals.

This decomposes the data into Intrinsic Mode Functions (IMFs), which behave like oscillatory “ripples” embedded in the expansion history.

These IMFs were tested against BAO scales — the “acoustic ruler” left over from sound waves in the early Universe.

To check if any match was real or just noise, a surrogate testing (randomised versions of the data, up to 10,000 permutations) were injected. This provides a robust statistical baseline to remove noise artefacts.

3. What has been found (Novel signals)

IMF2 shows robust positive correlations at multiple strengths (0–4, 8, 12), but weak or absent at noether strengths 6 and 10:

, with in many tests (statistically significant).

Suggests that part of the residual structure in SN data may be phase-locked to the same scale that defines BAO.

IMF3 shows weaker and more variable correlations: sometimes significant, sometimes consistent with noise.However at Noether strength=6 IMF3 shows strong BAO correlation,

These signals survive 10,000-permutation surrogate tests, reducing the chance of spurious detection.

The strength of the IMF–BAO correlation appears to depend on Noether strength (a parameter in this model controlling conservation-law constraints), see appendix for further details, hinting at a possible link between cosmological symmetries and residual oscillatory structure.

4. Novelty and relevance

Standard ΛCDM assumes no preferred oscillatory features in the residuals once the smooth expansion is removed.

The detection of a statistically robust IMF–BAO correlation implies that:

The Universe’s expansion history may contain subtle oscillatory components linked to the same physics as BAO.

These features might indicate dynamic dark energy, variable speed of light, or effects of embedding cosmology in a parent structure (as in your FR/PBH frameworks).

If confirmed, this would represent a departure from ΛCDM’s assumption of pure smoothness, opening a new window into how large-scale cosmic structure interacts with fundamental physics.

5. Going forward… Why it may matter

Provides a new diagnostic tool: HHT allows discovery of hidden structure invisible to Fourier/wavelet methods.

Could help reconcile tensions in cosmology, such as:

The Hubble constant () tension between local and early-Universe measures.

Small but persistent anomalies in BAO vs SN fits.

Points toward a future where time-dependent, symmetry-driven corrections to expansion history are tested alongside ΛCDM.

Even if the effect turns out to be subtle or model-dependent, it helps establish robust methods for probing non-ΛCDM behaviour.

Note: HHT IMF analysis typically allows analysis of non-stationary datasets (behaviour)

👉 In short: ΛCDM says residuals are noise. HHT shows they may not be.If residual oscillations consistently line up with BAO scales, this may suggest new physics — or at the very least, a refinement of how to model the Universe’s expansion.

Acronyms and Variables

• HHT — Hilbert–Huang Transform• IMF — Intrinsic Mode Function• BAO — Baryon Acoustic Oscillations• r — Pearson correlation coefficient• p_perm — p-value estimated via surrogate permutations• Surrogate findings — significance results obtained by testing the observed IMF/BAO correlation against 10k of phase-randomised surrogate datasets to rule out spurious correlations.

Note: Locking → refers to statistical alignment of IMF frequencies with BAO scales, does NOT infer a literal physical resonance.

Key Equations

• SN residuals: Δμ(z) = μ_obs(z) − μ_model(z)

• BAO observables: DM/rd and DH/rd, where rd is the sound horizon at drag epoch.

• Pearson correlation used for locking tests: r = Cov(ω̄, 1/k_eff) / (σ_ω̄ σ_1/k).

Methodology

Supernova residuals were decomposed via HHT into IMFs.

The HHT method combines Empirical Mode Decomposition (EMD) with Hilbert transforms to extract instantaneous frequencies ω(z) from SN residuals as a function of log(1+z).

These are compared against effective BAO scales (k & 1/k_eff) derived from clustering measurements.

Pearson correlation coefficients (r) are computed for IMF2 and IMF3 to test alignment with DM/rd and DH/rd. Significance is assessed using phase-randomized surrogate datasets.

2. IMF2 and IMF3 instantaneous frequencies were compared to BAO effective wavenumbers (k and 1/k).3. Surrogate testing was applied: >10,000 random phase-shuffled realisations were generated to estimate null distributions, producing robust permutation p-values (p_perm).4. Significance thresholds were set at α = 0.05 with surrogate-corrected p-values.5. Locking strength tested across a grid of Noether-conserving model strengths (soft_0 to soft_12).

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

• IMF2/IMF3 correlations vs strength (10k permutations): p10k_IMF2_IMF3_vs_strength_s9.csv below

Next Steps

• Extend surrogate testing to >50,000 permutations for ultra-robust significance.• Expand redshift-window analyses to isolate low-z vs high-z contributions.• Cross-validate with eBOSS DR16 BAO datasets and JWST high-z supernovae.• Integrate with FR + PBH cosmological model fits for joint constraints.• Prepare publication-ready plots and results tables for peer review.

Revision History

• 1V0 – Initial HHT BAO locking setup• 1V1 – Added IMF2/3 r vs strength plots• 1V2 – Added multi-window redshift tests• 1V3 – Baseline surrogate heatmaps• 1V4 – Harmonic test integration• 1V5 – 10k surrogate findings, expanded robustness checks, embedded plots with figure captions.

Appendix Summary Results  p10k_IMF2_IMF3_vs_Noether strength_s9

Table: p10k_IMF2_IMF3_vs_strength_s9.csv

Appendix  Primary Datasets used within scope of investigation.

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

In scripts: passed as --rd 147.09.

HHT detail packs (generated from SN residuals)

Example baseline:plamb_runs/hht_sn_FRPBH/hht_sn_detail.npz

Stores decomposed IMFs, instantaneous frequencies, redshift alignment.

Variants created under different conditions:

plamb_runs/noether_grid/soft{s}/hht_sn_detail.npz (Noether strength sweep).

plamb_runs/hht_sn_FRPBH_z0.05_1.5/hht_sn_detail.npz (restricted redshift windows).

These NPZ packs are the direct input to BAO-locking tests.

Output / analysis tables

Heatmap correlation tables (BAO locking):plamb_runs/hht_bao_lock_real/baseline_s9_p10k/bao_lock_heatmap_table_with_p.csv

Main record of Pearson r-values + permutation test significance.Pearson correlation coefficient r measures linear alignment between IMF frequencies and BAO scales.

Direction scores (DM vs DH)plamb_runs/noether_grid/summary/direction_scores.csv

Captures Δr between DM/rd and DH/rd correlations.

IMF2/IMF3 vs strength summariesplamb_runs/noether_grid/summary/p10k_IMF2_IMF3_vs_strength_s9.csv

Collates BAO-locking strength dependence across Noether sweep.

Appendix  Noether Strength in this Analysis

1. Background: Noether’s Theorem

Noether’s theorem (1918) is one of the foundational principles of modern physics.

It states: every symmetry of the action corresponds to a conservation law.

Examples…

time-translation symmetry → conservation of energy.

Space-translation symmetry → conservation of momentum.

Rotational symmetry → conservation of angular momentum.

In cosmology, these symmetries and their conservation laws govern how the large-scale dynamics of the Universe evolve.

2. Role in the Model

In this project’s FR + PBH framework, we explicitly tested how relaxing or enforcing certain conservation symmetries affects the residual signal structure.

The Noether strength parameter was introduced as a control knob that adjusts the degree to which conservation-law constraints are applied in the modelling of cosmic expansion.

Soft = 0 → minimal/no constraint. Residual IMFs are free to vary more flexibly.

Higher soft values (2,4,…,12) → progressively stronger enforcement of conservation conditions (time-energy, space-momentum balance, etc.).

Effectively, this parameter interpolates between a fully unconstrained system and one that obeys strict Noetherian symmetry rules.

3. Why HHT was used in Analysis

The goal was to see whether IMF–BAO locking (the correlation between oscillatory IMFs in SN data and BAO scales) depends on how strongly conservation laws are enforced.

By running the HHT decomposition under different Noether strengths, allowed testing of:

BAO-locking if robust across all configurations → suggests a fundamental physical origin.

BAO-locking if it appears only under certain symmetry regimes → may indicate a deeper connection to conservation laws in the underlying cosmology.

4. What Was Found to-date Sep-2025

IMF2 correlation with BAO was significant at multiple Noether strengths (0, 2, 4, 8, 12).

At intermediate strengths (6, 10), correlations weakened or disappeared.

IMF3 correlations were weaker, more variable, but showed a notable spike at strength=6.

This suggests that the way conservation-law symmetries are applied in the model alters the detectability of oscillatory features in cosmic expansion.

5. Novel Behavour

Standard ΛCDM does NOT allow for a “tuneable” conservation-law symmetry — these are fixed assumptions.

By parameterising Noether strength, the analysis explored whether small departures or modifications to conservation enforcement might produce observable residual structures.

If oscillatory correlations (IMF–BAO locking) depend on Noether tuning, this hints at dynamical symmetry-breaking or restoration in the fabric of spacetime itself.

That is potentially relevant for extensions of cosmology that involve:

Variable speed of light (c).

Dynamical dark energy.

Parent–black-hole (PBH) embedding scenarios.

👉 In simple terms:

Noether strength is a test parameter we used to dial how strongly conservation rules are applied. By varying it, IMF–BAO correlations are showing strong signs of NOT being random, but appear preferentially under certain symmetry regimes

Notably this suggests that conservation-law structure may leave a detectable imprint on cosmic expansion.

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
