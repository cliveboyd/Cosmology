# HHT_Locking_Results_Summary_v1

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\HHT_Locking_Results_Summary_v1.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2025-10-20T10:30:00Z`

## Extracted Text

Hilbert–Huang Transform (HHT) Locking Analysis – Results Summary

1. Purpose and Investigation Scope

The HHT locking analysis investigates whether supernova (SN Ia) residuals,—defined as deviations between observed and model-predicted distance moduli,—contain intrinsic quasi-periodic structure that correlates (‘locks’) with Baryon Acoustic Oscillation (BAO) scale evolution as a function of redshift.

This analysis aims to detect non-random, scale-coherent oscillations potentially arising from variable-speed-of-light or accretion-induced space-time fabric (STF) perturbations in the Full Relativity (FR) and Kerr–de Sitter cosmological frameworks.

2. Methodology and Testing Procedure

Key methodological components:

• Data Inputs: Pantheon+SH0ES supernova data, eBOSS DR16 BAO scales, and optional CMB/H(z) datasets.• Signal Extraction: Empirical Mode Decomposition (EMD), EEMD, and CEEMDAN were used to extract intrinsic mode functions (IMFs) from Δμ(z).• Hilbert Transform: Each IMF’s instantaneous frequency and amplitude were derived, defining energy spectra E_i(z) and frequency ridges f_i(z).• Surrogate Testing: 10,000 phase-randomized surrogates were used to form null distributions of correlation metrics.• Locking Analysis: Examined whether IMF frequency ridges clustered around BAO-effective frequencies, implying phase or frequency resonance.

3. Results Found to Date

• IMF2 Locking – Core Finding: IMF2 shows strong, statistically significant correlation (p ≈ 0.005–0.01) with BAO scale evolution, stable across seeds and methods. This suggests a genuine physical coupling between SN residual microstructure and the BAO acoustic scale.

• IMF3 Behaviour: IMF3 displays intermittent phase alignment and amplitude correlations with integrated entropy or AGN density proxies—possibly indicating entropy perturbations or accretion-driven oscillations.

• Surrogate Validation: Random-phase ensembles confirm IMF2 locking cannot arise from stochastic signals at the 99.5% confidence level.

• Noether Strength Effects: Parameter sweeps (0–12) reveal strongest locking around strengths 4–6, with excessive constraint (>10) dampening IMF activity. This indicates optimal partial symmetry breaking conditions for natural oscillatory behavior.

Noether Strength Effects – Detailed Interpretation

Conceptual Basis

In the Full Relativity / Kerr–de Sitter framework, Noether Strength (Nₛ) is a tuneable parameter that quantifies how strongly Noether conservation constraints (energy–momentum–angular-momentum symmetry) are enforced within the numerical model generating the supernova residuals.

Mathematically, the model assumes a generalized Lagrangian for the space–time fabric (STF) that includes a symmetry-coupling term:

where

is the base cosmological Lagrangian (matter + Λ + accretion terms),

encodes conservation coupling, i.e. how tightly the field obeys energy–momentum conservation at local curvature boundaries,

defines the Noether Strength, effectively the rigidity of conservation enforcement.

Thus:

Low Nₛ (≈ 0–2) → weak symmetry enforcement → highly plastic STF, energy flow not tightly conserved → noisy, incoherent IMFs.

Intermediate Nₛ (≈ 4–6) → balanced enforcement → partial symmetry breaking → emergence of coherent oscillatory modes (IMF₂, IMF₃).

High Nₛ (≈ 10–12) → over-constrained STF → suppression of adaptive oscillations → IMF energy flattens.

Physical Interpretation

At Nₛ ≈ 4–6, energy conservation is relaxed enough to allow dynamic curvature exchange between neighbouring redshift shells, while still constrained sufficiently to maintain phase coherence. This condition is equivalent to soft symmetry breaking, the same mechanism that allows spontaneous oscillations in many field theories (analogous to Goldstone or phonon modes in condensed-matter physics).

Empirical Signature

The HHT results show:

IMF₂ energy peaks sharply at Nₛ ≈ 5.

Locking correlation maximizes in the same range (0.09–0.10, p ≈ 0.005).

Phase dispersion across redshift bins is minimal (Δφ < 0.2 rad).

As Nₛ increases past 10, both IMF₂ amplitude and locking strength decay roughly exponentially, confirming that oscillations depend on partial—not total—symmetry preservation.

Cosmological Implication

This behaviour suggests that the observable Universe may sit near a dynamically critical Noether regime, where symmetry is neither fully enforced (ΛCDM-like rigidity) nor entirely broken (chaotic FR limit).Such a regime naturally supports persistent standing-wave patterns in the expansion field—precisely what IMF₂ reveals as BAO-scale locking.

4. Interpretation and Implications

Evidence of structured residuals — SN Ia residuals are not fully stochastic; they contain reproducible oscillatory components aligned with known cosmic scales.2. IMF2–BAO phase coherence — Suggests coupling between distance–modulus fluctuations and sound-horizon geometry, possibly due to variable-c or STF resonance effects.3. Entropy–Energy Connection — IMF3 amplitudes may trace energy feedback from Primordial (Parent) Black Hole PBH entropy growth or Active Galactic Nucleus AGN-scale accretion events.4. Consistency Across Seeds/Methods — Persistence across EMD/EEMD/CEEMDAN confirms intrinsic signal origin.5. Statistical Confidence — p ≈ 0.005 under 10,000 surrogates indicate ~3σ significance, motivating expanded surrogate testing.

Next Steps

• Extend surrogate testing to 50k–100k runs for higher robustness.• Isolate low-z vs high-z behaviour for mode evolution tracking.• Cross-validate with CMB and JWST high-z SN datasets.• Integrate with Bayesian cosmological fits to estimate locking strength as an inferred coupling parameter.

## Extracted Tables

### Table 1

Noether Strength (Nₛ) | STF Behaviour | Observed IMF Pattern | Interpretation
0 – 2 | Symmetry almost free; local metric unconstrained | Disordered IMFs, broad spectral energy | Chaotic regime; resembles unregulated entropy field
4 – 6 | Partial symmetry breaking, mild curvature-energy exchange | Narrow, coherent IMF₂ ridge aligned with BAO | Resonant zone — STF oscillates naturally under conservation–violation balance
8 – 10 | Symmetry re-tightening | IMF amplitudes reduce, phase drifts appear | Damped resonance; system regains rigidity
> 10 | Fully enforced symmetry | IMFs collapse toward baseline | Over-constraint suppresses oscillatory degrees of freedom

### Table 2

IMF | Behaviour | Locking Target | p-value Range | Interpretation
1 | High-frequency noise | — | — | Noise / numerical residue
2 | Stable oscillation | BAO scales | 0.005–0.01 | Real physical locking
3 | Intermittent, entropy-linked | AGN / PBH activity | 0.05–0.2 | Possible entropy resonance
4+ | Long-trend modes | — | > 0.5 | Baseline curvature
