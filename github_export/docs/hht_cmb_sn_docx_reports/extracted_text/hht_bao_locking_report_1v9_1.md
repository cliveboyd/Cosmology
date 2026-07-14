# HHT_BAO_Locking_Report_1V9 (1)

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\HHT_BAO_Locking_Report_1V9 (1).docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2025-10-12T02:47:00Z`

## Extracted Text

Cosmology HHT–BAO Locking Report — 1V9 (Concise)

Executive Summary

IMF2 is uniquely BAO‑linked: median locking r ≈ 0.973 at (pad=0.15, ngrid=384); q ≤ 0.0015 across methods.

IMF2 energy significance aligns: p ≈ 0.0015 (median across seeds).

Sliding‑τ localization: 17/18 windows significant (q≤0.10); peaks near τ≈0.62–0.82 (z≈0.85–1.26).

Method invariance: EMD/EEMD/CEEMDAN agree on top cells; strong ensembles match baseline.

New in 1V9

Distance‑gated specificity: Δr(IMF1) ≈ 0; Δr(IMF3) small (|Δr|≲5×10⁻³); IMF2 unaffected.

Method‑invariance overlays (CEEMDAN/EEMD/EMD) show negligible gating impact on IMF1/3.

Tiny‑δ windowed scans around k₂: stable IMF2 locking for δ ∈ {0, 0.01, 0.02, 0.05}.

Recommended Defaults

Decomposition: CEEMDAN (EEMD/EMD equivalent for top cells).

Grid: pad = 0.15, ngrid = 384; Seeds = {3, 7, 11, 23}.

Metrics: IMF2 locking r and energy ratio; Surrogates ns≈1000; BH‑FDR over grids/windows.

Consensus Best Cells (Energy + Locking)

Core Definitions / Equations

τ = log(1 + z)

Analytic signal: a(τ) = x(τ) + i·H{x(τ)} (Hilbert transform H{})

Instantaneous phase: φ(τ) = unwrap(arg[a(τ)])

Instantaneous wavenumber: k(τ) = (1 / 2π) · dφ/dτ

Locking r = corr( 1/|k_IMF2|_DH , 1/|k_IMF2|_DM ) over τ‑support

IMF2 energy ratio: E2 / ΣY², with surrogates via phase randomization (p‑value)

FDR: Benjamini–Hochberg over cells/windows

Acronyms & Variables

BAO — Baryon Acoustic Oscillation

HHT — Hilbert–Huang Transform; IMF — Intrinsic Mode Function

EMD/EEMD/CEEMDAN — Empirical Mode Decomposition variants

DH/rd, DM/rd — Radial/Transverse BAO distance over sound horizon

τ — log(1+z); z — redshift

k(τ) — instantaneous wavenumber; 1/|k| — instantaneous wavelength (τ units)

r — Pearson correlation of DH vs DM 1/|k| within a window

δ — distance gate: |k_imf − k₂| / |k₂| ≥ δ

ns — # surrogates; q_BH — BH‑FDR adjusted p‑value

Figures (1V9 Additions)

## Extracted Tables

### Table 1

Method | pad | ngrid | IMF2 energy p (med) | Locking r (med)
EMD | 0.15 | 384 | 0.0015 | 0.973
EEMD | 0.15 | 384 | 0.0015 | 0.973
CEEMDAN | 0.15 | 384 | 0.0015 | 0.973
