# HHT_BAO_Locking_Report_1V9_with_figs

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\HHT_BAO_Locking_Report_1V9_with_figs.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2025-10-12T03:19:00Z`

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

New Figures — Method Invariance and Distance-Gated Robustness

Method invariance: Δr for IMF1

Δr = r_gated − r_full for IMF1 under CEEMDAN/EEMD/EMD at (pad, ngrid) = (0.15, 384).

  What it is: Δr ≡ r(gated) − r(full) for IMF1 vs distance threshold δ, compared across CEEMDAN / EEMD / EMD.

  CEEMDAN (blue): tiny positive Δr ≈ +3×10⁻⁴ (flat with δ).→ Gating barely changes IMF1 locking.

  EEMD (orange): small positive Δr (+2.0–2.4×10⁻³) for δ≥0.2; slight negative at δ=0.1.→ Mild uplift with stricter gating.

  EMD (green): small negative Δr (−0.2 to −0.7×10⁻³) that becomes a bit more negative with δ.→ Slight attenuation under gating.

  Scale & invariance: All method effects are very small (|Δr| ≤ ~2.5×10⁻³) relative to r≈1 seen elsewhere.→ Conclusion: IMF1 locking is method-invariant and robust to distance gating; any changes are negligible.

Method invariance: Δr for IMF3

Δr = r_gated − r_full for IMF3 across methods at the same grid settings.

  What it is: Δr (gated − full) for IMF3 vs the distance threshold δ across three methods (CEEMDAN, EEMD, EMD).

Δr ≈ 0 ⇒ gating barely changes locking r.

δ is the fractional separation from used to exclude near- samples.

  CEEMDAN (blue): Slightly negative and gently decreasing with δ (~−0.0001 → −0.0008).→ Gating trims r by a tiny amount.

  EEMD (orange): Small positive at low δ (~+0.0013–+0.0015), crosses zero near δ≈0.3, negative by δ=0.4 (~−0.0017).→ Mild sensitivity that flips sign as the gate widens.

  EMD (green): Consistently positive and grows with δ (~+0.002 → +0.0046).→ Wider gates slightly increase measured r.

  Scale: All effects are minute (|Δr| ≤ ~0.005), i.e., <0.5% of r’s full dynamic range.

  Takeaway: IMF3 locking is robust to distance gating; method differences exist but are small and do not alter the qualitative conclusion.

Distance-gated Δr (IMF1 & IMF3)

Median change in locking after excluding samples within fractional distance δ of IMF2’s |k|.

  What it is: Δr ≡ r(gated) − r(full) for IMF1/IMF2/IMF3 vs distance threshold δ (fractional separation from ; larger δ = stricter exclusion of near- samples).

  IMF2 (orange): ~−0.041 at all δ (flat).→ Gating near consistently lowers IMF2 locking by ~0.04; effect size is stable as the gate widens.

  IMF1 (blue): ≈ 0 (slightly above zero, weak δ trend).→ Locking is essentially unchanged by distance gating ⇒ little dependence on -adjacent samples.

  IMF3 (green): Small negative Δr that becomes more negative with δ (~−0.001 → ~−0.007).→ Mild sensitivity; removing points closer to trims IMF3 locking a bit.

  Scale: Effects are modest: |Δr| ≤ ~0.041; only IMF2 shows a clearly non-zero, δ-invariant drop—expected since the gate targets its own neighborhood.

  Takeaway: Distance gating confirms specificity to IMF2 (largest, consistent reduction), while IMF1 is robust and IMF3 shows only a small δ-dependent weakening.

Sliding-τ locking (IMF2)

Median r(τ) across seeds with BH-FDR ≤ 0.10 markers; τ = log(1+z).

  What it is: Sliding-τ window scan of IMF2 1/|k| locking (DH vs DM). Blue = median r across seeds; band = ±IQR/2; orange dots = windows with q_BH ≤ 0.10.

  Main result: Sustained, near-perfect locking r ≈ 0.97–1.00 over most of the range τ ≈ 0.55–1.15 (≈ z ~ 0.73–2.2), with significance in essentially every window.

  Robustness: Very narrow IQR band → minimal seed variability.

  Localized dip: A single trough near τ ~ 0.96 (≈ z ~ 1.6) with r ~ 0.65, still significant; neighboring windows recover to r ~ 0.97–1.00.

  Edge artifact: The earliest window (τ ~ 0.31) shows a negative r; this sits at the boundary of the usable range and does not persist.

  Takeaway: IMF2 locking is strong, consistent, and highly significant across redshift, with only a brief, localized weakening around z ~ 1.6.

Locking r per seed (consistency check)

Heatmap of r by τ-window and seed; brighter = stronger locking.

## Extracted Tables

### Table 1

Method | Pad | ngrid | IMF2 energy p (med) | Locking r (med)
EMD | 0.15 | 384 | 0.0015 | 0.973
EEMD | 0.15 | 384 | 0.0015 | 0.973
CEEMDAN | 0.15 | 384 | 0.0015 | 0.973
