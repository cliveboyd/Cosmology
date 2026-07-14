# HHT_BAO_Locking_Report_2V0

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\HHT_BAO_Locking_Report_2V0.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2025-10-12T06:00:00Z`

## Extracted Text

Cosmology HHT–BAO Locking Report — 2V0

Concise update with latest distance-gated and method-invariance results, plus a glossary of equations, acronyms, and variables.

Executive Summary

IMF2 uniquely tracks the BAO scale: median locking r ≈ 0.973 at (pad=0.15, ngrid=384); q ≤ 0.0015 across methods.

Energy significance aligns with locking: IMF2 energy p ≈ 0.0015 (median across seeds).

Sliding‑τ localization: 17/18 windows significant (q≤0.10); peaks near τ≈0.62–0.82 (z≈0.85–1.26).

Method invariance: EMD/EEMD/CEEMDAN agree on best cells; strong ensembles match baseline.

Distance‑gated specificity: IMF2 shows the expected δ‑invariant reduction; IMF1 ≈ 0 change; IMF3 small |Δr| (≤5×10⁻³).

Recommended Defaults

Decomposition: CEEMDAN (EEMD/EMD equivalent for top cells).

Grid: pad = 0.15, ngrid = 384; Seeds = {3, 7, 11, 23}.

Metrics: IMF2 locking r and energy ratio; Surrogates ns≈1000; BH‑FDR over grids/windows.

Core Definitions & Equations

τ = log(1 + z)

Analytic signal: a(τ) = x(τ) + i·H{x(τ)} (Hilbert transform H{})

Instantaneous phase: φ(τ) = unwrap(arg[a(τ)])

Instantaneous wavenumber: k(τ) = (1 / 2π) · dφ/dτ

Instantaneous wavelength: 1/|k(τ)| (in τ units)

Locking: r = corr( 1/|k_IMF2|_DH , 1/|k_IMF2|_DM ) within a τ‑window or full τ‑support

IMF2 energy ratio: E2 / ΣY²; p via phase‑randomized surrogates

FDR: Benjamini–Hochberg (q_BH) over grids/windows

Distance gate: exclude samples with |k_imf − k₂|/|k₂| < δ (report Δr = r_gated − r_full)

Acronyms, Symbols, and Definitions

• HHT — Hilbert–Huang Transform; EMD/EEMD/CEEMDAN are its decomposition variants.• IMF — Intrinsic Mode Function. IMF2 denotes the second mode; IMF1/IMF3 are neighboring modes.• k(τ) — instantaneous angular frequency from the analytic signal phase; 1/|k| is the instantaneous wavelength proxy.• DH / DM — data subsets (e.g., halves or disjoint selections) compared for locking.• τ = log(1+z) — uniformizing variable; sliding windows are defined in τ.• r — Pearson correlation between DH and DM time series of 1/|k| within a window (locking coefficient).• Δr — difference after applying a distance gate: (r_gated − r_full).• Distance gate — keep samples with |k_imf − k2|/|k2| ≥ δ (excludes points too close to the IMF2 frequency).• BH-FDR (q) — Benjamini–Hochberg false discovery rate control applied across τ-windows.

Sliding-τ Locking (IMF2)

Figure: Median locking coefficient r between DH and DM 1/|k| (IMF2) in sliding τ-windows (width=0.25, step=0.05). Orange markers indicate windows passing BH-FDR q≤0.10.

Takeaways:• r≈1 from τ≈0.55–1.15 (z≈0.75–2.2), indicating near-perfect wavelength locking over most of the range.• Local dips appear near τ≈0.47 and τ≈0.96; still significant after FDR.• Shaded band (not visible in this export) corresponds to ±IQR/2 across seeds; dispersion is small.

Locking r per Seed (Consistency Check)

Figure: Color map of locking r by window (x-axis, τ) and seed (y-axis).

Takeaways:• High r is consistent across all seeds in most windows (yellow band), supporting robustness.• The same τ regions that dip in the sliding-τ curve also show slightly lower r here, matching the summary.

Distance-Gated Δr (Specificity Test)

Figure: Change in r (Δr = gated − full) as points increasingly close to the IMF2 instantaneous frequency (k2) are excluded: δ ∈ {0.1, 0.2, 0.3, 0.4}.

Interpretation:• IMF2: Δr≈−0.042 and nearly flat vs δ ⇒ removing near-k2 samples lowers r substantially and uniformly, consistent with locking being concentrated around the BAO-scale frequency.• IMF1: Δr≈0 (insensitive to gating) ⇒ little evidence that IMF1 contributes to the locking.• IMF3: small negative trend with δ (few×10⁻³), indicating weak proximity-driven contribution compared to IMF2.

Method Invariance: Δr (IMF1)

Figure: Distance-gated Δr for IMF1 across CEEMDAN, EEMD, and EMD at (pad, ngrid)=(0.15, 384).

Key points:• All methods yield |Δr|≲2–3×10⁻³.• Signs differ slightly (EEMD > 0, EMD < 0, CEEMDAN ≈ 0), but magnitudes are tiny ⇒ method choice does not materially change IMF1 specificity.

Method Invariance: Δr (IMF3)

Figure: Distance-gated Δr for IMF3 across CEEMDAN, EEMD, and EMD at (pad, ngrid)=(0.15, 384).

Key points:• CEEMDAN: small negative Δr (≈−10⁻³).• EEMD: Δr transitions from ≈+1.3×10⁻³ to ≈−1.7×10⁻³ as δ increases.• EMD: small positive Δr that grows with δ (≈+2×10⁻³ to +4.5×10⁻³).Overall magnitudes remain O(10⁻³), much smaller than IMF2’s effect.
