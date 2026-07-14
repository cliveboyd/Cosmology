# SNIA_BAO_HHT_Dechirp_Report_0V2

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\SNIA_BAO_HHT_Dechirp_Report_0V2.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2025-10-22T08:10:00Z`

## Extracted Text

SN Ia Residuals, HHT Locking, and BAO‑driven De‑chirp

Version 0V2  •  Date: October 22, 2025

Executive Summary

This update removes redundant figures, adds clear per‑figure descriptions, and keeps the prior layout.

Key findings are unchanged: IMF2–IMF3 instantaneous‑frequency ridges show negative locking to BAO proxies in L = ln(1+z).

A BAO‑informed de‑chirp U(L), built from a mix of transverse and radial BAO rulers (best α≈0.375), substantially narrows the IMF2–IMF3 frequency distributions (≈40–60% tightening).

Phase‑randomized surrogates maintain high significance for the focus‑IMF energy (p≈1×10⁻⁴).

Acronyms & Symbols

SN Ia: Type‑Ia supernovae.BAO: Baryon Acoustic Oscillation.HHT: Hilbert–Huang Transform (EMD/EEMD/CEEMDAN + Hilbert analytics).IMF: Intrinsic Mode Function.L: ln(1+z), natural log of one plus redshift.Δμ(L): SN distance‑modulus residual vs L.WL: Sliding‑window half‑width in L used for locking.DM_over_rd, DH_over_rd: Transverse/radial BAO rulers normalized by the sound horizon.U(L): De‑chirped coordinate from BAO proxy b(L): U(L) ∝ ∫|b(L)| dL, scaled to [0,1].MAD: Median Absolute Deviation.

Key Relations

Instantaneous frequency:  f(L) = (1 / 2π) · dφ/dL.

Local lock:  ρ(L₀) = corr_zscore( f[L₀±WL], b[L₀±WL] ), then Lock = ⟨ρ(L₀)⟩ over BAO‑supported L.

De‑chirp coordinate:  U(L) = [∫_{Lmin}^{L} |b(ℓ)| dℓ] / [∫_{Lmin}^{Lmax} |b(ℓ)| dℓ],  U∈[0,1].

Tightening:  R = MAD(f_U) / MAD(f_L),  with percent tightening = 100×(1−R).

Scripts Used

• imf23_probe.py — HHT/CEEMDAN runs, locking metrics, CSV dumps (imf_inst_freq*.csv).

• helpers_signal.py — df/dL bootstrap, U(L) construction, uniform‑U resampling.

• report_plots.py — histograms, α‑sweep, tightening plots.

• regen_freq_vs_bao_plots.py — refreshed IMF‑vs‑BAO scatterplots with distinct coloring.

Figures

Figure 1. IMF energy density vs L = ln(1+z)

Scatter of summed IMF energy (focus band) versus L.

The color scale maps the provided meta field (if present), allowing visual association between energy concentration and any external diagnostic.

Local peaks mark L regions where the HHT extracts higher oscillatory power from Δμ(L).

Figure 2. Standing‑wave probe on Δμ(L)

Power spectrum of Δμ(L) on a uniform L grid.

The dominance at low Fourier modes indicates slow modulation consistent with a chirp-like component; sharper features would point to quasi‑standing waves.

Figure 3. IMF1 instantaneous frequency vs BAO proxy

Scatter of IMF1 instantaneous frequency f(L) against the BAO proxy evaluated on the same L grid.

Each point represents a co‑located (b(L), f(L)).

The color field (if provided) is the chosen meta quantity.

IMF1 typically captures high‑frequency/edge content; weak trends here are expected.

Figure 4. IMF2 instantaneous frequency vs BAO proxy

Scatter of IMF2 instantaneous frequency against the BAO proxy.

A negative slope indicates that larger BAO scale correlates with lower instantaneous frequency, i.e., anti‑locking consistent with a BAO‑modulated chirp.

Figure 5. IMF3 instantaneous frequency vs BAO proxy

As for IMF2 but one mode lower. IMF3 exhibits a clearer negative association in this dataset, strengthening the interpretation that IMF2–IMF3 carry the BAO‑related phase modulation.

Figure 5xx. IMF4 instantaneous frequency vs BAO proxy

Figure 6. De‑chirp tightening for IMF2

Histogram of IMF2 instantaneous frequency before (L‑space) and after BAO‑informed de‑chirp (U‑space).

The U‑space distribution narrows substantially; the annotated MAD ratio quantifies the tightening.

Figure 7. De‑chirp tightening for IMF3

As for IMF2 but for IMF3. IMF3 shows an even stronger narrowing, indicating that the U(L) transform removes a significant fraction of the BAO‑driven non‑stationarity.

Figure 8. α‑mix sweep for de‑chirp driver

De‑chirp efficacy (MAD_U/MAD_L) versus α in b_α = α·DH_over_rd + (1–α)·DM_over_rd.

Lower is better.

IMF2–IMF3 prefer α≈0.375, suggesting an optimal linear combination of the two BAO proxies.

Figure 9. Percent tightening across IMFs (α ≈ 0.375)

Bar chart of tightening = 100×(1–MAD_U/MAD_L) across IMFs using the preferred α.

IMF2–IMF3 tighten by ~40–60%, with other modes showing varying behavior.

Relation to Prior HHT Locking Report

This 0V2 version streamlines the figure set (duplicates removed) and augments captions with what each panel demonstrates regarding SN–BAO coupling and de‑chirp efficacy.

The statistical inferences and qualitative trends match the earlier report.

Conclusion

IMF2–IMF3 show robust, negative locking to BAO proxies in L, and BAO‑informed U(L) de‑chirp substantially stabilizes their ridge frequencies.

The evidence points toward a BAO‑modulated chirp in SN Ia residuals.

Conjecture

An effective phase driver is well‑approximated by a linear mix of transverse and radial BAO rulers with α≈0.35–0.40.

Independent catalogs and pipelines should reproduce the same α‑preference if the coupling is physical.
