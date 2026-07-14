# SNIA_BAO_HHT_Dechirp_Report_0V1

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\SNIA_BAO_HHT_Dechirp_Report_0V1.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2025-10-22T10:13:00Z`

## Extracted Text

SN Ia Residuals, HHT Locking, and BAO‑driven De‑chirp

Version 0V1Date: October 22, 2025

Executive Summary

This report details an extension to the previous HHT Locking Reports with the addition of a BAO‑informed de‑chirp U(L), quantitative tightening metrics on IMF instantaneous frequency ridges, and a window‑width (WL) sensitivity study for the locking statistic. IMF2–IMF3 show consistent, negative locking with BAO proxies that strengthens as WL grows from 0.03 to 0.05 in L = ln(1+z). De‑chirping along a mixed BAO driver (α≈0.375 weighting of DH_over_rd relative to DM_over_rd) reduces the spread (MAD) of instantaneous frequencies by ≈40–60%. Phase‑randomized surrogates continue to yield high significance (p ≈ 1×10⁻⁴) for the focus‑IMF energy. Together these results support a picture in which a cosmological BAO‑related ruler modulates the phase of Δμ(L), producing a chirped component carried by IMF2–IMF3.

Acronyms & Symbols

SN Ia: Type‑Ia supernovae.

BAO: Baryon Acoustic Oscillation.

HHT: Hilbert–Huang Transform (EMD/EEMD/CEEMDAN + Hilbert analytics).

IMF: Intrinsic Mode Function extracted by EMD‑family algorithms.

L: ln(1+z) (natural log of one plus redshift).

Δμ(L): SN distance‑modulus residual vs L.

WL: Sliding‑window half‑width in L.

DM_over_rd, DH_over_rd: Transverse and radial BAO rulers normalized by sound horizon.

U(L): De‑chirped coordinate built from a BAO proxy b(L): U(L) ∝ ∫|b(L)| dL, rescaled to [0,1].

MAD: Median Absolute Deviation (robust spread).

Methods (brief)

Data are interpolated to a uniform grid in L. CEEMDAN extracts IMFs; the analytic signal a(L)·e^{iφ(L)} yields instantaneous amplitude a(L) and phase φ(L). Instantaneous frequency is f(L) = (1/2π) dφ/dL. For a BAO proxy b(L), we compute a local z‑scored correlation between f(L) and b(L) within windows of width 2·WL; the lock score is the mean of these local correlations over the BAO‑supported L‑range. Energy significance is assessed with phase‑randomized surrogates preserving the amplitude spectrum of Δμ(L).

Key Equations (text form)

Instantaneous frequency: f(L) = (1/2π) · dφ/dL.

Lock (local) at center L₀:  ρ(L₀) = corr_zscore( f[L₀±WL], b[L₀±WL] );  Lock = ⟨ρ(L₀)⟩ₗ.

De‑chirp coordinate:  U(L) = [∫_{Lmin}^{L} |b(ℓ)| dℓ] / [∫_{Lmin}^{Lmax} |b(ℓ)| dℓ],   U ∈ [0,1].

Tightening metric:  R = MAD(f_U) / MAD(f_L);  lower is better.

Code & Scripts Used

• imf23_probe.py       — CEEMDAN/HHT, lock metrics, CSVs (imf_inst_freq*.csv) and figures.

• helpers_signal.py       — df/dL spline bootstrap, U(L) builder, and uniform‑U resampling.

• report_plots.py       — convenience plotting script for report figures.

• Spyder/Notebook cells   — pipelines for U‑space IMFs and α‑mix sweeps.

Observed script hashes from runs: d74ea8e…, 5dd39def…, fc63a7f….

Results

1) IMF Energy vs L and Standing‑wave Probe

Figure: IMF energy density vs L = ln(1+z). Peaks near L≈0.8–0.9 indicate concentrated mode energy.

Figure: FFT of Δμ(L) on uniform L grid. Low‑frequency dominance consistent with slow modulation.

2) Instantaneous Frequency vs BAO Proxy

Figure: IMF1 f(L) vs BAO proxy (scatter). Broad, weak structure; not primary science band.

Figure: IMF2 f(L) vs BAO proxy. Negative trend indicates anti‑correlation (higher BAO → lower f).

Figure: IMF3 f(L) vs BAO proxy. Clearer negative trend with localized features.

Figure: IMF4 f(L) vs BAO proxy. Sign behavior differs from IMF2–3, perhaps secondary mechanism.

3) De‑chirp: Frequency Distribution Tightening

Figure 5 — IMF2 IF histograms before/after de‑chirp (U built with α=0.375). MAD_U/MAD_L≈0.57 shows clear tightening in U‑space (orange) vs. L‑space (blue).

Figure 6 — IMF3 IF histograms before/after de‑chirp (α=0.375). Dispersion halves relative to L‑space: MAD_U/MAD_L≈0.44.

Figure: Tightening across IMFs at α=0.375; IMF4 tightens most, but note sign differences.Several modes exhibit >50% tightening (IMF4 ≈83%, IMF5 ≈68%).

Numerical Highlights

Interpretation: The U(L) reparameterization informed by BAO distances removes a significant fraction of frequency drift in the SN residual IMFs, especially for IMF2–3, without invoking a strong linear correspondence between IF and the BAO proxy itself. This is consistent with a geometric (phase‑relabeling) explanation.

4) α‑mix Sweep (0→DM_over_rd, 1→DH_over_rd)

Figure: De‑chirp efficacy vs α for multiple IMFs. IMF2–IMF3 prefer α≈0.375.IMF2–3 attain minima near α≈0.375 (i.e., ~62.5% DM/rd + ~37.5% DH/rd).

5) Lock vs Window Width (WL) in L

Figure: Absolute lock strength |score| vs WL. IMF2–IMF3 grow monotonically from WL=0.03→0.05.

Figure: Signed lock score vs WL. IMF2–IMF3 stay negative; IMF4 flips sign around WL≈0.04.

Table: WL sweep summary (scores; p≈1e−4 for all).

Relation to Previous HHT Locking Report

The prior report established robust IMF2–IMF3 energy and BAO locking. This 0V1 update adds explicit chirp diagnostics (df/dL bootstrap), U‑space de‑chirping, α‑mix preference near 0.375, and WL sensitivity. All results strengthen the original interpretation.

Sequence of Tests to Date (incl. non‑conclusive / void controls)

• HHT on Δμ(L) via CEEMDAN; IMF energy + surrogate energy significance (p≈1e−4).

• Locking vs BAO proxies (DM_over_rd, DH_over_rd); IMF2–3 negative sign.

• Chirp slope df/dL via spline‑bootstrap; sign‑consistent fraction ≳0.8 and mean negative slope.

• De‑chirp U(L) using BAO proxies; MAD_U/MAD_L tightening for IMF2–3 ≈ 0.4–0.6.

• α‑mix sweep; best α≈0.375 for IMF2–3.

• WL sweep of lock strength; monotonic increase for IMF2–3; IMF4 sign flip near WL≈0.04.

• Non‑conclusive/void: beyond phase‑randomized surrogates, dedicated random/flat BAO drivers not yet fully executed; preliminary trials did not tighten.

Conclusion

The BAO‑informed de‑chirp U(L) substantially tightens IMF2–3 instantaneous‑frequency distributions of SN Ia residuals. The α≈0.375 optimum suggests a specific linear combination of transverse and radial BAO proxies that best linearizes the ridge. Consistently negative locking scores (stronger with larger WL) and p≈10⁻⁴ surrogate significance indicate the effect is unlikely to arise from chance.

Conjecture

IMF2–3 in Δμ(L) likely encode a slowly varying, BAO‑scale‑modulated chirp. Constructing U from bα partially removes this modulation, revealing narrower, more stationary frequency content. If confirmed across independent SN/BAO catalogs, this would be consistent with a genuine SN–LSS coupling at BAO‑effective scales rather than a processing artifact.

## Extracted Tables

### Table 1

Quantity | Value | Notes
Best α for IMF2–3 | α≈0.375 | α=0→DM/rd, α=1→DH/rd; mixture favors DM.
IMF2 MAD ratio | ≈0.54–0.57 | From histograms and sweep.
IMF3 MAD ratio | ≈0.44 | Strongest de‑chirp among mid‑order IMFs.
Energy‑weighted IMF2+3 | ≈0.51 | Mean tightening for the science‑relevant band.
Locking scores (example) | small/negative | Sliding‑window correlation metric; geometry effect dominates.

### Table 2

WL | stat_obs | p_value | IMF1 | IMF2 | IMF3 | IMF4
0.030000 | 48.969698 | 0.000100 | -0.050913 | -0.068552 | -0.048556 | 0.061210
0.035000 | 48.107214 | 0.000100 | -0.044921 | -0.071515 | -0.066169 | 0.076311
0.040000 | 50.622517 | 0.000100 | -0.039372 | -0.083078 | -0.090185 | -0.072064
0.045000 | 50.787299 | 0.000100 | -0.045422 | -0.105136 | -0.110862 | -0.065604
0.050000 | 50.438822 | 0.000100 | -0.053482 | -0.114618 | -0.129916 | -0.025598
