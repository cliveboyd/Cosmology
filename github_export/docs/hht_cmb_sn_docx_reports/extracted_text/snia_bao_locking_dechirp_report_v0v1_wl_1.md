# SNIA_BAO_Locking_Dechirp_Report_v0V1_WL (1)

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\SNIA_BAO_Locking_Dechirp_Report_v0V1_WL (1).docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2013-12-23T23:15:00Z`

## Extracted Text

SN Ia – BAO Locking & De‑chirp Analysis (v0V1)

Executive summary

We analyze Type Ia supernova (SN Ia) Hubble‑diagram residuals Δμ on a uniform L=log(1+z) grid with the Hilbert–Huang Transform (HHT; CEEMDAN). By constructing a de‑chirped coordinate U(L) from a BAO proxy bα = (1−α)·(DM/rd) + α·(DH/rd), we find that α ≈ 0.375 minimizes the scatter (MAD) of IMF2–3 instantaneous frequencies when re‑analyzed in U. For this configuration, IMF2 and IMF3 show MAD_U/MAD_L ≈ 0.57 and 0.44, respectively, consistent with the bar‑chart tightening across IMFs. Surrogate tests yield p ≈ 10⁻⁴ for the focus‑IMF energy statistic. A window‑width (WL) sweep of the locking metric shows increasingly negative scores for IMF2–3 with larger WL; by a simple |IMF2,IMF3| figure‑of‑merit, WL ≈ 0.05 performs best.

Acronyms & symbols

HHT: Hilbert–Huang Transform.  CEEMDAN: Complete Ensemble Empirical Mode Decomposition with Adaptive Noise.  IMF: Intrinsic Mode Function.  L: log(1+z).  U(L): de‑chirped coordinate derived from a BAO proxy.  BAO: Baryon Acoustic Oscillation.  DM/rd, DH/rd: transverse and radial BAO‑scale proxies normalized by the sound horizon rd.  MAD: median absolute deviation.  WL: sliding window width in L used in the locking score.

Key definitions

Instantaneous frequency (per sample) for IMF i:  f_i(L) = (1/2π) · d/dL arg(Hilbert(IMF_i(L))).

BAO mix:  b_α(L) = (1−α)·(DM/rd)(L) + α·(DH/rd)(L).

De‑chirp map:  U(L) = CDF(|b_α(L)|) = (∫|b_α(L′)| dL′) / (∫|b_α(L′)| dL′)_range  ∈ [0,1].

Tightening per IMF:  100×(1 − MAD_U/MAD_L)  [%].

Core figures

Fig. 1 — IMF2 instantaneous frequency histograms in L (blue) vs U (orange; α=0.375). The U‑space distribution narrows; reported ratio MAD_U/MAD_L ≈ 0.569.

Fig. 2 — IMF3 instantaneous frequency histograms in L (blue) vs U (orange; α=0.375). Pronounced tightening; MAD_U/MAD_L ≈ 0.439.

Fig. 3 — α sweep of MAD_U/MAD_L. IMF2–3 achieve their minimum near α≈0.375 (vertical line).

Fig. 4 — Percent tightening across IMFs at α=0.375. IMF4 shows the largest tightening, with IMF2–3 also significantly reduced.

Fig. 5 — IMF energy density (sum over focus IMFs) vs L. Peaks indicate localized energy concentrations in log(1+z).

Fig. 6 — Standing‑wave probe on Δμ(L). The near‑zero‑frequency dominance is consistent with a chirp rather than a strong stationary line.

Fig. 7 — IMF1 instantaneous frequency vs BAO proxy on the L grid.

Fig. 8 — IMF2 instantaneous frequency vs BAO proxy on the L grid.

Fig. 9 — IMF3 instantaneous frequency vs BAO proxy on the L grid.

WL sweep for locking metric

Lock scores (signed correlation‑like statistic) become more negative for IMF2–3 with increasing WL, supporting a robust anti‑correlation between IMF frequency and the BAO proxy when averaged over a broader L range. By a simple |IMF2,IMF3| figure‑of‑merit, WL≈0.050 is optimal in this sweep. All runs show p≈10⁻⁴ under phase‑randomized surrogates.

Fig. 10 — Lock *strength* (|score|) vs WL for IMFs 1–4.

Fig. 11 — Signed lock score vs WL for IMFs 1–4.

Methods & scripts used

• imf23_probe.py (v0V1): loads SN/BAO, runs CEEMDAN, computes Hilbert analytics, lock metrics, surrogates, optional de‑chirp in U.

• report_plots.py: builds publication figures (histograms L vs U, α‑sweep, tightening bar).

• helpers_signal.py: interpolation utilities, de‑chirp mapping helpers, and robust MAD/bootstrapping utilities used earlier.

• Notebooks/console harness: parameter sweeps (α, WL), dataset selection, and summary CSV generation.

Conclusion

The BAO‑informed de‑chirp U(L) substantially tightens IMF2–3 instantaneous‑frequency distributions of SN Ia residuals. The α≈0.375 optimum suggests a specific linear combination of transverse and radial BAO proxies that best linearizes the ridge. Consistently negative locking scores (stronger with larger WL) and p≈10⁻⁴ surrogate significance indicate the effect is unlikely to arise from chance.

Conjecture

IMF2–3 in Δμ(L) likely encode a slowly varying, BAO‑scale‑modulated chirp. Constructing U from bα partially removes this modulation, revealing narrower, more stationary frequency content. If confirmed across independent SN/BAO catalogs, this would be consistent with a genuine SN–LSS coupling at BAO‑effective scales rather than a processing artifact.

## Extracted Tables

### Table 1

WL | stat_obs | p_value | IMF1 | IMF2 | IMF3 | IMF4
0.030 | 48.970 | 0.000 | -0.050913 | -0.068552 | -0.048556 | 0.061210
0.035 | 48.107 | 0.000 | -0.044921 | -0.071515 | -0.066169 | 0.076311
0.040 | 50.623 | 0.000 | -0.039372 | -0.083078 | -0.090185 | -0.072064
0.045 | 50.787 | 0.000 | -0.045422 | -0.105 | -0.111 | -0.065604
0.050 | 50.439 | 0.000 | -0.053482 | -0.115 | -0.130 | -0.025598
