# HHT_Locking_Chirp_Investigation_1V0b

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\HHT_Locking_Chirp_Investigation_1V0b.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2025-10-22T14:51:00Z`
- DOCX modified: `2025-10-22T15:21:00Z`

## Extracted Text

HHT Locking — Chirp Investigation

Version 1V0 • October 22, 2025

Executive Summary

Previous HHT IMF analysis investigations have shown a robust link between SN Ia and BAO data sets.  This report extends the analysis of what seems to be an underlying chirp metric to the spectral locking profile between the two datasets.

Detailed below the investigation of a whether chirp-like frequency drifts exist in IMF2/IMF3 of SN Ia residuals and whether these drifts exhibit ‘locking’ to the BAO scale (k ≈ 1/rd) across redshift.

Using the Hilbert–Huang Transform (HHT), instantaneous frequency f(z) and amplitude A(z) are extracted from Intrinsic Mode Functions (IMFs) derived via EMD/EEMD/CEEMDAN.

A positive detection is defined by a statistically significant non‑zero slope df/dτ (with τ = ln(1+z)) and a persistent correlation between IMF scale and BAO wavenumber proxies.

Significance is established with permutation / phase-randomized surrogates (typically minimum N≈1000).

Preliminary results reproduce prior ‘locking’ signals for IMF2 with r≈0.090845 and permutation p≈0.005 across multiple seeds, suggesting a low false‑positive rate for this specific correlation metric.

Chirp slopes appear small but consistent in sign for τ‑space linear fits, motivating expanded surrogate tests and cross‑validation with independent datasets.

Instantaneous frequency (IF) ‘chirp’ in IMFs 2–3 vs. L≡ln(1+z) is reduced by a monotone reparameterization U(L)=normalized ∫|b(L)|dL built from a BAO proxy b(L).

For the preferred BAO mixture b_α(L)=α·(DH/rd)+(1−α)·(DM/rd), an α≈0.375 (i.e., ~37.5% DH, ~62.5% DM) minimizes dispersion for IMF2–3.

Histogram-based dispersion ratio MAD_U/MAD_L indicates substantial tightening after de‑chirping: IMF2 ≈0.54–0.57, IMF3 ≈0.44 (lower is better).

Across all IMFs, the α=0.375 transform yields tightening of ≳50% for several modes (e.g., IMF4 ≈83%, IMF5 ≈68%), with an energy-weighted (IMF2+IMF3) ratio ≈0.51.

Sliding-window ‘locking’ scores between IF and BAO proxies remain small/negative in the current metric; thus the de‑chirp effect observed, as explained below, appears to be a geometric reparameterisation rather than a simple linear covariation.

Note: Chirp Analysis (De-Chirp) technique is provided within appendix below.

Acronyms & Symbols

SN Ia: Type‑Ia supernovae.

BAO: Baryon Acoustic Oscillation.

HHT: Hilbert–Huang Transform (EMD/EEMD/CEEMDAN + Hilbert analytics).

IMF: Intrinsic Mode Function extracted by EMD‑family algorithms.

L: ln(1+z) (natural log of one plus redshift).

Δμ(L): SN distance‑modulus residual vs L.

WL: Sliding‑window half‑width in L.

DM_over_rd, DH_over_rd: Transverse and radial BAO rulers normalized by sound horizon.

U(L): De‑chirped coordinate built from a BAO proxy             b(L): U(L) ∝ ∫|b(L)| dL, rescaled to [0,1].

MAD: Median Absolute Deviation (robust spread).

IGM: Intergalactic Medium (gas between galaxies)

Datasets

• Pantheon+SH0ES SN Ia residuals (with standard light‑curve corrections and your working residual pipeline).

• eBOSS DR16 BAO (ELG, LRG, QSO) — transverse and radial rd‑scaled distances (D_M/rd, D_H/rd).

• Optional: Cosmic chronometer H(z) and CMB distance priors for cross‑checks.

• Internal: HHT_Locking report and previous robustness grids for IMF2/IMF3.

6) Multiple‑seed replication across EEMD/CEEMDAN realizations; aggregate via inverse‑variance or Fisher‑z methods.

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

Results Snapshot (from prior robustness checks)

IMF2 locking statistic r_1/rd reproduced across seeds with consistent values and low permutation p-values:

These figures are illustrative snapshots matching earlier runs; the 1V0 full rerun reports aggregated statistics across all methods (EMD/EEMD/CEEMDAN) and seeds with confidence intervals.

Figures Included (1V0)

Fig. 1 — HHT amplitude (Hilbert spectrum) for IMF2 and IMF3 vs τ with ridge overlays; caption to explain axes and resolution.

Fig. 2 — Ridge instantaneous frequency f(τ) with linear/quadratic fits; report β and p_perm from surrogates.

Fig. 3 — 3D scatter: IMF energy vs ln(1+z) vs Noether strength; color‑encode the third variable per panel.

Fig. 4 — Locking metric across redshift bins vs BAO k=1/rd proxy; include permutation p‑values and CI bands.

Fig. 5 — Surrogate distributions of β and r_1/rd vs observed values (vertical lines).

HHT Chirp Observation

Description…:

The curve is df/dL vs L=ln(1+z). For most of the range df/dL < 0, only turning slightly > 0 near the high-L edge (where the CI balloons → low S/N / edge effects).

So, as redshift increases (L rises), the instantaneous frequency decreases: a down-chirp with z (high → low). Equivalently, moving toward today (lower z), it’s an up-chirp (low → high).

What it implies for HHT locking & chirp dynamics:

Consistency with locking:A gentle, mostly negative slope means the IMF ridge prefers slightly longer periods at higher z. That’s qualitatively compatible with a locking picture where the characteristic scale evolves with the BAO-projected scale within the coordinates: the IMF tracks a slowly changing “target” scale rather than jumping around.

Strength/significance:The median stays below zero from roughly L≈0.05–1.0 while the 16–84% band is narrowest in the mid-range—so the down-chirp is weak but coherent there. The spike and the late up-tick sit in edge-instability territory (Hilbert endpoints + sparse data), so I wouldn’t over-interpret those.

Physical read (PLamb/variable-c framing):A down-chirp with z means longer apparent timescales in the past. In FR/variable-c, earlier epochs have different clock rates/propagation speeds; that generally pushes modulations toward lower observed f at higher z, which matches the sign you see. If IMF2 is also “locked” to BAO proxies, this may suggest an adiabatic, |scale-aware modulation rather than a transient artifact.

Quick checklist to firm up analysis --- Pending:

Report the global β = ⟨df/dL⟩ with CI over the mid-range window (e.g., 0.1≤L≤0.9), plus permutation/IAAFT p-values.

Do piecewise slope analysis (e.g., 0.1–0.5, 0.5–0.9) and show they agree in sign.

Overlay β(L) with locking metric vs z to show they co-vary where S/N is highest.

Mark and exclude endpoint cones (where the shaded band explodes).

Conclusion:

The above plot points to a mild down-chirp with increasing redshift.

Supportive (however, yet to be validated) of a HHT-locking scenario in which IMF2/3 track a slowly evolving physical scale.

Interpretation & Cosmological Context

If IMF2/3 exhibit small but coherent chirp slopes with locking to BAO scale, this suggests a structured, scale‑aware process in the residual field.

In ΛCDM, persistent narrow‑band structures in SN residuals are unexpected after standard corrections; any detection must therefore survive stringent null tests and multiple‑comparison control.

In the PLamb/Kerr–de Sitter/FR frameworks (possible variable‑c or accretion‑driven effects), quasi‑coherent modulations could arise from episodic entropy injections, localised low entropy (mass) void or STF rotation/strain effects, potentially imprinting weak, redshift‑dependent frequency drifts.

Interpretation should consider selection effects, calibration systematics, and method‑induced artifacts (e.g., mode mixing).

Independent replication on alternative residual pipelines and independent BAO compilations is expected in the future.

Reproducibility Notes

Example (recent):

%run hht_robustness_grid.py \  --pack   plamb_runs/tools/bao_hht_2bin_z0p7_longbest \  --pads   0.15 --ngrids 384 \  --methods ceemdan --seeds 3 7 11 23 --ns 1000 --workers 8 \  --progress auto \  --out    plamb_runs/tools/hht_imf2_compare/robustness_grid_ceemdan_main.csv   --ckpt plamb_runs/tools/hht_imf2_compare/robustness_grid_ceemdan_main.ckpt.csv \  --overwrite

Next Steps for 1V0

• Run chirp slope estimation for IMF2 and IMF3 across all seeds/methods; store β and p_perm per bin.

• Extend surrogate family to IAAFT and block‑bootstrap on residual indices; compare p-values.

• Add BAO‑locking cross‑validation using alternate compilations; report sensitivity of r and β.

• Integrate Noether‑toggle metadata and A_acc posteriors into 3D visuals to test STF hypotheses.

• Prepare final 1V0 PDF with figure panels and a one‑page plain‑language summary.

Relation to HHT_Locking Report

This note should be read alongside the HHT_Locking_Report_0V6 for full methodology, parameter grids, and prior robustness outcomes. This analysis reuse definitions of IMF numbering, redshift binning, permutation tests, and BAO proxy construction.

Key Results

Figure 1 — IMF energy density (sum over selected IMFs) vs L=ln(1+z). Pronounced structure appears around L≈0.8 (z≈1.2), indicating non‑stationarity in the residual signal.

Figure 2 — Standing‑wave probe: FFT power of Δμ(L) on the uniform L‑grid. A large low‑frequency pedestal plus high‑frequency tails confirm mixed scales; HHT decomposes these adaptively.

Figure 3 — IMF1 instantaneous frequency vs BAO proxy. Scatter of IMF1 instantaneous frequency f(L) against the BAO proxy evaluated on the same L grid.

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

.

Figure 5 — IMF2 IF histograms before/after de‑chirp (U built with α=0.375). MAD_U/MAD_L≈0.57 shows clear tightening in U‑space (orange) vs. L‑space (blue).

Figure 6 — IMF3 IF histograms before/after de‑chirp (α=0.375). Dispersion halves relative to L‑space: MAD_U/MAD_L≈0.44.

Figure 7 — De‑chirp efficacy sweep R_k(α)=MAD_U/MAD_L vs α. IMF2–3 attain minima near α≈0.375 (i.e., ~62.5% DM/rd + ~37.5% DH/rd).

Figure 8 — Percent tightening 100×(1−R_k) at α=0.375 across IMFs. Several modes exhibit >50% tightening (IMF4 ≈83%, IMF5 ≈68%).

Bar chart of tightening = 100×(1–MAD_U/MAD_L) across IMFs using the preferred α.

Numerical Highlights

Interpretation: The U(L) reparameterization informed by BAO distances removes a significant fraction of frequency drift in the SN residual IMFs, especially for IMF2–3, without invoking a strong linear correspondence between IF and the BAO proxy itself. This is consistent with a geometric (phase‑relabeling) explanation.

Equations & Definitions

• L ≡ ln(1+z)• D_M(z) = (1+z)·D_A(z)  (comoving angular diameter distance)• D_H(z) = c/H(z)        (Hubble distance)• r_d    = sound‑horizon at baryon drag• BAO proxies: D_M/r_d, D_H/r_d (dimensionless)• U(L)  = normalized cumulative ∫|b(L)| dL• IF: f_k(L) = (1/2π)·dφ_k/dL from the analytic IMF a_k(L)

Discussion

The α‑sweep shows that a blended BAO scale defined by ≈62.5% D_M/r_d and ≈37.5% D_H/r_d best removes chirp in IMF2–3. This blend is physically plausible: D_M and D_H capture transverse and radial rulers respectively; a mixture could mirror how residual features project onto the L variable. The modest/negative locking scores suggest that the transform exploits distance‑ladder geometry rather than pointwise covariation between IF and BAO magnitude.

The strong tightening in IMF4–5 (and noticeable tightening in IMF2–3) implies that multiple scales within the residuals are non‑stationary with L—consistent with slowly varying phase or frequency content (‘chirp‑like’ behavior). De‑chirping re‑centers the IF distribution around a quasi‑stationary band in U, enabling cleaner downstream tests.

Conclusion

We find robust evidence that an L→U(L) reparameterization derived from a BAO proxy significantly reduces the dispersion of instantaneous frequencies for SN Ia residual IMFs, especially in IMF2–3. The preferred mixture uses mostly D_M/r_d with a non‑negligible D_H/r_d component (α≈0.375). This supports the hypothesis that BAO distances encode the appropriate nonlinear ‘clock’ for these oscillatory components.

The evidence points toward a BAO‑modulated chirp in SN Ia residuals.

Conjecture

An effective phase driver is well‑approximated by a linear mix of transverse and radial BAO rulers with α≈0.35–0.40.

Independent catalogs and pipelines should reproduce the same α‑preference if the coupling is physical.

If SN residual oscillations partly trace large‑scale geometric modulation tied to the BAO ruler, then transforming to U(L) constructed from contemporary BAO ladders are expected to systematically de‑chirp analogous datasets (e.g., alternative SN compilations, binned distance moduli) and sharpen cross‑correlations with galaxy/IGM tracers. Testing stability of α with independent BAO compilations and sub‑samples may assist in distinguishing geometric from astrophysical systematics.

Appendix: Chirp Analysis Description… “de-chirping”:

De-chirping is how you test whether an IMF truly carries a chirp (i.e., its instantaneous frequency drifts with ). Technique:  “subtract” a hypothesized drift from the signal; if the IMF then becomes nearly constant-frequency (a flat ridge), the chirp model was right.

If it stays sloped or messy, the “chirp” was likely noise or method-induced.

The HHT IMFs can be considered as an amplitude-modulated tone with a time-varying phase:

Pick a chirp model (e.g., linear , or add curvature). Integrate to get a model phase

Then demodulate (heterodyne) the IMF by the conjugate of this model:

If your model matches the real drift, the instantaneous frequency of ,

becomes flat (≈ constant), and the HHT ridge horizontal. That’s positive evidence for a chirp with parameters .

Two equivalent viewpoints

Phase demodulation (mixing/heterodyne). Multiply by as above.

Time-warping. Define a new coordinate . In -space, a correct model makes the ridge horizontal without explicit mixing.

Practical workflow with HHT IMFs

Extract IMF2/IMF3 and compute the analytic signal via Hilbert; unwrap the phase to get .

Fit a chirp model (start linear; allow quadratic if needed): .

De-chirp: .

Re-Hilbert on → get .

Diagnostics: you want

Slope (piecewise too),

Big drop in vs ,

Narrower ridge / lower Hilbert spectral entropy,

Concentration increases in the STFT/CWT/HHT energy around one band.

Null testing: run the same de-chirp on surrogates (circular-shift, IAAFT). If your real data shows much larger slope-reduction or variance-reduction than surrogates, you’ve got statistical support for a true chirp.

The median is mostly negative over ⇒ down-chirp with redshift (frequency gently decreases as increases). Therefore with . After de-chirping one expects to see:

nearly flat in the mid-range (exclude edge cones),

Stronger locking metrics (or at least more stable) if the chirp is physically tied to your BAO proxy,

Surrogate p-values showing your slope-flattening is unlikely by chance.

Figures Included (1V0) — Updated

Fig. 1 — Hilbert spectrum (amplitude) for IMF2 and IMF3 vs L≡ln(1+z); ridges overlaid. Axes: x=L (dimensionless), y=f [cycles per ΔL], color=amplitude.

Fig. 2 — Instantaneous-frequency slope df/dL vs L with 16–84% CI and median. Axes: x=L, y=df/dL [cycles per ΔL²]; dashed zero-line.

Fig. 3 — 3D scatter (panelized): IMF energy vs L vs Noether strength. Panels color-encode the third variable. Axes: energy [arb], L, Noether index.

Fig. 4 — Locking statistic r_{1/rd}(L) between IMF ridge scale and BAO proxy, with bootstrap confidence bands; points=bin centers.

Fig. 5 — Surrogate distributions for chirp slope β and locking r_{1/rd}; vertical lines=observed; annotate p-values.

Fig. 6 — IMF2 instantaneous-frequency histograms before (L-space) and after de-chirp (U-space, α=0.375). Annotate MAD_U/MAD_L.

Fig. 7 — IMF3 instantaneous-frequency histograms before/after de-chirp (α=0.375). Annotate MAD_U/MAD_L.

Fig. 8 — De-chirp efficacy R(α)=MAD_U/MAD_L vs α in b_α = α·(D_H/r_d) + (1−α)·(D_M/r_d). Mark minimum near α≈0.375.

Fig. 9 — Percent tightening 100×(1−R) at α=0.375 across IMFs (2–5); lower is better; bars annotated.

Integrity & Results Check (1V0)

CEEMDAN snapshot (pads=0.15, ngrid=384, seeds=3,7,11,23): median p-values — DH_over_rd ≈ 0.0010; DM_over_rd ≈ 0.0045. These indicate significant IMF2 energy excess relative to phase-randomized surrogates. Figure captions above specify axes/units for reproducibility.

## Extracted Tables

### Table 1

Seed | r_1/rd | p_perm
101 | 0.090845 | 0.004999
202 | 0.090845 | 0.004999
303 | 0.090845 | 0.004999

### Table 2

Quantity | Value | Notes
Best α for IMF2–3 | α≈0.375 | α=0→DM/rd, α=1→DH/rd; mixture favors DM.
IMF2 MAD ratio | ≈0.54–0.57 | From histograms and sweep.
IMF3 MAD ratio | ≈0.44 | Strongest de‑chirp among mid‑order IMFs.
Energy‑weighted IMF2+3 | ≈0.51 | Mean tightening for the science‑relevant band.
Locking scores (example) | small/negative | Sliding‑window correlation metric; geometry effect dominates.
