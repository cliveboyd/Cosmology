# SNIA_HHT_BAO_dechirp_report

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\SNIA_HHT_BAO_dechirp_report.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2025-10-22T09:45:00Z`

## Extracted Text

SN Ia Residuals, HHT IMFs,

and BAO-Driven De‑chirping

Report timestamp: 2025-10-21T15:07:53.952723Z

This report summarizes the current state of the Hilbert–Huang Transform (HHT) analysis of Type Ia supernova (SN Ia) Hubble-residuals, their apparent ‘locking’ to BAO distance scales, and a new de‑chirping transform U(L) constructed from a BAO proxy.

Executive Summary

Instantaneous frequency (IF) ‘chirp’ in IMFs 2–3 vs. L≡ln(1+z) is reduced by a monotone reparameterization U(L)=normalized ∫|b(L)|dL built from a BAO proxy b(L).

For the preferred BAO mixture b_α(L)=α·(DH/rd)+(1−α)·(DM/rd), an α≈0.375 (i.e., ~37.5% DH, ~62.5% DM) minimizes dispersion for IMF2–3.

Histogram-based dispersion ratio MAD_U/MAD_L indicates substantial tightening after de‑chirping: IMF2 ≈0.54–0.57, IMF3 ≈0.44 (lower is better).

Across all IMFs, the α=0.375 transform yields tightening of ≳50% for several modes (e.g., IMF4 ≈83%, IMF5 ≈68%), with an energy-weighted (IMF2+IMF3) ratio ≈0.51.

Sliding-window ‘locking’ scores between IF and BAO proxies remain small/negative in the current metric; thus the de‑chirp effect appears to be a geometric reparameterization rather than a simple linear covariation.

Note: Chirp Analysis De-Chirp technique is provided within appendix below.

Acronyms

BAO — Baryon Acoustic OscillationCEEMDAN — Complete Ensemble Empirical Mode Decomposition with Adaptive NoiseEMD/EEMD — (Ensemble) Empirical Mode DecompositionHHT — Hilbert–Huang TransformIF — Instantaneous FrequencyIMF — Intrinsic Mode FunctionMAD — Median Absolute DeviationSN Ia — Type Ia SupernovaU(L) — BAO‑driven monotone reparameterization of L=ln(1+z)

Data & Preprocessing

• SN sample: Pantheon+SH0ES residuals file (columns zCMB, MU_SH0ES, MU_SH0ES_ERR_DIAG, DMU).

• BAO proxies: dimensionless distances DM/rd and DH/rd compiled into a single redshift ladder.

• Independent variable: L≡ln(1+z). Data are interpolated to a uniform L‑grid (nbins=2000).

• HHT: CEEMDAN (trials=100, noise_width=0.2, seed=11), extracting up to 6 IMFs and the trend.

Methods

Hilbert–Huang Transform (HHT).

Given the residual series y(L), we compute Intrinsic Mode Functions (IMFs) {c_k(L)} via CEEMDAN. The analytic signal a_k(L)=c_k(L)+i·H{c_k(L)} provides instantaneous amplitude A_k(L)=|a_k(L)| and phase φ_k(L)=arg a_k(L); the instantaneous frequency (IF) is f_k(L)= (1/2π)·dφ_k/dL (cycles per unit ΔL).

BAO‑driven de‑chirp U(L).

Let b(L) denote a BAO proxy on the L‑grid (e.g., DM/rd or DH/rd, or a mixture b_α(L)=α·(DH/rd)+(1−α)·(DM/rd)). We define a monotone reparameterization by the cumulative magnitude of b(L):    U(L) = [∫^{L}_{L_min} |b(ℓ)| dℓ − min] / range,    with U∈[0,1].This compresses/expands regions of L according to the BAO geometry and tends to ‘straighten’ chirps in f_k(L).

Dispersion metric.

For each IMF k we compare the dispersion of IF in L‑space vs. U‑space using the robust median absolute deviation (MAD):    MAD(x) = 1.4826 · median_i |x_i − median(x)|,    R_k(α) = MAD(f_k^U) / MAD(f_k^L).Values R_k<1 indicate that the IF distribution narrows after the U‑transform, consistent with de‑chirping.

Locking score (qualitative).

We also use a sliding-window, z‑restricted correlation between standardized IF and standardized b(L); low/negative scores suggest weak linear covariation despite the de‑chirp benefit.

HHT Chirp Observation

Description…:

The curve is df/dL vs L=ln(1+z). For most of the range df/dL < 0, only turning slightly > 0 near the high-L edge (where the CI balloons → low S/N / edge effects).

So, as redshift increases (L rises), the instantaneous frequency decreases: a down-chirp with z (high → low). Equivalently, moving toward today (lower z), it’s an up-chirp (low → high).

What it implies for HHT locking & chirp dynamics:

Consistency with locking:A gentle, mostly negative slope means the IMF ridge prefers slightly longer periods at higher z. That’s qualitatively compatible with a locking picture where the characteristic scale evolves with the BAO-projected scale within the coordinates: the IMF tracks a slowly changing “target” scale rather than jumping around.

Strength/significance:The median stays below zero from roughly L≈0.05–1.0 while the 16–84% band is narrowest in the mid-range—so the down-chirp is weak but coherent there. The spike and the late up-tick sit in edge-instability territory (Hilbert endpoints + sparse data), so I wouldn’t over-interpret those.

Physical read (PLamb/variable-c framing):A down-chirp with z means longer apparent timescales in the past. In FR/variable-c, earlier epochs have different clock rates/propagation speeds; that generally pushes modulations toward lower observed f at higher z, which matches the sign you see. If IMF2 is also “locked” to BAO proxies, this may suggest an adiabatic, scale-aware modulation rather than a transient artifact.

Quick checklist to firm up analysis --- Pending ToDo:

Report the global β = ⟨df/dL⟩ with CI over the mid-range window (e.g., 0.1≤L≤0.9), plus permutation/IAAFT p-values.

Do piecewise slope analysis (e.g., 0.1–0.5, 0.5–0.9) and show they agree in sign.

Overlay β(L) with locking metric vs z to show they co-vary where S/N is highest.

Mark and exclude endpoint cones (where the shaded band explodes).

Conclusion:

The above plot points to a mild down-chirp with increasing redshift.

Supportive (but not yet to be validated) of an HHT-locking scenario in which IMF2/3 track a slowly evolving physical scale.

Key Results

Figure 1 — IMF energy density (sum over selected IMFs) vs L=ln(1+z). Pronounced structure appears around L≈0.8 (z≈1.2), indicating non‑stationarity in the residual signal.

Figure 2 — Standing‑wave probe: FFT power of Δμ(L) on the uniform L‑grid. A large low‑frequency pedestal plus high‑frequency tails confirm mixed scales; HHT decomposes these adaptively.

Figure 3 — IMF1 instantaneous frequency vs BAO proxy. IF scatters with mild trends across the BAO ladder; note the broadening at large proxy values.

Figure 4 — IMF3 instantaneous frequency vs BAO proxy. Higher‑order IMF displays narrower IF with localized excursions, consistent with chirp segments.

Figure 5 — IMF2 IF histograms before/after de‑chirp (U built with α=0.375). MAD_U/MAD_L≈0.57 shows clear tightening in U‑space (orange) vs. L‑space (blue).

Figure 6 — IMF3 IF histograms before/after de‑chirp (α=0.375). Dispersion halves relative to L‑space: MAD_U/MAD_L≈0.44.

Figure 7 — De‑chirp efficacy sweep R_k(α)=MAD_U/MAD_L vs α. IMF2–3 attain minima near α≈0.375 (i.e., ~62.5% DM/rd + ~37.5% DH/rd).

Figure 8 — Percent tightening 100×(1−R_k) at α=0.375 across IMFs. Several modes exhibit >50% tightening (IMF4 ≈83%, IMF5 ≈68%).

Numerical Highlights

Interpretation: The U(L) reparameterization informed by BAO distances removes a significant fraction of frequency drift in the SN residual IMFs, especially for IMF2–3, without invoking a strong linear correspondence between IF and the BAO proxy itself. This is consistent with a geometric (phase‑relabeling) explanation.

Equations & Definitions

• L ≡ ln(1+z)• D_M(z) = (1+z)·D_A(z)  (comoving angular diameter distance)• D_H(z) = c/H(z)        (Hubble distance)• r_d    = sound‑horizon at baryon drag• BAO proxies: D_M/r_d, D_H/r_d (dimensionless)• U(L)  = normalized cumulative ∫|b(L)| dL• IF: f_k(L) = (1/2π)·dφ_k/dL from the analytic IMF a_k(L)

Discussion

The α‑sweep shows that a blended BAO scale defined by ≈62.5% D_M/r_d and ≈37.5% D_H/r_d best removes chirp in IMF2–3. This blend is physically plausible: D_M and D_H capture transverse and radial rulers respectively; a mixture could mirror how residual features project onto the L variable. The modest/negative locking scores suggest that the transform exploits distance‑ladder geometry rather than pointwise covariation between IF and BAO magnitude.

The strong tightening in IMF4–5 (and noticeable tightening in IMF2–3) implies that multiple scales within the residuals are non‑stationary with L—consistent with slowly varying phase or frequency content (‘chirp‑like’ behavior). De‑chirping re‑centers the IF distribution around a quasi‑stationary band in U, enabling cleaner downstream tests.

Conclusion

We find robust evidence that an L→U(L) reparameterization derived from a BAO proxy significantly reduces the dispersion of instantaneous frequencies for SN Ia residual IMFs, especially in IMF2–3. The preferred mixture uses mostly D_M/r_d with a non‑negligible D_H/r_d component (α≈0.375). This supports the hypothesis that BAO distances encode the appropriate nonlinear ‘clock’ for these oscillatory components.

Conjecture

If SN residual oscillations partly trace large‑scale geometric modulation tied to the BAO ruler, then transforming to U(L) constructed from contemporary BAO ladders should systematically de‑chirp analogous datasets (e.g., alternative SN compilations, binned distance moduli) and sharpen cross‑correlations with galaxy/IGM tracers. Testing stability of α with independent BAO compilations and sub‑samples will help distinguish geometric from astrophysical systematics.

Chirp Analysis Description… “de-chirping”:

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

Concentration increase in the STFT/CWT/HHT energy around one band.

Null testing: run the same de-chirp on surrogates (circular-shift, IAAFT). If your real data shows much larger slope-reduction or variance-reduction than surrogates, you’ve got statistical support for a true chirp.

How this answers your plot

Your median is mostly negative over ⇒ down-chirp with redshift (frequency gently decreases as increases). So try with . After de-chirping you should see:

nearly flat in the mid-range (exclude edge cones),

Stronger locking metrics (or at least more stable) if the chirp is physically tied to your BAO proxy,

Surrogate p-values showing your slope-flattening is unlikely by chance.

## Extracted Tables

### Table 1

Quantity | Value | Notes
Best α for IMF2–3 | α≈0.375 | α=0→DM/rd, α=1→DH/rd; mixture favors DM.
IMF2 MAD ratio | ≈0.54–0.57 | From histograms and sweep.
IMF3 MAD ratio | ≈0.44 | Strongest de‑chirp among mid‑order IMFs.
Energy‑weighted IMF2+3 | ≈0.51 | Mean tightening for the science‑relevant band.
Locking scores (example) | small/negative | Sliding‑window correlation metric; geometry effect dominates.
