# Quaia Redshift Report Update Dec02

- Source: `Quaia Redshift Report Update Dec02.docx`
- Source size: 31401 bytes
- Source modified: 2025-12-02T21:07:00
- Extracted: 2026-07-14
- Word count estimate: 2122

## Extracted Text
Quaia Redshift-Dipole Analysis

FR / Mittal-Masked Update (0V2)

Prepared by: Clive Boyd & Peter Lamb (working notes)
Date: December 2025

1. Overview

This document summarises the current state of the Quaia redshift-dipole analysis, with the pipeline now aligned explicitly to the Full Relativity (FR) framework and updated to include a Mittal-style Galactic Centre (GC) masking strategy.

The goal is to test whether there exists a statistically robust component of the Quaia ⟨z⟩ dipole that:

is aligned with the CMB dipole,

has a kinematic amplitude consistent with FR expectations, and

is stable under changes in redshift range, Galactic latitude cuts, and sky jackknife selections.

2. Pipeline and Scripts

The current analysis pipeline is executed in IPython/Spyder using the following sequence of scripts:

# 1. Build main <z> maps & slices

%run /Users/boyde/.spyder-py3/quaia_build_zmean_maps.py

# 2. FR-style dipole fits (0V2, with CMB projection)

%run /Users/boyde/.spyder-py3/quaia_dipole_fit.py

# 3. Remove N-trend and re-fit dipoles per slice

%run /Users/boyde/.spyder-py3/quaia_dipole_residual_z_vs_counts_slices.py

# 4. Shuffle significance per slice (5000 shuffles, amplitude distributions)

%run /Users/boyde/.spyder-py3/quaia_dipole_shuffle_significance_slices.py

# 5. Build |b| > 10/20/30 maps (with optional GC hole)

%run /Users/boyde/.spyder-py3/quaia_build_zmean_maps_bcut.py

# 6. FR-style bcut dipole grid (0V2)

%run /Users/boyde/.spyder-py3/quaia_dipole_bcut_grid.py

# 7. Detailed |b|>20 dipole summary per slice

%run /Users/boyde/.spyder-py3/quaia_dipole_fit_bcut.py

# 8. Jackknife tests (hemispheres & quadrants, |b|>20)

%run /Users/boyde/.spyder-py3/quaia_dipole_jackknife_bcut20.py

# 9. Skymaps (mean & residual) for full + each z-slice, |b|>20

%run /Users/boyde/.spyder-py3/quaia_plot_maps.py

A “next generation” FR + Mittal-style run focuses on the masked maps and is obtained by rerunning the b-cut and dipole steps after enabling the GC hole in the builder:

%run /Users/boyde/.spyder-py3/quaia_build_zmean_maps_bcut.py

%run /Users/boyde/.spyder-py3/quaia_dipole_bcut_grid.py

%run /Users/boyde/.spyder-py3/quaia_dipole_fit_bcut.py

%run /Users/boyde/.spyder-py3/quaia_dipole_jackknife_bcut20.py

%run /Users/boyde/.spyder-py3/quaia_plot_maps.py

3. Redshift Binning and FR Context

The analysis now adopts a refined set of redshift bins:

Z_BINS = [

[0.10, 0.30],

[0.30, 0.50],

[0.50, 0.75],

[0.75, 1.00],

[1.00, 1.50],

[1.50, 2.50] ]

These bins:

Separate the very low-z regime (0.10–0.30), where peculiar motions and spectroscopic redshifts dominate,

From higher-z regimes where survey selection and large-scale structure become more important.

Within FR, the expectation is that the redshift dipole is primarily a kinematic Doppler effect arising from our motion through an otherwise statistically stationary galaxy background. 

This leads to a cos θ pattern with an amplitude ∝ v/c, which should be approximately independent of redshift once observational systematics are controlled.

4. Mittal-style Galactic Centre Mask

Following the analysis of Mittal et al., an additional mask has been introduced to remove a circular region centred near the Galactic Centre, in addition to the standard Galactic latitude cuts.

The updated builder script quaia_build_zmean_maps_bcut.py (0V2) now applies:

|b| > b_cut with b_cut ∈ {10°, 20°, 30°}

A circular “GC hole” with parameters:

(l₀, b₀) = (0°, 20°)

R = 40°

so that only sources with angular separation from (l₀, b₀) greater than R are retained.

This mask reduces contamination from Galactic dust, stars, and known systematics concentrated toward the Galactic Centre, while still leaving a large and reasonably uniform extragalactic footprint.

The updated runs reported here use the GC hole in combination with |b| cuts to construct new ⟨z⟩ maps and dipole fits.

5. Results with GC Hole and |b| Cuts

5.1 Full Sample, |b| > 20°, GC Hole

With |b| > 20° and the GC hole, the full-sample ⟨z⟩ dipole is measured as:

|b| = 5.52 × 10⁻³

(l, b) = (173.2°, −0.9°)

Separation from CMB dipole ≈ 91.2°

Compared to earlier runs without the GC hole (where the amplitude was ≈ 1.0 × 10⁻² and the direction was closer to (110°, 18°)), the masked run shows:

A reduction in amplitude by roughly a factor of ~2

A rotation of the dipole direction to be almost orthogonal to the CMB dipole

This behaviour is consistent with the idea that part of the previous full-sample signal was driven or biased by structure in/near the Galactic Centre region now removed by the mask.

5.2 Low-z Slice (0.10 ≤ z < 0.30, |b| > 20°)

The FR-critical low redshift slice retains a strong dipole even after applying the GC hole.

For |b| > 20°, with the Mittal-style mask, we obtain:

N_good ≈ 199 pixels

|b| ≈ 8.06 × 10⁻³

(l, b) ≈ (319.1°, 49.4°)

Separation from CMB dipole ≈ 35.5°

In earlier runs without the GC hole, the same slice yielded:

|b| ≈ 9.0 × 10⁻³

(l, b) ≈ (325.5°, 51.1°)

Separation from CMB ≈ 39°

Thus:

The direction of the low-z ⟨z⟩ dipole is remarkably stable under the new masking; it remains in the same quadrant and latitude band, within ~35–40° of the CMB dipole.

The amplitude decreases slightly but remains strong.

This robustness is important if the signal is to be interpreted as having a genuine kinematic component in the FR sense.

5.3 Higher-z Slices

For higher redshift bins (0.30–0.50, 0.50–0.75, 0.75–1.00, 1.00–1.50, 1.50–2.50) with |b| > 20° and the GC hole:

The measured dipole directions vary substantially and are often far from the CMB direction (typically 100–160° separation).

The amplitudes are at the level of ∼(1.4–2.0) × 10⁻³ and show no clean convergence toward a single CMB-aligned direction.

This pattern suggests that at higher redshifts the ⟨z⟩ dipole is dominated by:

Survey footprint effects, and

Large-scale structure,

rather than a single coherent FR-style kinematic dipole.

The FR-relevant information is therefore likely concentrated in the low-z bin, where angular sampling is still good and systematics are better controlled.

5.4 Jackknife Tests (|b| > 20°, GC Hole)

The jackknife analysis in quaia_dipole_jackknife_bcut20.py removes entire hemispheres or longitude quadrants and refits the dipole.

For the full sample, amplitudes can vary by factors of 2–4 and directions swing significantly, indicating that the global dipole is not a simple FR-like, sky-uniform signal.

For the low-z slice (0.10–0.30):

Jackknife amplitudes fluctuate strongly (as expected, given the small number of good pixels),

But many of the jackknife directions remain clustered in the same broad region as the full-slice solution and stay in the CMB-adjacent quadrant.

This provides additional, though still noisy, support for a persistent low-z component that is not entirely driven by one single patch of sky.

6. Relation to the FR Target Structure

The updated pipeline now answers a sharply FR-framed question:

Is there a statistically robust component of the Quaia redshift dipole that is aligned with the CMB dipole, has approximately the expected kinematic amplitude, and is independent of redshift once observational systematics are removed?

Key points linking the current analysis to FR are:

The ⟨z⟩ field is treated as an empirical sky map, without imposing ΛCDM structure, distances, or expansion-based modelling.

The core fitter (quaia_dipole_fit.py 0V2) now explicitly projects the measured dipole onto the CMB direction, yielding amp_par and frac_par per slice.

The refined redshift bins and |b| cuts, together with the new GC hole, allow us to identify a regime (especially low z and |b| > 20–30°) where the dipole direction is CMB-adjacent and relatively stable.

The shuffle significance tests, N-trend residual fits, and jackknife resampling all work together to rule out trivial explanations based purely on survey geometry or number-count gradients.

Preliminary results suggest that the cleanest FR-like behaviour appears in the low-z bin 0.10 ≤ z < 0.30, where the ⟨z⟩ dipole direction remains close to the CMB dipole even after strong masking, while higher-z slices become increasingly dominated by survey and structure effects.

Further work will quantify the CMB-aligned component amp_par for this slice and compare it directly to the FR-predicted value (≈ 0.00123 in the cos θ = 1 direction).

7. Next Steps

Planned follow-up work includes:

Extracting and tabulating the CMB-aligned component amp_par and its uncertainty for the low-z slice, and comparing directly to the FR target (~0.00123).

Systematically exploring different GC-hole parameters (centre and radius) to verify that the low-z CMB-adjacent behaviour is not an artefact of a particular masking choice.

Increasing the number of random shuffles (beyond 5000) for the most critical bins, if higher precision on very small p-values is required.

Cross-checking Quaia results with independent radio or optical number-count datasets where comparable masking strategies can be implemented.

Integrating these dipole constraints into the broader FR cosmological framework, alongside SN, BAO, and CMB fits already developed in the PLamb MCMC scripts.

x1. Low-z CMB-aligned component across various holes

Focusing on the 0.10 ≤ z < 0.30, |b|>20° slice, your runs give:

(l₀=0°, b₀=40°, R=45°)
amp_par ≈ 7.09×10⁻³, frac_par ≈ +0.83, sep ≈ 34°

(0°, 35°, 45°)
amp_par ≈ 7.41×10⁻³, frac_par ≈ +0.86, sep ≈ 30°

(0°, 25°, 45°)
amp_par ≈ 6.15×10⁻³, frac_par ≈ +0.75, sep ≈ 41°

(0°, 25°, 40°)
amp_par ≈ 6.15×10⁻³, frac_par ≈ +0.77, sep ≈ 40°

(0°, 15°, 40°)
amp_par ≈ 7.08×10⁻³, frac_par ≈ +0.84, sep ≈ 33°

(0°, 40°, 42°)
amp_par ≈ 7.13×10⁻³, frac_par ≈ +0.90, sep ≈ 26°

Plus the earlier ones you ran:

(0°, 20°, 40°) → amp_par ≈ 6.56×10⁻³, sep ≈ 35.5°

(0°, 20°, 30°) → amp_par ≈ 7.07×10⁻³, sep ≈ 36.3°

So across all these fairly different GC holes:

amp_par is always ~ (6.1 – 7.4)×10⁻³

i.e. scatter ≲ ~10–15% around ≈ 6.8–7.0×10⁻³

frac_par ≈ +0.75–0.90, so the dipole vector is strongly aligned with the CMB direction in this slice.

The angle to the CMB dipole lives between ≈ 25° and 41°, typically around 30–35°.

That’s very strong evidence that the low-z CMB-aligned component is not an artefact of any single GC-hole choice in this sensible range of centres/radii.

x2. Comparison to the FR “target” amplitude

FR-style kinematic expectation for a pure Doppler dipole using v≈370 km/s is:

ΔzzFR∼vc≈1.23×10-3.

Your measured CMB-aligned component in the Quaia <z> dipole is roughly:

amp∥∼(6.1-7.4)×10-3,

so:

amp∥FR target∼5–6.

So the story you can tell in the write-up is:

Direction: strongly CMB-aligned at low z and robust to mask variations.

Amplitude: consistently several times larger than the naïve FR kinematic prediction, even when we aggressively vary GC-hole centre and radius.

(We still need σ_par printed to quantify this properly in σ-units, but qualitatively the mismatch is big.)

x3. What I’d do right now with these results

Lock in some labelled summaries

For a few representative GC-holes (say three):

A: (0°,  20°,  40°)

B: (0°,  25°,  40°)

C: (0°,  40°,  42°)

In IPython:

!cp skymap2/quaia_outputs/quaia_dipole_summary_bcut20.txt \

skymap2/quaia_outputs/quaia_dipole_summary_bcut20_GChole_C.txt

(After each run, with the appropriate label.)

Add σ_par to your printout

In quaia_dipole_fit_bcut.py, where you have:

amp, amp_par, sigma_par, frac_par = project_dipole_on_cmb_with_cov(b_vec, cov)

expand the print:

print(f"  component along CMB dir: amp_par = {amp_par:.4e} "

f"(sigma_par ≈ {sigma_par:.4e}, frac_par = {frac_par:+.3f})")

Then re-run a couple of key masks to get actual uncertainties on the ~7×10⁻³.

Prepare a small table for the docx

Something like:

GC hole (l₀,b₀,R) | N_good_pix | amp_par (×10⁻³) | frac_par | sep_CMB

(0, 20, 40) | 199 | 6.56 | +0.81 | 35.5°

(0, 25, 40) | 197 | 6.15 | +0.77 | 39.9°

(0, 40, 42) | 189 | 7.13 | +0.90 | 25.9°

(0, 15, 40) | 201 | 7.08 | +0.84 | 32.7°

Then one line of text:

“Across a range of reasonable Galactic Centre masks, the low-z (0.10≤z<0.30), |b|>20° Quaia <z> dipole remains strongly aligned with the CMB dipole with a 
CMB-parallel component ≈(6–7)×10⁻³, i.e. ≳5× the simple FR kinematic expectation (0.00123).”

x4. Next directions (once this is locked)

After you’ve got the σ_par printout and a table:

Fold this into the updated docx section on “CMB-aligned component & mask robustness”.

Then we can:

either move to a quick mock-based calibration: what does a v/c=1.23×10⁻3 kinematic signal actually look like in this <z> dipole estimator?

or jump to the amp_par vs z plot (including error bars) to show the sign flip and strong anti-alignment at high z.
