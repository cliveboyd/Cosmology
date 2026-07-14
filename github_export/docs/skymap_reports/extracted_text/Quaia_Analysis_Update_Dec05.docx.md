# Quaia Analysis Update Dec05

- Source: `Quaia Analysis Update Dec05.docx`
- Source size: 357448 bytes
- Source modified: 2025-12-26T18:50:18
- Extracted: 2026-07-14
- Word count estimate: 4953

## Extracted Text
1. Overview

This report documents analysis of large-scale anisotropy in the mean redshift field ⟨z⟩ of the Quaia quasar catalogue (G < 20). The goals are:

Measure the dipole of ⟨z⟩ in different Galactic latitude cuts (|b| > b_cut) and redshift slices.

Compare the measured dipoles with the CMB kinematic dipole to look for parallel or anti-parallel alignment.

Use rotation Monte Carlo (MC) to estimate the significance of any observed alignment.

Recently: increase the minimum number of objects per pixel from 3 to 5
(N_MIN_PER_PIX = 5) to improve robustness, then re-run the full rotation MC suite.

There is emerging evidence for a redshift-modulated alignment pattern:

mid-z bins (roughly 0.5–1.0) tend to show CMB-parallel dipoles;

high-z bins (1.5–2.5) tend to show CMB-anti-parallel dipoles;

one mid-z bin (0.50–0.75) has an unusually large dipole amplitude, independent of CMB alignment.

At present, the statistical significance is modest (typically ≤ 2σ), but the pattern is persistent across masks and the new N≥5 maps, and is therefore worthy of further investigation.

Note:	This report has a number of report updates that may look a bit disjointed due to the addition of forthcoming results from various programs and analysis options.

Re-Indexing has not been addressed

At this stage HHT has not been actioned. 

Planning on enhancing z slices to 0.2 in order to produce more bins and possibly modify min pixel count to be adaptive.

2. Data, maps and masks

Catalogue: Quaia G20.0 (755,850 objects).

Magnitude cut: G < 20.

HEALPix resolution: NSIDE = 64 (npix = 49,152).

Galactic latitude cuts

∣b∣>bcut,bcut=10∘,15∘,20∘,25∘,30∘,35∘.

Additional Galactic centre (GC) hole

centre (l, b) = (0°, 20°), radius 40°.

Minimum per pixel (new setting)

Nmin=NMIN_PER_PIX=5

i.e. only pixels containing at least 5 quasars contribute to the ⟨z⟩ dipole fits.

Redshift slices currently used

Non-overlapping slices:

0.10 ≤ z < 0.30

0.30 ≤ z < 0.50

0.50 ≤ z < 0.75

0.75 ≤ z < 1.00

1.00 ≤ z < 1.50

1.50 ≤ z < 2.50

Each slice is tagged in the code as:

z0p10_0p30,

z0p30_0p50,

z0p50_0p75,

z0p75_1p00,

z1p00_1p50,

z1p50_2p50.

3. Dipole model and definitions

Let ⟨z⟩(n̂) be the mean redshift in direction n̂ on the sky. 

We fit a monopole + dipole model

⟨z⟩(n)=z0+A⋅n,

where

z0is the monopole,

Ais the dipole vector with amplitude

A=∣A∣.

The best-fit direction is expressed in Galactic coordinates (l, b).

The CMB dipole direction is denoted dCMB.

We project the ⟨z⟩ dipole onto this direction:

A∥=A⋅dCMB,

and define the CMB-parallel fraction

f∥=A∥A.

Interpretation:

f∥≈+1: dipole almost parallel to CMB dipole.

f∥≈0: roughly orthogonal.

f∥≈-1: almost perfectly anti-parallel.

We also measure the angular separation between the ⟨z⟩ dipole direction and the CMB dipole:

θsep=arccos⁡ ⁣AdCMB,A=A∣A∣.

4. Rotation Monte Carlo (MC) method

For each bcutz-slice combination:

Take the observed ⟨z⟩ map (with mask and N ≥ 5 per-pixel cut).

Generate nmock=2000 random rotations of the entire map on the sphere
(preserving the sky pattern but scrambling its orientation relative to the CMB dipole).

For each rotated map, recompute the dipole:

amplitude Ai,

CMB-parallel component A∥i,

fraction f∥i.

Compute empirical p-values:

for the amplitude

p(A)=#{i:Ai≥Aobs}Nmock,

for alignment

pA∥=#{i:∣A∥i∣≥∣A∥,obs∣}Nmock,

p(∣f∥∣)=#{i:∣f∥i∣≥∣f∥,obs∣}Nmock.

These are effectively two-sided tail probabilities relative to the rotation ensemble.

Rough σ-equivalents:

p≈0.32→ ~1σ

p≈0.10→ ~1.6σ

p≈0.05→ ~1.96σ

No single bin currently exceeds ≈ 2σ, but some are in the 1.5–2σ range, especially at high z and for certain latitude cuts.

5. Effect of raising N_MIN_PER_PIX from 3 to 5

We rebuilt the ⟨z⟩ maps and dipole grid with N_MIN_PER_PIX = 5, then re-ran the full rotation MC suite.

Key consequences

Low-z slice 0.10–0.30:

Ngoodpixels with N ≥ 5 collapses to only 1–3 pixels across all |b| cuts.

This slice is not usable for robust dipole estimation with N ≥ 5.

Any apparent “strong” low-z signal from earlier N ≥ 3 runs must therefore be treated as fragile and possibly dominated by a few noisy pixels.

Mid- and high-z slices (z ≥ 0.30):

Still have hundreds to tens of thousands of good pixels.

Dipole measurements remain stable and well-constrained.

Summary

The N ≥ 5 cut substantially improves robustness at the cost of discarding the shallow 0.10–0.30 slice.

All current interpretations therefore focus on z ≥ 0.30 and on the new N ≥ 5 rotation MC results.

6. Current dipole grid: qualitative summary (N ≥ 5)

From quaia_dipole_bcut_grid.txt (N_MIN_PER_PIX = 5, NSIDE 64):

6.1 Full sample (all z, varying |b| cuts)

Amplitudes A≈(4.6–5.7)×10-3 across bcut=10–35∘.

Dipole directions are roughly in the Galactic anti-centre hemisphere, with latitudes around b ≈ 0° to −16°.

Angular separations from the CMB dipole are large,
θsep≈85–100∘, i.e. approximately orthogonal to the CMB dipole.

Corresponding fractions f∥are small: −0.18 ≲ f∥≲ +0.06.

So the global ⟨z⟩ dipole does not align strongly with the CMB dipole.

6.2 Redshift dependence (example: |b| > 20° cut)

Below are the notable trends for |b| > 20° (numbers approximate):

0.30 ≤ z < 0.50

A≈3.4×10-3, Ngood≈180.

Dipole direction roughly at (l ≈ 87°, b ≈ 42°).

Separation from CMB: θsep≈90∘.

f∥≈0(nearly orthogonal).

0.50 ≤ z < 0.75

A≈4.5×10-3.

Direction (l ≈ 128°, b ≈ 43°).

θsep≈82∘(still roughly orthogonal, but now in the same quadrant as CMB).

f∥≈+0.15(mildly CMB-parallel).

0.75 ≤ z < 1.00

A≈1.1×10-3.

Direction (l ≈ 114°, b ≈ 69°).

θsep≈61∘.

f∥≈+0.48→ moderately CMB-parallel.

1.00 ≤ z < 1.50

A≈2.0×10-3.

Direction (l ≈ 60°, b ≈ 13°).

θsep≈116∘.

f∥≈-0.43→ moderately anti-parallel.

1.50 ≤ z < 2.50

A≈2.1×10-3.

Direction (l ≈ 90°, b ≈ −25°).

θsep≈157∘.

f∥≈-0.92→ strongly CMB-anti-parallel.

This is the key pattern:

Mid-z bins (~0.5–1.0) tend to align parallel to the CMB dipole, whereas 
high-z bins (≥ 1.5, especially 1.5–2.5) tend to align anti-parallel to the CMB.

The same qualitative behaviour persists across bcut=15∘, 20∘, 25∘, 30∘, 35∘:

low-to-mid z: f∥>0(CMB-parallel),

high z: f∥<0(almost anti-parallel).

7. Selected rotation MC results (N ≥ 5 runs)

The latest rotation MC run (N_MIN_PER_PIX = 5, n_mock = 2000, seed = 42) gives updated p-values that are fully consistent with the current pixel cut.

7.1 Slice 0.50 ≤ z < 0.75: strong dipole amplitude, weak CMB link

This bin shows the clearest amplitude anomaly, but no strong preference for CMB alignment.

|b|>10°:

Aobs≈3.82×10-3,

f∥,obs≈+0.40.

p(A)=0.001→ very rare amplitude in rotations.

p(∣f∥∣)≈0.61→ alignment not unusual.

|b|>15°: p(A)=0.045.

|b|>20°: p(A)=0.009.

|b|>25°: p(A)=0.002.

|b|>30°: p(A)=0.083.

|b|>35°: p(A)=0.049.

Across these cuts, the dipole amplitude is consistently high relative to the rotated skies (p(A) typically 0.001–0.05), while p(∣f∥∣)stays large (~0.65–0.70).

This suggests a real, strong ⟨z⟩ dipole in the 0.50–0.75 bin, but not one that is specifically tied to the CMB dipole direction.

7.2 Slice 0.75 ≤ z < 1.00: CMB-parallel tendency at high |b|

Amplitudes in this slice are modest, but at high latitude cuts the dipole becomes strongly CMB-parallel:

|b|>30°:

f∥,obs≈+0.90.

p(∣f∥∣)≈0.095(~1.7σ).

|b|>35°:

f∥,obs≈+0.94.

p(∣f∥∣)≈0.064(~1.9σ).

Amplitude p-values in this slice are large (p(A) ≈ 0.18–1.0), i.e. unremarkable in size.

The interest here is purely directional: a marginal (~1.7–1.9σ) preference for CMB-parallel alignment at 0.75–1.00 when restricting to high latitudes.

7.3 Slice 1.50 ≤ z < 2.50: persistent CMB-anti-parallel signal

This slice shows the most consistent CMB-anti-parallel behaviour:

|b|>10°:

f∥,obs≈-0.90,

p(∣f∥∣)=0.10.

|b|>15°:

f∥,obs≈-0.95,

p(∣f∥∣)=0.045(~2σ).

|b|>20°:

f∥,obs≈-0.92,

p(∣f∥∣)=0.082.

|b|>25°:

f∥,obs≈-0.85,

p(∣f∥∣)=0.141.

|b|>30°:

f∥,obs≈-0.97,

p(∣f∥∣)=0.023(most extreme case, ~2.0–2.3σ).

|b|>35°:

f∥,obs≈-0.84,

p(∣f∥∣)=0.145.

Amplitude p(A) values here are mostly unremarkable (0.35–0.74, except p(A)=0.038 for |b|>30°).

The primary anomaly in this slice is directional: a high-z dipole consistently anti-parallel to the CMB, with the strongest case (|b|>30°) at p(|f_\parallel|) ≈ 0.023.

7.4 Full-sample bins

For the full-z sample at each |b| cut:

p(A)≈0.26–0.98,

p(∣f∥∣)≈0.81–0.97.

The full-sample dipole amplitude and orientation are completely consistent with the rotation ensemble, i.e. no single all-z kinematic-like anomaly.

8. Possible redshift-modulated alignment (“modulation artefact”)

The current results, especially from quaia_dipole_bcut_grid.txt and the updated N ≥ 5 rotation MC, suggest a systematic evolution of f∥ with redshift:

f∥ transitions from positive (CMB-parallel) at mid-z to strongly negative 
(CMB-anti-parallel) at high-z.

This pattern is present across multiple |b| cuts, indicating that it is not simply a local mask artefact.

The new MC confirms that several of these bins have p-values in the 0.02–0.10 range (roughly 1.6–2.3σ), particularly:

amplitude in 0.50–0.75,

CMB-parallel alignment in 0.75–1.00 at |b|≥30–35°,

CMB-anti-parallel alignment in 1.50–2.50 at |b|≥15°, especially |b|>30°.

There are at least three broad possibilities:

Astrophysical / cosmological modulation
Large-scale structure and peculiar motions could produce a redshift-dependent dipole pattern that changes phase relative to the CMB rest frame.

More exotic interpretations could involve evolving anisotropy or foreground structures.

Survey / selection systematics
Redshift-dependent completeness, extinction residuals, or calibration drifts across the sky could imprint a spurious dipole whose direction varies with z.

Masking and variable depth near the Galactic plane could contribute.

Statistical fluctuations
Given multiple bins and cuts, some 1.5–2σ features are expected by chance.

At this stage, the main novel observation is:

A redshift-modulated sign flip in the CMB-parallel fraction f∥(z), with mid-z slices mildly parallel to the CMB dipole and high-z slices strongly anti-parallel, plus one mid-z bin with an anomalously large dipole amplitude.

9. Proposed strategy: modulation analysis with sliding z-windows and HHT

To better understand whether this is a real modulation or just a few noisy bins, we can move beyond disjoint slices and treat f∥(or related quantities) as a function of redshift.

9.1 Sliding redshift windows

Instead of fixed non-overlapping bins, define overlapping windows, e.g. width Δz ≈ 0.4–0.5, step δz ≈ 0.1–0.2, such as:

0.10–0.50, 0.20–0.60, 0.30–0.70, …, 1.50–2.50.

For each window:

Build ⟨z⟩ map with N ≥ 5.

Fit dipole and compute A,A∥,f∥.

Optionally record θsepand the raw dipole components in Cartesian form.

This yields a quasi-continuous function of redshift:

f∥(zmid), A(zmid), A∥(zmid).

9.2 HHT / EMD-based modulation analysis

Once we have f∥(zmid)sampled at many z_mid values, we can treat it as a 1-D “time series” in z and apply Hilbert–Huang Transform (HHT) methods already used in the HHT/BAO work:

Empirical Mode Decomposition (EMD)
Decompose f∥(z)into Intrinsic Mode Functions (IMFs) plus a trend.

Identify whether there is a dominant mode representing a slow modulation 
(e.g. transition from + to −).

Hilbert spectral analysis
For each IMF, compute instantaneous amplitude and phase vs z.

Look specifically for phase-coherent modes across b_cuts or different catalogue subsets (e.g. NGC vs SGC).

Surrogate tests
Generate surrogates by randomly scrambling z or by rotating sky positions within each z window, then recompute the f∥(z)trajectory.

Apply the same HHT pipeline to surrogates.

Quantify how often a modulation as strong/coherent as the observed one appears by chance.

Cross-comparisons
Compare the modulation pattern between different |b| cuts and between NGC/SGC hemispheres.

A truly cosmological modulation should be robust across reasonable mask changes.

This would let us move from “a few 1.5–2.3σ bins” to a global significance test of the entire f∥(z)pattern.

10. Current status summary

Maps and dipole grid rebuilt with N_MIN_PER_PIX = 5 (NSIDE 64, G < 20).

Low-z 0.10–0.30 slice effectively unusable (1–3 good pixels only). 
We now treat any apparent low-z alignment as unreliable.

Mid- and high-z slices remain robust, with hundreds to tens of thousands of good pixels.

Dipole grid shows a clear redshift trend in CMB-parallel fraction:

mid-z (0.5–1.0): f∥>0, roughly CMB-parallel;

high-z (≥ 1.5): f∥≈-0.8to −0.95, strongly anti-parallel.

Updated rotation MC (N ≥ 5, n_mock = 2000):

Full-sample amplitudes and orientations are unremarkable (p(A) ≈ 0.26–0.98, p(|f_\parallel|) ≈ 0.81–0.97).

The 0.50–0.75 slice shows a strong amplitude anomaly with p(A) typically in the 0.001–0.05 range across |b| cuts, but with no special CMB alignment.

0.75–1.00 shows CMB-parallel alignment at high |b|, with 
p(|f_\parallel|) ≈ 0.06–0.10 (≈1.7–1.9σ).

1.50–2.50 shows CMB-anti-parallel alignment, most strongly for 
|b|>30° where p(|f_\parallel|) ≈ 0.023 (~2.0–2.3σ).

Overall, several bins now sit in the p ≈ 0.02–0.10 range (roughly 1.6–2.3σ), 
but no single result yet exceeds ≈2.5σ.

Next major step: implement a sliding-window, HHT-based modulation pipeline for f∥(z)to test whether the apparent sign flip and modulation with redshift is statistically meaningful or a chance/systematic artefact.

11. Acronyms & notation (for quick reference)

CMB 		– Cosmic Microwave Background.

GC 		– Galactic Centre.

HEALPix 	– Hierarchical Equal Area isoLatitude Pixelization.

NSIDE 	– HEALPix resolution parameter.

MC 		– Monte Carlo.

EMD 		– Empirical Mode Decomposition.

HHT 		– Hilbert–Huang Transform.

IMF 		– Intrinsic Mode Function (in HHT/EMD).

N_MIN_PER_PIX (N_min) 
                         – minimum number of objects per pixel required to include that 
                            pixel in dipole fits.

A 		– dipole amplitude of ⟨z⟩.

A_∥ 		– component of the dipole parallel to the CMB dipole.

f_∥ 		– CMB-parallel fraction of the dipole, f∥=A∥/A.

θ_sep 		– angular separation between the ⟨z⟩ dipole and the CMB dipole.

1. What the shuffle MC with N ≥ 5 is saying

Big picture: nothing looks crazy under pixel–shuffle, even in the bins that looked a bit hot in the rotation MC.

Key points from your log (N_MIN_PER_PIX = 5, n_mock = 2000):

Low-z (0.10–0.30)

|b|>10° and 15°: p(A) ≈ 0.33–0.34, p(|A∥|) ≈ 0.49–0.50.

But this slice only has N_good = 3 pixels, so it’s not trustworthy anyway.

0.30–0.50 (all b-cuts)

p(A) ≈ 0.64, 0.80, 0.87, 0.66, 0.98, 0.98 → completely unremarkable.

No sign of an amplitude excess; any alignment/anti-alignment we saw in rotation MC is purely an orientation effect.

0.50–0.75 (the band where rotation MC gave low p(A) for some b-cuts)

p(A) ≈ 0.26, 0.27, 0.17, 0.17, 0.34, 0.46 across b_cuts = 10–35°.

So under shuffle, the amplitude is entirely consistent with isotropic noise + mask.

Conclusion: the “low p(A)” from rotation MC here was not an intrinsic ⟨z⟩–dipole excess, just how that fixed pattern happens to sit inside the rotating mask.

1.00–1.50 (this is where something mildly interesting appears)

p(A) by b_cut:

|b|>10°: 0.116

|b|>15°: 0.048

|b|>20°: 0.137

|b|>25°: 0.053

|b|>30°: 0.083

|b|>35°: 0.030

So there’s a consistent trend of slightly low p(A) for 1.0–1.5, especially at higher b_cuts (0.048, 0.053, 0.030 → ~1.7–2.1σ at best, before trials factors).

p(|A∥|) are all ~0.14–0.89 → no strong CMB–direction anomaly here; it’s an amplitude effect.

1.50–2.50 (where rotation MC showed strong anti-parallel f∥)

p(A) ≈ 0.56, 0.51, 0.41, 0.34, 0.78, 0.77 across b_cuts.

p(|A∥|) ≈ 0.20, 0.13, 0.12, 0.11, 0.32, 0.38.

So shuffle sees nothing special in amplitude; and A∥ along the CMB axis is also not particularly extreme.

The high-z “anti-parallel” behaviour therefore looks like a pure orientation effect relative to the CMB, not extra dipole power.

So the main update is:

Only the 1.0–1.5 slice shows a mild, b_cut-robust excess in dipole amplitude under pixel–shuffle (p(A) ≈ 0.03–0.05 at high |b|). 

Everything else, including the 0.5–0.75 band and the high-z 1.5–2.5 anti-alignment, is consistent with an isotropic ⟨z⟩ field once you account for the mask.

2. Text you can paste into the report (replace / extend sections)

You can drop this in as a replacement for your current Section 7 and tweak Section 10.

7. Rotation vs shuffle MC results (N ≥ 5 maps)

With the ⟨z⟩ maps rebuilt at NSIDE = 64 and N_MIN_PER_PIX = 5, we have now run both:

the rotation MC (quaia_mc_bcut_rotate.py), and

the pixel–shuffle MC (quaia_mc_bcut_shuffle.py, n_mock = 2000, seed = 52).

These probe complementary nulls:

Rotation MC keeps the ⟨z⟩ pattern fixed on the sky and randomly rotates it relative to the CMB dipole → tests whether the observed dipole is unusually aligned or anti-aligned with the CMB given the mask.

Shuffle MC keeps the mask and pixel sampling fixed, but randomly reassigns redshifts, thus erasing any real large-scale structure → tests whether the dipole amplitude itself is larger than expected from shot noise + mask.

All quoted p-values below correspond to the new N ≥ 5 configuration.

7.1 Low-z, 0.10 ≤ z < 0.30

For all b_cuts (|b| > 10°, 15°, …, 35°) we have N_good ≈ 1–3 pixels only.

Shuffle MC gives p(A) ≈ 0.33–0.34 and p(|A∥|) ≈ 0.49–0.50, i.e. nothing unusual under the null, but in practice this slice is not usable and we discard it from interpretation.

7.2 0.30 ≤ z < 0.50

Across all b_cuts:

Shuffle MC: p(A) ≈ 0.64–0.98, p(|A∥|) ≈ 0.75–0.94.

Thus there is no evidence for an anomalous ⟨z⟩ dipole amplitude in this band. 

Any interesting behaviour seen in rotation MC (e.g. apparent anti-alignment at some b_cuts) is a pure orientation effect relative to the CMB, not an excess of dipole power.

7.3 0.50 ≤ z < 0.75

This band previously showed low p(A) in the rotation MC for some b_cuts 
(particularly |b| > 10°–25°), suggesting a potentially large dipole amplitude.  However:

Shuffle MC now gives p(A) ≈ 0.26 (|b|>10°), 0.27 (15°), 0.17 (20°), 0.17 (25°), 0.34 (30°), 0.46 (35°).

All corresponding p(|A∥|) are ≈ 0.40–0.72.

So under the isotropy+noise null, this amplitude is entirely consistent with chance. 

The low p(A) from rotation MC was therefore driven by how the fixed, real dipole pattern sits within the rotated mask, rather than by a genuinely “too big” ⟨z⟩ dipole.

7.4 0.75 ≤ z < 1.00

Shuffle MC: p(A) ≈ 0.77–0.89 and p(|A∥|) ≈ 0.61–0.66 for |b|>10°–25°, and still high at larger b_cuts.

So there is no amplitude anomaly. 

The moderate positive f∥ seen in the dipole grid and rotation MC (CMB-parallel tendency) is once again an orientation feature only, with amplitude well within expectations.

7.5 1.00 ≤ z < 1.50 – mild amplitude excess

This slice is the only one where shuffle MC suggests a mildly enhanced ⟨z⟩–dipole amplitude:

Shuffle p(A) by b_cut:

|b|>10°: 0.116

|b|>15°: 0.048

|b|>20°: 0.137

|b|>25°: 0.053

|b|>30°: 0.083

|b|>35°: 0.030

The trend is that for higher latitude cuts (|b| ≥ 15°–35°), p(A) lies in the 
range ≈ 0.03–0.05, i.e. roughly 1.7–2.1σ before accounting for multiple trials.

At the same time:

p(|A∥|) is modest (≈ 0.14–0.89), and

f∥ ≈ −0.3 to −0.45 in the dipole grid,

so this is best described as a slightly stronger-than-expected ⟨z⟩ dipole amplitude, with only a mild preference for anti-alignment in CMB coordinates.

7.6 1.50 ≤ z < 2.50 – anti-parallel orientation but normal amplitude

For the highest redshift band:

Shuffle MC: p(A) ≈ 0.56, 0.51, 0.41, 0.34, 0.78, 0.77 across b_cuts 10°–35°.

p(|A∥|) ≈ 0.20, 0.13, 0.12, 0.11, 0.32, 0.38.

Thus:

Amplitude: perfectly compatible with isotropy + noise; no evidence for extra dipole power.

Orientation: rotation MC and dipole grid still show strongly negative f∥ (≈ −0.85 to −0.95) and large θ_sep (≈ 150–170°), i.e. the dipole is strongly anti-parallel to the CMB, but that anti-alignment is not accompanied by an anomalously large amplitude.

In short, the high-z anisotropy is a directional effect only, at ≲ 2σ significance per bin.

Updated bullet for Section 10 (you can replace the two bullets about MC there):

– Rotation MC (n_mock = 2000, N ≥ 5): confirms the redshift-dependent orientation pattern in f∥(z): mid-z dipoles tend to be CMB-parallel; high-z dipoles (1.5–2.5) tend to be strongly anti-parallel. Individual bins reach at most ≈ 1.5–2σ in |f∥|.

– Shuffle MC (n_mock = 2000, N ≥ 5): shows that, under an isotropic ⟨z⟩ field plus noise, the dipole amplitudes are mostly unremarkable, 

with the exception of a mild, b_cut-robust excess in the 1.0–1.5 slice (p(A) ≈ 0.03–0.05 at high |b|). 

The 0.5–0.75 band and the high-z 1.5–2.5 band, in particular, are fully consistent with the null in amplitude, so their interesting behaviour is purely in orientation relative to the CMB.

3. How to call the shuffle plotting script

Your last error is just because quaia_plot_mc_bcut.py needs arguments. For example, to plot the shuffle MC for:

|b|>30°,

0.50 ≤ z < 0.75,

you’d do:

%run /Users/boyde/.spyder-py3/quaia_plot_mc_bcut.py \

--bcut 30 --zlo 0.50 --zhi 0.75

and similarly for other bands, e.g.:

# High-z anti-parallel band

%run /Users/boyde/.spyder-py3/quaia_plot_mc_bcut.py \

--bcut 30 --zlo 1.50 --zhi 2.50

# The mildly anomalous amplitude band

%run /Users/boyde/.spyder-py3/quaia_plot_mc_bcut.py \

--bcut 35 --zlo 1.00 --zhi 1.50

(If the script has a --save option you can tack that on the end as well.)

If you’d like, next step we can sketch the sliding-window f∥(z) + HHT pipeline and decide concrete z-steps and window widths based on the N_good counts you’re seeing.

Update after running… quaia_mc_bcut_shuffle.py

1. What the new shuffle plots show

You plotted three interesting cases from the pixel-shuffle MC (N_MIN_PER_PIX = 5):

|b| > 30°, 0.50 ≤ z < 0.75 (z ≃ 0.62)

Amplitude A: observed value sits toward the upper side of the distribution but well inside the bulk → consistent with the printed p(A) ≈ 0.34.

CMB-parallel component A_∥: observed value is slightly positive, again well within the central Gaussian region → p(|A_∥|) ≈ 0.49.

Fraction f_∥: histogram is ~flat, red line at f_∥ ≈ +0.36 is nowhere near an extreme tail.

➜ No evidence for unusual CMB alignment in this slice once you randomise pixel redshifts.

|b| > 30°, 1.50 ≤ z < 2.50 (z ≃ 2.0)

Amplitude A: observed near the mode of the distribution → p(A) ≈ 0.78.

A_∥: observed is fairly negative but not in the tails → p(|A_∥|) ≈ 0.32.

f_∥: again essentially flat between −1 and +1, with the observed f_∥ ≈ −0.97 on the extreme left edge, but that is by construction: in the shuffled ensemble f_∥ is uniform, so any given f_∥ value (even near −1) has probability ≈ a few percent in 2000 mocks. That matches the p(|A_∥|) value.

➜ Shuffling wipes out the “strong anti-parallel” look – it’s not an outlier relative to a uniform f_∥ distribution.

|b| > 35°, 1.00 ≤ z < 1.50 (z ≃ 1.25)

Amplitude A: the observed vertical line is clearly on the high side of the distribution → p(A) ≈ 0.03 (as in the summary file). That’s ~2σ in amplitude only.

f_∥: distribution is again flat and the observed f_∥ ≈ −0.04 sits near zero → p(|A_∥|) ≈ 0.89.

➜ Here the interesting thing is amplitude, not direction; alignment with the CMB is entirely ordinary.

A general pattern in all f_∥ histograms:

In the shuffle MC, f_∥ is very close to uniform on [−1, +1], as expected when you randomise the redshifts per pixel and destroy any coherent sky pattern.

The observed f_∥ (red line) is never at an extreme tail for the shuffle ensemble, consistent with the p(|A_∥|) values.

2. How this updates the written report

Shuffle MC with N_MIN_PER_PIX = 5
To test whether the apparent redshift-dependent alignment could arise purely from the per-pixel redshift distribution and mask geometry, we ran a complementary “pixel shuffle” Monte Carlo (quaia_mc_bcut_shuffle.py) with N_mock = 2000 and N_MIN_PER_PIX = 5.
In each mock, the list of pixel ⟨z⟩ values is randomly permuted while holding the sky positions (and N_pix) fixed. 

This destroys any coherent large-scale structure signal but preserves the overall distribution of ⟨z⟩ and the mask.

Results:

For all b-cuts and redshift slices, the shuffle f_∥ distributions are essentially uniform on [−1, +1], as expected when the dipole direction is random. 

The observed f_∥ values lie comfortably within this uniform ensemble, with typical p(|A_∥|) ~ 0.1–0.9.

In particular, the high-z bin 1.50 ≤ z < 2.50, |b| > 30°, which showed a suggestive CMB anti-alignment in the rotation MC (p(|A_∥|) ≃ 0.02–0.03), is not unusual in the shuffle MC (p(|A_∥|) ≃ 0.32).

One bin, 1.00 ≤ z < 1.50, |b| > 35°, shows a moderately large amplitude relative to the shuffled ensemble (p(A) ≃ 0.03, ≈ 2σ), but its alignment with the CMB is completely ordinary (p(|A_∥|) ≃ 0.89).

Interpretation:

The shuffle MC confirms that the directional anomalies seen in the rotation tests (e.g. high-z anti-alignment) are not simply an artefact of the per-pixel ⟨z⟩ distribution or mask; when we destroy the sky pattern by shuffling pixels, any preference relative to the CMB disappears and f_∥ becomes uniform.

The modest p(A) ≃ 0.03 in one mid-z, high-b bin suggests a slightly high dipole amplitude there, but with no special CMB orientation.

Overall, the redshift-modulated sign flip in f_∥(z) remains a feature of the real sky orientation (rotation MC) rather than of the local pixel statistics (shuffle MC), but its significance is still at the ≲ 2σ level.

These results strengthen the case for a more global, sliding-window analysis of f_∥(z) (and A(z)), combined with HHT/EMD and surrogate tests, to determine whether the apparent modulation is statistically meaningful or consistent with chance alignment plus cosmic variance.

7.4 Extra MC check: ∣b∣>35∘,  zmid≃1.25

For the ∣b∣>35∘,  1.00≤z<1.50slice we now have:

Rotate-MC (sky rotations, N_mock = 2000)

Observed amplitude: Aobs≈3.4×10-3.

The rotation ensemble for Ais roughly Gaussian with mean ∼1.7×10-3; 
Aobslies on the high-amplitude tail (see “A MC distribution
 ∣b∣>35,z∼1.25” plot), but still within the cloud of realizations.

The CMB-parallel component A∥is very close to zero, and the fraction f∥=A∥/A≈-0.04. 

The corresponding histograms (“A∥ MC distribution” and “f∥ MC distribution”) show that both A∥and f∥sit essentially at the centre of their rotation distributions. 

There is no evidence for either parallel or anti-parallel alignment in this bin – only a moderately large total dipole amplitude.

Shuffle-MC (pixel shuffling, N_mock = 2000, N_{\min}=5)

For the same bin, the pixel-shuffle experiment gives
p(A)≈0.03, p(∣A∥∣)≈0.89.

This means that, given the actual sky mask and N(z) structure in this slice, a dipole as large as Aobsoccurs only in ≈3% of random shufflings of the pixel values, but its projection along the CMB dipole is completely unremarkable.

Interpretation.

The ∣b∣>35∘,  zmid≃1.25slice hosts a somewhat high ⟨z⟩ dipole amplitude, but this excess is not preferentially aligned with the CMB dipole. 

Combined with the shuffle-MC result, this points to a genuine large-scale pattern in the ⟨z⟩ field (i.e. real structure plus mask geometry) rather than a CMB-related kinematic effect. 

In other words, this bin contributes to the overall redshift-modulated amplitude pattern, but does not strengthen the case for a CMB-parallel / anti-parallel modulation beyond what is already seen in the high-z anti-parallel bins.
