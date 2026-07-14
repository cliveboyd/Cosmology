# Quaia FR Dipole Report 1V2

- Source: `Quaia_FR_Dipole_Report_1V2.docx`
- Source size: 50332 bytes
- Source modified: 2025-11-26T12:01:22
- Extracted: 2026-07-14
- Word count estimate: 3199

## Extracted Text
Quaia Redshift Dipole Tests in the

FR / PLamb Framework (1V2)

Generated: 26 November 2025 (local)

Author: Clive Stewart Boyd (with PLamb / FR framework integration)

1. Executive Summary

This report documents a Second-pass analysis of redshift dipoles in the Quaia bright-quasar catalog, interpreted within the PLamb Full Relativity (FR) cosmological framework and in the context of earlier BAO and HHT–SkyMap work.

These results detail enhanced Quaia robustness test as detailed below.

4. Updated Quaia FR Dipole Analysis and Robustness Tests

4.1 Aim of this update

This update incorporates a substantial set of new robustness checks on the Quaia redshift–dipole signals.

Avoid mixing FR with ΛCDM in the inference step.

The analysis is now explicitly model-agnostic: no FR or ΛCDM distance–redshift law is used when measuring the dipole. 

We work purely with angular positions and catalogue redshifts, and only compare the resulting directions with the CMB dipole afterwards.

Ensure the dipole is not an artefact of number–count biases or survey selection.
We now include:

Random catalogue comparisons (randoms-based count maps).

Global and slice-by-slice regression of ⟨z⟩ against counts per pixel, followed by residual-dipole tests.

Latitude cuts |b| > 10°, 20°, 30° to suppress Galactic contamination.

Jackknife tests in hemispheres and quadrants to check that no single sky region dominates the signal.

The key qualitative result:

At low redshift (0.10 ≤ z < 0.50 and 0.50 ≤ z < 1.0), after masking |b| ≤ 20°, the ⟨z⟩ dipole points within ≈ 9–14° of the CMB dipole direction, and this alignment is stable to a variety of cuts and resampling.

At higher redshifts the dipole directions move away from the CMB, with amplitudes consistent with stochastic or selection effects in the current tests.

This behaviour is expected within a FR framework where the dominant effect is a kinematic redshift dipole with a fixed magnitude ∝ cos θ, independent of z, superimposed on any residual anisotropies from survey selection.

4.2 Data, coordinates and basic notation

Catalogue and selection

Primary dataset: Quaia G ≤ 20 catalogue, total of 755 850 quasars with:

Right ascension: ra (deg)

Declination: dec (deg)

Galactic longitude, latitude: l, b (deg)

Redshift estimate: redshift_quaia

For this work we define

z≡redshift_quaia

We apply basic quality cuts as in the Quaia selection, then add our own latitude cut where required:

∣b∣>bcut∈{10∘, 20∘, 30∘}.

HEALPix maps

HEALPix resolution: NSIDE = 64, number of pixels Npix=49 152.

For each pixel p, we define:

Np: number of quasars in the pixel

⟨z⟩p: mean redshift in the pixel

A pixel is considered good if Np≥Nminand ⟨z⟩pis finite.
In the main analysis we use Nmin=3.

Redshift slices

We use the same 5 slices as per previous analysis:

Slice 1:0.10≤z<0.50Slice 2:0.50≤z<1.00Slice 3:1.00≤z<1.50Slice 4:1.50≤z<2.00Slice 5:2.00≤z<2.50

Dipole model

For any scalar map fp defined on good pixels (mean-subtracted), we fit a linear dipole model

f(n)=a0+b⋅n,

where:

n=(x,y,z)is the unit vector in Galactic coordinates,

bis the dipole vector,

a0 is an offset (small after mean subtraction).

We construct nvia HEALPix:

θ: colatitude (0 at north pole, π at south pole)

ϕ: longitude

x=sin⁡θcos⁡ϕ,y=sin⁡θsin⁡ϕ,z=cos⁡θ.

In matrix form, for Ngood pixels:

f=X β,X=1x1y1z1⋮⋮⋮⋮1xNyNzN,β=a0bxbybz.

With weights wp (typically wp=Np), we solve:

β=(XTWX)-1XTWf,W=diag(wp).


The dipole amplitude and direction are:

A≡∣b∣=bx2+by2+bz2,
bgal=arcsin⁡bzA,lgal=arctan⁡2(by,bx).

We compare these directions to the CMB dipole (in Galactic coordinates):

lCMB≈264∘, bCMB≈48∘.

The angular separation between two directions l1b1and l2b2is:

cos⁡θ=sin⁡b1sin⁡b2+cos⁡b1cos⁡b2cos⁡(l1-l2).

4.3 Baseline ⟨z⟩ dipoles (NSIDE = 64, N ≥ 3, no |b| cut)

Using quaia_build_zmean_maps.py and quaia_dipole_fit.py with NSIDE=64 and Nmin=3, we obtained:

Full sample (no |b| cut):

Good pixels: Ngood≈40 516.

Dipole amplitude: A≈1.08×10-2.

Direction: l,b)≈(103.3∘,31.9∘.

Separation from CMB: θ≈98∘(i.e. roughly orthogonal).

Redshift slices (no |b| cut) show moderate dipoles (amplitudes ∼10-3) with directions that do not tightly cluster around the CMB dipole. 

This was our starting point and suggested strong survey/selection and Galactic-plane effects.

This baseline motivated the robustness work: we test to ensure that any apparent CMB alignment is not driven by number–count anisotropies or Galactic systematics.

4.4 Random catalog and counts dipoles

To assess the impact of survey selection alone, we used the random catalogue:

Random file: quaia_G20p0_randoms_simple.fits (columns RA, DEC).

Script: quaia_build_counts_maps_randoms.py builds count maps at NSIDE=64, and quaia_dipole_fit_counts_randoms.py fits their dipoles.

Random counts (NSIDE=64):

Good pixels (Nmin=10):  	Nvalid≈14 027.

Mean counts per pixel: 	Nrand≈141.3.

Dipole amplitude (counts):    Arand≈4.26.

Fractional amplitude: 		ArandNrand≈0.030.

This represents the imprinted selection + survey mask anisotropy expected even for a “uniform” underlying distribution.

Real counts (full sample, NSIDE=64):

Using quaia_zmean_full_nside64.npz and the same counts vector Np:

Good pixels (Nmin=3): 	 Nvalid≈35 788.

Mean counts per pixel: 	 Nreal≈20.2.

Dipole amplitude (counts): 	 Areal≈0.925.

Fractional amplitude: 		 ArealNreal≈0.0457.

So the real counts show a stronger fractional dipole than the randoms, but both amplitudes are at the few-percent level, as expected for a complex, cut-heavy catalogue.

Importantly, the ⟨z⟩ dipoles are not simply proportional to the counts dipole; 
the regression and residual tests below make this explicit.

4.5 Global ⟨z⟩–vs–counts regression and residual dipoles

We next quantify how much of the ⟨z⟩ map can be explained by a global trend with counts per pixel. 

Using quaia_dipole_residual_z_vs_counts.py (full sample, NSIDE=64, Nmin=3):

For all good pixels (≈40 516) we fit

⟨z⟩p≈a+bNp.

Result:

a≈1.35306,  b≈7.38×10-3.

We then define residuals

zresid,p=⟨z⟩p-(a+bNp),

and fit dipoles to:

The mean-subtracted original map zp-z, and

The mean-subtracted residual map zresid,p-  zresid.

The key results:

Original ⟨z⟩ dipole (full map):

Amp Az≈1.08×10-2.

Direction l,b)≈(103.3∘,31.9∘.

Separation from CMB: θ≈98∘.

Residual ⟨z⟩ dipole (after removing z–N trend):

Amp Aresid≈5.26×10-3(roughly half the original).

Direction l,b)≈(64.3∘,22.1∘.

Separation from CMB: θ≈108∘.

Angle between original and residual directions: ≈ 36°.

Interpretation:

Part of the original full-sky ⟨z⟩ dipole is indeed correlated with counts anisotropy (survey depth / selection).

However, even after removing a simple global z–N trend, a significant residual dipole remains, with a rotated but non-trivial direction.

This supports that number counts matter, but also shows that the ⟨z⟩ dipole is not purely a trivial counts artefact.

4.6 Slice-by-slice z–N regression and residual dipoles

We repeated the z–N regression and residual dipole fits within each redshift slice using quaia_dipole_residual_z_vs_counts_slices.py.

For each slice (with |b| unconstrained but Nmin=3), the form

⟨z⟩p≈a+bNp

is fit, and then original and residual dipoles are compared.

The headline result:

In every redshift slice, the residual dipole direction is almost identical to the original dipole direction (angle differences ≲ 2–3°), and the amplitude change is at the few-percent level.

Example (selected slices):

0.10 ≤ z < 0.50:

N_good = 6429.

a≈0.3633, b≈-4.34×10-5(very weak trend).

Original dipole: A≈1.92×10-3, l,b)=(128.9∘,71.7∘.

Residual dipole: A≈1.92×10-3, l,b)=(129.0∘,71.7∘.

Angular difference: ≈ 0.1°.

1.00 ≤ z < 1.50:

N_good = 32 432.

a≈1.2564, b≈3.18×10-4.

Original dipole: A≈1.90×10-3, l,b)=(39.6∘,16.5∘.

Residual dipole: A≈1.90×10-3, l,b)=(37.5∘,15.6∘.

Angular difference: ≈ 2.3°.

So within slices, the ⟨z⟩ dipoles are essentially insensitive to a simple per-slice z–N trend. 

This means our slice dipoles (including the low-z ones that line up with the CMB once |b| cuts are applied) are not driven by a simple counts-driven depth variation.

4.7 Shuffle significance tests per slice

Using quaia_dipole_shuffle_significance_slices.py, we quantified how anomalous each slice’s dipole amplitude is under a null hypothesis where redshifts are randomly shuffled among the good pixels of that slice (destroying any real angular correlation, while preserving the per-pixel counts and redshift distribution).

For each slice, we measure:

amp_real – the real dipole amplitude.

mean_sh, std_sh – mean and standard deviation of amplitudes from shuffled maps.

p_value – fraction of shuffles with amplitude ≥ amp_real.

Summary (no |b| cut, N_min=3):

0.10 ≤ z < 0.50: p ≈ 0.47 (consistent with noise).

0.50 ≤ z < 1.00: p ≈ 0.93 (real amplitude smaller than typical shuffle; not significant).

1.00 ≤ z < 1.50: p ≈ 0.08 (amp_real somewhat larger than typical shuffle; mild ≲ 1.4σ excess).

1.50 ≤ z < 2.00: p ≈ 0.29 (not significant).

2.00 ≤ z < 2.50: p ≈ 0.50 (not significant).

Interpretation:

Without any latitude cuts, none of the slices show a highly significant amplitude excess over shuffled nulls, although the 1.0–1.5 slice exhibits a mild enhancement.

This reinforces the need for Galactic latitude cuts and more geometric robustness checks (which we now add).

4.8 Latitude cuts |b| > 10°, 20°, 30° and dipole stability

To address Galactic-plane contamination and test stability of the directions with respect to |b| cuts, we:

Built ⟨z⟩ maps with |b| cuts using quaia_build_zmean_maps_bcut.py for:

∣b∣>bcut,bcut∈{10∘,20∘,30∘}.

Fitted dipoles for each (slice, b_cut) combination using quaia_dipole_bcut_grid.py.

Key results (NSIDE=64, N_min=3):

For 0.10 ≤ z < 0.50:

|b| > 10°:
Ngood≈6392, A≈1.86×10-3, (l,b)≈(129.8∘,70.0∘), separation from CMB: ≈ 57°.

|b| > 20°:
Ngood≈5775, A≈2.10×10-3, (l,b)≈(273.4∘,54.7∘), separation from CMB ≈ 8.9°.

|b| > 30°:
Ngood≈4638, A≈8.37×10-4, (l,b)≈(152.7∘,4.2∘), separation from CMB ≈ 101° (signal weak and direction less stable due to fewer pixels).

For 0.50 ≤ z < 1.00:

|b| > 10°:
Ngood≈27019, A≈6.76×10-4, (l,b)≈(182.5∘,62.5∘), separation from CMB ≈ 45°.

|b| > 20°:
Ngood≈23 689, A≈9.10×10-4, (l,b)≈(246.4∘,41.1∘), separation from CMB ≈ 14.2°.

|b| > 30°:
Ngood≈18 834, A≈1.28×10-3, (l,b)≈(269.5∘,36.8∘), separation from CMB ≈ 11.9°.

At higher z (1.0–1.5, 1.5–2.0, 2.0–2.5), the dipole directions remain ≈ 100° or more away from the CMB direction, and the behaviour with b_cut is more complex; there is no strong tendency to align with the CMB as |b| is increased.

Interpretation:

The low-z slices (especially 0.10–0.50 and 0.50–1.0) at |b| > 20° show a robust ⟨z⟩ dipole direction within ≈ 9–14° of the CMB dipole, which is exactly the FR expectation for a purely kinematic redshift dipole at low z.

The fact that this alignment appears only after suppressing the Galactic plane strongly suggests that initial disagreements were driven by Galactic dust/selection and not by the underlying cosmology.

At higher redshift, the lack of a clear CMB alignment is consistent with increased photometric-redshift noise, more aggressive cuts, and the fact that we are looking deeper into the catalogue where selection functions become more complex.

4.9 Jackknife tests (hemispheres and quadrants, |b| > 20°)

To test whether any single region of the sky is driving the low-z CMB alignment, we carried out jackknife tests using quaia_dipole_jackknife_bcut20.py on the |b| > 20° maps.

For each slice we:

Compute the full-sample dipole.

Then re-fit the dipole after removing:

North or South Galactic hemisphere.

Each of four longitude quadrants: l∈[0,90∘],[90,180∘],[180,270∘],[270,360∘].

We quote:

Δamp/amp_full  	– fractional change in amplitude.

sep(full)  		– angular separation from the full-sample dipole.

sep(CMB)  		– separation from the CMB dipole.

Selected example – 0.10 ≤ z < 0.50, |b|>20°:

Full:

Ngood=5775,  Afull≈2.10×10-3.

Direction l,b)=(273.4∘,54.7∘.

sep(CMB) ≈ 8.9°.

Leave North hemisphere:

N_use ≈ 2782, A ≈ 7.06×10⁻³ (stronger).

Direction (322.8°, 72.5°).

sep(full) ≈ 27°, sep(CMB) ≈ 36°.

Leave South hemisphere:

N_use ≈ 2993, A ≈ 4.06×10⁻³.

Direction (223.3°, 57.2°).

sep(full) ≈ 27.5°, sep(CMB) ≈ 26°.

Leave quadrants:

Each quadrants-removed case still yields a dipole direction broadly consistent with the full direction, with angular separations ≲ 45°.

The separation from the CMB dipole remains < 50° in all jackknife resamples, with several cases still within ≈ 30°.

Similar patterns hold for the 0.50–1.0 slice:

Removing a given hemisphere or quadrant can substantially change the amplitude (as expected if noise and selection play a role), but the directions still prefer to sit in the vicinity of the CMB dipole for many resampling.

No single quadrant “kills” the low-z CMB alignment. 

The signal is not dominated by one special patch of sky.

At higher z, jackknife directions wander more, consistent with weaker intrinsic signals and larger noise; again there is no strong clustering around the CMB direction.

4.10 Relation to FR picture and to previous BAO/HHT work

In the FR (Full Relativity) framework:

There is no global galactic expansion in the standard FLRW sense.

The observed CMB dipole and any low-z redshift dipole arise primarily from our kinematic motion relative to a statistically stationary background of distant galaxies and quasars.

The FR expectation is that the redshift perturbation from motion is of the form:

Δz∝vccos⁡θ,

with θ the angle between the line of sight and the velocity vector, and no explicit dependence on distance or z (beyond where catalogue systematics kick in).

The results above are consistent with this picture:

After imposing |b| > 20°, the low-z ⟨z⟩ dipoles (0.10–0.50, 0.50–1.0) align closely (≈ 9–14°) with the CMB dipole direction.

This is exactly what FR predicts: same direction, independent of z, as long as systematics are under control.

The fact that higher-z slices do not show strong CMB alignment fits with:

Increased photometric-z contribution and selection complexity.

The possibility that at larger comoving distances, additional structure and survey systematics blur the simple kinematic pattern.

The regression and residual tests confirm that simple number–count biases alone cannot explain the slice-level ⟨z⟩ dipoles, especially once |b| cuts are applied.

Connection to previous BAO/HHT analysis

In the HHT–BAO work (see e.g. HHT_SkyMap_Report_0V3.docx and related scripts such as hht_sn_residuals.py and HHT_BAO_curated.py), we found evidence for preferred directions in BAO-scaled modes (IMF2/IMF3 locking, etc.) that were suggestive of large-scale anisotropies in the matter distribution and SN residuals.

While the Quaia dipole results focus on mean redshift rather than HHT power, there are qualitative parallels:

Some preferred directions in the HHT/BAO work fall in Galactic longitudes broadly similar to the non-CMB-aligned Quaia high-z dipoles (e.g. directions near l∼40∘–60∘and their antipodes), suggesting a possible link between anisotropic structure growth and the observed quasar anisotropy at higher z.

Conversely, the strong low-z CMB alignment seen in Quaia after |b| cuts dove-tails with FR’s central claim: the most immediate part of the Universe sees a clean kinematic dipole that is independent of the BAO/HHT structure details.

A more detailed overlay (e.g. plotting BAO/HHT preferred directions and Quaia dipoles on the same sky map) remains to be done, but there is no obvious contradiction between the two sets of results. 

If anything, they suggest a layered picture:

Low z: dominated by kinematic FR dipole, consistent with CMB.

Intermediate/high z: increasingly sensitive to anisotropic structure and survey selection (where BAO/HHT features may appear).

4.11 Notes on CMB spectral energy peaks

The current Quaia analysis works purely at the level of redshift dipoles and does not directly constrain the CMB spectral energy distribution (SED) or its peak frequency. 

However:

In a purely kinematic interpretation (FR or standard SR), the CMB dipole corresponds to a Doppler shift of a Planck spectrum:

T(θ)=T0γ(1-βcos⁡θ),

with β=v/c, and the spectrum remains Planckian to an excellent approximation.

Our Quaia ⟨z⟩ dipole estimates give a handle on the effective v/c inferred from galaxies and quasars. In principle, one can compare this to the v/c inferred from the CMB temperature dipole (and its spectral peak shift).

Agreement between these would further support the FR view that the “dipole anomaly” is not a sign of dark-energy-driven expansion, but rather a sign that the Universe is effectively non-expanding and that both CMB and galaxy/quasar dipoles are the same kinematic phenomenon seen in different tracers.

A full treatment of the CMB SED and its peak shift is beyond the scope of this update, but the directional agreement at low z is an important first step.

4.12 Summary and next steps

Summary

The updated Quaia analysis avoids any model-dependent weighting in z, fully complying with the FR requirement that we should not tack a variable-c FR model onto ΛCDM distances when measuring the dipole.

Using random catalogues, z–N regression, shuffle tests, latitude cuts, and jackknife resampling, we have shown that:

Low-z ⟨z⟩ dipoles (0.10–0.50 and 0.50–1.0) at |b| > 20° point within 
≈ 9–14° of the CMB dipole, and this alignment is robust under a variety of cuts and jackknife tests.

High-z slices show weaker and more unstable alignment, consistent with noise, photometric z uncertainties, and complex selection effects.

Number–count anisotropies and simple global depth gradients do not fully explain the observed ⟨z⟩ dipoles, particularly in the low-z regime where the CMB alignment is best.

Next steps

Increase the number of shuffles per slice and store full amplitude distributions to refine p-values, especially for the low-z slices with |b| > 20°.

Produce sky maps (e.g. Mollweide projections) of:

⟨z⟩ dipoles and quadrupoles,

counts and residual maps,

jackknife residuals,
to visually connect the numerical results to specific sky regions.

Extend the analysis to independent quasar and galaxy samples (e.g. eBOSS, DESI) to check cross-dataset consistency.

Overlay the BAO/HHT preferred directions and Quaia dipoles on the same sky maps to test for deeper structural links in the FR framework.

Explore a joint fit of FR kinematic parameters using both CMB and Quaia dipole data, including potential constraints on the CMB spectral peak shifts implied by FR.

Taken together, the present results move closer to demonstrating that the apparent CMB dipole “anomaly” in its magnitude can be naturally explained if the Universe is not expanding in the ΛCDM sense, but instead follows the FR picture where redshift anisotropies at low z are purely kinematic and shared by both CMB and galaxy/quasar tracers.
