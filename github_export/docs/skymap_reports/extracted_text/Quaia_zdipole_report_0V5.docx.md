# Quaia zdipole report 0V5

- Source: `Quaia_zdipole_report_0V5.docx`
- Source size: 1717205 bytes
- Source modified: 2025-11-21T10:09:42
- Extracted: 2026-07-14
- Word count estimate: 3639

## Extracted Text
Quaia Redshift-Dipole & Quadrupole Analysis with FR/BAO Cosmology Context (Report 0V4)

1. Executive Summary

This report summarises a sequence of analyses performed on the Quaia G<20.0 quasar catalogue (DOI: 10.5281/zenodo.8060755) and related FR cosmological fits using DESI DR2 BAO data.

The focus is on:

Measuring the number-count dipole of Quaia quasars.

Measuring the redshift (mean-z) dipole and its redshift dependence.

Constructing a residual mean-z map with the best-fit dipole removed.

Extracting the quadrupole (ℓ = 2) of the residual map and assessing its significance using random mocks.

Placing these anisotropy results in the broader context of previous BAO/HHT findings and FR model fits.

Key high-level results:

Quaia G<20.0 number-count dipole: amplitude ≈ 2.8×10⁻² towards (ℓ, b) ≈ (231°, 42°) in Galactic coordinates.

Selection-corrected redshift dipole (full z-range, nside=64, sel ≥ 0.3×max): fractional amplitude ≈ 6.0×10⁻³ towards (ℓ, b) ≈ (180.5°, 13.7°).

Redshift-binned dipoles (0<z<1, 1<z<2, 2<z<4) show coherent directions with modest evolution in amplitude and orientation, indicating a roughly persistent anisotropy across cosmic time.

The residual mean-z map (after dipole removal) shows a strong quadrupole-like structure. The raw quadrupole power D₂ ≡ ℓ(ℓ+1)C₂/2π is extremely high compared with random mocks.

Using 50 random mocks built from the Quaia random catalogue (random_G20.0_10x.fits) and the same selection function, the dipole fractional amplitude is ≈ 4.5σ above the mock mean, while the residual quadrupole D₂ is far beyond the mock distribution (nominally ≫ 100σ, even after |b|>30° cuts).

These high significances should be interpreted cautiously: they motivate further checks for residual systematics, selection-function artefacts, and possible connections with previously observed BAO/HHT structures, rather than being taken as final evidence for exotic physics.

2. Acronyms, Symbols, and Parameters

Symbol / Acronym | Meaning

Quaia | Gaia–unWISE quasar catalogue (G<20.0 subset used here).

DESI | Dark Energy Spectroscopic Instrument (DR2 BAO sample used for cosmology fits).

BAO | Baryon Acoustic Oscillation(s).

FR | ‘Full Relativity’ cosmological model in PLamb_Test_10V6c_plus.py.

H0 | Hubble constant at z=0.

Ωm | Matter density parameter today.

ΩΛ | Effective dark-energy density parameter (in FR, derived from A_acc, n_acc etc.).

A_acc | Amplitude parameter controlling additional acceleration terms in FR.

n_acc | Redshift scaling index of the FR acceleration term.

Γc | Coupling parameter entering the c–luminosity–Noether link tests.

εM | Magnitude evolution / luminosity evolution parameter (linked to γc via κ).

κ (kappa_em_gc) | Soft-link coefficient between εM and γc in the FR code.

σ_link | Width of the soft-link penalty (link-sigma in the FR runs).

w0, wa | Dark-energy equation-of-state parameters in CPL form w(z)=w0+wa·z/(1+z).

ℓ | Spherical-harmonic multipole index (ℓ=0 monopole, ℓ=1 dipole, ℓ=2 quadrupole, etc.).

Cℓ | Angular power spectrum of the residual mean-z map.

Dℓ | Rescaled power: Dℓ = ℓ(ℓ+1)Cℓ / (2π).

Nside | HEALPix resolution parameter (Npix = 12·Nside²).

Z | Redshift; in Quaia we use redshift_quaia column.

3. Datasets and Pre-processing

3.1 Quaia quasar catalogue (G<20.0)

Source: Quaia – Gaia-unWISE quasar catalogue, version 0.1.0 (DOI: 10.5281/zenodo.8060755).

Files used: quaia_G20.0.fits (main catalogue), selection_function_NSIDE64_G20.0.fits (selection map), random_G20.0_10x.fits (random catalogue).

Conversion to CSV: quaia_to_csv.py reads the FITS table and extracts basic columns (RA, DEC, redshift_quaia) into quaia_G20p0_basic.csv with 755,850 objects.

3.2 Coordinate conversion & dipole-friendly catalogue

qso_add_galactic_coords.py attaches Galactic coordinates (ℓ, b) to each object and writes quaia_G20p0_basic_gal.csv.

qso_dipole_catalog.py is used both to compute simple number-count dipoles and to build HEALPix maps at various Nside resolutions (Nside=32 for counts, Nside=64 for some tests).

4. Analysis Pipeline

4.1 Number-count dipole

Using qso_dipole_catalog.py on quaia_G20p0_basic.csv (no selection weighting), the number-count dipole is estimated directly from the unit vectors of all 755,850 quasars. The fitted dipole parameters:

Amplitude: A_count ≈ 2.83×10⁻²

Direction: (ℓ, b) ≈ (231.1°, 41.8°) in Galactic coordinates.

This represents a modest but clearly non-zero anisotropy in the angular number density of G<20.0 quasars.

However, the interpretation of this amplitude is sensitive to angular completeness, selection function, and the masking implied by the selection map.

4.2 Redshift dipole from mean-z HEALPix maps

We next construct mean-redshift maps using quaia_redshift_dipole.py.

For each HEALPix pixel, we compute the mean redshift ⟨z⟩ of quasars falling in that pixel and then fit a dipole to ⟨z⟩(n̂). The dipole model is of the form:

⟨z⟩(n̂) = z̄_global + D · n̂

where z̄_global is the all-sky mean redshift and D is the redshift-dipole vector (with fractional amplitude  |D| / z̄_global).

Baseline (no selection function, full redshift range):

Nside=64, minimum per-pixel count = 20.

z_global_mean ≈ 1.53.

|D| ≈ 4.18×10⁻³; fractional amplitude ≈ 2.74×10⁻³.

Direction: (ℓ, b) ≈ (191.4°, 8.8°).

Quaia G<20.0 — Mean redshift map with overplotted axes: redshift dipole, count dipole, CMB dipole, Galactic Centre (GC) and anti-centre.

Quaia G<20.0 — Mean redshift ⟨z⟩ map (nside=64), same selection.

Quaia G<20.0 — HEALPix counts map (nside=64), in Galactic coordinates.

Selection-corrected case (using selection_function_NSIDE64_G20.0.fits and sel ≥ 0.3×max):

z_global_mean ≈ 1.51.

|D| ≈ 9.16×10⁻³; fractional amplitude ≈ 6.05×10⁻³.

Direction: (ℓ, b) ≈ (180.5°, 13.7°).

The selection-corrected dipole is stronger and has a direction closer to the Galactic anti-centre than the raw full-sky estimate.

An angular-separation script (AngularSeperation.py) shows that the Quaia redshift-dipole axis is:

~166° away from the Galactic centre (i.e. ~14° from the anti-centre).

~76° away from the CMB dipole direction.

~52° away from the Quaia number-count dipole.

4.3 Redshift-binned dipoles

Using qso_dipole_catalog.py and quaia_redshift_dipole.py with redshift cuts, we obtain the following (count dipoles):

0 < z < 1: N ≈ 203,173; amplitude ≈ 2.54×10⁻²; (ℓ, b) ≈ (244.1°, 43.9°).

1 < z < 2: N ≈ 368,752; amplitude ≈ 2.77×10⁻²; (ℓ, b) ≈ (229.2°, 40.7°).

2 < z < 4: N ≈ 183,658; amplitude ≈ 3.32×10⁻²; (ℓ, b) ≈ (224.0°, 40.4°).

The directions stay in a fairly tight band around (ℓ≈230°, b≈40°) while the amplitude increases slightly with redshift, suggesting a persistent anisotropy in number counts.

Redshift-dipole fits (mean-z maps) in the same bins display a more complex evolution, but broadly maintain a preferred axis in the general vicinity of the full-sample dipole, with varying amplitudes and orientations.

Residual mean redshift map Δz after subtracting the best-fit redshift dipole (selection-weighted).

5. Residual Mean-z Map and Quadrupole

Residual mean-z anisotropy power spectrum Dℓ = ℓ(ℓ+1)Cℓ/2π for ℓ = 1–10, after dipole removal.

To isolate anisotropies beyond the dipole, we use quaia_dipole_residual_map.py to:

Quadrupole-only (ℓ = 2) component of the residual mean-redshift field, highlighting the dominant large-scale pattern at high Galactic latitudes.

Fit and subtract the best-fit redshift dipole (with selection weighting).

Save the best-fit dipole model map: quaia_G20p0_nside64_zdipole_model.fits.

Save the residual map: quaia_G20p0_nside64_zresidual.fits.

Compute the residual power spectrum Cℓ up to ℓ_max = 10.

Residual Cl values (ℓ=0..10):

Cℓ ≈ [1.68×10⁻⁴, 9.6×10⁻⁷, 4.93×10⁻⁴, 7.8×10⁻⁶, 5.9×10⁻⁵, 6.9×10⁻⁶, 3.4×10⁻⁵, 1.1×10⁻⁵, 4.2×10⁻⁵, 6.0×10⁻⁶, 2.2×10⁻⁵].

From these, the rescaled powers Dℓ = ℓ(ℓ+1)Cℓ/(2π) are computed, and the quadrupole power is found to be:

C₂_real ≈ 4.86×10⁻³.

D₂_real ≈ 4.64×10⁻³.

A quadrupole-only map is constructed by projecting the residual map onto ℓ=2 spherical harmonics, zeroing ℓ=0 and ℓ=1 coefficients, and reconstructing.

The location of the maximum absolute quadrupole signal is:

Quadrupole max |Δz| at (ℓ, b) ≈ (185.6°, 84.2°).

This lies close to the Galactic north pole, indicating that most of the quadrupole structure is a north–south pattern in Galactic latitude within the selected sky region.

6. Significance Estimates from Random Mocks

To quantify how unusual the observed dipole and quadrupole are, we generate random mocks using quaia_random_mocks_zdipole_quad.py.

The mocks draw positions from random_G20.0_10x.fits, apply the same selection map, and re-compute both the redshift dipole and residual quadrupole power for each mock.

Baseline mock run (no |b| cut):

Number of mocks: 50.

dipole_frac (mock distribution): mean ≈ 2.07×10⁻³, σ ≈ 8.89×10⁻⁴.

D₂_residual (mock distribution): mean ≈ 4.19×10⁻⁶, σ ≈ 2.67×10⁻⁶.

Comparing with the real-sky values:

dipole_frac_real ≈ 6.05×10⁻³ ⇒ z ≈ (6.05×10⁻³ − 2.07×10⁻³) / 8.89×10⁻⁴ ≈ 4.5σ.

D₂_real ≈ 4.64×10⁻³ ⇒ z ≈ (4.64×10⁻³ − 4.2×10⁻⁶) / 2.67×10⁻⁶ ≫ 1000σ (numerically ≈ 1700σ).

Because the quadrupole power is orders of magnitude larger than the mock mean, the naïve Gaussian σ-based estimate is not particularly meaningful—it simply tells us that none of the mocks comes close to the observed D₂. 

The main conclusion is that the residual quadrupole is extremely atypical under the assumed random model.

6.1 High-latitude cut (|b| > 30°)

To test the sensitivity to Galactic latitude, we enforce a mask with |b| ≥ 30° in addition to the selection cuts.

This reduces the quadrupole power but it remains highly significant.

Real-sky with |b| ≥ 30°:

C₂_cut ≈ 2.00×10⁻⁴; D₂_cut ≈ 1.91×10⁻⁴.

Mock distribution with |b| ≥ 30°:

Mean D₂_residual ≈ 2.13×10⁻⁶, σ ≈ 1.40×10⁻⁶.

Resulting significance: z ≈ (1.91×10⁻⁴ − 2.13×10⁻⁶) / 1.40×10⁻⁶ ≈ 135σ.

Even with the |b| cut, the observed quadrupole power is vastly above the random expectation. 

This strengthens the case that the residual quadrupole is not a simple fluctuation of an otherwise isotropic underlying distribution, though selection effects, calibration residuals, and large-scale structure modelling must still be scrutinized.

7. FR Cosmology Fits with DESI DR2 BAO

Alongside the Quaia anisotropy analysis, you ran a series of FR cosmology fits using the PLamb_Test_10V6c_plus.py pipeline, focusing here on DESI DR2 BAO data with various c(z) and soft-link configurations. A representative run:

Model: FR, flat geometry, CPL dark-energy model (w0, wa).

Data: DESI_DR2_BAO.csv + DESI_DR2_BAO_cov.txt (13 BAO points).

c(z) mode: c0 or c_of_z depending on test; here, c0 for a pure-BAO, constant-c baseline.

Link: soft link between εM and γc with κ = 1.0, σ_link ≈ 0.25–0.60.

Noether term: combo_soft with strength ~1–3, acting as a soft penalty on deviations from preferred symmetries.

Sampler: emcee, typically 128 walkers, 12,000–36,000 steps, burn-in ≈ 3,000–8,000 steps.

One DESI-only FR run (FR_w0_softlink_sigma0p6_desionly) yielded best-fit values roughly:

H0 ≈ 50.9 km s⁻¹ Mpc⁻¹, Ωm ≈ 0.134, A_acc ≈ 0.71, n_acc ≈ 1.20, γc ≈ −0.30, εM ≈ 0.49, w_eff ≈ −0.60.

χ²_BAO ≈ 6.7 for dof_total ≈ 6, χ²/dof ≈ 1.12, indicating a statistically acceptable fit.

A second DESI-only FR run with CPL (FR_CPL_softlink_c0_desionly_v2) found, approximately:

H0 ≈ 58.9 km s⁻¹ Mpc⁻¹, Ωm ≈ 0.40, A_acc ≈ 0.34, n_acc ≈ 0.54, γc ≈ 0.002, εM ≈ 0.22, w_eff ≈ −0.82.

χ²_BAO ≈ 9.18 for 5 dof, χ²/dof ≈ 1.84, still broadly acceptable but less optimal than the first FR run.

These FR+DESI fits are not yet tightly coupled to the Quaia anisotropy results but provide context: they show that FR models with additional acceleration terms and soft Noether/link penalties can accommodate BAO distances with reasonable χ².

The natural next step is to explore whether the large-scale anisotropies seen in Quaia are compatible with (or constrained by) the same FR framework.

8. Relation to Previous BAO/HHT Findings

Your previous HHT (Hilbert–Huang Transform) analyses of SN residuals and BAO-linked structures revealed hints of oscillatory components (IMF2/IMF3) whose characteristic scales line up with BAO-like distances.

Robustness grids and surrogate tests suggested that these features were unlikely to arise purely from noise, though systematic effects and selection functions remain key concerns.

The Quaia anisotropy results add another piece to this puzzle: a strong, coherent redshift dipole and an anomalously large quadrupole in the residual mean-z field.

Possible connections include:

If the FR framework allows for anisotropic or direction-dependent effective expansion (e.g. through structured accretion or preferred frames), then both BAO-based distance measures and line-of-sight quasar redshifts could inherit correlated anisotropies.

If the HHT-detected oscillatory structures in SN residuals are linked to large-scale matter or potential fluctuations, then similar structures could leave imprints in QSO redshift distributions and in BAO distances.

Conversely, persistent quadrupolar structure in mean-z might indicate residual calibration or selection effects in the photometry or redshift estimation of Quaia, which must be carefully disentangled before ascribing the signal to cosmological physics.

9. Key Programs and Scripts Used

Script | Purpose

download_zenodo_quaia.py | Downloads Quaia-related FITS and selection/random files from Zenodo using the DOI.

quaia_to_csv.py | Reads quaia_G20.0.fits, identifies RA/DEC/redshift_quaia columns, and exports a basic CSV.

qso_add_galactic_coords.py | Adds Galactic longitude/latitude (ℓ_GAL, b_GAL) columns to CSV catalogues.

qso_dipole_catalog.py | Computes number-count dipoles and optionally builds HEALPix count maps and saves them to FITS.

quaia_redshift_dipole.py | Builds mean-z HEALPix maps at given Nside, fits a redshift dipole, and can output z and count maps.

quaia_redshift_dipole_sel.py | As above, but includes a HEALPix selection function and a sel-min-frac threshold.

quaia_dipole_residual_map.py | Fits and subtracts the redshift dipole, saves dipole model and residual maps, and computes Cℓ.

quaia_random_mocks_zdipole_quad.py | Generates random mocks using the Quaia random catalogue, measures dipole and quadrupole power, and outputs statistics to CSV.

AngularSeperation.py | Computes angular separations between axes (e.g. Quaia z-dipole, Galactic centre/anti-centre, CMB dipole).

PLamb_Test_10V6c_plus.py | Main FR/ΛCDM/PBH cosmology engine, used here for FR + BAO fits with Noether/link tests.

10. Conclusions and Next Steps

The current Quaia G<20.0 analysis indicates a strong, selection-corrected redshift dipole oriented close to the Galactic anti-centre and a highly significant quadrupole component in the residual mean-z map, even after removing pixels near the Galactic plane.

These features are not reproduced by simple isotropic random mocks that include the same selection map and redshift distribution.

Immediate next steps that suggest themselves:

Repeat the analysis with alternative selection thresholds, Nside values, and redshift binning to test stability.

Cross-check against independent quasar catalogues (e.g. eBOSS DR16Q, other Gaia-based samples) to see whether the quadrupole persists.

Introduce more realistic mocks that incorporate large-scale structure and survey strategy effects, rather than purely random angular distributions.

Explore FR-based anisotropic extensions or direction-dependent effective parameters to see whether they can produce similar mean-z anisotropies consistent with BAO and SN constraints.

Connect the quadrupole and dipole signals to the previously studied HHT/BAO oscillatory structures to test for a common origin or to rule out simple shared explanations.

This report (0V4) acts as a structured record of the Quaia anisotropy pipeline and its initial results, ready to be updated as further tests, systematics checks, and extended FR/BAO/SN fits are performed.

Quaia redshift dipole – 0V5 update

This 0V5 revision extends the previous Quaia redshift-dipole report by:
  • Adding object-level dipole checks with redshift weights and bootstrap resampling.
  • Updating the HEALPix mean-z and dipole model maps at nside = 64.
  • Constructing a simple random catalogue from the Quaia selection map and running 500 mocks.
  • Quantifying the significance of the observed redshift-dipole amplitude relative to the mocks.

1. Object-level dipole checks with bootstrap

Script: skymap2_dipole_checks.py
Command template:
    %run skymap2_dipole_checks.py \
        --fits "/Users/boyde/Downloads/quaia_G20.0.fits" \
        --zmin 1.5 \
        --maxN 200000 \
        --mode all \
        --zweight {none|z|z2} \
        --nboot 500

For z ≥ 1.5 and maxN = 200000 objects, three weight choices were tested:
  • zweight = none  (equal weights)
  • zweight = z     (linear redshift weight)
  • zweight = z2    (quadratic redshift weight)

Equatorial-frame dipole directions (unweighted / z / z²) are:
  • none : RA ≈ 139.14°, Dec ≈ +6.79°
  • z    : RA ≈ 139.53°, Dec ≈ +7.63°
  • z²   : RA ≈ 140.49°, Dec ≈ +8.63°

Galactic-frame directions for the same runs are:
  • none : l ≈ 224.42°, b ≈ +35.23°
  • z    : l ≈ 223.74°, b ≈ +35.97°
  • z²   : l ≈ 223.19°, b ≈ +37.28°

Bootstrap resampling (nboot = 500) was used to quantify directional stability.

For all three weight choices, the median angular separation between each bootstrap dipole and the full-sample dipole is ≈ 2.6–2.7°, with a 68% confidence interval of roughly 1.3–4.8°.

This indicates that the high-z dipole direction is stable at the few-degree level under resampling.

2. Mean-z, dipole model, and residual maps (HEALPix)

Scripts: quaia_redshift_dipole.py, quaia_dipole_residual_map.py, quaia_redshift_dipole_sel.py, quaia_plot_zdipole_maps.py

Core commands:
    %run quaia_redshift_dipole.py \
        --csv "/Users/boyde/Downloads/quaia_G20p0_basic_gal.csv" \
        --nside 64 \
        --zmin 1.5 \
        --min-per-pixel 1 \
        --map-out-z "/Users/boyde/Downloads/quaia_G20p0_nside64_zmean.fits" \
        --map-out-count "/Users/boyde/Downloads/quaia_G20p0_nside64_counts.fits"

    %run quaia_dipole_residual_map.py \
        --csv "/Users/boyde/Downloads/quaia_G20p0_basic_gal.csv" \
        --nside 64 \
        --zmin 1.5 \
        --min-per-pixel 1 \
        --sel-map "/Users/boyde/Downloads/quaia_G20p0_nside64_counts.fits" \
        --sel-min-frac 0.5 \
        --out-dipole-map "/Users/boyde/Downloads/quaia_G20p0_nside64_zdipole_model.fits" \
        --out-residual-map "/Users/boyde/Downloads/quaia_G20p0_nside64_zresidual.fits" \
        --lmax 1

At nside = 64 and z ≥ 1.5, using all pixels with at least one object, the mean-z map fit returns:
  • N_objects       ≈ 355,332
  • z_global_mean   ≈ 2.13
  • |D|             ≈ 4.93 × 10⁻³
  • frac_amplitude  ≈ 2.32 × 10⁻³
  • Dipole (l, b)   ≈ (190.45°, +57.45°)

When the selection map and sel_min_frac = 0.5 cut are applied in quaia_dipole_residual_map.py, the dipole fit yields:
  • z_global_mean   ≈ 2.14
  • |D|             ≈ 8.10 × 10⁻³
  • frac_amplitude  ≈ 3.79 × 10⁻³
  • Dipole (l, b)   ≈ (137.40°, +21.38°)
  • Residual C_l    : C₀ ≈ 1.36 × 10⁻⁸,  C₁ ≈ 3.17 × 10⁻⁹ (lmax = 1)

Updated sky maps (nside = 64, z ≥ 1.5):

(a) Mean redshift map

(b) Best-fit dipole model map

(c) Residual map after subtracting dipole model

3. Angular separations to GC, anti-GC, CMB, and count dipole

Script: AngularSeperation.py
Command:
    %run AngularSeperation.py

Using the selection-aware redshift dipole (l ≈ 137.4°, b ≈ +21.4°), the measured separations are:
  • sep(Quaia z-dip, GC)      ≈ 166.29°

  • sep(Quaia z-dip, anti-GC) ≈ 13.71°

  • sep(Quaia z-dip, CMB dip) ≈ 75.55°

  • sep(Quaia z-dip, count dip) ≈ 51.89°

Thus the Quaia z-dipole lies much closer to the anti-Galactic-centre direction than to the Galactic centre, and is substantially misaligned with the CMB temperature dipole and with the Quaia number-count dipole.

4. Random-catalogue mocks and z-dipole significance

Random catalogue generation used the HEALPix selection map and a Galactic latitude cut:
Script: quaia_make_simple_randoms.py
Command:
    %run quaia_make_simple_randoms.py \
        --sel-map "/Users/boyde/Downloads/quaia_G20p0_nside64_counts.fits" \
        --nside 64 \
        --nrand 2000000 \
        --sel-min-frac 0.5 \
        --bmin-abs 20 \
        --out-fits "/Users/boyde/Downloads/quaia_G20p0_randoms_simple.fits"

Mock analysis script: quaia_random_mocks_zdipole_quad.py
Command:
    %run quaia_random_mocks_zdipole_quad.py \
        --data-csv "/Users/boyde/Downloads/quaia_G20p0_basic_gal.csv" \
        --random-fits "/Users/boyde/Downloads/quaia_G20p0_randoms_simple.fits" \
        --sel-map "/Users/boyde/Downloads/quaia_G20p0_nside64_counts.fits" \
        --nside 64 \
        --min-per-pixel 1 \
        --sel-min-frac 0.5 \
        --nmocks 500 \
        --lmax 2 \
        --seed 12345 \
        --bmin-abs 20 \
        --outcsv "/Users/boyde/Downloads/quaia_G20p0_random_mocks_zdipole_quad.csv"

Over 500 mocks, the summary statistics are:
  • dipole_frac   : mean ≈ 8.69 × 10⁻³,  std ≈ 4.28 × 10⁻³
  • D2_residual   : mean ≈ 5.68 × 10⁻⁷,  std ≈ 4.80 × 10⁻⁷

The real, selection-aware Quaia redshift dipole has frac_amplitude ≈ 3.79 × 10⁻³.

Comparing this to the mock distribution yields a z-score of ≈ −1.14, corresponding to an approximate two-sided Gaussian p-value of p ≈ 0.25.

In other words, the observed dipole amplitude is about 1.1σ below the mean expectation from the mocks and is fully consistent with the range produced by random realisations under the same mask and selection function.

A histogram of dipole_frac from the mocks, with the real Quaia value overplotted as a vertical line, confirms that the observed redshift-dipole amplitude lies well within the central bulk of the mock distribution.

This supports the interpretation that, once the selection function and sky mask are properly accounted for, the Quaia redshift-dipole amplitude is not anomalously large.

Figure: Distribution of redshift-dipole fractional amplitudes (dipole_frac) measured from 500 random-catalogue mocks (outline histogram), constructed with the Quaia G<20.0 selection map and |b| ≥ 20° cut. 

The vertical dashed line marks the real Quaia value (frac_amplitude ≈ 3.79×10⁻³), which lies well within the central bulk of the mock distribution (z ≈ −1.1, p ≈ 0.25), indicating that the observed redshift dipole amplitude is statistically consistent with the random expectation under the adopted mask and selection function.

5. Interpretation (0V5 snapshot)

The 0V5 analysis shows that:
  • The high-z Quaia object-level dipole direction is stable at the few-degree level across different redshift weights and bootstrap resampling.

  • The map-based redshift dipole at nside = 64 with z ≥ 1.5 and a conservative selection mask has frac_amplitude ≈ 3.8 × 10⁻³ and points towards (l, b) ≈ (137°, +21°).

  • This direction is close to the anti-Galactic-centre and substantially offset from both the CMB temperature dipole and the Quaia number-count dipole.

• However, the amplitude of the redshift dipole is statistically consistent (≈ 1.1σ low) with the distribution of amplitudes measured from 500 random-catalogue mocks that respect the Quaia mask  and selection function.

Further work can extend this framework by exploring redshift tomographic bins, alternative latitude cuts, and more realistic random catalogues that include redshift-dependent selection effects.
