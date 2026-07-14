# Quaia Dipole Quadrupole Report 0V1 v2

- Source: `Quaia_Dipole_Quadrupole_Report_0V1_v2.docx`
- Source size: 42987 bytes
- Source modified: 2025-11-17T15:05:18
- Extracted: 2026-07-14
- Word count estimate: 2248

## Extracted Text
Quaia Redshift Dipole and Quadrupole Analysis (0V1)

Quaia Redshift Dipole and Quadrupole Analysis (0V1)

Author: Clive Stewart Boyd (analysis notebooks and scripts)

Collaborator context: Peter Lamb (quasar anisotropy hypothesis)

Date: 17 November 2025

1. Executive Summary

This report summarises the current Quaia-based analysis of large-scale anisotropy

in the mean quasar redshift field. Using the Gaia–unWISE Quaia catalogue (G < 20),

we construct HEALPix maps of the mean redshift ⟨z⟩, fit a dipole, and then study

the residual field and its quadrupole (ℓ = 2) component. A suite of random mocks

based on the Quaia random catalogue and selection function is used to estimate the

expected noise level and to assess the significance of the observed signals.

Key points:

The mean redshift field shows a clear dipole aligned close to the Galactic

anti-centre, with fractional amplitude Δz / z ≈ 6.0 × 10⁻³, corresponding to

≈ 4.5σ above selection + shot-noise mocks.

After subtracting the best-fit dipole, the residual mean-redshift field displays

a very strong quadrupole whose symmetry axis is almost aligned with the North

Galactic Pole. The quadrupole power at ℓ = 2 remains roughly two orders of

magnitude above the mock expectation even when restricting to |b| > 30° and

applying the official selection function.

The combination of a strong dipole and a strong quadrupole indicates that the

distance–redshift relation traced by quasars is not isotropic on the sky, at

least under the assumptions of the mock model. Whether this is cosmological

(e.g. anisotropic metric / H(z) / variable-c) or due to subtle survey

systematics requires further work.

These findings are qualitatively consistent with earlier PLamb/FR analyses of

BAO data, where anisotropic DH/DM ratios (with DH/DM ≈ 0.35 along a preferred

axis) suggested differential radial vs transverse expansion. Quaia provides an

independent quasar-based probe of similar large-scale structure in the STF

(space-time fabric).

2. Acronyms and Notation

Quaia    : Gaia–unWISE quasar catalogue (version 0.1.0, G < 20.0 used here).

HEALPix  : Hierarchical Equal Area isoLatitude Pixelisation of the sphere.

NSIDE    : HEALPix resolution parameter (npix = 12 × NSIDE²).

ℓ, m     : Spherical harmonic multipole indices.

Cℓ       : Angular power spectrum coefficients of a field on the sphere.

Dℓ       : ℓ(ℓ+1) Cℓ / (2π), used for visualising angular power.

⟨z⟩      : Mean redshift in a given sky pixel.

Dipole   : ℓ = 1 component; one-lobed pattern (one hemisphere high, one low).

Quadrupole: ℓ = 2 component; two-lobed / ring vs poles pattern.

(l, b)   : Galactic longitude and latitude (degrees).

GC       : Galactic Centre (l = 0°, b = 0°).

anti-GC  : Galactic anti-centre (l = 180°, b = 0°).

BAO      : Baryon Acoustic Oscillations.

DH       : BAO radial distance measure (∝ 1 / H(z)).

DM       : BAO transverse comoving distance measure.

DH/DM    : Ratio of radial to transverse BAO distances.

FR       : Full Relativity cosmological framework in PLamb analysis.

STF      : Space-Time Fabric (PLamb terminology for large-scale metric behaviour).

3. Datasets Used

3.1 Quaia quasar catalogue (G < 20)

Primary quasar catalogue:

quaia_G20.0.fits  (main quasar sample with G < 20.0)

random_G20.0_10x.fits  (random angular positions at 10× density)

selection_function_NSIDE64_G20.0.fits  (NSIDE = 64 selection function map)

For the present analysis the essential columns are:

ra              : Right Ascension (ICRS, degrees)

dec             : Declination (ICRS, degrees)

redshift_quaia  : photometric redshift estimate (used as z)

3.2 Reduced CSV catalogue

To simplify the analysis pipeline, a reduced CSV was created:

/Users/boyde/Downloads/quaia_G20p0_basic.csv

Columns:

RA  : Right Ascension (degrees, ICRS)

DEC : Declination (degrees, ICRS)

Z   : redshift_quaia

Galactic coordinates were added via qso_add_galactic_coords.py, giving:

/Users/boyde/Downloads/quaia_G20p0_basic_gal.csv

Additional columns:

L_GAL : Galactic longitude (degrees)

B_GAL : Galactic latitude (degrees)

3.3 Random catalogue and selection function

The Quaia random catalogue and selection function are used to:

1) Define a selection-aware mask (sel ≥ sel_min_frac × max(sel), typically 0.3).

2) Build mocks by drawing random positions from random_G20.0_10x.fits and

assigning redshifts by resampling the real Quaia redshift distribution.

4. Methodology

4.1 HEALPix mean-redshift maps

Using NSIDE = 64, each quasar is converted from (RA, DEC) to Galactic (l, b)

and then to a HEALPix pixel p. For each pixel we accumulate:

N_p   : number of quasars in pixel p

Σz_p  : sum of redshifts in pixel p

⟨z⟩_p : Σz_p / N_p  (if N_p > 0; otherwise pixel is set to UNSEEN)

This yields:

counts_map  N_p

mean_z_map  ⟨z⟩_p

4.2 Redshift dipole model

The first-order anisotropy model for the mean-z field is:

⟨z( n̂ )⟩ ≈ z₀ + D · n̂

where

z₀   : global mean redshift

D    : 3D dipole vector (units of redshift)

n̂   : unit vector in the direction (l, b)

For pixel p we define

δz_p = ⟨z⟩_p − z₀ ≈ D · n̂_p,

and solve for D by weighted least squares, with weights ∝ N_p.

4.3 Pixel selection

To ensure a robust fit, pixels are required to satisfy:

N_p ≥ N_min (commonly N_min = 10)

selection function sel_p ≥ sel_min_frac × max(sel); typically sel_min_frac = 0.3

This excludes sparsely populated pixels and regions with poor completeness.

4.4 Residual map and quadrupole

Once D and z₀ are determined, we construct a dipole model map:

z_model( n̂ ) = z₀ + D · n̂

The residual map is

z_res( n̂ ) = ⟨z( n̂ )⟩ − z_model( n̂ )

We then compute the angular power spectrum Cℓ of z_res( n̂ ) using healpy.anafast,

optionally with a |b| cut (e.g. |b| ≥ 30°) in addition to the selection mask.

From Cℓ we define

Dℓ = ℓ (ℓ + 1) Cℓ / (2π)

and examine particularly ℓ = 2 (quadrupole). A strong peak at ℓ = 2 implies a

strong quadrupole in the residual mean-z field.

To visualise the quadrupole pattern, we can project the residual field onto

spherical harmonics up to ℓ = 2, set all modes except ℓ = 2 to zero, and

reconstruct the corresponding sky map. This isolates the quadrupole component

and its axis.

4.5 Random mocks and significance

Random mocks are generated as follows:

1) Draw N_real random positions from random_G20.0_10x.fits (with replacement if needed).

2) Draw N_real redshift values by resampling the real Quaia redshift array (with replacement).

3) Build mean-z HEALPix maps, fit the dipole, and compute the residual map and Cℓ

exactly as for the real data.

From each mock we record:

dipole_frac  = |D| / z₀

D2_residual  = Dℓ at ℓ = 2 (quadrupole of residual map)

The distribution of these quantities across nmocks provides the expected level

of fluctuations from the mask, selection, and shot noise. Comparing the real

values to this distribution yields an approximate “σ-level” significance under

this null model.

5. Results

5.1 Global redshift dipole

For Quaia G < 20, NSIDE = 64, sel ≥ 0.3 × max(sel), and N_p ≥ 10, the dipole fit yields:

N_objects       ≈ 755 850

z_global_mean   ≈ 1.5135

|D|             ≈ 9.16 × 10⁻³

dipole_frac     ≈ 6.05 × 10⁻³

l_dip           ≈ 180.5° (Galactic longitude)

b_dip           ≈ 13.7° (Galactic latitude)

Thus the mean quasar redshift is about 0.6% higher in the direction

(l, b) ≈ (180.5°, 13.7°) (near the Galactic anti-centre) than in the opposite

direction.

5.2 Dipole significance from mocks

A suite of 50 mocks with the same NSIDE, selection cut, and N_p threshold yields:

⟨dipole_frac⟩_mocks ≈ 2.07 × 10⁻³

σ(dipole_frac)_mocks ≈ 8.89 × 10⁻⁴

The observed value then corresponds to

z_dip ≈ (6.05 × 10⁻³ − 2.07 × 10⁻³) / (8.89 × 10⁻⁴) ≈ 4.5 σ

Within the assumptions of the mock model, the redshift dipole is therefore a

highly significant departure from isotropy.

5.3 Residual quadrupole (full mask)

After dipole subtraction, the residual mean-z map shows strong power at ℓ = 2:

C2_real ≈ 4.86 × 10⁻³

D2_real ≈ 4.64 × 10⁻³

The mock suite gives

⟨D2_residual⟩_mocks ≈ 4.19 × 10⁻⁶

σ(D2_residual)_mocks ≈ 2.67 × 10⁻⁶

Formally,

z_D2 ≫ 1000 σ

This reflects the fact that the mock model produces a residual quadrupole

roughly two orders of magnitude smaller than the real one. The exact σ value is

not to be taken literally, but the qualitative conclusion is that the residual

quadrupole in Quaia is far larger than what these simple mocks can generate.

5.4 Residual quadrupole with |b| > 30°

To test whether the quadrupole is dominated by the Galactic plane, we impose

|b| ≥ 30° in both the real map and mocks. For the real residual map:

C2_real(|b|>30°) ≈ 2.00 × 10⁻⁴

D2_real(|b|>30°) ≈ 1.91 × 10⁻⁴

The mocks (same |b| cut) give:

⟨D2_residual⟩_mocks ≈ 2.13 × 10⁻⁶

σ(D2_residual)_mocks ≈ 1.40 × 10⁻⁶

So

z_D2(|b|>30°) ≈ (1.91 × 10⁻⁴ − 2.13 × 10⁻⁶) / (1.40 × 10⁻⁶) ≈ 135 σ

Even outside the Galactic plane, the quadrupole is still roughly 100× higher

than expected from these mocks.

5.5 Quadrupole axis

Using the ℓ = 2 component of the residual map, the direction of maximum |Δz|

is found to be:

(l_quad, b_quad) ≈ (185.6°, 84.2°)

This is very close to the North Galactic Pole, so the quadrupole pattern is

roughly a ring vs poles structure with a symmetry axis almost aligned with

the Galactic Z-axis.

Approximate angular separations:

sep(quad, z-dip)      ~ 70°

sep(quad, count dip)  ~ 44° (using the count dipole from Quaia counts)

sep(quad, CMB dipole) ~ 41°

sep(quad, GC)         ~ 96°

sep(quad, anti-GC)    ~ 84°

The quadrupole axis is therefore distinct from the dipole axis yet still part

of the broader zoo of large-scale preferred directions.

6. Relation to BAO DH/DM Anisotropy

Previous PLamb/FR/STF work has used BAO measurements to probe anisotropic

expansion, especially through the ratio DH/DM of radial to transverse BAO

distances. In some directions and redshift ranges, analyses have suggested

DH/DM ≈ 0.35, indicating suppressed radial expansion relative to transverse

modes along a preferred axis.

Qualitatively, the Quaia results fit into this picture as follows:

The redshift dipole in Quaia is a direct signature of directional

modulation in the distance–redshift relation.

The quadrupole indicates that the anisotropy is not purely dipolar: there

is a strong two-lobed component, likely encoding more complex angular

structure in the STF or in survey systematics.

Together with BAO DH/DM deviations and HHT-based IMF locking to BAO scales

in SN residuals, the Quaia anisotropies strengthen the case for a common

underlying large-scale structure (or symmetry breaking) that affects

distance measures across multiple tracers.

From a purely phenomenological point of view, Quaia provides an independent

quasar-based confirmation that the sky is not perfectly isotropic in redshift

space under simple ΛCDM + isotropic H(z) assumptions. Whether this eventually

maps onto a rotating STF, variable-c, or some more conventional set of

anisotropic or inhomogeneous models remains an open question.

7. Scripts Used (overview and sample calls)

quaia_to_csv.py

Extracts RA, DEC, redshift_quaia from quaia_G20.0.fits into a compact CSV.

Example:

%run /Users/boyde/.spyder-py3/quaia_to_csv.py \

--fits /Users/boyde/.spyder-py3/quaia_G20.0.fits \

--out  /Users/boyde/Downloads/quaia_G20p0_basic.csv

qso_add_galactic_coords.py

Adds Galactic coordinates to a CSV with RA, DEC, Z.

Example:

%run /Users/boyde/.spyder-py3/qso_add_galactic_coords.py \

--csv "/Users/boyde/Downloads/quaia_G20p0_basic.csv"

quaia_redshift_dipole.py

Builds a mean-z map and fits a redshift dipole (no selection map).

Example:

%run /Users/boyde/.spyder-py3/quaia_redshift_dipole.py \

--csv  "/Users/boyde/Downloads/quaia_G20p0_basic.csv" \

--nside 64 \

--min-per-pixel 10

quaia_redshift_dipole_sel.py

As above, but includes selection_function_NSIDE64_G20.0.fits and a

sel_min_frac threshold to exclude low-completeness pixels.

Example:

%run /Users/boyde/.spyder-py3/quaia_redshift_dipole_sel.py \

--csv          "/Users/boyde/Downloads/quaia_G20p0_basic.csv" \

--nside        64 \

--min-per-pixel 10 \

--sel-map      "/Users/boyde/.spyder-py3/selection_function_NSIDE64_G20.0.fits" \

--sel-min-frac 0.3

quaia_dipole_residual_map.py

Fits the redshift dipole, builds the dipole model map, subtracts it to get

the residual mean-z map, and computes the residual Cℓ spectrum up to lmax.

Example:

%run /Users/boyde/.spyder-py3/quaia_dipole_residual_map.py \

--csv              "/Users/boyde/Downloads/quaia_G20p0_basic.csv" \

--nside            64 \

--min-per-pixel    10 \

--sel-map          "/Users/boyde/.spyder-py3/selection_function_NSIDE64_G20.0.fits" \

--sel-min-frac     0.3 \

--out-dipole-map   "/Users/boyde/Downloads/quaia_G20p0_nside64_zdipole_model.fits" \

--out-residual-map "/Users/boyde/Downloads/quaia_G20p0_nside64_zresidual.fits" \

--lmax             10

quaia_random_mocks_zdipole_quad.py

Generates mock catalogues based on the Quaia random file, fits the redshift

dipole, and computes residual quadrupole power D2_residual for each mock,

with optional |b| > b_min_abs cut.

Example:

%run /Users/boyde/.spyder-py3/quaia_random_mocks_zdipole_quad.py \

--data-csv      "/Users/boyde/Downloads/quaia_G20p0_basic.csv" \

--random-fits   "/Users/boyde/.spyder-py3/random_G20.0_10x.fits" \

--sel-map       "/Users/boyde/.spyder-py3/selection_function_NSIDE64_G20.0.fits" \

--nside         64 \

--min-per-pixel 10 \

--sel-min-frac  0.3 \

--nmocks        50 \

--lmax          10 \

--bmin-abs      30.0 \

--seed          42 \

--outcsv        "/Users/boyde/Downloads/quaia_G20p0_random_mocks_b30_stats.csv"

8. Conclusions and Next Steps

The Quaia G < 20 catalogue exhibits a statistically robust redshift dipole

and an even stronger quadrupole in the residual mean-redshift field after

dipole removal. Under simple random + selection null models, both signals are

strong outliers, especially the quadrupole, which remains significant even

outside the Galactic plane.

From a PLamb/FR/STF perspective this provides another piece of evidence that

the large-scale STF may not be isotropic, and that radial vs transverse

distance measures (BAO DH/DM, quasar ⟨z⟩, SN residuals) all show related

anisotropic structure. From a survey-analysis perspective it motivates the

construction of more sophisticated mocks and systematics models to test

whether any combination of calibration, extinction, and completeness effects

could reproduce the observed low-ℓ structure.

Further work could include:

Extending the analysis to Quaia G < 20.5 using the corresponding random and

selection files.

Incorporating more realistic systematics into the mock pipeline.

Jointly fitting BAO, SN, and Quaia anisotropy within the FR/STF framework

to test variable-c or rotating-STF hypotheses quantitatively.
