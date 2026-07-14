# Quaia FR Quaia BAO Report 0V2

- Source: `Quaia_FR_Quaia_BAO_Report_0V2.rtf`
- Source size: 73595 bytes
- Source modified: 2025-11-18T10:20:46
- Extracted: 2026-07-14
- Word count estimate: 3758

## Extracted Text
* 02020603050405020304 Times New Roman; * 02040503050406030204 Cambria Math;

* 020b0004020202020204 Aptos; * 02020603050405020304 Times New Roman;

* 02020603050405020304 Times New Roman; * 020b0004020202020204 Aptos Display;

* 02020603050405020304 Times New Roman; * 02020603050405020304 Times New Roman;

* 02020603050405020304 Times New Roman; * 020b0004020202020204 Aptos;

* 02020603050405020304 Times New Roman; Times New Roman CE; Times New Roman Cyr;

Times New Roman Greek; Times New Roman Tur; Times New Roman (Hebrew); Times New Roman (Arabic);

Times New Roman Baltic; Times New Roman (Vietnamese); Cambria Math CE; Cambria Math Cyr;

Cambria Math Greek; Cambria Math Tur; Cambria Math Baltic; Cambria Math (Vietnamese);

Aptos CE; Aptos Cyr; Aptos Greek; Aptos Tur; Aptos Baltic;

Aptos (Vietnamese); Times New Roman CE; Times New Roman Cyr;

Times New Roman Greek; Times New Roman Tur; Times New Roman (Hebrew);

Times New Roman (Arabic); Times New Roman Baltic; Times New Roman (Vietnamese);

Times New Roman CE; Times New Roman Cyr; Times New Roman Greek;

Times New Roman Tur; Times New Roman (Hebrew); Times New Roman (Arabic);

Times New Roman Baltic; Times New Roman (Vietnamese); Aptos Display CE;

Aptos Display Cyr; Aptos Display Greek; Aptos Display Tur;

Aptos Display Baltic; Aptos Display (Vietnamese); Times New Roman CE;

Times New Roman Cyr; Times New Roman Greek; Times New Roman Tur;

Times New Roman (Hebrew); Times New Roman (Arabic); Times New Roman Baltic;

Times New Roman (Vietnamese); Times New Roman CE; Times New Roman Cyr;

Times New Roman Greek; Times New Roman Tur; Times New Roman (Hebrew);

Times New Roman (Arabic); Times New Roman Baltic; Times New Roman (Vietnamese);

Times New Roman CE; Times New Roman Cyr; Times New Roman Greek;

Times New Roman Tur; Times New Roman (Hebrew); Times New Roman (Arabic);

Times New Roman Baltic; Times New Roman (Vietnamese); Aptos CE;

Aptos Cyr; Aptos Greek; Aptos Tur; Aptos Baltic;

Aptos (Vietnamese); Times New Roman CE; Times New Roman Cyr;

Times New Roman Greek; Times New Roman Tur; Times New Roman (Hebrew);

Times New Roman (Arabic); Times New Roman Baltic; Times New Roman (Vietnamese);

; ; ; ; ; ; ; ; ; ; ; ;

; ; ; ; ; ; ; * *

Normal; *

Default Paragraph Font; *

Normal Table; *

Clive Boyd * http://schemas.microsoft.com/office/word/2

003/wordml

* 2450 * . *

. * . * ) * ( ) *

( ) * ( ) * ( ) *

( )

Quaia Redshift Anisotropy & FR Cosmology Fits Updated Results (0V2)

EXECUTIVE SUMMARY

This document summarises the updated analysis of large-scale anisotropy in the Quaia G<20.0 quasar catalogue and its connection to FR (Full Relativity) cosmological fits using DESI DR2 BAO and Pantheon+SH0ES data. It extends the previous 0V1 report by:

Including the full command and processing pipeline used to obtain the Quaia G<20 selection.

Adding the selection corrected redshift-dipole fits and redshift-binned dipoles.

Describing the residual quadrupole in the mean-redshift map and its significance versus randomised mocks, including a |b|>30 mask variant.

Summarising new FR cosmology fits to DESI DR2 BAO alone, with different c(z) coupling choices and Noether-link settings.

The key result is that the Quaia mean-redshift field exhibits a statistically extreme large-scale quadrupole (even after removing a best-fit dipole and masking the Galactic plane), while BAO-only FR fits point to non-standard combinations of H0,

m and a n additional acceleration sector with moderate per degree of freedom.

1. ACRONYMS, SYMBOLS, AND BASIC DEFINITIONS

GC: Galactic Centre (l = 0 , b = 0 ).

anti-GC: Galactic anti-centre (l = 180 , b = 0 ).

CMB dipole: Direction of the CMB temperature dipole (l 264 , b 48 ).

Quaia: Gaia unWISE quasar catalogue (version 0.1.0, DOI 10.5281/zenodo.8060755).

DESI DR2 BAO: Baryon Acoustic Oscillation distance measurements from DESI Data Release

2.

FR: Full Relativity cosmological framework with an additional acceleration term A_acc and exponent n_acc.

H0: Hubble constant at z = 0 (km s

Mpc

).

m, ,

k: Present-day matter, dark-energy and curvature density parameters.

A_acc, n_acc: Amplitude and exponent of the extra acceleration sector in FR.

c,

M: FR parameters controlling effective c(z) behaviour and luminosity/magnitude evolution, linked via the soft Noether constraint.

w0, wa: CPL dark-energy equation-of-state parameters.

l, b: Galactic longitude and latitude.

z: Residual mean redshift after subtracting the fitted dipole model.

C , D : Angular power spectrum and its rescaled form D

= ( +1)C / 2 .

D2: Quadrupole power at = 2, expressed as D2 = 2 3 C2 / 2

.

-significance: (measured mean(mock)) / std(mock) based on random mocks.

2. DATASETS AND PRE-PROCESSING

2.1 Quaia G<20.0 Sample

The Quaia catalogue (DOI 10.5281/zenodo.8060755) provides a Gaia

unWISE quasar sample with improved photometric redshifts. For this analysis, the G<20.0 subset is used, with the following main columns: source_id, redshift_quaia, ra, dec, l, b, and photo metric magnitudes.

Key conversion steps:

1. Download the full Quaia suite via download_zenodo_quaia.py.

2. Extract the G<20.0 science catalogue (quaia_G20.0.fits).

3. Convert to a compact CSV via quaia_to_csv.py, selecting RA, DEC, and redshift_quaia:

O utput: quaia_G20p0_basic.csv (columns: RA, DEC, Z).

4. Optionally add Galactic coordinates (L_GAL, B_GAL) via qso_add_galactic_coords.py.

Basic statistics for the Quaia G<20.0 sample:

N 755,850 objects.

RA [0 , 360 ), DEC [ 90

, 90 ].

Z [0.084, 4.537] with median z 1.45 1.50.

2.2 Random Catalogue and Selection Function

The Quaia release includes random_G20.0_10x.fits and selection_function_NSIDE64_G20.0.fits.

The random catalogue approximates the angular and magnitude selection; the selection map encodes the relative completeness per HEALPix pixel at nside = 64.

These are used to:

Build randomised mock realisations with the same angular selection as the data.

Define an analysis mask by requiring sel f max(sel) with f [0.2, 0.5], typically f = 0.3.

2.3 DESI DR2 BAO and Pantheon+SH0ES

The FR cosmology runs use:

DESI_DR2_BAO.csv and DESI_DR2_BAO_cov.txt (13 BAO data points).

Optionally Pantheon+SH0ES SN Ia (Pantheon+SH0ES.dat with the full STAT+SYS covariance) for joint FR fits. One such joint run was attempted but halted when emcee encountered infinite parameter values for aggressive link-strength settings.

3. QUAIA COUNTS AND REDSHIFT DIPOLES

3.1 Number-Count Dipole

The initial analysis computed a simple number-count dipole in Galactic coordinates by summing unit vectors for all quasars.

Global dipole (all 0.08 z 4.5):

N = 755,850, amplitude 2.83 10 .

Direction (l, b) (231.1 , 41.8 ).

Redshift-binned number-count dipoles (z bins: [0,1], [1,2], [2,4]) show broadly similar orientations with modest variation in amplitude, reinforcing a mild number-count anisotropy.

3.2 Mean-Redshift Dipole (No Selection Function)

A HEALPix map of the me an redshift z was constructed at nside = 64, requiring a minimum per-pixel count.

Full redshift range, no selection function:

z 1.53, |D| 4.18 10

, fractional amplitude 2.74 10 .

(l, b) (191.4 , 8.8 ).

3.3 Selection-Function Corrected Redshift Dipole

Including the Quaia selection map at nside = 64, with sel 0.3 max(sel) and 10 objects per pixel:

z 1.5135, |D| 9.16 10

, fractional amplitude 6.05 10 .

Direction (l, b) (180.5 , 13.7 ).

Varying sel-min-frac between 0.2 and 0.5 leaves the dipole direction broadly stable (l 175 185 ).

3.4 Redshift-Binned Dipoles with Selection Function

Selection-function corrected redshift dipoles in three z bins:

0.0 z < 1.0:

z 0.667, |D| 3.0 10 , frac amplitude

4.5 10 ,

(l, b) (172 , 10 ).

1.0 z < 2.0:

z 1.48, |D| 2.6 10 , frac amplitude

1.7 10 ,

(l, b) (166 , 26 ).

2.0 z < 4.0:

z 2.50, |D| 1.07 10 , frac amplitude

4.3 10 ,

(l, b) (102 , 33 ).

The low- and mid-z bins maintain an anti-centre like direction; the high-z bin rotates significantly.

4. RESIDUAL QUADRUPOLE IN MEAN-z AND RANDOM MOCKS

4.1 Dipole Removal and Residual Map

Using quaia_dipole_residual_map.py, the best-fit redshift dipole was subtracted to form a residual map

z(n ) = z (n ) z_dipole_model(n ), with hp.UNSEEN in unobserved pixels.

The angular power spectrum C of z up to max = 10 gives, for the full masked sky:

C 4.86 10 , D

4.64 10 .

4.2 Quadrupole-Only Map and Axis

A pure quadrupole map was constructed via map2alm up to =2, zeroing =0,1 and reconstructing via alm2map.

Max | z| in the quadrupole map :

(l, b) (185.6 , 84.2 ) (near the Galactic north pole).

4.3 Random Mocks and Significance (No |b| Cut)

quaia_random_mocks_zdipole_quad.py built 50 random realisations using random_G20.0_10x and the same selection mask.

Mock ensemble:

dipole_frac (mock): mean 2.07 10 , std 8.89

10 .

D _residual (mock): mean 4.19 10 , std

2.67 10 .

Real map (same mask, no |b| cut):

dipole_frac_real 6.05 10

4.5 above the mock mean.

D _real 4.64 10 formally

1.7 10 above the mock mean.

This huge D significance shows that the residual quadrupole is not explained by the random catalogue alone.

4.4 |b|>30 Cut and Revised Quadrupole Significa nce

Imposing an additional |b| 30 cut:

Real map:

C _cut 2.00 10 , D _cut

1.91 10 .

50 mocks with identical |b|-cut:

D _residual (mock, b 30 ): mean 2.13

10 , std 1.40 10 .

z_D2_b30 135 .

Thus even outside the Galactic plane, the residual quadrupole is extraordinarily large.

5. UPDATED FR COSMOLOGY FITS WITH DESI DR2 BAO

5.1 FR + DESI DR2 BAO, c_of_z, Soft Link

Command (schematic):

PLamb_Test_10V6c_plus.py --out plamb_runs/desi_tests/FR_w0_softlink_sigma0p6_desionly --model FR --flat --de-model cpl --bao DESI_DR2_BAO.csv --bao-cov DESI_DR2_BAO_cov.txt --rd 147.09 --bao-c-mode c_of_z --c-model log1pz --link-

soft --kappa-em-gc 1.0 --link-sigma 0.60 --noether combo_soft --noether-strength 1.0 ...

Active parameters: H0, m, A_acc, n_acc, c, M, w0 with wa fixed = 0.

Best-fit (approx):

H 0 50.9 km s

Mpc

.

m 0.134,

0.866 (flat FR).

A_acc 0.71, n_acc 1.20.

c 0.30, M 0.49.

w_eff 0.60.

Goodness-of-fit:

_BAO 6.73 for dof 6

/dof 1.12.

Link + Noether penalty 4.0 in .

The chain is mildly under-sampled (steps_post/ _max 25.6) but indicates that DESI BAO alone prefer a low-H0, moderate-

m FR solution with non-trivial extra acceleration.

5.2 FR + DESI DR2 BAO, c0 BAO Mode, Stronger Noether Link

Command (schemat ic):

PLamb_Test_10V6c_plus.py --out plamb_runs/desi_tests/FR_CPL_softlink_c0_desionly_v2 --model FR --flat --de-model cpl --bao DESI_DR2_BAO.csv --bao-cov DESI_DR2_BAO_cov.txt --rd 147.09 --bao-c-mode c0 --c-model log1pz --link-soft -

-kappa-em-gc 1.0 --link-sigma 0.25 --noether combo_soft --noether-strength 1.5 ...

Best-fit summary:

H0 58.9 km s

Mpc

.

m 0.40,

0.60.

A_acc 0.34, n_acc 0.54.

c 0.0024, M 0.223.

w_ eff 0.82.

Goodness-of-fit:

_BAO 9.18 for dof 5

/dof 1.84.

Link residual M c 0.22

_link 0.78; Noether 1.04; total extra 1.82.

This run pushes towards higher H0 and m with a more CDM-like effective w but a worse /dof than the c_of_z run.

5.3 Joint SN+BAO FR Run and Instability

A joint FR run including Pantheon+SH0ES SN + DESI DR2 BAO in a full-covariance setup, with strong link-soft and Noether-strength, halted when emcee reported infinite parameter values, signalling numerical or prior issues. This suggests that the combined d

ata plus strong symmetry-breaking priors tightly constrain the allowed FR parameter region and ma y require re-tuned priors or reparametrisation.

6. RELATION TO PREVIOUS BAO/HHT RESULTS

Previous HHT_Locking and BAO-residual work identified characteristic BAO scales and potential quasi-periodic residuals in SN Ia data. The current Quaia analysis adds a complementary sky-map perspective:

Instead of SN residuals, it examines the mean-redshift field of a deep quasar survey.

The extremely large quadrupole in z , after dipole removal and |b|>30

masking, hints at either a major systematic or a real large-scale anisotropy.

If physical, such a quadrupole could be linked to the same underlying acceleration or c(z) physics probed in FR/HHT analyses, potentially tied to large-scale structure or symmetry-breaking effects.

7. SCRIPTS AND PIPELINES USED

Main scripts:

download_zenodo_quaia.py download Quaia FITS and auxiliary files.

quaia_to_csv.py convert quaia_G20.0.fits to RA/DEC/Z CSV.

qso_add_galactic_coords.py add Galactic (l, b) to CSV.

qso_dipole_catalog.py co mpute number-count dipoles and count maps.

quaia_redshift_dipole.py build mean-z maps and fit redshift dipoles.

quaia_redshift_dipole_sel.py as above, but with selection map masking.

quaia_dipole_residual_map.py produce dipole model and residual mean-z maps; compute C .

quaia_random_mocks_zdipole_quad.py generate random mocks and dipole/quadrupole statistics.

AngularSeperation.py compute angular separations between axes (Quaia, GC, CMB, etc.).

PLamb_Test_10V6c_p lus.py FR cosmology engine for SN+BAO with flexible c(z), (z), and Noether options.

8. CONCLUSIONS AND NEXT STEPS

1. The Quaia G<20.0 catalogue exhibits a robust mean-redshift dipole aligned roughly with the Galactic anti-centre but at positive latitude, with fractional amplitude 6

10 under selection masking.

2. After dipole removal, the residual mean-z field shows an extraordinarily large quadrupole. D is orders of magnitude above random-mock expectations, even for |b|>30

, suggesting either a strong unmodelled systematic or a genuine cosmic anisotropy.

3. FR BAO-only fits using DESI DR2 favour non-standard combinations of H0, m, and A_acc with reasonable

/dof, and are sensitive to the assumed c(z) model and Noether-link strength.

4. Fully joint FR fits with SN+BAO and strong symmetry-breaking priors are numerically delicate, indicating highly constrained or unstable directions in parameter space.

5. Combining Quaia anisotropy, BAO distance constraints, and SN/HHT residual structure offers a promising multi-probe framework to test FR-style extensions and cosmic anisotropy hypotheses.

Planned next steps:

Extend the mock analysis (more mocks, alternative selection-function treatments).

Explore redshift-dependent masks and cross-correlation of Quaia z with large-scale structure tracers.

Refi ne FR joint fits (SN+BAO+Planck) with more robust priors and lower-dimensional parameterisations.

Compare orientations of Quaia quadrupole/dipole axes with BAO- and SN-derived anisotropy axes from previous HHT work.

* 504b030414000600080000002100e9de0fbfff0000001c020000130000005b436f6e74656e745f54797065735d2e786d6cac91cb4ec3301045f748fc83e52d4a

9cb2400825e982c78ec7a27cc0c8992416c9d8b2a755fbf74cd25442a820166c2cd933f79e3be372bd1f07b5c3989ca74aaff2422b24eb1b475da5df374fd9ad

5689811a183c61a50f98f4babebc2837878049899a52a57be670674cb23d8e90721f90a4d2fa3802cb35762680fd800ecd7551dc18eb899138e3c943d7e503b6

b01d583deee5f99824e290b4ba3f364eac4a430883b3c092d4eca8f946c916422ecab927f52ea42b89a1cd59c254f919b0e85e6535d135a8de20f20b8c12c3b0

0c895fcf6720192de6bf3b9e89ecdbd6596cbcdd8eb28e7c365ecc4ec1ff1460f53fe813d3cc7f5b7f020000ffff0300504b030414000600080000002100a5d6

a7e7c0000000360100000b0000005f72656c732f2e72656c73848fcf6ac3300c87ef85bd83d17d51d2c31825762fa590432fa37d00e1287f68221bdb1bebdb4f

c7060abb0884a4eff7a93dfeae8bf9e194e720169aaa06c3e2433fcb68e1763dbf7f82c985a4a725085b787086a37bdbb55fbc50d1a33ccd311ba548b6309512

0f88d94fbc52ae4264d1c910d24a45db3462247fa791715fd71f989e19e0364cd3f51652d73760ae8fa8c9ffb3c330cc9e4fc17faf2ce545046e37944c69e462

a1a82fe353bd90a865aad41ed0b5b8f9d6fd010000ffff0300504b0304140006000800000021006b799616830000008a0000001c0000007468656d652f746865

6d652f7468656d654d616e616765722e786d6c0ccc4d0ac3201040e17da17790d93763bb284562b2cbaebbf600439c1a41c7a0d29fdbd7e5e38337cedf14d59b

4b0d592c9c070d8a65cd2e88b7f07c2ca71ba8da481cc52c6ce1c715e6e97818c9b48d13df49c873517d23d59085adb5dd20d6b52bd521ef2cdd5eb9246a3d8b

4757e8d3f729e245eb2b260a0238fd010000ffff0300504b030414000600080000002100d0557692fa0700000d220000160000007468656d652f7468656d652f

7468656d65312e786d6cec5a4b8f1bb911be07c87f68f45d5677eb3db0bcd0d3b3f68c3db064077ba4244a4d0fbb2934a99911160b04de532e0b2cb00972c802

b9e5100459200b64914b7e8c011bc9e647a448b65aa444791e30022398994b37f555f16355b1aa9add0f3fbb4aa87781334e58daf6c30781efe174ca66245db4

fd97e361a9e97b5ca07486284b71db5f63ee7ff6e897bf78888e448c13ec817cca8f50db8f85581e95cb7c0ac3883f604b9cc26f73962548c06db628cf327409

7a135a8e82a05e4e10497d2f4509a87d3e9f9329f6c652a5ff68a37c40e136155c0e4c693692aab125a1b0b3f35022f89af768e65d20daf6619e19bb1ce32be1

7b1471013fb4fd40fdf9e5470fcbe82817a2e280ac2137547fb95c2e303b8fd49cd962524c1a0ca266352cf42b0015fbb84153fe17fa14004da7b052cdc5d419

d6ea4133cab106485f3a74b71a61c5c61bfa2b7b9cc356bd1b552dfd0aa4f557f7f0c1b035e8d72cbc02697c6d0fdf09a26eab62e11548e3eb7bf8eaa0d38806

165e81624ad2f37d74bdd16cd673740199337aec84b7eaf5a0d1cfe15b144443115d728a394bc5a1584bd06b960d0120811409927a62bdc473348528ee2c05e3

5e9ff025456bdf5ba29471180ea23084d0ab0651f1af2c8e8e3032a4252f60c2f786241f8f4f33b2146dff0968f50dc8bb9f7e7afbe6c7b76ffefef6ebafdfbe

f9ab774216b1d0aa2cb963942e4cb99ffff4ed7fbeffb5f7efbffdf1e7ef7eebc67313fffe2fbf79ff8f7f7e483d6cb5ad29defdee87f73ffef0eef7dffcebcf

df39b477323431e1639260ee3dc397de0b96c00295296cfe7892dd4e621c23624a74d205472992b338f40f446ca19fad11450e5c17db767c9541aa71011faf5e

5b844771b612c4a1f1699c58c053c66897654e2b3c957319661eafd2857bf26c65e25e2074e19abb8752cbcb83d512722c71a9ecc5d8a27946512ad002a75878

f237768eb163755f1062d9f5944c33c6d95c785f10af8b88d3246332b1a2692b744c12f0cbda4510fc6dd9e6f495d765d4b5ea3ebeb091b0371075901f636a99

f1315a0994b8548e51424d839f2011bb488ed6d9d4c40db8004f2f3065de60863977c93ccf60bd86d39f22c86e4eb79fd27562233341ce5d3a4f106326b2cfce

7b314a962eec88a4b189fd9c9f438822ef8c0917fc94d93b44de831f507ad0ddaf08b6dc7d7d36780959cea4b40d10f9cb2a73f8f2316656fc8ed6748eb02bd5

74b2c44ab19d8c38a3a3bb5a58a17d823145976886b1f7f27307832e5b5a36df927e12435639c6aec07a82ec5895f729e6d02bc9e6663f4f9e106e85ec082fd8

013ea7eb9dc4b3466982b2439a9f81d74d9b0f26196c460785e7747a6e029f11e801215e9c4679ce418711dc07b59ec5c82a60f29ebbe3759d59febbc91e837d

f9daa271837d0932f8d63290d84d990fda668ca835c13660c6887827ae740b2296fbb722b2b82ab195536e6e6fdaad1ba03bb29a9e84a4d77440ffbbce07fa8b

777ff8de11821fa7db712bb652d52dfb9c43a9e478a7bb3984dbed697a2c9b914fbfa5e9a3557a86a18aece7abfb8ee6bea3f1ffef3b9a43fbf9be8f39d46ddc

f7313ef417f77d4c7eb4f271fa986deb025d8d3c5ed0c73cead0273978e63327948ec49ae213ae8e7d383ccdcc863028e5d479272ece0097315cca32071358b8

4586948c9731f12b22e2518c96703614fa52c982e7aa17dc5b320e47466ad8a95be2e92a3965337dd4a9ce96025d593912dbf1a006874e7a1c8ea98446d71bf9

a0e4a7ce5381af62bb50c7ac1b0252f636248cc96c12150789c666f01a12f2d4ece3b068395834a5fa8dabf64c01d40aafc0e3b6070fe96dbf569584e08c9c4f

a1359f493f69576fbcab9cf9313d7dc8985604c0b1a25e091cca179e6e49ae07972757a743ed069eb64828a7e8b0b24928cba8068fc7f0109c47a71cbd098ddb

fabab575a9454f9a42cd07f1bda5d1687e88c55d7d0d72bbb981a666a6a0a977097b3c824de77b53b46cfb73383386cb6409c1c3e52317a20b78f1321599def1

77492dcb8c8b3ee2b1b6b8ca3ada3f091138f32849dabe5c7fe1079aaa24a2c9b560eb7eaae422b9e13e3572e075dbcb783ec75361fadd189196d6b790e275b2

70feaac4ef0e96926c05ee1ec5b34b6f4257d90b0421566b84d2bb33c2e1d541a85d3d23f02eacc864dbf8dba94c79f6375f46a918d2e3882e63949714339b6b

b82a28051d7557d8c0b8cbd70c06354c9257c2c9425658d3a856392d6a97e670b0ec5e2f242d6764cd6dd1b4d28a2c9bee3466cdb0a9033bb6bc5b9537586d4c

0c49cd2cf13a77efe6dcd626d9ed340a4599008317f6bb5bed37a86d27b3a849c6fb795826ed7cd42e1e9b055e43ed2655c248fbf58dda1dbb1545c2391d0cde

a9f483dc6ed4c2d07cd3582a4bab97e6e67b6d36790dc9a30f6dee8aea37dd34853b19957c799629df4ed86c9d5f52ae138df6b96c4a2592a62ff0dc23b3abb6

1fb93a47fdb635ccbb01859662b2781582ce6ecf16ccf152546fd84258077811547a53dac285849a197aef42589d28ba688bab0d65d9ab035e9990eb55836973

4bc1d5be15e1743c43d0db8e5467a7732fd0be12797e812b6f9591b6ff6550eb547b51ad570a9ab541a95aa906a566ad5329756ab54a38a88541bf1b7d05f444

9c8435fdd5c3105e02d175feed831adffbfe21d9bce77a30654999a9ef1bcacafbeafb87303afcfd033812684583b01a75a25ea9d70feba56ad4af979a8d4aa7

d48beafda80345bb3eec7ce57b170a1c76fbfde1b01695ea3dc055834eadd4e9567aa57a73d08d86e1a0da0f009c979f2b788a019b6d6c01978ad7a3ff020000

ffff0300504b0304140006000800000021000dd1909fb60000001b010000270000007468656d652f7468656d652f5f72656c732f7468656d654d616e61676572

2e786d6c2e72656c73848f4d0ac2301484f78277086f6fd3ba109126dd88d0add40384e4350d363f2451eced0dae2c082e8761be9969bb979dc9136332de3168

aa1a083ae995719ac16db8ec8e4052164e89d93b64b060828e6f37ed1567914b284d262452282e3198720e274a939cd08a54f980ae38a38f56e422a3a641c8bb

d048f7757da0f19b017cc524bd62107bd5001996509affb3fd381a89672f1f165dfe514173d9850528a2c6cce0239baa4c04ca5bbabac4df000000ffff030050

4b01022d0014000600080000002100e9de0fbfff0000001c0200001300000000000000000000000000000000005b436f6e74656e745f54797065735d2e786d6c

504b01022d0014000600080000002100a5d6a7e7c0000000360100000b00000000000000000000000000300100005f72656c732f2e72656c73504b01022d0014

0006000800000021006b799616830000008a0000001c00000000000000000000000000190200007468656d652f7468656d652f7468656d654d616e616765722e

786d6c504b01022d0014000600080000002100d0557692fa0700000d2200001600000000000000000000000000d60200007468656d652f7468656d652f746865

6d65312e786d6c504b01022d00140006000800000021000dd1909fb60000001b0100002700000000000000000000000000040b00007468656d652f7468656d652f5f72656c732f7468656d654d616e616765722e786d6c2e72656c73504b050600000000050005005d010000ff0b00000000

* 3c3f786d6c2076657273696f6e3d22312e302220656e636f64696e673d225554462d3822207374616e64616c6f6e653d22796573223f3e0d0a3c613a636c724d

617020786d6c6e733a613d22687474703a2f2f736368656d61732e6f70656e786d6c666f726d6174732e6f72672f64726177696e676d6c2f323030362f6d6169

6e22206267313d226c743122207478313d22646b3122206267323d226c743222207478323d22646b322220616363656e74313d22616363656e74312220616363

656e74323d22616363656e74322220616363656e74333d22616363656e74332220616363656e74343d22616363656e74342220616363656e74353d22616363656e74352220616363656e74363d22616363656e74362220686c696e6b3d22686c696e6b2220666f6c486c696e6b3d22666f6c486c696e6b222f3e

* Normal; heading 1;

heading 2; heading 3; heading 4;

heading 5; heading 6; heading 7;

heading 8; heading 9; index 1;

index 2; index 3; index 4; index 5;

index 6; index 7; index 8; index 9;

toc 1; toc 2; toc 3;

toc 4; toc 5; toc 6;

toc 7; toc 8; toc 9; Normal Indent;

footnote text; annotation text; header; footer;

index heading; caption; table of figures;

envelope address; envelope return; footnote reference; annotation reference;

line number; page number; endnote reference; endnote text;

table of authorities; macro; toa heading; List;

List Bullet; List Number; List 2; List 3;

List 4; List 5; List Bullet 2; List Bullet 3;

List Bullet 4; List Bullet 5; List Number 2; List Number 3;

List Number 4; List Number 5; Title; Closing;

Signature; Default Paragraph Font; Body Text; Body Text Indent;

List Continue; List Continue 2; List Continue 3; List Continue 4;

List Continue 5; Message Header; Subtitle; Salutation;

Date; Body Text First Indent; Body Text First Indent 2; Note Heading;

Body Text 2; Body Text 3; Body Text Indent 2; Body Text Indent 3;

Block Text; Hyperlink; FollowedHyperlink; Strong;

Emphasis; Document Map; Plain Text; E-mail Signature;

HTML Top of Form; HTML Bottom of Form; Normal (Web); HTML Acronym;

HTML Address; HTML Cite; HTML Code; HTML Definition;

HTML Keyboard; HTML Preformatted; HTML Sample; HTML Typewriter;

HTML Variable; Normal Table; annotation subject; No List;

Outline List 1; Outline List 2; Outline List 3; Table Simple 1;

Table Simple 2; Table Simple 3; Table Classic 1; Table Classic 2;

Table Classic 3; Table Classic 4; Table Colorful 1; Table Colorful 2;

Table Colorful 3; Table Columns 1; Table Columns 2; Table Columns 3;

Table Columns 4; Table Columns 5; Table Grid 1; Table Grid 2;

Table Grid 3; Table Grid 4; Table Grid 5; Table Grid 6;

Table Grid 7; Table Grid 8; Table List 1; Table List 2;

Table List 3; Table List 4; Table List 5; Table List 6;

Table List 7; Table List 8; Table 3D effects 1; Table 3D effects 2;

Table 3D effects 3; Table Contemporary; Table Elegant; Table Professional;

Table Subtle 1; Table Subtle 2; Table Web 1; Table Web 2;

Table Web 3; Balloon Text; Table Grid; Table Theme; Placeholder Text;

No Spacing; Light Shading; Light List; Light Grid; Medium Shading 1; Medium Shading 2;

Medium List 1; Medium List 2; Medium Grid 1; Medium Grid 2; Medium Grid 3; Dark List;

Colorful Shading; Colorful List; Colorful Grid; Light Shading Accent 1; Light List Accent 1;

Light Grid Accent 1; Medium Shading 1 Accent 1; Medium Shading 2 Accent 1; Medium List 1 Accent 1; Revision;

List Paragraph; Quote; Intense Quote; Medium List 2 Accent 1; Medium Grid 1 Accent 1;

Medium Grid 2 Accent 1; Medium Grid 3 Accent 1; Dark List Accent 1; Colorful Shading Accent 1; Colorful List Accent 1;

Colorful Grid Accent 1; Light Shading Accent 2; Light List Accent 2; Light Grid Accent 2; Medium Shading 1 Accent 2;

Medium Shading 2 Accent 2; Medium List 1 Accent 2; Medium List 2 Accent 2; Medium Grid 1 Accent 2; Medium Grid 2 Accent 2;

Medium Grid 3 Accent 2; Dark List Accent 2; Colorful Shading Accent 2; Colorful List Accent 2; Colorful Grid Accent 2;

Light Shading Accent 3; Light List Accent 3; Light Grid Accent 3; Medium Shading 1 Accent 3; Medium Shading 2 Accent 3;

Medium List 1 Accent 3; Medium List 2 Accent 3; Medium Grid 1 Accent 3; Medium Grid 2 Accent 3; Medium Grid 3 Accent 3;

Dark List Accent 3; Colorful Shading Accent 3; Colorful List Accent 3; Colorful Grid Accent 3; Light Shading Accent 4;

Light List Accent 4; Light Grid Accent 4; Medium Shading 1 Accent 4; Medium Shading 2 Accent 4; Medium List 1 Accent 4;

Medium List 2 Accent 4; Medium Grid 1 Accent 4; Medium Grid 2 Accent 4; Medium Grid 3 Accent 4; Dark List Accent 4;

Colorful Shading Accent 4; Colorful List Accent 4; Colorful Grid Accent 4; Light Shading Accent 5; Light List Accent 5;

Light Grid Accent 5; Medium Shading 1 Accent 5; Medium Shading 2 Accent 5; Medium List 1 Accent 5; Medium List 2 Accent 5;

Medium Grid 1 Accent 5; Medium Grid 2 Accent 5; Medium Grid 3 Accent 5; Dark List Accent 5; Colorful Shading Accent 5;

Colorful List Accent 5; Colorful Grid Accent 5; Light Shading Accent 6; Light List Accent 6; Light Grid Accent 6;

Medium Shading 1 Accent 6; Medium Shading 2 Accent 6; Medium List 1 Accent 6; Medium List 2 Accent 6;

Medium Grid 1 Accent 6; Medium Grid 2 Accent 6; Medium Grid 3 Accent 6; Dark List Accent 6; Colorful Shading Accent 6;

Colorful List Accent 6; Colorful Grid Accent 6; Subtle Emphasis; Intense Emphasis;

Subtle Reference; Intense Reference; Book Title; Bibliography;

TOC Heading; Plain Table 1; Plain Table 2; Plain Table 3; Plain Table 4;

Plain Table 5; Grid Table Light; Grid Table 1 Light; Grid Table 2; Grid Table 3; Grid Table 4;

Grid Table 5 Dark; Grid Table 6 Colorful; Grid Table 7 Colorful; Grid Table 1 Light Accent 1; Grid Table 2 Accent 1;

Grid Table 3 Accent 1; Grid Table 4 Accent 1; Grid Table 5 Dark Accent 1; Grid Table 6 Colorful Accent 1;

Grid Table 7 Colorful Accent 1; Grid Table 1 Light Accent 2; Grid Table 2 Accent 2; Grid Table 3 Accent 2;

Grid Table 4 Accent 2; Grid Table 5 Dark Accent 2; Grid Table 6 Colorful Accent 2; Grid Table 7 Colorful Accent 2;

Grid Table 1 Light Accent 3; Grid Table 2 Accent 3; Grid Table 3 Accent 3; Grid Table 4 Accent 3;

Grid Table 5 Dark Accent 3; Grid Table 6 Colorful Accent 3; Grid Table 7 Colorful Accent 3; Grid Table 1 Light Accent 4;

Grid Table 2 Accent 4; Grid Table 3 Accent 4; Grid Table 4 Accent 4; Grid Table 5 Dark Accent 4;

Grid Table 6 Colorful Accent 4; Grid Table 7 Colorful Accent 4; Grid Table 1 Light Accent 5; Grid Table 2 Accent 5;

Grid Table 3 Accent 5; Grid Table 4 Accent 5; Grid Table 5 Dark Accent 5; Grid Table 6 Colorful Accent 5;

Grid Table 7 Colorful Accent 5; Grid Table 1 Light Accent 6; Grid Table 2 Accent 6; Grid Table 3 Accent 6;

Grid Table 4 Accent 6; Grid Table 5 Dark Accent 6; Grid Table 6 Colorful Accent 6; Grid Table 7 Colorful Accent 6;

List Table 1 Light; List Table 2; List Table 3; List Table 4; List Table 5 Dark;

List Table 6 Colorful; List Table 7 Colorful; List Table 1 Light Accent 1; List Table 2 Accent 1; List Table 3 Accent 1;

List Table 4 Accent 1; List Table 5 Dark Accent 1; List Table 6 Colorful Accent 1; List Table 7 Colorful Accent 1;

List Table 1 Light Accent 2; List Table 2 Accent 2; List Table 3 Accent 2; List Table 4 Accent 2;

List Table 5 Dark Accent 2; List Table 6 Colorful Accent 2; List Table 7 Colorful Accent 2; List Table 1 Light Accent 3;

List Table 2 Accent 3; List Table 3 Accent 3; List Table 4 Accent 3; List Table 5 Dark Accent 3;

List Table 6 Colorful Accent 3; List Table 7 Colorful Accent 3; List Table 1 Light Accent 4; List Table 2 Accent 4;

List Table 3 Accent 4; List Table 4 Accent 4; List Table 5 Dark Accent 4; List Table 6 Colorful Accent 4;

List Table 7 Colorful Accent 4; List Table 1 Light Accent 5; List Table 2 Accent 5; List Table 3 Accent 5;

List Table 4 Accent 5; List Table 5 Dark Accent 5; List Table 6 Colorful Accent 5; List Table 7 Colorful Accent 5;

List Table 1 Light Accent 6; List Table 2 Accent 6; List Table 3 Accent 6; List Table 4 Accent 6;

List Table 5 Dark Accent 6; List Table 6 Colorful Accent 6; List Table 7 Colorful Accent 6; Mention;

Smart Hyperlink; Hashtag; Unresolved Mention; Smart Link; *
