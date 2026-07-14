#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  2 00:56:14 2025

@author: boyde

TODO...

# 1. Build main <z> maps & slices
%run /Users/boyde/.spyder-py3/quaia_build_zmean_maps.py

# 2. FR-style dipole fits (now 0V2, with CMB projection)
%run /Users/boyde/.spyder-py3/quaia_dipole_fit.py

# 3. Remove N-trend and re-fit dipoles per slice
%run /Users/boyde/.spyder-py3/quaia_dipole_residual_z_vs_counts_slices.py

# 4. Shuffle significance per slice (5000 shuffles, amp distributions)
%run /Users/boyde/.spyder-py3/quaia_dipole_shuffle_significance_slices.py

# 5. Build |b| > 10/20/30 maps
%run /Users/boyde/.spyder-py3/quaia_build_zmean_maps_bcut.py

# 6. FR-style bcut dipole grid (now 0V2, with CMB projection if you wired it in)
%run /Users/boyde/.spyder-py3/quaia_dipole_bcut_grid.py

# 7. Detailed |b|>20 dipole summary per slice
%run /Users/boyde/.spyder-py3/quaia_dipole_fit_bcut.py

# 8. Jackknife tests (hemispheres & quadrants, |b|>20)
%run /Users/boyde/.spyder-py3/quaia_dipole_jackknife_bcut20.py

# 9. Skymaps (mean & residual) for full + each z-slice, |b|>20
%run /Users/boyde/.spyder-py3/quaia_plot_maps.py


Your “next gen” FR/Mittal-style run is:
%run /Users/boyde/.spyder-py3/quaia_build_zmean_maps_bcut.py
%run /Users/boyde/.spyder-py3/quaia_dipole_bcut_grid.py
%run /Users/boyde/.spyder-py3/quaia_dipole_fit_bcut.py
%run /Users/boyde/.spyder-py3/quaia_dipole_jackknife_bcut20.py
%run /Users/boyde/.spyder-py3/quaia_plot_maps.py


DESCRIPTION...
Summary of Recent Quaia Pipeline Changes and Relation to the FR Target Structure

Over the last iteration we have rebuilt the Quaia analysis pipeline so that it 
aligns much more cleanly with the Full Relativity (FR) framework, 
rather than “ΛCDM + variable-c”. 

The key changes are:

1. Finer redshift binning consistent with FR goals

What changed

The redshift binning has been refined to:

𝑍_𝐵𝐼𝑁𝑆 =[ 0.10, 0.30
         0.30, 0.50
         0.50, 0.75
         0.75, 1.00
         1.00, 1.50
         1.50, 2.50]  
	
Z_BINS=​  0.10, 0.30, 0.50, 0.75, 1.00, 1.50
	​
Corresponding maps were rebuilt with:

quaia_build_zmean_maps.py (full + all slices),

quaia_build_zmean_maps_bcut.py (for |b| > 10°, 20°, 30°).


Why this matters for FR

FR predicts that the redshift dipole is a kinematic Doppler effect from our 
motion through a stationary galaxy background, with a cos θ pattern of fixed 
amplitude (∝ v/c) that does not depend on z once systematics are removed.

The finer binning allows us to:

Separate the very low-z regime (0.10–0.30) where peculiar motions and 
spectroscopic redshifts dominate.

Track how the measured dipole direction and amplitude evolve from 0.10 to 2.50, 
and test whether the component aligned with the CMB dipole tends toward a 
constant value at low z, as FR requires.

2. Core dipole fitter now FR-aware (projection onto CMB direction)

What changed

quaia_dipole_fit.py has been upgraded (0V2) to:

Fit a linear dipole in each map:


4. Galactic latitude cuts and jackknife tests

What changed

|b| cuts grid:

quaia_build_zmean_maps_bcut.py and quaia_dipole_bcut_grid.py now build and 
analyse:

|b| > 10°, 20°, 30°,

For full sample and all new z-bins.

For each (bcut, z-bin), we record:

N_good pixels,

Dipole amplitude and (l, b),

Separation from the CMB dipole.

Jackknife robustness:

quaia_dipole_jackknife_bcut20.py runs jackknife tests for |b| > 20°:

Remove the Northern or Southern Galactic hemisphere.

Remove each of four longitude quadrants 0–90°, 90–180°, 180–270°, 270–360° in turn.

Refit the dipole each time and compare amplitude and direction to the full-sample 
result.

Why this matters for FR

FR is claiming a global, kinematic effect, not a local patch driven by one 
weird region (e.g. near the Galactic centre, or a single survey stripe).

The |b| cuts and jackknifes show:

How stable the CMB-aligned, low-z dipole is when:

We tighten Galactic cuts to reduce dust and stellar contamination.

We remove large chunks of the sky.

Whether any one hemisphere or quadrant is unduly dominating the signal.

If the FR-consistent component (amp_par along CMB) persists across:

Multiple |b| cuts,

Jackknife subsamples,

then it is much harder to explain the result away as a survey artefact.


5. Skymap visualisation for FR interpretation

What changed

Added Option USE_GC_HOLE=True into quaia_dipole_bcut_grid.py

quaia_plot_maps.py now generates for |b| > 20°:

Mean <z> maps,

Residual maps (after subtraction of the fitted dipole),

For:

Full sample,

Each of the new redshift bins.

These plots are written as PNGs in skymap2/quaia_outputs/ and can be inserted 
directly into the report.

Why this matters for FR

The maps let us visually check:

That the large-scale dipole structure resembles a smooth cos θ pattern,

That after subtracting the best-fit dipole the residual maps look statistically 
consistent with noise (no obvious residual large-scale gradients),

Whether the low-z maps qualitatively resemble the FR expectation of a kinematic 
flow, and how that appearance changes with redshift.

How this all fits into the FR target structure

Putting it together, the recent changes move the Quaia analysis from a 
generic dipole-hunting exercise towards a purpose-built FR test:


No ΛCDM structure is imposed:

⟨z⟩ as an empirical field on the sky.

No FR-incompatible z-weighting, distances, or expansion-based modelling appear 
in the pipeline.

We now explicitly measure the component of the Quaia dipole along the C
MB dipole (amp_par) in each redshift bin and for different |b| cuts.

The refined binning and |b| grid allow us to search for a regime 
(especially low z and |b| > 20–30°) where:

The dipole direction aligns with the CMB,

The CMB-aligned amplitude approaches the FR-expected value 
(≈ 0.00123 in the cos θ = 1 direction),

This component is stable against Galactic cuts and jackknife resampling.

The shuffle and residual tests guarantee that any FR-like signal we see 
is not trivially explained by:

Uneven sky coverage,

Varying number counts,

Or a peculiar choice of mask.

In short, the pipeline is now structured to answer a sharply FR-framed question:

Is there a statistically robust component of the Quaia redshift dipole that is 
aligned with the CMB dipole, of approximately the expected kinematic amplitude, 
and independent of redshift once observational systematics are controlled?


"""

