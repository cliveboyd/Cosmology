#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 08:45:30 2025

@author: boyde

TL;DR ordered list

If you want a one-liner sequence to follow in a fresh session (assuming config is correct and maps don’t exist yet):

quaia_build_zmean_maps.py

quaia_dipole_fit.py

quaia_build_counts_maps_randoms.py

quaia_dipole_fit_counts_randoms.py (+ the small comparison block you used)

quaia_dipole_residual_z_vs_counts.py

quaia_dipole_residual_z_vs_counts_slices.py

quaia_dipole_shuffle_significance_slices.py

quaia_build_zmean_maps_bcut.py

quaia_dipole_bcut_grid.py

quaia_dipole_jackknife_bcut20.py


SEQUENCE...
You said:
Explain jackknife tests
You said:
ok lets enhance the following Within Quaia (short-term, robustness)

Run more shuffles per slice, and store full amplitude distributions.

Do a small grid of |b| cuts (10°, 20°, 30°) and show that the low-z CMB alignment is stable.

Add jackknife tests (hemispheres / quadrants) to show no single region dominates.
You said:
ok these files have been created or updated  ... show me the order they should be run

RESULTS...
Python 3.11.14 | packaged by conda-forge | (main, Oct 22 2025, 22:56:31) [Clang 19.1.7 ]
Type "copyright", "credits" or "license" for more information.

IPython 9.6.0 -- An enhanced Interactive Python. Type '?' for help.

%run /Users/boyde/.spyder-py3/quaia_dipole_fit.py
[quaia] dipole fits from <z> maps

[slice] full sample
  N_good_pix = 40516
  dipole amplitude |b| = 1.0797e-02
  dipole direction (l,b) = (103.3°, 31.9°)
  a0 (offset, after mean-subtract) = 1.9194e-02

[slice] 0.10 ≤ z < 0.50
  N_good_pix = 6429
  dipole amplitude |b| = 1.9195e-03
  dipole direction (l,b) = (128.9°, 71.7°)
  a0 (offset, after mean-subtract) = -6.4787e-05

[slice] 0.50 ≤ z < 1.00
  N_good_pix = 27430
  dipole amplitude |b| = 5.6372e-04
  dipole direction (l,b) = (181.5°, 35.9°)
  a0 (offset, after mean-subtract) = 1.4902e-05

[slice] 1.00 ≤ z < 1.50
  N_good_pix = 32432
  dipole amplitude |b| = 1.9033e-03
  dipole direction (l,b) = (39.6°, 16.5°)
  a0 (offset, after mean-subtract) = 2.8965e-04

[slice] 1.50 ≤ z < 2.00
  N_good_pix = 29812
  dipole amplitude |b| = 1.4459e-03
  dipole direction (l,b) = (56.6°, 18.2°)
  a0 (offset, after mean-subtract) = -3.0052e-05

[slice] 2.00 ≤ z < 2.50
  N_good_pix = 19507
  dipole amplitude |b| = 1.5105e-03
  dipole direction (l,b) = (83.3°, -45.7°)
  a0 (offset, after mean-subtract) = 7.0554e-05

[save] wrote dipole summary -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_dipole_summary.txt
[quaia] done.

%run /Users/boyde/.spyder-py3/quaia_build_counts_maps_randoms.py
[quaia-rand] building COUNT maps from randoms
[load-rand] loading randoms from /Users/boyde/Downloads/quaia_G20p0_randoms_simple.fits
[build-rand] NSIDE=64  npix=49152
[build-rand] pixels with N>0: 18264 / 49152
[build-rand] pixels with N>=3: 16661 / 49152
[save-rand] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/randoms/quaia_random_counts_nside64.npz
[quaia-rand] done.

%run /Users/boyde/.spyder-py3/quaia_dipole_fit_counts_randoms.py
[quaia-rand] dipole fit from COUNT map (randoms)
  N_valid_pix = 16661
  mean N per valid pixel = 119.9
  dipole amplitude |b| = 4.5459e+00 (counts)
  dipole direction (l,b) = (156.1°,  64.8°)
  separation from CMB dipole ≈  54.2°
  a0 (offset term) = -2.7281e-01

[save-rand] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/randoms/quaia_randoms_counts_dipole.txt
[quaia-rand] done.

%run /Users/boyde/.spyder-py3/quaia_build_counts_maps_randoms.py
[quaia-rand] building COUNT maps from randoms
[load-rand] loading randoms from /Users/boyde/Downloads/quaia_G20p0_randoms_simple.fits
[build-rand] NSIDE=64  npix=49152
[build-rand] pixels with N>0: 18264 / 49152
[build-rand] pixels with N>=3: 16661 / 49152
[save-rand] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/randoms/quaia_random_counts_nside64.npz
[quaia-rand] done.

%run /Users/boyde/.spyder-py3/quaia_dipole_fit_counts_randoms.py
[quaia-rand] dipole fit from COUNT map (randoms)
  N_valid_pix = 16661
  mean N per valid pixel = 119.9
  dipole amplitude |b| = 4.5459e+00 (counts)
  dipole direction (l,b) = (156.1°,  64.8°)
  separation from CMB dipole ≈  54.2°
  a0 (offset term) = -2.7281e-01

[save-rand] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/randoms/quaia_randoms_counts_dipole.txt
[quaia-rand] done.

import numpy as np
from pathlib import Path
from quaia_config import OUT_DIR, NSIDE
from quaia_dipole_fit_counts_randoms import fit_dipole_from_counts, ang_sep, L_CMB, B_CMB

real_path = OUT_DIR / f"quaia_zmean_full_nside{NSIDE}.npz"
data_real = np.load(real_path, allow_pickle=False)

N_real = data_real["N"]

fit_real = fit_dipole_from_counts(N_real)
amp_r = fit_real["amp"]
l_r   = fit_real["l"]
b_r   = fit_real["b"]

theta_r = ang_sep(l_r, b_r, L_CMB, B_CMB)
sep_r_deg = np.rad2deg(theta_r)

print(f"[real counts] N_valid = {fit_real['N_valid']}")
print(f"[real counts] mean N = {fit_real['N_mean']:.1f}")
print(f"[real counts] amp = {amp_r:.4e}")
print(f"[real counts] (l,b) = ({np.rad2deg(l_r):.1f}°, {np.rad2deg(b_r):.1f}°)")
print(f"[real counts] separation from CMB dipole ≈ {sep_r_deg:.1f}°")

# Fractional amplitudes vs mean:
N_mean_rand = 141.3    # your printed value from console
amp_rand    = 4.2630   # from randoms console

frac_rand = amp_rand / N_mean_rand
print(f"[randoms] fractional amp ≈ {frac_rand:.4f}")

frac_real = amp_r / fit_real['N_mean']
print(f"[real counts] fractional amp ≈ {frac_real:.4f}")
[real counts] N_valid = 40516
[real counts] mean N = 18.6
[real counts] amp = 1.1683e+00
[real counts] (l,b) = (141.4°, 30.8°)
[real counts] separation from CMB dipole ≈ 86.0°
[randoms] fractional amp ≈ 0.0302
[real counts] fractional amp ≈ 0.0628

%run /Users/boyde/.spyder-py3/quaia_dipole_residual_z_vs_counts.py
[resid] good pixels: 40516 / 49152
[resid] z ≈ a + b*N with a=1.35306, b=7.37685e-03

[original <z> dipole]
  amp = 1.080e-02
  (l,b) = (103.31°, 31.89°)
  sep to CMB = 98.3°

[residual <z> dipole after removing N-trend]
  amp = 5.263e-03
  (l,b) = (64.25°, 22.14°)
  sep to CMB = 107.7°
  angle between original and residual = 35.9°

%run /Users/boyde/.spyder-py3/quaia_dipole_residual_z_vs_counts_slices.py
[resid-slices] Quaia <z> residual dipoles per redshift bin

[slice z0p10_0p50] 0.10 ≤ z < 0.50
  good pixels = 6429 / 49152
  z ≈ a + b*N with a=0.36328, b=-4.33896e-05
  [orig <z> dipole]
    amp = 1.919e-03
    (l,b) = (128.85°, 71.69°)
    sep to CMB = 56.2°
  [resid <z> dipole]
    amp = 1.921e-03
    (l,b) = (129.00°, 71.71°)
    sep to CMB = 56.1°
    angle(orig,resid) = 0.1°

[slice z0p50_1p00] 0.50 ≤ z < 1.00
  good pixels = 27430 / 49152
  z ≈ a + b*N with a=0.77362, b=4.74357e-05
  [orig <z> dipole]
    amp = 5.637e-04
    (l,b) = (181.49°, 35.94°)
    sep to CMB = 59.5°
  [resid <z> dipole]
    amp = 5.599e-04
    (l,b) = (181.83°, 35.93°)
    sep to CMB = 59.3°
    angle(orig,resid) = 0.3°

[slice z1p00_1p50] 1.00 ≤ z < 1.50
  good pixels = 32432 / 49152
  z ≈ a + b*N with a=1.25641, b=3.18491e-04
  [orig <z> dipole]
    amp = 1.903e-03
    (l,b) = (39.63°, 16.48°)
    sep to CMB = 104.3°
  [resid <z> dipole]
    amp = 1.896e-03
    (l,b) = (37.46°, 15.60°)
    sep to CMB = 104.1°
    angle(orig,resid) = 2.3°

[slice z1p50_2p00] 1.50 ≤ z < 2.00
  good pixels = 29812 / 49152
  z ≈ a + b*N with a=1.73460, b=-3.33252e-05
  [orig <z> dipole]
    amp = 1.446e-03
    (l,b) = (56.56°, 18.17°)
    sep to CMB = 109.4°
  [resid <z> dipole]
    amp = 1.447e-03
    (l,b) = (56.73°, 18.27°)
    sep to CMB = 109.4°
    angle(orig,resid) = 0.2°

[slice z2p00_2p50] 2.00 ≤ z < 2.50
  good pixels = 19507 / 49152
  z ≈ a + b*N with a=2.23596, b=1.90844e-04
  [orig <z> dipole]
    amp = 1.510e-03
    (l,b) = (83.33°, -45.66°)
    sep to CMB = 177.6°
  [resid <z> dipole]
    amp = 1.505e-03
    (l,b) = (83.33°, -46.16°)
    sep to CMB = 178.1°
    angle(orig,resid) = 0.5°

[resid-slices] wrote summary -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_dipole_residual_slices.txt

%run /Users/boyde/.spyder-py3/quaia_dipole_shuffle_significance_slices.py
[shuffle] Quaia <z> dipole significance per redshift slice

[shuffle] slice full  full sample
  good pixels = 40516 / 49152
  [real] amp = 1.080e-02, (l,b)=(103.31°, 31.89°)
         sep to CMB = 98.3°
    shuffle 5000/5000
  [shuffles] mean amp = 2.937e-03, std = 1.242e-03
             p ≈ 0.0000 (fraction with amp >= real)
  [save] amplitude distribution -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_dipole_shuffle_full_amps.npz

[shuffle] slice z0p10_0p50  0.10 ≤ z < 0.50
  good pixels = 6429 / 49152
  [real] amp = 1.919e-03, (l,b)=(128.85°, 71.69°)
         sep to CMB = 56.2°
    shuffle 5000/5000
  [shuffles] mean amp = 1.866e-03, std = 8.042e-04
             p ≈ 0.4396 (fraction with amp >= real)
  [save] amplitude distribution -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_dipole_shuffle_z0p10_0p50_amps.npz

[shuffle] slice z0p50_1p00  0.50 ≤ z < 1.00
  good pixels = 27430 / 49152
  [real] amp = 5.637e-04, (l,b)=(181.49°, 35.94°)
         sep to CMB = 59.5°
    shuffle 5000/5000
  [shuffles] mean amp = 1.274e-03, std = 5.465e-04
             p ≈ 0.9204 (fraction with amp >= real)
  [save] amplitude distribution -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_dipole_shuffle_z0p50_1p00_amps.npz

[shuffle] slice z1p00_1p50  1.00 ≤ z < 1.50
  good pixels = 32432 / 49152
  [real] amp = 1.903e-03, (l,b)=(39.63°, 16.48°)
         sep to CMB = 104.3°
    shuffle 5000/5000
  [shuffles] mean amp = 1.098e-03, std = 4.683e-04
             p ≈ 0.0592 (fraction with amp >= real)
  [save] amplitude distribution -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_dipole_shuffle_z1p00_1p50_amps.npz

[shuffle] slice z1p50_2p00  1.50 ≤ z < 2.00
  good pixels = 29812 / 49152
  [real] amp = 1.446e-03, (l,b)=(56.56°, 18.17°)
         sep to CMB = 109.4°
    shuffle 5000/5000
  [shuffles] mean amp = 1.200e-03, std = 5.096e-04
             p ≈ 0.2968 (fraction with amp >= real)
  [save] amplitude distribution -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_dipole_shuffle_z1p50_2p00_amps.npz

[shuffle] slice z2p00_2p50  2.00 ≤ z < 2.50
  good pixels = 19507 / 49152
  [real] amp = 1.510e-03, (l,b)=(83.33°, -45.66°)
         sep to CMB = 177.6°
    shuffle 5000/5000
  [shuffles] mean amp = 1.651e-03, std = 7.171e-04
             p ≈ 0.5420 (fraction with amp >= real)
  [save] amplitude distribution -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_dipole_shuffle_z2p00_2p50_amps.npz

[shuffle] wrote summary -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_dipole_shuffle_slices.txt

%run /Users/boyde/.spyder-py3/quaia_build_zmean_maps_bcut.py
[quaia-bcut] building <z> maps for BCUT grid: [10.0, 20.0, 30.0]
[bcut-build] loading Quaia from /Users/boyde/Downloads/quaia_G20.0.fits
[bcut-build] b_cut = 10.0 deg -> kept 744834 / 755850 objects
[bcut-build] NSIDE=64  npix=49152
[bcut-build] tag=full, |b| > 10.0 deg
               pixels with N>0      : 40149 / 49152
               pixels with N>=3: 39147 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_full_bcut10_nside64.npz
[bcut-build] tag=z0p10_0p50, |b| > 10.0 deg
               pixels with N>0      : 28464 / 49152
               pixels with N>=3: 6392 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_z0p10_0p50_bcut10_nside64.npz
[bcut-build] tag=z0p50_1p00, |b| > 10.0 deg
               pixels with N>0      : 38101 / 49152
               pixels with N>=3: 27019 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_z0p50_1p00_bcut10_nside64.npz
[bcut-build] tag=z1p00_1p50, |b| > 10.0 deg
               pixels with N>0      : 38598 / 49152
               pixels with N>=3: 32040 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_z1p00_1p50_bcut10_nside64.npz
[bcut-build] tag=z1p50_2p00, |b| > 10.0 deg
               pixels with N>0      : 37923 / 49152
               pixels with N>=3: 29526 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_z1p50_2p00_bcut10_nside64.npz
[bcut-build] tag=z2p00_2p50, |b| > 10.0 deg
               pixels with N>0      : 35289 / 49152
               pixels with N>=3: 19416 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_z2p00_2p50_bcut10_nside64.npz
[bcut-build] loading Quaia from /Users/boyde/Downloads/quaia_G20.0.fits
[bcut-build] b_cut = 20.0 deg -> kept 660631 / 755850 objects
[bcut-build] NSIDE=64  npix=49152
[bcut-build] tag=full, |b| > 20.0 deg
               pixels with N>0      : 32614 / 49152
               pixels with N>=3: 32388 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_full_bcut20_nside64.npz
[bcut-build] tag=z0p10_0p50, |b| > 20.0 deg
               pixels with N>0      : 24212 / 49152
               pixels with N>=3: 5775 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_z0p10_0p50_bcut20_nside64.npz
[bcut-build] tag=z0p50_1p00, |b| > 20.0 deg
               pixels with N>0      : 31555 / 49152
               pixels with N>=3: 23689 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_z0p50_1p00_bcut20_nside64.npz
[bcut-build] tag=z1p00_1p50, |b| > 20.0 deg
               pixels with N>0      : 32037 / 49152
               pixels with N>=3: 28058 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_z1p00_1p50_bcut20_nside64.npz
[bcut-build] tag=z1p50_2p00, |b| > 20.0 deg
               pixels with N>0      : 31784 / 49152
               pixels with N>=3: 26264 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_z1p50_2p00_bcut20_nside64.npz
[bcut-build] tag=z2p00_2p50, |b| > 20.0 deg
               pixels with N>0      : 30130 / 49152
               pixels with N>=3: 17849 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_z2p00_2p50_bcut20_nside64.npz
[bcut-build] loading Quaia from /Users/boyde/Downloads/quaia_G20.0.fits
[bcut-build] b_cut = 30.0 deg -> kept 530364 / 755850 objects
[bcut-build] NSIDE=64  npix=49152
[bcut-build] tag=full, |b| > 30.0 deg
               pixels with N>0      : 24889 / 49152
               pixels with N>=3: 24778 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_full_bcut30_nside64.npz
[bcut-build] tag=z0p10_0p50, |b| > 30.0 deg
               pixels with N>0      : 18775 / 49152
               pixels with N>=3: 4638 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_z0p10_0p50_bcut30_nside64.npz
[bcut-build] tag=z0p50_1p00, |b| > 30.0 deg
               pixels with N>0      : 24230 / 49152
               pixels with N>=3: 18834 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_z0p50_1p00_bcut30_nside64.npz
[bcut-build] tag=z1p00_1p50, |b| > 30.0 deg
               pixels with N>0      : 24592 / 49152
               pixels with N>=3: 22073 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_z1p00_1p50_bcut30_nside64.npz
[bcut-build] tag=z1p50_2p00, |b| > 30.0 deg
               pixels with N>0      : 24476 / 49152
               pixels with N>=3: 20860 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_z1p50_2p00_bcut30_nside64.npz
[bcut-build] tag=z2p00_2p50, |b| > 30.0 deg
               pixels with N>0      : 23435 / 49152
               pixels with N>=3: 14614 / 49152
[save-bcut] wrote /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_zmean_z2p00_2p50_bcut30_nside64.npz
[quaia-bcut] done.

%run /Users/boyde/.spyder-py3/quaia_dipole_bcut_grid.py
[bcut-grid] Quaia <z> dipoles for |b| cuts: [10.0, 20.0, 30.0]
[bcut-grid] tag=full, |b|>10.0: N_good=39147, amp=1.190e-02, (l,b)=(105.0°, 24.3°), sepCMB=105.3°
[bcut-grid] tag=z0p10_0p50, |b|>10.0: N_good=6392, amp=1.857e-03, (l,b)=(129.8°, 70.0°), sepCMB=57.4°
[bcut-grid] tag=z0p50_1p00, |b|>10.0: N_good=27019, amp=6.755e-04, (l,b)=(182.5°, 62.5°), sepCMB=45.2°
[bcut-grid] tag=z1p00_1p50, |b|>10.0: N_good=32040, amp=1.941e-03, (l,b)=(45.9°, 14.0°), sepCMB=109.3°
[bcut-grid] tag=z1p50_2p00, |b|>10.0: N_good=29526, amp=1.384e-03, (l,b)=(52.3°, 20.4°), sepCMB=105.9°
[bcut-grid] tag=z2p00_2p50, |b|>10.0: N_good=19416, amp=1.525e-03, (l,b)=(84.6°, -51.9°), sepCMB=176.1°
[bcut-grid] tag=full, |b|>20.0: N_good=32388, amp=8.962e-03, (l,b)=(110.9°, 18.5°), sepCMB=109.2°
[bcut-grid] tag=z0p10_0p50, |b|>20.0: N_good=5775, amp=2.103e-03, (l,b)=(273.4°, 54.7°), sepCMB=8.9°
[bcut-grid] tag=z0p50_1p00, |b|>20.0: N_good=23689, amp=9.104e-04, (l,b)=(246.4°, 41.1°), sepCMB=14.2°
[bcut-grid] tag=z1p00_1p50, |b|>20.0: N_good=28058, amp=1.662e-03, (l,b)=(34.9°, 12.3°), sepCMB=105.6°
[bcut-grid] tag=z1p50_2p00, |b|>20.0: N_good=26264, amp=1.175e-03, (l,b)=(59.3°, 10.4°), sepCMB=117.6°
[bcut-grid] tag=z2p00_2p50, |b|>20.0: N_good=17849, amp=1.592e-03, (l,b)=(81.5°, -71.6°), sepCMB=156.4°
[bcut-grid] tag=full, |b|>30.0: N_good=24778, amp=5.187e-03, (l,b)=(123.3°, 14.1°), sepCMB=108.7°
[bcut-grid] tag=z0p10_0p50, |b|>30.0: N_good=4638, amp=8.373e-04, (l,b)=(152.7°, 4.2°), sepCMB=100.8°
[bcut-grid] tag=z0p50_1p00, |b|>30.0: N_good=18834, amp=1.275e-03, (l,b)=(269.5°, 36.8°), sepCMB=11.9°
[bcut-grid] tag=z1p00_1p50, |b|>30.0: N_good=22073, amp=2.321e-03, (l,b)=(41.8°, 21.9°), sepCMB=100.6°
[bcut-grid] tag=z1p50_2p00, |b|>30.0: N_good=20860, amp=9.684e-04, (l,b)=(80.1°, -19.9°), sepCMB=151.7°
[bcut-grid] tag=z2p00_2p50, |b|>30.0: N_good=14614, amp=1.126e-03, (l,b)=(329.1°, -20.0°), sepCMB=89.4°
[bcut-grid] wrote summary -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_dipole_bcut_grid.txt

%run /Users/boyde/.spyder-py3/quaia_dipole_jackknife_bcut20.py
[JK-bcut20] Quaia <z> jackknife tests with |b| > 20 deg

[JK-bcut20] slice full  full sample
  [full] N_good=32388, amp=8.962e-03, (l,b)=(110.9°, 18.5°), sep(CMB)=109.2°
  [JK leave_N_hemisphere] N_use=16120, amp=2.554e-02 (Δamp/amp_full=185.03%), (l,b)=(60.8°, 12.1°), sep(full)=48.7°, sep(CMB)=116.4°
  [JK leave_S_hemisphere] N_use=16268, amp=2.061e-02 (Δamp/amp_full=129.93%), (l,b)=(197.5°, 13.6°), sep(full)=82.5°, sep(CMB)=64.3°
  [JK leave_Q1_0-90] N_use=23388, amp=1.493e-02 (Δamp/amp_full=66.60%), (l,b)=(116.2°, 54.7°), sep(full)=36.4°, sep(CMB)=73.8°
  [JK leave_Q2_90-180] N_use=25221, amp=1.426e-02 (Δamp/amp_full=59.06%), (l,b)=(113.0°, 5.1°), sep(full)=13.5°, sep(CMB)=121.1°
  [JK leave_Q3_180-270] N_use=23328, amp=1.081e-02 (Δamp/amp_full=20.63%), (l,b)=(94.4°, -31.2°), sep(full)=52.2°, sep(CMB)=161.5°
  [JK leave_Q4_270-360] N_use=25227, amp=5.837e-03 (Δamp/amp_full=-34.87%), (l,b)=(122.6°, 38.9°), sep(full)=22.8°, sep(CMB)=86.6°

[JK-bcut20] slice z0p10_0p50  0.10 ≤ z < 0.50
  [full] N_good=5775, amp=2.103e-03, (l,b)=(273.4°, 54.7°), sep(CMB)=8.9°
  [JK leave_N_hemisphere] N_use=2782, amp=7.060e-03 (Δamp/amp_full=235.77%), (l,b)=(322.8°, 72.5°), sep(full)=27.0°, sep(CMB)=35.6°
  [JK leave_S_hemisphere] N_use=2993, amp=4.056e-03 (Δamp/amp_full=92.89%), (l,b)=(223.3°, 57.2°), sep(full)=27.5°, sep(CMB)=25.9°
  [JK leave_Q1_0-90] N_use=4236, amp=3.235e-03 (Δamp/amp_full=53.87%), (l,b)=(300.8°, 72.1°), sep(full)=20.9°, sep(CMB)=29.3°
  [JK leave_Q2_90-180] N_use=4424, amp=2.501e-03 (Δamp/amp_full=18.92%), (l,b)=(288.6°, 38.1°), sep(full)=19.5°, sep(CMB)=20.4°
  [JK leave_Q3_180-270] N_use=4068, amp=2.834e-03 (Δamp/amp_full=34.77%), (l,b)=(259.5°, 20.3°), sep(full)=36.0°, sep(CMB)=28.0°
  [JK leave_Q4_270-360] N_use=4597, amp=1.635e-03 (Δamp/amp_full=-22.24%), (l,b)=(185.5°, 60.9°), sep(full)=43.7°, sep(CMB)=44.4°

[JK-bcut20] slice z0p50_1p00  0.50 ≤ z < 1.00
  [full] N_good=23689, amp=9.104e-04, (l,b)=(246.4°, 41.1°), sep(CMB)=14.2°
  [JK leave_N_hemisphere] N_use=11762, amp=3.996e-03 (Δamp/amp_full=338.91%), (l,b)=(293.1°, 62.9°), sep(full)=34.8°, sep(CMB)=21.9°
  [JK leave_S_hemisphere] N_use=11927, amp=1.375e-03 (Δamp/amp_full=50.99%), (l,b)=(170.0°, 12.6°), sep(full)=71.6°, sep(CMB)=83.3°
  [JK leave_Q1_0-90] N_use=17206, amp=2.392e-03 (Δamp/amp_full=162.79%), (l,b)=(53.7°, 13.7°), sep(full)=123.9°, sep(CMB)=112.7°
  [JK leave_Q2_90-180] N_use=18223, amp=2.971e-03 (Δamp/amp_full=226.31%), (l,b)=(295.7°, 7.7°), sep(full)=54.9°, sep(CMB)=48.4°
  [JK leave_Q3_180-270] N_use=16996, amp=3.336e-03 (Δamp/amp_full=266.46%), (l,b)=(243.0°, 11.7°), sep(full)=29.6°, sep(CMB)=40.4°
  [JK leave_Q4_270-360] N_use=18642, amp=2.486e-03 (Δamp/amp_full=173.07%), (l,b)=(147.9°, 15.8°), sep(full)=85.9°, sep(CMB)=94.6°

[JK-bcut20] slice z1p00_1p50  1.00 ≤ z < 1.50
  [full] N_good=28058, amp=1.662e-03, (l,b)=(34.9°, 12.3°), sep(CMB)=105.6°
  [JK leave_N_hemisphere] N_use=13918, amp=1.883e-03 (Δamp/amp_full=13.29%), (l,b)=(85.5°, 10.3°), sep(full)=49.6°, sep(CMB)=121.6°
  [JK leave_S_hemisphere] N_use=14140, amp=2.566e-03 (Δamp/amp_full=54.35%), (l,b)=(1.8°, 0.8°), sep(full)=34.8°, sep(CMB)=94.6°
  [JK leave_Q1_0-90] N_use=20452, amp=9.541e-04 (Δamp/amp_full=-42.60%), (l,b)=(289.8°, 29.7°), sep(full)=96.6°, sep(CMB)=26.9°
  [JK leave_Q2_90-180] N_use=21655, amp=2.322e-03 (Δamp/amp_full=39.70%), (l,b)=(72.5°, 21.3°), sep(full)=37.0°, sep(CMB)=109.9°
  [JK leave_Q3_180-270] N_use=20105, amp=2.591e-03 (Δamp/amp_full=55.88%), (l,b)=(42.3°, 8.5°), sep(full)=8.3°, sep(CMB)=112.6°
  [JK leave_Q4_270-360] N_use=21962, amp=2.029e-03 (Δamp/amp_full=22.06%), (l,b)=(13.3°, -2.4°), sep(full)=26.0°, sep(CMB)=104.6°

[JK-bcut20] slice z1p50_2p00  1.50 ≤ z < 2.00
  [full] N_good=26264, amp=1.175e-03, (l,b)=(59.3°, 10.4°), sep(CMB)=117.6°
  [JK leave_N_hemisphere] N_use=13049, amp=1.637e-03 (Δamp/amp_full=39.30%), (l,b)=(97.3°, 72.7°), sep(full)=66.3°, sep(CMB)=58.9°
  [JK leave_S_hemisphere] N_use=13215, amp=3.432e-03 (Δamp/amp_full=191.96%), (l,b)=(47.7°, -52.7°), sep(full)=63.9°, sep(CMB)=156.6°
  [JK leave_Q1_0-90] N_use=19116, amp=1.884e-03 (Δamp/amp_full=60.28%), (l,b)=(73.7°, -3.5°), sep(full)=20.0°, sep(CMB)=134.7°
  [JK leave_Q2_90-180] N_use=20231, amp=1.624e-03 (Δamp/amp_full=38.20%), (l,b)=(37.1°, 22.8°), sep(full)=24.5°, sep(CMB)=97.7°
  [JK leave_Q3_180-270] N_use=18912, amp=1.457e-03 (Δamp/amp_full=23.97%), (l,b)=(39.6°, 45.0°), sep(full)=38.5°, sep(CMB)=79.2°
  [JK leave_Q4_270-360] N_use=20533, amp=8.790e-04 (Δamp/amp_full=-25.22%), (l,b)=(58.2°, -63.6°), sep(full)=74.0°, sep(CMB)=159.0°

[JK-bcut20] slice z2p00_2p50  2.00 ≤ z < 2.50
  [full] N_good=17849, amp=1.592e-03, (l,b)=(81.5°, -71.6°), sep(CMB)=156.4°
  [JK leave_N_hemisphere] N_use=8834, amp=1.991e-03 (Δamp/amp_full=25.04%), (l,b)=(238.1°, 33.7°), sep(full)=140.1°, sep(CMB)=24.1°
  [JK leave_S_hemisphere] N_use=9015, amp=3.978e-03 (Δamp/amp_full=149.85%), (l,b)=(61.9°, -50.1°), sep(full)=23.3°, sep(CMB)=165.4°
  [JK leave_Q1_0-90] N_use=12897, amp=3.601e-03 (Δamp/amp_full=126.18%), (l,b)=(109.2°, -70.5°), sep(full)=9.0°, sep(CMB)=154.5°
  [JK leave_Q2_90-180] N_use=13669, amp=1.085e-03 (Δamp/amp_full=-31.88%), (l,b)=(88.5°, -49.2°), sep(full)=22.6°, sep(CMB)=176.8°
  [JK leave_Q3_180-270] N_use=12804, amp=1.481e-03 (Δamp/amp_full=-7.00%), (l,b)=(40.1°, -17.4°), sep(full)=59.4°, sep(CMB)=133.0°
  [JK leave_Q4_270-360] N_use=14177, amp=1.507e-03 (Δamp/amp_full=-5.35%), (l,b)=(92.2°, -76.4°), sep(full)=5.6°, sep(CMB)=151.4°

[JK-bcut20] wrote summary -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_dipole_jackknife_bcut20.txt

%run /Users/boyde/.spyder-py3/quaia_dipole_fit_bcut.py
[bcut] Quaia <z> dipole fits with |b| > 20.0 deg

[bcut slice] full sample
  N_good_pix = 32388
  dipole amplitude |b| = 8.9620e-03
  dipole direction (l,b) = (110.9°, 18.5°)
  a0 (offset, after mean-subtract) = 5.9903e-03
  separation from CMB dipole ≈ 109.2°

[bcut slice] 0.10 ≤ z < 0.50
  N_good_pix = 5775
  dipole amplitude |b| = 2.1027e-03
  dipole direction (l,b) = (273.4°, 54.7°)
  a0 (offset, after mean-subtract) = -3.5965e-05
  separation from CMB dipole ≈ 8.9°

[bcut slice] 0.50 ≤ z < 1.00
  N_good_pix = 23689
  dipole amplitude |b| = 9.1037e-04
  dipole direction (l,b) = (246.4°, 41.1°)
  a0 (offset, after mean-subtract) = -9.7878e-05
  separation from CMB dipole ≈ 14.2°

[bcut slice] 1.00 ≤ z < 1.50
  N_good_pix = 28058
  dipole amplitude |b| = 1.6622e-03
  dipole direction (l,b) = (34.9°, 12.3°)
  a0 (offset, after mean-subtract) = 9.7139e-05
  separation from CMB dipole ≈ 105.6°

[bcut slice] 1.50 ≤ z < 2.00
  N_good_pix = 26264
  dipole amplitude |b| = 1.1754e-03
  dipole direction (l,b) = (59.3°, 10.4°)
  a0 (offset, after mean-subtract) = -1.1248e-04
  separation from CMB dipole ≈ 117.6°

[bcut slice] 2.00 ≤ z < 2.50
  N_good_pix = 17849
  dipole amplitude |b| = 1.5923e-03
  dipole direction (l,b) = (81.5°, -71.6°)
  a0 (offset, after mean-subtract) = 1.6109e-04
  separation from CMB dipole ≈ 156.4°

[bcut] wrote summary -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_dipole_summary_bcut20.txt

"""

