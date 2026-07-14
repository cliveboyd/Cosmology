# quaia outputs 23N0V25

- Source: `quaia outputs 23N0V25.docx`
- Source size: 7235645 bytes
- Source modified: 2025-11-23T09:01:10
- Extracted: 2026-07-14
- Word count estimate: 670

## Extracted Text
Summary of output plots generated via:

%run /Users/boyde/.spyder-py3/quaia_plot_dipole_quadrupole_maps.py

quaia_plot_dipole_quadrupole_maps.py

#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""

Created on Sun Nov 23 08:36:40 2025

@author:  clive.boyd@hotmail.com

"""

#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""

quaia_plot_dipole_quadrupole_maps.py  (0V1)

Make HEALPix maps for:

- full ⟨z⟩ field

- dipole-only component

- quadrupole-only component

for:

- full sample

- the standard redshift slices

Works with:

quaia_zmean_full_nside{NSIDE}.npz

quaia_zmean_zXXXX_YYYY_nside{NSIDE}.npz

Note: Optional toggle a b-cut flavour… edit TAGS and FILE_FMT.

"""

import numpy as np

import healpy as hp

import matplotlib.pyplot as plt

from pathlib import Path

from quaia_config import OUT_DIR, NSIDE, Z_BINS

# ---------------------------------------------------------------------

# Config: which tags to plot

# ---------------------------------------------------------------------

# This matches your existing naming convention at nside=64:

#   full

#   z0p10_0p50, z0p50_1p00, z1p00_1p50, z1p50_2p00, z2p00_2p50

#

# For the |b|>20° set, just change FILE_FMT below to "..._bcut20_nside{NSIDE}.npz"

# ---------------------------------------------------------------------

def slice_tag(z_lo, z_hi):

return f"z{z_lo:.2f}_{z_hi:.2f}".replace(".", "p")

TAGS = ["full"] + [slice_tag(z_lo, z_hi) for (z_lo, z_hi) in Z_BINS]

# Base filename pattern (no b-cut):

FILE_FMT = "quaia_zmean_{tag}_nside{nside}.npz"

# If you want to use the b-cut maps instead, change to:

# FILE_FMT = "quaia_zmean_{tag}_bcut20_nside{nside}.npz"

# Output PNG prefix

PNG_PREFIX = "quaia_map"

# ---------------------------------------------------------------------

# Helpers

# ---------------------------------------------------------------------

def load_zmean_map(tag):

"""

Load N and zmean HEALPix maps for a given tag.

Returns (N_map, zmean_map).

"""

path = OUT_DIR / FILE_FMT.format(tag=tag, nside=NSIDE)

data = np.load(path)

N_map   = data["N"]

zmean   = data["zmean"]

return N_map, zmean

def build_good_map(zmean):

"""

Build a full-sky HEALPix map, with bad pixels set to UNSEEN.

"""

npix = hp.nside2npix(NSIDE)

m = np.full(npix, hp.UNSEEN, dtype=float)

good = np.isfinite(zmean)

m[good] = zmean[good]

return m

def decompose_dipole_quadrupole(m):

"""

Given a HEALPix map m (monopole+dipole+quadrupole+...), compute:

- dipole-only map

- quadrupole-only map

using spherical harmonic decomposition up to l=2.

"""

# Replace UNSEEN with 0 for the harmonic transform

m_work = m.copy()

unseen_mask = (m_work == hp.UNSEEN)

m_work[unseen_mask] = 0.0

# Get alm up to l=2

alm = hp.map2alm(m_work, lmax=2)

# Index mapping for lmax=2:

# index(l,m) = l*(l+1)//2 + m

# l=0: m=0 -> index 0

# l=1: m=0,1 -> indices 1,2

# l=2: m=0,1,2 -> indices 3,4,5

alm_dip = np.zeros_like(alm)

alm_quad = np.zeros_like(alm)

# dipole terms l=1 -> indices 1 and 2

alm_dip[1:3] = alm[1:3]

# quadrupole terms l=2 -> indices 3,4,5

alm_quad[3:6] = alm[3:6]

map_dip  = hp.alm2map(alm_dip, NSIDE, verbose=False)

map_quad = hp.alm2map(alm_quad, NSIDE, verbose=False)

# Restore UNSEEN where original map was unseen

map_dip[unseen_mask] = hp.UNSEEN

map_quad[unseen_mask] = hp.UNSEEN

return map_dip, map_quad

def plot_moll(m, title, outpath, min_val=None, max_val=None, unit=None):

"""

Standard Mollweide projection of a HEALPix map, saved to PNG.

"""

plt.figure(figsize=(8, 5))

hp.mollview(

m,

title=title,

unit=unit,

min=min_val,

max=max_val,

hold=True,

)

hp.graticule()

plt.savefig(outpath, dpi=200, bbox_inches="tight")

plt.close()

def main():

print("[plot] Quaia ⟨z⟩ / dipole / quadrupole maps")

OUT_DIR.mkdir(parents=True, exist_ok=True)

for tag in TAGS:

print(f"[plot] tag={tag}")

N_map, zmean = load_zmean_map(tag)

# Full ⟨z⟩ map (mean-subtracted for visual contrast)

m_z = build_good_map(zmean)

good = (m_z != hp.UNSEEN)

z_mean_val = np.mean(m_z[good])

m_z_zero = m_z.copy()

m_z_zero[good] = m_z[good] - z_mean_val

# Dipole / quadrupole decomposition

map_dip, map_quad = decompose_dipole_quadrupole(m_z_zero)

# Auto scale dipole/quadrupole around 0 with symmetric limits

for name, m_comp in [("zmean", m_z_zero),

("dipole", map_dip),

("quadrupole", map_quad)]:

good_comp = (m_comp != hp.UNSEEN)

if not np.any(good_comp):

print(f"  [warn] no good pixels for {name} in tag={tag}")

continue

v = np.abs(m_comp[good_comp]).max()

vmin, vmax = -v, v

fname = f"{PNG_PREFIX}_{name}_{tag}_nside{NSIDE}.png"

outpath = OUT_DIR / fname

if name == "zmean":

title = f"Quaia ⟨z⟩ (mean-subtracted), {tag}"

unit = "Δz"

elif name == "dipole":

title = f"Quaia ⟨z⟩ dipole component, {tag}"

unit = "Δz (dipole)"

else:

title = f"Quaia ⟨z⟩ quadrupole component, {tag}"

unit = "Δz (quadrupole)"

print(f"  [save] {outpath}")

plot_moll(m_comp, title, outpath, min_val=vmin, max_val=vmax, unit=unit)

print("[plot] done.")

if __name__ == "__main__":

main()
