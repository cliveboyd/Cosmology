#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_build_zmean_maps.py  (0V1)

Build N(z) and <z> HEALPix maps from the Quaia catalogue,
using config in quaia_config.py.

- For each HEALPix pixel:
    N_pix   = number of objects
    z_mean  = mean redshift_quaia in that pixel

- Also builds the same maps in redshift bins defined in Z_BINS.

Outputs go into OUT_DIR from quaia_config.
"""

from pathlib         import Path

import numpy         as      np
from   astropy.table import Table
import healpy        as     hp

from quaia_config    import (
                             QUAIA_FITS,
                             OUT_DIR,
                             NSIDE,
                             N_MIN_PER_PIX,
                             COL_RA,
                             COL_DEC,
                             COL_Z,
                             Z_BINS,
                             Z_BINS_MASK
                            )


def build_maps_for_sample(mask_name, z_min=None, z_max=None):
    """
    Build N and <z> maps, optionally restricted to a redshift slice.

    Parameters
    ----------
    mask_name : str
        Label appended to output filenames (e.g. "full", "z0p5_1p0").
    z_min, z_max : float or None
        Redshift bounds. If both None, use full sample.
    """
    print(f"[build] loading Quaia from {QUAIA_FITS}")
    t   = Table.read(QUAIA_FITS)

    ra  = np.array(t[COL_RA], dtype=float)
    dec = np.array(t[COL_DEC], dtype=float)
    z   = np.array(t[COL_Z], dtype=float)

    # Apply z-cut if requested
    if (z_min is not None) or (z_max is not None):
        m = np.ones_like(z, dtype=bool)
        
        if z_min is not None:
            m &= (z >= z_min)
        
        if z_max is not None:
            m &= (z < z_max)
        ra, dec, z = ra[m], dec[m], z[m]
        print(f"[build] redshift slice {z_min}–{z_max}: kept {m.sum()} / {len(m)} objects")

    # (ra, dec) -> HEALPix pixel index
    theta    = np.deg2rad(90.0 - dec)  # colatitude
    phi      = np.deg2rad(ra)          # longitude
    pix      = hp.ang2pix(NSIDE, theta, phi, nest=False)

    npix     = hp.nside2npix(NSIDE)
    print(f"[build] NSIDE={NSIDE}  npix={npix}")

    # Accumulate counts and z-sums per pixel
    N_map    = np.zeros(npix, dtype=np.int32)
    sumz_map = np.zeros(npix, dtype=float)

    for p, zz in zip(pix, z):
        N_map[p]    += 1
        sumz_map[p] += zz

    # Mean z per pixel where we have at least one object
    zmean_map       = np.full(npix, np.nan, dtype=float)
    good            = N_map > 0
    zmean_map[good] = sumz_map[good] / N_map[good]

    # Diagnostics on occupancy
    n_nonzero = np.count_nonzero(N_map)
    print(f"[build] pixels with N>0: {n_nonzero} / {npix}")
    print(f"[build] pixels with N>={N_MIN_PER_PIX}: {np.count_nonzero(N_map >= N_MIN_PER_PIX)} / {npix}")

    # Save maps as simple .npz bundle
    out = OUT_DIR / f"quaia_zmean_{mask_name}_nside{NSIDE}.npz"
    np.savez(
             out,
             N     = N_map,
             zmean = zmean_map,
             NSIDE = NSIDE,
             z_min = z_min,
             z_max = z_max,
            )
    print(f"[save] wrote {out}")


def main():
    print("[quaia] building <z> maps")

    # 1) Full sample
    build_maps_for_sample(mask_name="full")

    # 2) Redshift slices
    for (z_lo, z_hi) in Z_BINS:
        tag = f"z{z_lo:.2f}_{z_hi:.2f}".replace(".", "p")
        build_maps_for_sample(mask_name=tag, z_min=float(z_lo), z_max=float(z_hi))

    print("[quaia] done.")


if __name__ == "__main__":
    main()
