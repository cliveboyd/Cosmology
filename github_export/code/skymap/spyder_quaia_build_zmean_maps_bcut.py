#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_build_zmean_maps_bcut.py  (0V2)

Build <z> HEALPix maps from Quaia for a grid of Galactic latitude cuts:

   |b| > b_cut   with b_cut in BCUT_LIST

For each b_cut and each z-slice in Z_BINS (plus full sample), write:
    quaia_zmean_{tag}_bcut{int(b_cut)}_nside{NSIDE}.npz

Fields:
    N      : counts per pixel
    zmean  : mean redshift per pixel (NaN if N < N_MIN_PER_PIX)
    NSIDE  : NSIDE value

Typical usage:
    %run /Users/boyde/.spyder-py3/quaia_build_zmean_maps_bcut.py
    
Optional ---> USE_GC_HOLE=True

Outputs:
    quaia_zmean_full_bcut10_nside64.npz, bcut20, bcut30
    plus all the z-slice versions.
"""

import numpy         as     np
from   pathlib       import Path
import healpy        as     hp
from   astropy.table import Table

from quaia_config import (
                          QUAIA_FITS,
                          OUT_DIR,
                          NSIDE,
                          N_MIN_PER_PIX,
                          COL_RA,
                          COL_DEC,
                          COL_Z,
                          Z_SLICES,
                          USE_GC_HOLE,
                          GC_HOLE_L0_DEG,
                          GC_HOLE_B0_DEG,
                          GC_HOLE_RADIUS_DEG,
                         )

# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------
BCUT_LIST          = [10.0, 15.0, 20.0, 25.0, 30.0, 35.0]
N_MIN_PER_PIX      = 5      # Was = 3

# Mittal-style extra mask around Galactic Centre
USE_GC_HOLE        = True      # set False to revert to plain |b|>b_cut

"""
GC_HOLE_L0_DEG     = 0.0       # centre longitude (Galactic)
GC_HOLE_B0_DEG     = 20.0      # centre latitude (Galactic)
GC_HOLE_RADIUS_DEG = 40.0      # radius of the hole in degrees
"""
GC_HOLE_L0_DEG     =  0.0      # centre longitude (Galactic)
GC_HOLE_B0_DEG     = 20.0      # centre latitude  (Galactic)
GC_HOLE_RADIUS_DEG = 40.0      # radius of the hole in degrees

def ang_sep_deg(l1_deg, b1_deg, l2_deg, b2_deg):
    """
    Great-circle separation (deg) between two Galactic directions.

    l1_deg, b1_deg can be arrays; l2_deg, b2_deg are usually scalars.
    """
    
    l1      = np.deg2rad(l1_deg)
    b1      = np.deg2rad(b1_deg)
    
    l2      = np.deg2rad(l2_deg)
    b2      = np.deg2rad(b2_deg)

    x1      = np.cos(b1)*np.cos(l1)
    y1      = np.cos(b1)*np.sin(l1); z1 = np.sin(b1)
    
    x2      = np.cos(b2)*np.cos(l2)
    y2      = np.cos(b2)*np.sin(l2); z2 = np.sin(b2)

    dot     = x1*x2 + y1*y2 + z1*z2
    dot     = np.clip(dot, -1.0, 1.0)
    
    return np.rad2deg(np.arccos(dot))


def load_quaia():
    print(f"[bcut-build] loading Quaia from {QUAIA_FITS}")
    t       = Table.read(QUAIA_FITS)
    
    # Adjust column names if needed.
    ra      = np.array(t["ra"],             dtype=float)
    dec     = np.array(t["dec"],            dtype=float)
    z       = np.array(t["redshift_quaia"], dtype=float)
    b       = np.array(t["b"],              dtype=float)    # Galactic latitude in deg
    l       = np.array(t["l"],              dtype=float)    # Galactic longitude in deg
    return ra, dec, z, b, l


def build_maps_for_bcut(b_cut_deg):
    """
    Build maps for full and all z-slices for a given |b| cut,
    optionally with a Mittal-style GC hole.
    """
    ra, dec, z, b, l = load_quaia()

    # 1) Standard |b| > b_cut
    mask_b           = np.abs(b) > b_cut_deg

    # 2) Optional GC hole: remove sources within GC_HOLE_RADIUS_DEG of (l0,b0)
    if USE_GC_HOLE:
        sep_gc_deg = ang_sep_deg(l, b, GC_HOLE_L0_DEG, GC_HOLE_B0_DEG)
        mask_gc    = sep_gc_deg > GC_HOLE_RADIUS_DEG
    else:
        mask_gc    = np.ones_like(mask_b, dtype=bool)

    # 3) Combined mask
    mask_all = mask_b & mask_gc

    ra_cut   = ra[mask_all]
    dec_cut  = dec[mask_all]
    z_cut    = z[mask_all]

    print(f"[bcut-build] b_cut = {b_cut_deg:.1f} deg -> kept {mask_all.sum()} / {len(z)} objects")
    if USE_GC_HOLE:
        print(f"             GC hole: centre=({GC_HOLE_L0_DEG:.1f}°, {GC_HOLE_B0_DEG:.1f}°), "
              f"radius={GC_HOLE_RADIUS_DEG:.1f}°")

    # --- HEALPix binning ---
    theta    = np.deg2rad(90.0 - dec_cut)  # colatitude
    phi      = np.deg2rad(ra_cut)          # longitude
    ipix     = hp.ang2pix(NSIDE, theta, phi)
    npix     = 12 * NSIDE**2

    print(f"[bcut-build] NSIDE={NSIDE}  npix={npix}")

    def accumulate_zmean(z_vals, ipix_vals, tag):
        """
        Accumulate counts and mean z in each pixel and save to NPZ.
        """
        N     = np.zeros(npix, dtype=np.int64)
        z_sum = np.zeros(npix, dtype=float)

        for zp, p in zip(z_vals, ipix_vals):
            if np.isfinite(zp):
                N[p]     += 1
                z_sum[p] += zp

        zmean       = np.full(npix, np.nan, dtype=float)
        good        = N >= N_MIN_PER_PIX
        zmean[good] = z_sum[good] / N[good]

        print(f"[bcut-build] tag={tag}, |b| > {b_cut_deg:.1f} deg")
        print(f"               pixels with N>0          : {(N>0).sum()} / {npix}")
        print(f"               pixels with N>={N_MIN_PER_PIX}: {good.sum()} / {npix}")

        out = OUT_DIR / f"quaia_zmean_{tag}_bcut{int(b_cut_deg)}_nside{NSIDE}.npz"
        np.savez(out, N=N, zmean=zmean, NSIDE=float(NSIDE))
        print(f"[save-bcut] wrote {out}")

    # 1) Full sample for this b_cut
    accumulate_zmean(z_cut, ipix, "full")

    # 2) Redshift slices (on the |b| and GC-cut sample)
    for (z_lo, z_hi, label, tag) in Z_SLICES:
        sel = (z_cut >= z_lo) & (z_cut < z_hi)
        n_sel = sel.sum()
        print(f"[bcut-build] slice {label}, |b| > {b_cut_deg:.1f}°: kept {n_sel} objects")
    
        if n_sel == 0:
            print(f"[bcut-build] slice {label}: no objects after cuts, skipping")
            continue
    
        accumulate_zmean(z_cut[sel], ipix[sel], tag)

def main():
    print("[quaia-bcut] building <z> maps for BCUT grid:", BCUT_LIST)
    for b_cut_deg in BCUT_LIST:
        build_maps_for_bcut(b_cut_deg)
        
    print("[quaia-bcut] done.")


if __name__ == "__main__":
    main()
