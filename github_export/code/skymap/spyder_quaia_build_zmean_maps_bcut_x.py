#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 22 22:46:23 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_build_zmean_maps_bcut.py  (0V1)

Build <z> HEALPix maps for Quaia with a Galactic latitude cut |b| > b_cut:

  - Full sample
  - Redshift slices defined in quaia_config.Z_BINS

Outputs:
  OUT_DIR/quaia_zmean_full_bcut{b_cut}_nside{NSIDE}.npz
  OUT_DIR/quaia_zmean_z{...}_bcut{b_cut}_nside{NSIDE}.npz
"""

import numpy as np
import healpy as hp
from   pathlib import Path
from   astropy.table import Table

from   quaia_config import (
                            QUAIA_FITS,
                            OUT_DIR,
                            NSIDE,
                            Z_BINS,
                            COL_RA,
                            COL_DEC,
                            COL_Z,
                           )

# --- set your cut here ---
B_CUT_DEG = 20.0   # e.g. keep |b| > 20 deg; bump to 30.0 if you want


def build_map(tag, mask, b_cut_deg):
    """Build <z> and counts map for a given boolean mask on the full table."""
    print(f"[build-bcut] tag={tag}, |b| > {b_cut_deg:.1f} deg")

    t         = Table.read(QUAIA_FITS)
    ra        = np.array(t[COL_RA], dtype=float)[mask]
    dec       = np.array(t[COL_DEC], dtype=float)[mask]
    z         = np.array(t[COL_Z], dtype=float)[mask]

    # HEALPix angles
    ra_rad    = np.deg2rad(ra)
    dec_rad   = np.deg2rad(dec)
    theta     = 0.5 * np.pi - dec_rad
    phi       = ra_rad

    npix      = hp.nside2npix(NSIDE)
    N_map     = np.zeros(npix, dtype=int)
    z_sum_map = np.zeros(npix, dtype=float)

    ipix = hp.ang2pix(NSIDE, theta, phi)

    for ix, zz in zip(ipix, z):
        N_map[ix]     += 1
        z_sum_map[ix] += zz

    zmean = np.full(npix, np.nan, dtype=float)
    mask_nonzero = N_map > 0
    zmean[mask_nonzero] = z_sum_map[mask_nonzero] / N_map[mask_nonzero]

    out_path = OUT_DIR / f"quaia_zmean_{tag}_bcut{int(b_cut_deg):02d}_nside{NSIDE}.npz"
    np.savez(out_path, N=N_map, zmean=zmean, NSIDE=float(NSIDE),
             b_cut_deg=float(b_cut_deg))
    print(f"[save-bcut] wrote {out_path}")


def main():
    print(f"[quaia-bcut] building <z> maps with |b| > {B_CUT_DEG:.1f} deg")
    t = Table.read(QUAIA_FITS)

    # Galactic latitude column in Quaia is 'b'
    b = np.array(t["b"], dtype=float)
    z_all = np.array(t[COL_Z], dtype=float)

    # --- full-sample mask ---
    mask_full = np.isfinite(z_all) & (np.abs(b) > B_CUT_DEG)
    print(f"[build-bcut] full: kept {mask_full.sum()} / {len(z_all)} objects")
    build_map("full", mask_full, B_CUT_DEG)

    # --- redshift slices ---
    for (z_lo, z_hi) in Z_BINS:
        tag = f"z{z_lo:.2f}_{z_hi:.2f}".replace(".", "p")
        z_mask = (z_all >= z_lo) & (z_all < z_hi)
        mask   = z_mask & (np.abs(b) > B_CUT_DEG) & np.isfinite(z_all)
        print(f"[build-bcut] slice {tag}: kept {mask.sum()} / {len(z_all)} objects")
        build_map(tag, mask, B_CUT_DEG)

    print("[quaia-bcut] done.")


if __name__ == "__main__":
    main()
