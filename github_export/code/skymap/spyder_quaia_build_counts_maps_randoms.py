#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 22 20:03:17 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_build_counts_maps_randoms.py  (0V1)

Build per-pixel COUNT maps from the Quaia randoms catalogue.

Outputs:
    OUT_DIR/randoms/quaia_random_counts_nside{NSIDE}.npz

with arrays:
    N      – counts per pixel
    NSIDE  – HEALPix NSIDE (float)
"""

from   pathlib       import Path
import numpy         as     np
import healpy        as     hp
from   astropy.table import Table

from   quaia_config  import NSIDE, N_MIN_PER_PIX, OUT_DIR

# Path to randoms file (you already confirmed this exists)
RANDOM_FITS = Path("/Users/boyde/Downloads/quaia_G20p0_randoms_simple.fits")

# Column names in the randoms catalogue
COL_RA_RAND  = "RA"
COL_DEC_RAND = "DEC"

OUT_DIR_RANDOM = OUT_DIR / "randoms"
OUT_DIR_RANDOM.mkdir(parents=True, exist_ok=True)


def main():
    print("[quaia-rand] building COUNT maps from randoms")

    print(f"[load-rand] loading randoms from {RANDOM_FITS}")
    t = Table.read(RANDOM_FITS)

    ra_deg  = np.array(t[COL_RA_RAND],  dtype=float)
    dec_deg = np.array(t[COL_DEC_RAND], dtype=float)

    npix = hp.nside2npix(NSIDE)

    # Convert to HEALPix angles
    theta = np.radians(90.0 - dec_deg)   # colatitude
    phi   = np.radians(ra_deg)           # longitude
    ipix  = hp.ang2pix(NSIDE, theta, phi)

    # Counts per pixel
    N = np.bincount(ipix, minlength=npix)

    print(f"[build-rand] NSIDE={NSIDE}  npix={npix}")
    print(f"[build-rand] pixels with N>0: {np.count_nonzero(N>0)} / {npix}")
    print(f"[build-rand] pixels with N>={N_MIN_PER_PIX}: {np.count_nonzero(N>=N_MIN_PER_PIX)} / {npix}")

    out_path = OUT_DIR_RANDOM / f"quaia_random_counts_nside{NSIDE}.npz"
    np.savez(out_path, N=N, NSIDE=float(NSIDE))

    print(f"[save-rand] wrote {out_path}")
    print("[quaia-rand] done.")


if __name__ == "__main__":
    main()
