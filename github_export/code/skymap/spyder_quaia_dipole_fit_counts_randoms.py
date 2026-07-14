#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 22 20:08:59 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_dipole_fit_counts_randoms.py  (0V1)

Fit a dipole to COUNT maps built from the Quaia randoms catalogue.

Reads:
    OUT_DIR/randoms/quaia_random_counts_nside{NSIDE}.npz

Fits:
    N - <N> = a0 + bx * nx + by * ny + bz * nz

using only pixels with N >= N_MIN_PER_PIX.

Reports amplitude, direction, and separation from the CMB dipole.
"""

import numpy   as     np
import healpy  as     hp
from   pathlib import Path

from   quaia_config import OUT_DIR, NSIDE, N_MIN_PER_PIX

OUT_DIR_RANDOM = OUT_DIR / "randoms"

# CMB dipole direction (galactic, radians)
L_CMB = np.deg2rad(264.0)
B_CMB = np.deg2rad(48.0)


def ang_sep(l1, b1, l2, b2):
    """Angular separation between two (l,b) directions in radians."""
    x1 = np.cos(b1) * np.cos(l1)
    y1 = np.cos(b1) * np.sin(l1)
    z1 = np.sin(b1)

    x2 = np.cos(b2) * np.cos(l2)
    y2 = np.cos(b2) * np.sin(l2)
    z2 = np.sin(b2)

    dot = x1*x2 + y1*y2 + z1*z2
    dot = np.clip(dot, -1.0, 1.0)
    return np.arccos(dot)


def fit_dipole_from_counts(N_map):
    """
    Fit a dipole to counts per pixel.

    We fit:
        N - <N> = a0 + bx * nx + by * ny + bz * nz

    where <N> is the mean over valid pixels.
    """
    npix = len(N_map)

    valid = (N_map >= N_MIN_PER_PIX)
    if not np.any(valid):
        raise RuntimeError("No valid pixels for dipole fit.")

    ipix = np.arange(npix)[valid]
    N    = N_map[valid]

    theta, phi = hp.pix2ang(NSIDE, ipix)
    nx = np.sin(theta) * np.cos(phi)
    ny = np.sin(theta) * np.sin(phi)
    nz = np.cos(theta)

    N_mean = np.mean(N)
    N_centered = N - N_mean

    X = np.column_stack([np.ones_like(nx), nx, ny, nz])
    coeffs, *_ = np.linalg.lstsq(X, N_centered, rcond=None)
    a0, bx, by, bz = coeffs

    amp = np.sqrt(bx**2 + by**2 + bz**2)

    if amp > 0:
        nx_d = bx / amp
        ny_d = by / amp
        nz_d = bz / amp
        b = np.arcsin(nz_d)
        l = np.arctan2(ny_d, nx_d)
        if l < 0:
            l += 2.0 * np.pi
    else:
        l = 0.0
        b = 0.0

    N_valid = np.count_nonzero(valid)
    return dict(a0=a0, amp=amp, l=l, b=b, N_mean=N_mean, N_valid=N_valid)


def main():
    print("[quaia-rand] dipole fit from COUNT map (randoms)")

    path = OUT_DIR_RANDOM / f"quaia_random_counts_nside{NSIDE}.npz"
    data = np.load(path, allow_pickle=False)

    N_map = data["N"]
    NSIDE_loaded = int(data["NSIDE"])
    assert NSIDE_loaded == NSIDE, f"NSIDE mismatch: {NSIDE_loaded} vs {NSIDE}"

    fit = fit_dipole_from_counts(N_map)
    amp = fit["amp"]
    l   = fit["l"]
    b   = fit["b"]
    a0  = fit["a0"]
    N_mean  = fit["N_mean"]
    N_valid = fit["N_valid"]

    theta = ang_sep(l, b, L_CMB, B_CMB)
    sep_deg = np.rad2deg(theta)

    print(f"  N_valid_pix = {N_valid}")
    print(f"  mean N per valid pixel = {N_mean:.1f}")
    print(f"  dipole amplitude |b| = {amp:.4e} (counts)")
    print(f"  dipole direction (l,b) = ({np.rad2deg(l):5.1f}°, {np.rad2deg(b):5.1f}°)")
    print(f"  separation from CMB dipole ≈ {sep_deg:5.1f}°")
    print(f"  a0 (offset term) = {a0:.4e}")

    out_path = OUT_DIR_RANDOM / "quaia_randoms_counts_dipole.txt"
    with out_path.open("w") as f:
        f.write("# N_valid  N_mean  amp_counts  l_deg  b_deg  sep_deg  a0\n")
        f.write(f"{N_valid:d}  {N_mean:.3f}  {amp:.6e}  "
                f"{np.rad2deg(l):.3f}  {np.rad2deg(b):.3f}  "
                f"{sep_deg:.3f}  {a0:.6e}\n")

    print(f"\n[save-rand] wrote {out_path}")
    print("[quaia-rand] done.")


if __name__ == "__main__":
    main()
