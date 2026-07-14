#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
quaia_dipole_residual_z_vs_counts_slices.py  (0V2)

For each redshift slice (Z_BINS from quaia_config):

- Load <z> and counts N per pixel from quaia_zmean_*.npz
- Regress zmean vs N: z ≈ a + b*N (within the slice)
- Build residual map z_resid = zmean - (a + b*N)
- Fit dipoles to:
    * original z (mean-subtracted)
    * residual z (mean-subtracted)
- Report amplitudes, directions, and angles to the CMB.
"""

import numpy   as     np
import healpy  as     hp
from   pathlib import Path

from   quaia_config     import OUT_DIR, NSIDE, Z_BINS
from   quaia_dipole_fit import fit_dipole_linear, dipole_summary_from_b

# CMB dipole direction in Galactic coordinates (degrees)
L_CMB_DEG = 264.0
B_CMB_DEG = 48.0


def ang_sep_deg(l1_deg, b1_deg, l2_deg, b2_deg):
    """Angular separation between two (l,b) directions, all in degrees."""
    l1 = np.deg2rad(l1_deg);    b1 = np.deg2rad(b1_deg)
    l2 = np.deg2rad(l2_deg);    b2 = np.deg2rad(b2_deg)

    x1 = np.cos(b1)*np.cos(l1); y1 = np.cos(b1)*np.sin(l1); z1 = np.sin(b1)
    x2 = np.cos(b2)*np.cos(l2); y2 = np.cos(b2)*np.sin(l2); z2 = np.sin(b2)

    dot = x1*x2 + y1*y2 + z1*z2
    dot = np.clip(dot, -1.0, 1.0)
    return np.rad2deg(np.arccos(dot))


def process_slice(z_lo, z_hi):
    tag   = f"z{z_lo:.2f}_{z_hi:.2f}".replace(".", "p")
    path  = OUT_DIR / f"quaia_zmean_{tag}_nside{NSIDE}.npz"
    data  = np.load(path)

    N_map = data["N"]
    zmean = data["zmean"]

    # same cut as for full sample: N>=3 and finite z
    good  = (N_map >= 3) & np.isfinite(zmean)
    pix   = np.where(good)[0]

    print(f"\n[slice {tag}] {z_lo:.2f} ≤ z < {z_hi:.2f}")
    print(f"  good pixels = {good.sum()} / {len(N_map)}")

    if good.sum() < 10:
        print("  [warn] too few pixels for a meaningful dipole; skipping.")
        return None

    N_good              = N_map[good].astype(float)
    z_good              = zmean[good].astype(float)

    # 1) z ≈ a + b*N within this slice
    A                   = np.vstack([np.ones_like(N_good), N_good]).T
    beta, *_            = np.linalg.lstsq(A, z_good, rcond=None)
    a, b                = beta
    print(f"  z ≈ a + b*N with a={a:.5f}, b={b:.5e}")

    z_model             = a + b * N_good
    z_resid             = z_good - z_model

    # 2) geometry for these pixels
    theta, phi          = hp.pix2ang(NSIDE, pix)
    x                   = np.sin(theta) * np.cos(phi)
    y                   = np.sin(theta) * np.sin(phi)
    z                   = np.cos(theta)
    # IMPORTANT: shape (N_good, 3), not (3, N_good)
    n_hat               = np.vstack([x, y, z]).T

    # 3) original <z> dipole (mean-subtracted)
    z_zero              = z_good - z_good.mean()
    a0_z, b_z, _        = fit_dipole_linear(z_zero, n_hat, w=N_good)
    amp_z, l_z, b_z_deg = dipole_summary_from_b(b_z)
    sep_z_cmb           = ang_sep_deg(l_z, b_z_deg, L_CMB_DEG, B_CMB_DEG)

    # 4) residual <z> dipole (mean-subtracted)
    z_resid_zero        = z_resid - z_resid.mean()
    a0_r, b_r, _        = fit_dipole_linear(z_resid_zero, n_hat, w=N_good)
    amp_r, l_r, b_r_deg = dipole_summary_from_b(b_r)
    sep_r_cmb           = ang_sep_deg(l_r, b_r_deg, L_CMB_DEG, B_CMB_DEG)
    sep_z_resid         = ang_sep_deg(l_z, b_z_deg, l_r, b_r_deg)

    print("  [orig <z> dipole]")
    print(f"    amp = {amp_z:.3e}")
    print(f"    (l,b) = ({l_z:.2f}°, {b_z_deg:.2f}°)")
    print(f"    sep to CMB = {sep_z_cmb:.1f}°")

    print("  [resid <z> dipole]")
    print(f"    amp = {amp_r:.3e}")
    print(f"    (l,b) = ({l_r:.2f}°, {b_r_deg:.2f}°)")
    print(f"    sep to CMB = {sep_r_cmb:.1f}°")
    print(f"    angle(orig,resid) = {sep_z_resid:.1f}°")

    return dict(
                tag         = tag, 
                z_lo        = z_lo, 
                z_hi        = z_hi,
                N_good      = int(good.sum()),
                a           = a, 
                b           = b,
                amp_z       = amp_z, l_z=l_z, b_z=b_z_deg, sep_z_cmb=sep_z_cmb,
                amp_r       = amp_r, l_r=l_r, b_r=b_r_deg, sep_r_cmb=sep_r_cmb,
                sep_z_resid = sep_z_resid,
               )


def main():
    print("[resid-slices] Quaia <z> residual dipoles per redshift bin")
    results = []
    for (z_lo, z_hi) in Z_BINS:
        res = process_slice(z_lo, z_hi)
        if res is not None:
            results.append(res)

    out_path = OUT_DIR / "quaia_dipole_residual_slices.txt"
    with out_path.open("w") as f:
        f.write("# tag  z_lo  z_hi  N_good  a  b  amp_z  l_z  b_z  sep_z_cmb  amp_r  l_r  b_r  sep_r_cmb  sep_z_resid\n")
        for r in results:
            f.write(
                f"{r['tag']:12s}  {r['z_lo']:.2f}  {r['z_hi']:.2f}  {r['N_good']:6d}  "
                f"{r['a']:.5f}  {r['b']:.5e}  "
                f"{r['amp_z']:.3e}  {r['l_z']:.2f}  {r['b_z']:.2f}  {r['sep_z_cmb']:.1f}  "
                f"{r['amp_r']:.3e}  {r['l_r']:.2f}  {r['b_r']:.2f}  {r['sep_r_cmb']:.1f}  {r['sep_z_resid']:.1f}\n"
            )
    print(f"\n[resid-slices] wrote summary -> {out_path}")


if __name__ == "__main__":
    main()
