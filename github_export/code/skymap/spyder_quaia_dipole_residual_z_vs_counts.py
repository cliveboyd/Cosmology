#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
quaia_dipole_residual_z_vs_counts.py  (0V2)

- Load full-sample <z> map (Quaia).
- Regress zmean vs counts N per pixel: zmean ≈ a + b * N.
- Build residual map z_resid = zmean - (a + b*N).
- Fit a dipole to z_resid and compare:
    - to original <z> dipole
    - to CMB dipole
"""

import  numpy           as     np
from   pathlib          import Path
import healpy           as     hp

from   quaia_config     import OUT_DIR, NSIDE
from   quaia_dipole_fit import fit_dipole_linear, dipole_summary_from_b

# CMB dipole direction in Galactic coordinates (degrees)
L_CMB_DEG = 264.0
B_CMB_DEG = 48.0


def ang_sep_deg(l1_deg, b1_deg, l2_deg, b2_deg):
    """Angular separation between two (l,b) directions, all in degrees."""
    l1 = np.deg2rad(l1_deg); b1 = np.deg2rad(b1_deg)
    l2 = np.deg2rad(l2_deg); b2 = np.deg2rad(b2_deg)

    x1 = np.cos(b1)*np.cos(l1); y1 = np.cos(b1)*np.sin(l1); z1 = np.sin(b1)
    x2 = np.cos(b2)*np.cos(l2); y2 = np.cos(b2)*np.sin(l2); z2 = np.sin(b2)

    dot = x1*x2 + y1*y2 + z1*z2
    dot = np.clip(dot, -1.0, 1.0)
    return np.rad2deg(np.arccos(dot))


def main():
    path      = OUT_DIR / f"quaia_zmean_full_nside{NSIDE}.npz"
    data      = np.load(path)

    N_map     = data["N"]
    zmean     = data["zmean"]

    # Mask: finite z, N >= 3 (match your latest choice)
    good      = (N_map >= 3) & np.isfinite(zmean)
    pix       = np.where(good)[0]

    print(f"[resid] good pixels: {good.sum()} / {len(N_map)}")

    N_good     = N_map[good].astype(float)
    z_good     = zmean[good].astype(float)

    # --- 1) Regress zmean on N: z ≈ a + b * N
    A          = np.vstack([np.ones_like(N_good), N_good]).T
    beta, *_   = np.linalg.lstsq(A, z_good, rcond=None)
    a, b       = beta
    print(f"[resid] z ≈ a + b*N with a={a:.5f}, b={b:.5e}")

    z_model    = a + b * N_good
    z_resid    = z_good - z_model

    # --- 2) Build n_hat for these pixels
    theta, phi = hp.pix2ang(NSIDE, pix)  # theta colatitude, phi longitude (radians)

    # Unit vectors for each pixel
    x          = np.sin(theta) * np.cos(phi)
    y          = np.sin(theta) * np.sin(phi)
    z          = np.cos(theta)
    
    n_hat      = np.vstack([x, y, z]).T  # shape (N_good, 3)


    # --- 3) Fit dipole to original <z> (after mean-subtraction)
    z_zero              = z_good - z_good.mean()
    a0_z, b_z, _        = fit_dipole_linear(z_zero, n_hat, w=N_good)
    amp_z, l_z, b_z_deg = dipole_summary_from_b(b_z)

    # --- 4) Fit dipole to residual <z> (after mean-subtraction)
    z_resid_zero        = z_resid - z_resid.mean()
    a0_r, b_r, _        = fit_dipole_linear(z_resid_zero, n_hat, w=N_good)
    amp_r, l_r, b_r_deg = dipole_summary_from_b(b_r)

    # --- 5) Compare angles
    sep_z_resid         = ang_sep_deg(l_z, b_z_deg, l_r, b_r_deg)
    sep_z_cmb           = ang_sep_deg(l_z, b_z_deg, L_CMB_DEG, B_CMB_DEG)
    sep_r_cmb           = ang_sep_deg(l_r, b_r_deg, L_CMB_DEG, B_CMB_DEG)

    print("\n[original <z> dipole]")
    print(f"  amp = {amp_z:.3e}")
    print(f"  (l,b) = ({l_z:.2f}°, {b_z_deg:.2f}°)")
    print(f"  sep to CMB = {sep_z_cmb:.1f}°")

    print("\n[residual <z> dipole after removing N-trend]")
    print(f"  amp = {amp_r:.3e}")
    print(f"  (l,b) = ({l_r:.2f}°, {b_r_deg:.2f}°)")
    print(f"  sep to CMB = {sep_r_cmb:.1f}°")
    print(f"  angle between original and residual = {sep_z_resid:.1f}°")


if __name__ == "__main__":
    main()
