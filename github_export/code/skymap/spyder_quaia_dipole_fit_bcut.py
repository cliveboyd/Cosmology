#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 22 22:47:57 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_dipole_fit_bcut.py  (0V1)

Dipole fits to Quaia <z> maps with Galactic latitude cut |b| > b_cut:

  - Full sample
  - Per redshift slice (Z_BINS)

Reads maps produced by quaia_build_zmean_maps_bcut.py.

Outputs a summary text file:
  OUT_DIR/quaia_dipole_summary_bcut{b_cut}.txt
"""

import numpy            as np
import healpy           as hp
from   pathlib          import Path

from   quaia_config     import OUT_DIR, NSIDE, Z_BINS
from   quaia_dipole_fit import fit_dipole_linear, dipole_summary_from_b

from   quaia_dipole_fit import (
                                fit_dipole_linear,
                                dipole_summary_from_b,
                                unit_vec_gal,
                                L_CMB_DEG,
                                B_CMB_DEG,
                                project_dipole_on_cmb_with_cov,   # new helper
                               )


# Must match the builder
B_CUT_DEG = 20.0

# CMB dipole direction (Galactic, degrees)
L_CMB_DEG = 264.0
B_CMB_DEG = 48.0

def unit_vec_gal(l_deg, b_deg):
    """Unit vector in Galactic coords (l,b in deg)."""
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    x = np.cos(b) * np.cos(l)
    y = np.cos(b) * np.sin(l)
    z = np.sin(b)
    return np.array([x, y, z])


def ang_sep_deg(l1_deg, b1_deg, l2_deg, b2_deg):
    """Angular separation between two (l,b) directions, all in degrees."""
    l1   = np.deg2rad(l1_deg); b1 = np.deg2rad(b1_deg)
    l2   = np.deg2rad(l2_deg); b2 = np.deg2rad(b2_deg)

    x1   = np.cos(b1)*np.cos(l1); y1 = np.cos(b1)*np.sin(l1); z1 = np.sin(b1)
    x2   = np.cos(b2)*np.cos(l2); y2 = np.cos(b2)*np.sin(l2); z2 = np.sin(b2)

    dot  = x1*x2 + y1*y2 + z1*z2
    dot  = np.clip(dot, -1.0, 1.0)
    return np.rad2deg(np.arccos(dot))


def load_map_bundle_bcut(tag):
    path = OUT_DIR / f"quaia_zmean_{tag}_bcut{int(B_CUT_DEG):02d}_nside{NSIDE}.npz"
    data = np.load(path)
    return data["N"], data["zmean"], float(data["NSIDE"])


def analyze_slice(tag, human_label):
    """
    Load b-cut maps for a given tag and fit dipole.
    Returns dict with results, or None if not enough pixels.
    """
    N_map, zmean_map, NSIDE_loaded = load_map_bundle_bcut(tag)
    assert int(NSIDE_loaded) == NSIDE, f"NSIDE mismatch in {tag}: {NSIDE_loaded} vs {NSIDE}"

    # Same cut as before
    good = (N_map >= 3) & np.isfinite(zmean_map)
    pix  = np.where(good)[0]
    if pix.size < 3:
        print(f"[bcut] slice '{human_label}' has N_good={pix.size} < 3; skipping.")
        return dict(
                    tag       = tag,
                    label     = human_label,
                    N_good    = pix.size,
                    amp       = np.nan,
                    l_deg     = np.nan,
                    b_deg     = np.nan,
                    a0        = np.nan,
                    sep_cmb   = np.nan,
                    amp_par   = np.nan,
                    sigma_par = np.nan,
                    frac_par  = np.nan,
                   )

    # Per-pixel counts and redshifts for good pixels
    counts     = N_map[good].astype(float)
    z_good     = zmean_map[good].astype(float)

    # Pixel directions
    nside_int  = int(NSIDE_loaded)
    theta, phi = hp.pix2ang(nside_int, pix)
    x          = np.sin(theta) * np.cos(phi)
    y          = np.sin(theta) * np.sin(phi)
    z          = np.cos(theta)
    n_hat      = np.stack([x, y, z], axis=1)  # (N_good, 3)

    # Mean-subtract ⟨z⟩
    f = z_good - z_good.mean()

    # Weighted by counts
    w = counts
    a0, b_vec, cov                    = fit_dipole_linear(f, n_hat, w=w)
    amp, l_deg, b_deg                 = dipole_summary_from_b(b_vec)
    amp, amp_par, sigma_par, frac_par = project_dipole_on_cmb_with_cov(b_vec, cov)
    sep_cmb                           = ang_sep_deg(l_deg, b_deg, L_CMB_DEG, B_CMB_DEG)

    print(f"\n[bcut slice] {human_label}")
    print(f"  N_good_pix = {pix.size}")
    print(f"  dipole amplitude |b|                 = {amp:.4e}")
    print(f"  dipole direction (l,b)               = ({l_deg:.1f}°, {b_deg:.1f}°)")
    print(f"  a0 (offset, after mean-subtract)     = {a0:.4e}")
    print(f"  component along CMB dir: amp_par     = {amp_par:.4e} "
          f"(sigma_par ≈ {sigma_par:.4e}, frac_par = {frac_par:+.3f})")
    print(f"  separation from CMB dipole           ≈ {sep_cmb:.1f}°")
    

    return dict(
                tag       = tag,
                label     = human_label,
                N_good    = pix.size,
                amp       = amp,
                l_deg     = l_deg,
                b_deg     = b_deg,
                a0        = a0,
                sep_cmb   = sep_cmb,
                amp_par   = amp_par,
                sigma_par = sigma_par,
                frac_par  = frac_par,
               )


def main():
    print(f"[bcut] Quaia <z> dipole fits with |b| > {B_CUT_DEG:.1f} deg")

    results = []

    # Full sample
    results.append(analyze_slice("full", "full sample"))

    # Redshift slices
    for (z_lo, z_hi) in Z_BINS:
        tag   = f"z{z_lo:.2f}_{z_hi:.2f}".replace(".", "p")
        label = f"{z_lo:.2f} ≤ z < {z_hi:.2f}"
        results.append(analyze_slice(tag, label))

    # Save summary table
    out_path = OUT_DIR / f"quaia_dipole_summary_bcut{int(B_CUT_DEG):02d}.txt"
    with out_path.open("w") as f:
        f.write("# tag  label  N_good  amp  l_deg  b_deg  a0  sep_cmb\n")
        for r in results:
            f.write(
                    f"{r['tag']:12s} "
                    f"{r['label']:28s} "
                    f"{r['N_good']:8d} "
                    f"{r['amp']:9.4e} "
                    f"{r['l_deg']:8.3f} "
                    f"{r['b_deg']:8.3f} "
                    f"{r['a0']:9.4e} "
                    f"{r['sep_cmb']:7.2f}\n"
                  )

    print(f"\n[bcut] wrote summary -> {out_path}")


if __name__ == "__main__":
    main()
