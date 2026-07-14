#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_dipole_jackknife_bcut20.py  (0V1)

Jackknife robustness tests for Quaia <z> dipoles with |b| > 20 deg.

For each redshift slice (including full sample):
  - Use maps quaia_zmean_{tag}_bcut20_nside{NSIDE}.npz
  - Compute the full-sample dipole
  - Perform jackknife runs by LEAVING OUT:
        * each Galactic hemisphere (b>0, b<0)
        * each of 4 longitude quadrants:
              Q1: 0<=l<90, Q2: 90<=l<180,
              Q3: 180<=l<270, Q4: 270<=l<360
  - Compare amplitudes and directions:
        * Δamp = (amp_jk - amp_full)/amp_full
        * angle(full, jk)
        * angle(jk, CMB)

%run /Users/boyde/.spyder-py3/quaia_dipole_jackknife_bcut20.py

Output:
  - Human-readable summary  -> quaia_dipole_jackknife_bcut20.txt
  


"""

import numpy            as     np
from   pathlib          import Path
import healpy           as     hp

from   quaia_config     import OUT_DIR, NSIDE, Z_BINS
from   quaia_dipole_fit import fit_dipole_linear, dipole_summary_from_b

N_MIN_PER_PIX = 3
BCUT          = 20.0  # deg; must match the builder

L_CMB_DEG     = 264.0
B_CMB_DEG     = 48.0


def ang_sep_deg(l1_deg, b1_deg, l2_deg, b2_deg):
    l1 = np.deg2rad(l1_deg);    b1 = np.deg2rad(b1_deg)
    l2 = np.deg2rad(l2_deg);    b2 = np.deg2rad(b2_deg)

    x1 = np.cos(b1)*np.cos(l1); y1 = np.cos(b1)*np.sin(l1); z1 = np.sin(b1)
    x2 = np.cos(b2)*np.cos(l2); y2 = np.cos(b2)*np.sin(l2); z2 = np.sin(b2)

    dot = x1*x2 + y1*y2 + z1*z2
    dot = np.clip(dot, -1.0, 1.0)
    return np.rad2deg(np.arccos(dot))


def load_map(tag):
    path = OUT_DIR / f"quaia_zmean_{tag}_bcut{int(BCUT)}_nside{NSIDE}.npz"
    data = np.load(path, allow_pickle=False)
    return data["N"], data["zmean"]


def compute_full_dipole(N_map, zmean):
    good              = (N_map >= N_MIN_PER_PIX) & np.isfinite(zmean)
    pix               = np.where(good)[0]
    N_good            = N_map[good].astype(float)
    z_good            = zmean[good].astype(float)

    theta, phi        = hp.pix2ang(NSIDE, pix)
    x                 = np.sin(theta)*np.cos(phi)
    y                 = np.sin(theta)*np.sin(phi)
    z                 = np.cos(theta)
    n_hat             = np.column_stack([x, y, z])

    f                 = z_good - z_good.mean()
    a0, b_vec, _      = fit_dipole_linear(f, n_hat, w=N_good)
    amp, l_deg, b_deg = dipole_summary_from_b(b_vec)
    sep_cmb           = ang_sep_deg(l_deg, b_deg, L_CMB_DEG, B_CMB_DEG)

    return dict(
        good_mask = good,
        pix       = pix,
        N_good    = N_good,
        z_good    = z_good,
        n_hat     = n_hat,
        amp       = amp,
        l_deg     = l_deg,
        b_deg     = b_deg,
        sep_cmb   = sep_cmb,
    )


def jackknife_regions(l_deg, b_deg):
    """
    Given per-pixel (l,b) in degrees, return a dict of region masks
    (on the full pixel grid) indicating which pixels belong to each region.
    Here we build hemisphere and quadrant masks.

    We assume l in [0, 360), b in [-90, 90].
    """
    # Hemispheres
    regions = {
        "leave_N_hemisphere": (b_deg >= 0.0),   # pixels *inside* this region will be dropped
        "leave_S_hemisphere": (b_deg <  0.0),
    }

    # Longitudinal quadrants, ignoring exact boundary issues (rare)
    # Q1: 0<=l<90, Q2: 90<=l<180, Q3: 180<=l<270, Q4: 270<=l<360
    Q1 = (l_deg >=   0.0) & (l_deg <  90.0)
    Q2 = (l_deg >=  90.0) & (l_deg < 180.0)
    Q3 = (l_deg >= 180.0) & (l_deg < 270.0)
    Q4 = (l_deg >= 270.0) & (l_deg < 360.0)

    regions.update({
        "leave_Q1_0-90"   : Q1,
        "leave_Q2_90-180" : Q2,
        "leave_Q3_180-270": Q3,
        "leave_Q4_270-360": Q4,
    })

    return regions


def fit_jackknife(N_map, zmean, full_info):
    """
    Perform jackknife fits by leaving out each region in turn.
    full_info: dict from compute_full_dipole
    """
    # Recompute directions on the FULL sky grid (for region definition)
    npix               = N_map.size
    pix_all            = np.arange(npix)
    theta_all, phi_all = hp.pix2ang(NSIDE, pix_all)
    l_all_deg          = np.rad2deg(phi_all)
    b_all_deg          = 90.0 - np.rad2deg(theta_all)

    reg_masks          = jackknife_regions(l_all_deg, b_all_deg)

    results = []

    for name, reg_mask in reg_masks.items():
        # Use pixels that are good AND NOT in the region being left out
        good_full = full_info["good_mask"]
        use       = good_full & (~reg_mask)

        pix_use   = np.where(use)[0]
        if pix_use.size < 3:
            print(f"  [JK {name}] too few pixels ({pix_use.size}); skipping")
            continue

        N_use      = N_map[use].astype(float)
        z_use      = zmean[use].astype(float)

        # Directions for these pixels
        theta, phi = hp.pix2ang(NSIDE, pix_use)
        x          = np.sin(theta)*np.cos(phi)
        y          = np.sin(theta)*np.sin(phi)
        z          = np.cos(theta)
        n_hat      = np.column_stack([x, y, z])
        f          = z_use - z_use.mean()
        
        a0_jk, b_jk, _         = fit_dipole_linear(f, n_hat, w=N_use)
        amp_jk, l_jk, b_jk_deg = dipole_summary_from_b(b_jk)
        sep_cmb_jk             = ang_sep_deg(l_jk, b_jk_deg, L_CMB_DEG, B_CMB_DEG)
        sep_full               = ang_sep_deg(l_jk, b_jk_deg, full_info["l_deg"], full_info["b_deg"])

        d_amp_frac             = (amp_jk - full_info["amp"]) / full_info["amp"]

        print(f"  [JK {name}] N_use={pix_use.size}, amp={amp_jk:.3e} "
              f"(Δamp/amp_full={d_amp_frac:.2%}), "
              f"(l,b)=({l_jk:.1f}°, {b_jk_deg:.1f}°), "
              f"sep(full)={sep_full:.1f}°, sep(CMB)={sep_cmb_jk:.1f}°")

        results.append(dict(
            region     = name,
            N_use      = int(pix_use.size),
            amp_jk     = amp_jk,
            l_jk       = l_jk,
            b_jk_deg   = b_jk_deg,
            sep_full   = sep_full,
            sep_cmb    = sep_cmb_jk,
            d_amp_frac = d_amp_frac,
        ))

    return results


def main():
    print("[JK-bcut20] Quaia <z> jackknife tests with |b| > 20 deg")
    out_lines = []
    out_lines.append("# tag  region  N_use  amp_full  amp_jk  d_amp_frac  "
                     "l_full  b_full  l_jk  b_jk  sep_full  sep_cmb_full  sep_cmb_jk")

    # Full sample + z-slices
    all_tags = [("full", "full sample")]
    for (z_lo, z_hi) in Z_BINS:
        tag   = f"z{z_lo:.2f}_{z_hi:.2f}".replace(".", "p")
        label = f"{z_lo:.2f} ≤ z < {z_hi:.2f}"
        all_tags.append((tag, label))

    for tag, label in all_tags:
        print(f"\n[JK-bcut20] slice {tag}  {label}")
        N_map, zmean = load_map(tag)

        full = compute_full_dipole(N_map, zmean)
        print(f"  [full] N_good={len(full['N_good'])}, "
              f"amp={full['amp']:.3e}, (l,b)=({full['l_deg']:.1f}°, {full['b_deg']:.1f}°), "
              f"sep(CMB)={full['sep_cmb']:.1f}°")

        jk_res = fit_jackknife(N_map, zmean, full)

        for r in jk_res:
            out_lines.append(
                f"{tag:12s}  {r['region']:16s}  {r['N_use']:6d}  "
                f"{full['amp']:.3e}  {r['amp_jk']:.3e}  {r['d_amp_frac']:.3f}  "
                f"{full['l_deg']:.2f}  {full['b_deg']:.2f}  "
                f"{r['l_jk']:.2f}  {r['b_jk_deg']:.2f}  "
                f"{r['sep_full']:.2f}  {full['sep_cmb']:.2f}  {r['sep_cmb']:.2f}"
            )

    out_txt = OUT_DIR / "quaia_dipole_jackknife_bcut20.txt"
    out_txt.write_text("\n".join(out_lines) + "\n")
    print(f"\n[JK-bcut20] wrote summary -> {out_txt}")


if __name__ == "__main__":
    main()
