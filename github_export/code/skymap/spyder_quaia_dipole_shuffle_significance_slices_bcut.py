#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 22 22:57:34 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_dipole_shuffle_significance_slices_bcut.py  (0V1)

Shuffle tests for Quaia <z> dipoles in redshift slices,
using maps with a Galactic latitude cut |b| > 20 deg.

- Uses the precomputed maps:
    quaia_zmean_full_bcut20_nside{NSIDE}.npz
    quaia_zmean_zXXXX_YYYY_bcut20_nside{NSIDE}.npz

- For each z-bin:
    * load map
    * select pixels with N >= 3 and finite zmean
    * fit dipole to original <z> (mean-subtracted)
    * shuffle z among these pixels N_shuffle times and refit dipole
    * estimate p-value = fraction of shuffles with amp >= real amp
    * write a summary text file in OUT_DIR
"""
 
import numpy            as     np
import healpy           as     hp
from   pathlib          import Path

from   quaia_config     import OUT_DIR, NSIDE, Z_BINS
from   quaia_dipole_fit import fit_dipole_linear, dipole_summary_from_b

# CMB dipole direction (Galactic) in degrees
L_CMB_DEG = 264.0
B_CMB_DEG = 48.0

# Number of shuffles per slice
N_SHUFFLE = 500

def ang_sep_deg(l1_deg, b1_deg, l2_deg, b2_deg):
    """Angular separation between two (l,b) directions, all in degrees."""
    l1  = np.deg2rad(l1_deg);    b1 = np.deg2rad(b1_deg)
    l2  = np.deg2rad(l2_deg);    b2 = np.deg2rad(b2_deg)

    x1  = np.cos(b1)*np.cos(l1); y1 = np.cos(b1)*np.sin(l1); z1 = np.sin(b1)
    x2  = np.cos(b2)*np.cos(l2); y2 = np.cos(b2)*np.sin(l2); z2 = np.sin(b2)

    dot = x1*x2 + y1*y2 + z1*z2
    dot = np.clip(dot, -1.0, 1.0)
    return np.rad2deg(np.arccos(dot))


def process_slice(z_lo, z_hi):
    """
    Shuffle-test for a single redshift slice (with bcut20 maps).
    Returns a dict of results.
    """
    tag          = f"z{z_lo:.2f}_{z_hi:.2f}".replace(".", "p")
    fname        = OUT_DIR / f"quaia_zmean_{tag}_bcut20_nside{NSIDE}.npz"

    data         = np.load(fname)
    N_map        = data["N"]
    zmean        = data["zmean"]
    NSIDE_loaded = int(data["NSIDE"])

    # Mask: N>=3, finite z
    good         = (N_map >= 3) & np.isfinite(zmean)
    pix          = np.where(good)[0]

    print(f"[shuffle-bcut] slice {tag}  {z_lo:.2f} ≤ z < {z_hi:.2f}")
    print(f"  good pixels = {good.sum()} / {len(N_map)}")
    if good.sum() < 3:
        print("  [warn] N_good < 3, skipping.")
        return None

    N_good     = N_map[good].astype(float)
    z_good     = zmean[good].astype(float)

    # Pixel directions -> unit vectors n_hat (N,3)
    theta, phi = hp.pix2ang(NSIDE_loaded, pix)
    x          = np.sin(theta) * np.cos(phi)
    y          = np.sin(theta) * np.sin(phi)
    z          = np.cos(theta)
    n_hat      = np.stack([x, y, z], axis=1)

    # Real <z> dipole (mean-subtracted)
    z_zero                       = z_good - z_good.mean()
    a0_real, b_real, _           = fit_dipole_linear(z_zero, n_hat, w=N_good)
    amp_real, l_real, b_real_deg = dipole_summary_from_b(b_real)
    sep_real_cmb                 = ang_sep_deg(l_real, b_real_deg, L_CMB_DEG, B_CMB_DEG)

    print(f"  [real] amp = {amp_real:.3e}, (l,b)=({l_real:.2f}°, {b_real_deg:.2f}°)")
    print(f"         sep to CMB = {sep_real_cmb:.1f}°")

    # Shuffle tests
    amps = np.empty(N_SHUFFLE, dtype=float)

    rng = np.random.default_rng()
    for i in range(N_SHUFFLE):
        z_sh         = rng.permutation(z_good)
        z_sh_zero    = z_sh - z_sh.mean()
        _, b_sh, _   = fit_dipole_linear(z_sh_zero, n_hat, w=N_good)
        amp_sh, _, _ = dipole_summary_from_b(b_sh)
        amps[i]      = amp_sh

    mean_sh = amps.mean()
    std_sh  = amps.std(ddof=1) if amps.size > 1 else 0.0
    p_val   = np.mean(amps >= amp_real)

    print(f"  [shuffles] mean amp = {mean_sh:.3e}, std = {std_sh:.3e}")
    print(f"             p ≈ {p_val:.4f} (fraction with amp >= real)")

    return dict(
        tag       = tag,
        z_lo      = z_lo,
        z_hi      = z_hi,
        N_good    = int(good.sum()),
        amp_real  = amp_real,
        mean_sh   = mean_sh,
        std_sh    = std_sh,
        p_value   = p_val,
        l_real    = l_real,
        b_real    = b_real_deg,
        sep_cmb   = sep_real_cmb,
    )


def main():
    print("[shuffle-bcut] Quaia <z> dipole significance per redshift slice (|b|>20°)")

    results = []
    for (z_lo, z_hi) in Z_BINS:
        res = process_slice(z_lo, z_hi)
        if res is not None:
            results.append(res)

    # Write summary table
    out_path = OUT_DIR / "quaia_dipole_shuffle_slices_bcut20.txt"
    with out_path.open("w") as f:
        f.write("# tag  z_lo  z_hi  N_good  amp_real  mean_sh  std_sh  p_value  l_real  b_real  sep_cmb\n")
        for r in results:
            f.write(
                f"{r['tag']}  {r['z_lo']:.2f}  {r['z_hi']:.2f}  {r['N_good']:6d} "
                f"{r['amp_real']:.3e} {r['mean_sh']:.3e} {r['std_sh']:.3e} "
                f"{r['p_value']:7.4f}  {r['l_real']:6.2f}  {r['b_real']:6.2f} "
                f"{r['sep_cmb']:6.2f}\n"
            )

    print(f"\n[shuffle-bcut] wrote summary -> {out_path}")


if __name__ == "__main__":
    main()
