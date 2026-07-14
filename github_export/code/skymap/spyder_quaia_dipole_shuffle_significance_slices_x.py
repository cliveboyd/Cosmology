#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 22 21:57:16 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_dipole_shuffle_significance_slices.py  (0V1)

For each Quaia redshift slice:
  - Load <z> and N HEALPix maps (already built by quaia_build_zmean_maps.py)
  - Compute the actual <z> dipole amplitude.
  - Shuffle z across pixels (keeping N fixed) many times and refit dipole.
  - Estimate a p-value: fraction of shuffles with amp >= real amp.

Uses:
  - quaia_config: OUT_DIR, NSIDE, Z_BINS
  - quaia_dipole_fit: fit_dipole_linear, dipole_summary_from_b
"""

import numpy          as np
import healpy         as hp
from   pathlib        import Path

from quaia_config     import OUT_DIR, NSIDE, Z_BINS
from quaia_dipole_fit import fit_dipole_linear, dipole_summary_from_b


def ang_sep_deg(l1_deg, b1_deg, l2_deg, b2_deg):
    """Angular separation between two (l,b) directions, all in degrees."""
    l1 = np.deg2rad(l1_deg); b1 = np.deg2rad(b1_deg)
    l2 = np.deg2rad(l2_deg); b2 = np.deg2rad(b2_deg)

    x1 = np.cos(b1)*np.cos(l1); y1 = np.cos(b1)*np.sin(l1); z1 = np.sin(b1)
    x2 = np.cos(b2)*np.cos(l2); y2 = np.cos(b2)*np.sin(l2); z2 = np.sin(b2)

    dot = x1*x2 + y1*y2 + z1*z2
    dot = np.clip(dot, -1.0, 1.0)
    return np.rad2deg(np.arccos(dot))


def process_slice(z_lo, z_hi, n_trials=500, seed=12345):
    """
    Run shuffle significance for one redshift slice.
    Returns dict with results, or None if not enough pixels.
    """
    tag  = f"z{z_lo:.2f}_{z_hi:.2f}".replace(".", "p")
    path = OUT_DIR / f"quaia_zmean_{tag}_nside{NSIDE}.npz"
    if not path.exists():
        print(f"[shuffle] missing map for slice {tag}: {path}")
        return None

    data  = np.load(path)
    N_map = data["N"]
    zmean = data["zmean"]

    # Same cut as your latest runs: N >= 3, finite z
    good  = (N_map >= 3) & np.isfinite(zmean)
    pix   = np.where(good)[0]

    if pix.size < 3:
        print(f"[shuffle] slice {tag} has only {pix.size} good pixels; skipping.")
        return None

    N_good = N_map[good].astype(float)
    z_good = zmean[good].astype(float)

    print(f"\n[shuffle] slice {tag}  {z_lo:.2f} ≤ z < {z_hi:.2f}")
    print(f"  good pixels = {pix.size} / {N_map.size}")

    # Pixel directions
    theta, phi = hp.pix2ang(NSIDE, pix)
    x          = np.sin(theta) * np.cos(phi)
    y          = np.sin(theta) * np.sin(phi)
    z          = np.cos(theta)
    n_hat      = np.stack([x, y, z], axis=1)  # (N_good, 3)

    # --- real dipole ---
    z_zero                       = z_good - z_good.mean()
    a0_real, b_real, _           = fit_dipole_linear(z_zero, n_hat, w=N_good)
    amp_real, l_real, b_real_deg = dipole_summary_from_b(b_real)

    print(f"  [real] amp = {amp_real:.3e}, (l,b)=({l_real:.2f}°, {b_real_deg:.2f}°)")

    # --- shuffles ---
    rng  = np.random.default_rng(seed)
    amps = np.empty(n_trials, dtype=float)

    for i in range(n_trials):
        z_shuffled     = rng.permutation(z_good)
        z_sh_zero      = z_shuffled - z_shuffled.mean()
        a0_sh, b_sh, _ = fit_dipole_linear(z_sh_zero, n_hat, w=N_good)
        amp_sh, _, _   = dipole_summary_from_b(b_sh)
        amps[i]        = amp_sh

    mean_sh = amps.mean()
    std_sh  = amps.std(ddof=1)

    # one-sided p-value: P(amp_sh >= amp_real)
    p = (1.0 + np.sum(amps >= amp_real)) / (n_trials + 1.0)

    print(f"  [shuffles] mean amp = {mean_sh:.3e}, std = {std_sh:.3e}")
    print(f"             p ≈ {p:.4f} (fraction of shuffles with amp >= real)")

    return dict(
                tag      = tag,
                z_lo     = z_lo,
                z_hi     = z_hi,
                N_good   = pix.size,
                amp_real = amp_real,
                mean_sh  = mean_sh,
                std_sh   = std_sh,
                p_value  = p,
                l_real   = l_real,
                b_real   = b_real_deg,
               )


def main():
    print("[shuffle] Quaia <z> dipole significance per redshift slice")
    results = []

    for (z_lo, z_hi) in Z_BINS:
        res = process_slice(z_lo, z_hi, n_trials=500, seed=12345)
        if res is not None:
            results.append(res)

    # Save summary table
    out_path = OUT_DIR / "quaia_dipole_shuffle_slices.txt"
    with out_path.open("w") as f:
        f.write("# tag  z_lo  z_hi  N_good  amp_real  mean_sh  std_sh  p_value  l_real  b_real\n")
        for r in results:
            f.write(
                f"{r['tag']:10s} "
                f"{r['z_lo']:5.2f} {r['z_hi']:5.2f} "
                f"{r['N_good']:6d} "
                f"{r['amp_real']:9.3e} "
                f"{r['mean_sh']:9.3e} "
                f"{r['std_sh']:9.3e} "
                f"{r['p_value']:7.4f} "
                f"{r['l_real']:7.2f} "
                f"{r['b_real']:7.2f}\n"
            )

    print(f"\n[shuffle] wrote summary -> {out_path}")


if __name__ == "__main__":
    main()
