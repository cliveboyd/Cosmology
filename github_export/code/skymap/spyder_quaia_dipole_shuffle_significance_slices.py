#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_dipole_shuffle_significance_slices.py  (0V2)

- For each redshift slice in Z_BINS:
  * Load the <z> HEALPix map from Quaia.
  * Compute the real dipole amplitude + direction.
  * Generate many shuffled realisations by shuffling z among pixels.
  * Fit dipole for each shuffle; build full amplitude distribution.
  * Compute p-value (fraction of shuffles with amp >= real amp).
  * Save:
      - text summary      -> quaia_dipole_shuffle_slices.txt
      - per-slice .npz    -> quaia_dipole_shuffle_{tag}_amps.npz
        containing the amplitude array and (optionally) directions.

Use NSHUFF to control how many shuffles you want (default 5000).

%run /Users/boyde/.spyder-py3/quaia_dipole_shuffle_significance_slices.py

Output:::
    Updated quaia_dipole_shuffle_slices.txt (summary),

    Per-slice NPZ files quaia_dipole_shuffle_{tag}_amps.npz containing the 
    full amp distributions for later plotting (histograms, CDFs, etc).
    
"""

import numpy as np
from pathlib import Path
import healpy as hp

from quaia_config import OUT_DIR, NSIDE, Z_BINS
from quaia_dipole_fit import fit_dipole_linear, dipole_summary_from_b

# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------
NSHUFF          = 5000        # increase/decrease if needed
N_MIN_PER_PIX   = 3           # match your current N cut
RNG_SEED_BASE   = 12345       # base seed for reproducibility

# CMB dipole (deg) – just for context in logs / optional diagnostics
L_CMB_DEG       = 264.0
B_CMB_DEG       = 48.0


def ang_sep_deg(l1_deg, b1_deg, l2_deg, b2_deg):
    """Angular separation between two (l,b) directions, all in degrees."""
    l1 = np.deg2rad(l1_deg); b1 = np.deg2rad(b1_deg)
    l2 = np.deg2rad(l2_deg); b2 = np.deg2rad(b2_deg)

    x1 = np.cos(b1)*np.cos(l1); y1 = np.cos(b1)*np.sin(l1); z1 = np.sin(b1)
    x2 = np.cos(b2)*np.cos(l2); y2 = np.cos(b2)*np.sin(l2); z2 = np.sin(b2)

    dot = x1*x2 + y1*y2 + z1*z2
    dot = np.clip(dot, -1.0, 1.0)
    return np.rad2deg(np.arccos(dot))


def load_slice_map(tag):
    """
    Load N and zmean arrays for a given slice tag.
    tag in {"full", "z0p10_0p50", ...}
    """
    path = OUT_DIR / f"quaia_zmean_{tag}_nside{NSIDE}.npz"
    data = np.load(path, allow_pickle=False)
    return data["N"], data["zmean"]


def analyze_slice(z_lo, z_hi, idx):
    """
    Do shuffle analysis for a single redshift slice.
    z_lo, z_hi: redshift bounds
    idx: index of slice in Z_BINS (for seed offset)
    """
    # Build tag and label to match your existing naming scheme
    if z_lo is None:
        tag   = "full"
        label = "full sample"
    else:
        tag   = f"z{z_lo:.2f}_{z_hi:.2f}".replace(".", "p")
        label = f"{z_lo:.2f} ≤ z < {z_hi:.2f}"

    # Load map
    N_map, zmean = load_slice_map(tag)

    # Good pixels
    good  = (N_map >= N_MIN_PER_PIX) & np.isfinite(zmean)
    pix   = np.where(good)[0]
    N_good = N_map[good].astype(float)
    z_good = zmean[good].astype(float)

    print(f"\n[shuffle] slice {tag}  {label}")
    print(f"  good pixels = {good.sum()} / {len(N_map)}")

    # Healpix directions
    theta, phi = hp.pix2ang(NSIDE, pix)  # radians
    x = np.sin(theta)*np.cos(phi)
    y = np.sin(theta)*np.sin(phi)
    z = np.cos(theta)
    n_hat = np.column_stack([x, y, z])   # (N_good,3)

    # --- Real dipole ---
    f_real         = z_good - z_good.mean()
    a0_real, b_real, _ = fit_dipole_linear(f_real, n_hat, w=N_good)
    amp_real, l_real, b_real_deg = dipole_summary_from_b(b_real)
    sep_cmb_real   = ang_sep_deg(l_real, b_real_deg, L_CMB_DEG, B_CMB_DEG)

    print(f"  [real] amp = {amp_real:.3e}, (l,b)=({l_real:.2f}°, {b_real_deg:.2f}°)")
    print(f"         sep to CMB = {sep_cmb_real:.1f}°")

    # --- Shuffles ---
    rng   = np.random.default_rng(RNG_SEED_BASE + idx)
    amps  = np.empty(NSHUFF, dtype=float)
    # (optionally store directions if you want later)
    # l_sh  = np.empty(NSHUFF, dtype=float)
    # b_sh  = np.empty(NSHUFF, dtype=float)

    for k in range(NSHUFF):
        z_shuff       = rng.permutation(z_good)
        f_sh          = z_shuff - z_shuff.mean()
        _, b_shuff, _ = fit_dipole_linear(f_sh, n_hat, w=N_good)
        amp_k, _, _   = dipole_summary_from_b(b_shuff)
        amps[k]       = amp_k

        if (k+1) % max(1, NSHUFF // 10) == 0:
            print(f"    shuffle {k+1}/{NSHUFF}", end="\r")

    print("")  # newline after progress

    mean_amp = amps.mean()
    std_amp  = amps.std(ddof=1)
    p_value  = np.mean(amps >= amp_real)

    print(f"  [shuffles] mean amp = {mean_amp:.3e}, std = {std_amp:.3e}")
    print(f"             p ≈ {p_value:.4f} (fraction with amp >= real)")

    # Save amplitude distribution for this slice
    out_npz = OUT_DIR / f"quaia_dipole_shuffle_{tag}_amps.npz"
    np.savez(
                out_npz,
                amps        = amps,
                amp_real    = amp_real,
                l_real      = l_real,
                b_real      = b_real_deg,
                z_lo        = z_lo if z_lo is not None else -1.0,
                z_hi        = z_hi if z_hi is not None else -1.0,
            )
    print(f"  [save] amplitude distribution -> {out_npz}")

    # Return for summary
    return dict(
                tag      = tag,
                z_lo     = z_lo if z_lo is not None else 0.0,
                z_hi     = z_hi if z_hi is not None else 0.0,
                N_good   = int(good.sum()),
                amp_real = amp_real,
                mean_sh  = mean_amp,
                std_sh   = std_amp,
                p_value  = p_value,
                l_real   = l_real,
                b_real   = b_real_deg,
               )


def main():
    print("[shuffle] Quaia <z> dipole significance per redshift slice")

    results = []

    # Full sample first
    results.append(analyze_slice(None, None, idx=0))

    # Redshift slices
    for i, (z_lo, z_hi) in enumerate(Z_BINS, start=1):
        results.append(analyze_slice(z_lo, z_hi, idx=i))

    # Write summary table
    out_txt = OUT_DIR / "quaia_dipole_shuffle_slices.txt"
    with out_txt.open("w") as f:
        f.write("# tag  z_lo  z_hi  N_good  amp_real  mean_sh  std_sh  p_value  l_real  b_real\n")
        for r in results:
            f.write(
                f"{r['tag']}  {r['z_lo']:.2f}  {r['z_hi']:.2f}  {r['N_good']:6d} "
                f"{r['amp_real']:.3e} {r['mean_sh']:.3e} {r['std_sh']:.3e} "
                f"{r['p_value']:.4f}  {r['l_real']:.2f}  {r['b_real']:.2f}\n"
            )

    print(f"\n[shuffle] wrote summary -> {out_txt}")


if __name__ == "__main__":
    main()
