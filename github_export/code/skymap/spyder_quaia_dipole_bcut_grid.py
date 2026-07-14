#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 23 20:41:20 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_dipole_bcut_grid.py  (0V2)

For a set of Galactic latitude cuts |b| > b_cut (10,20,30 deg), and for each
redshift slice (plus full sample), fit the <z> dipole and measure its
alignment with the CMB dipole.

%run /Users/boyde/.spyder-py3/quaia_dipole_bcut_grid.py

Outputs:
  - Text summary: quaia_dipole_bcut_grid.txt in OUT_DIR
"""

import numpy            as     np
from   pathlib          import Path
import healpy           as     hp

from   quaia_config     import OUT_DIR, NSIDE, Z_BINS
from   quaia_dipole_fit import fit_dipole_linear, dipole_summary_from_b

BCUT_LIST      = [10.0, 15.0, 20.0, 25.0, 30.0, 35.0]
N_MIN_PER_PIX  = 3

L_CMB_DEG      = 264.0
B_CMB_DEG      = 48.0


def ang_sep_deg(l1_deg, b1_deg, l2_deg, b2_deg):
    l1 = np.deg2rad(l1_deg); b1 = np.deg2rad(b1_deg)
    l2 = np.deg2rad(l2_deg); b2 = np.deg2rad(b2_deg)

    x1 = np.cos(b1)*np.cos(l1); y1 = np.cos(b1)*np.sin(l1); z1 = np.sin(b1)
    x2 = np.cos(b2)*np.cos(l2); y2 = np.cos(b2)*np.sin(l2); z2 = np.sin(b2)

    dot = x1*x2 + y1*y2 + z1*z2
    dot = np.clip(dot, -1.0, 1.0)
    return np.rad2deg(np.arccos(dot))





def load_map(tag, bcut):
    path = OUT_DIR / f"quaia_zmean_{tag}_bcut{int(bcut)}_nside{NSIDE}.npz"
    data = np.load(path, allow_pickle=False)
    return data["N"], data["zmean"]


def analyze_slice(tag, label, bcut):
    N_map, zmean = load_map(tag, bcut)

    good         = (N_map >= N_MIN_PER_PIX) & np.isfinite(zmean)
    pix          = np.where(good)[0]
    N_good       = N_map[good].astype(float)
    z_good       = zmean[good].astype(float)

    if pix.size < 3:
        print(f"[bcut-grid] tag={tag}, |b|>{bcut:.1f}: too few good pixels ({pix.size})")
        return None

    theta, phi        = hp.pix2ang(NSIDE, pix)
    x                 = np.sin(theta) * np.cos(phi)
    y                 = np.sin(theta) * np.sin(phi)
    z                 = np.cos(theta)
    n_hat             = np.column_stack([x, y, z])

    f                 = z_good - z_good.mean()
    a0, b_vec, _      = fit_dipole_linear(f, n_hat, w=N_good)
    amp, l_deg, b_deg = dipole_summary_from_b(b_vec)

    # Project onto CMB direction
    amp_chk, amp_par, frac_par = project_dipole_on_cmb(b_vec)
    # Optional sanity check: |b_vec| from projection should match amp
    # (you can comment this out later if you like)
    if not np.allclose(amp_chk, amp, rtol=1e-6, atol=1e-9):
        print(f"[warn] amp mismatch: dipole_summary={amp:.3e}, proj_norm={amp_chk:.3e}")

    sep_cmb = ang_sep_deg(l_deg, b_deg, L_CMB_DEG, B_CMB_DEG)

    print(
          f"[bcut-grid] tag={tag}, |b|>{bcut:.1f}: "
          f"N_good={pix.size}, amp={amp:.3e}, "
          f"(l,b)=({l_deg:.1f}°, {b_deg:.1f}°), sepCMB={sep_cmb:.1f}°, "
          f"amp_par={amp_par:.3e}, frac_par={frac_par:+.3f}"
         )

    return dict(
        bcut     = bcut,
        tag      = tag,
        label    = label,
        N_good   = int(pix.size),
        amp      = amp,
        l_deg    = l_deg,
        b_deg    = b_deg,
        sep_cmb  = sep_cmb,
        amp_par  = amp_par,
        frac_par = frac_par,
    )



def unit_vec_gal(l_deg, b_deg):
    """
    Unit vector in Galactic coords (l,b in deg).
    """
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    x = np.cos(b) * np.cos(l)
    y = np.cos(b) * np.sin(l)
    z = np.sin(b)
    return np.array([x, y, z])


def project_dipole_on_cmb(b_vec):
    """
    Given dipole vector b_vec (3,), return:
    - amp = |b_vec|
    - amp_par = component of b_vec along CMB direction (scalar)
    - frac_par = amp_par / amp  (can be +/-)
    """
    b_vec = np.asarray(b_vec, dtype=float)
    amp   = np.linalg.norm(b_vec)
    if amp == 0.0:
        return 0.0, 0.0, 0.0

    n_cmb    = unit_vec_gal(L_CMB_DEG, B_CMB_DEG)
    amp_par  = float(np.dot(b_vec, n_cmb))      # signed component
    frac_par = amp_par / amp
    
    return amp, amp_par, frac_par



def main():
    print("[bcut-grid] Quaia <z> dipoles for |b| cuts:", BCUT_LIST)
    results = []

    # Full sample + slices for each bcut
    for bcut in BCUT_LIST:
        # Full
        res = analyze_slice("full", "full sample", bcut)
        if res is not None:
            results.append(res)

        # Slices
        for (z_lo, z_hi) in Z_BINS:
            tag   = f"z{z_lo:.2f}_{z_hi:.2f}".replace(".", "p")
            label = f"{z_lo:.2f} ≤ z < {z_hi:.2f}"
            res   = analyze_slice(tag, label, bcut)
            if res is not None:
                results.append(res)

    # Save summary
    out_txt = OUT_DIR / "quaia_dipole_bcut_grid.txt"
    with out_txt.open("w") as f:
        f.write("# bcut  tag  N_good  amp  l_deg  b_deg  sep_cmb  amp_par  frac_par\n")
        for r in results:
            f.write(
                    f"{r['bcut']:4.1f}  {r['tag']:12s}  {r['N_good']:6d}  "
                    f"{r['amp']:.3e}  {r['l_deg']:.2f}  {r['b_deg']:.2f}  "
                    f"{r['sep_cmb']:.2f}  {r['amp_par']:.3e}  {r['frac_par']:+.3f}\n"
                   )


    print(f"[bcut-grid] wrote summary -> {out_txt}")


if __name__ == "__main__":
    main()
