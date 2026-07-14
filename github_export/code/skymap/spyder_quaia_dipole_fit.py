#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_dipole_fit.py  (0V2)

Fit a simple dipole model to the Quaia <z> maps built by
quaia_build_zmean_maps.py, for full sample and each redshift slice.

Model:
    f(n) = a0 + b · n

where:
    f(n) is the per-pixel quantity (here we use <z>_pix),
    n is the unit vector on the sphere.

We solve for a0 and the 3-component vector b via weighted least squares,
using only pixels with N >= N_MIN_PER_PIX.

Outputs:
  - Text summary of dipole amplitude and (l,b) direction per slice.
  
  
Updated Run Example:
    %run /Users/boyde/.spyder-py3/quaia_dipole_fit.py

"""

from   pathlib import Path

import numpy   as np
import healpy  as hp

from quaia_config import (
                          OUT_DIR,
                          NSIDE,
                          N_MIN_PER_PIX,
                          Z_BINS,
                        )


# CMB dipole direction (Planck) in Galactic coords
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
    """
    Angular separation between two Galactic directions (l1,b1) and (l2,b2), in degrees.
    """
    l1 = np.deg2rad(l1_deg); b1 = np.deg2rad(b1_deg)
    l2 = np.deg2rad(l2_deg); b2 = np.deg2rad(b2_deg)

    x1 = np.cos(b1) * np.cos(l1)
    y1 = np.cos(b1) * np.sin(l1)
    z1 = np.sin(b1)

    x2 = np.cos(b2) * np.cos(l2)
    y2 = np.cos(b2) * np.sin(l2)
    z2 = np.sin(b2)

    dot = x1 * x2 + y1 * y2 + z1 * z2
    dot = np.clip(dot, -1.0, 1.0)

    return np.rad2deg(np.arccos(dot))


def load_map_bundle(tag):
    path = OUT_DIR / f"quaia_zmean_{tag}_nside{NSIDE}.npz"
    # allow_pickle=True needed because z_min/z_max ended up as object arrays
    data = np.load(path, allow_pickle=True)
    return (
            data["N"],
            data["zmean"],
            float(data["NSIDE"]),
            data.get("z_min", None),
            data.get("z_max", None),
        )


def unit_vectors_for_nside(nside):
    """
    Return unit vectors n_hat[pix] for all pixels in RING ordering.
    """
    npix = hp.nside2npix(nside)
    pix  = np.arange(npix)
    theta, phi = hp.pix2ang(nside, pix, nest=False)
    # healpy uses colatitude theta, longitude phi
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    n_hat = np.vstack([x, y, z]).T  # shape (npix, 3)
    return n_hat


def fit_dipole_linear(f, n_hat, w=None):
    """
    Fit f(n) = a0 + b · n via linear least squares.

    Parameters
    ----------
    f : (N,) array
        Scalar field values per pixel (e.g. <z>).
    n_hat : (N,3) array
        Unit vectors for each pixel.
    w : (N,) array or None
        Optional weights; if None, all equal.

    Returns
    -------
    a0 : float
    b  : (3,) array
    cov : (4,4) array
        Covariance matrix of (a0, b_x, b_y, b_z).
    """
    f = np.asarray(f, dtype=float)
    n_hat = np.asarray(n_hat, dtype=float)
    N = f.size

    # Design matrix: [1, n_x, n_y, n_z]
    X = np.empty((N, 4), dtype=float)
    X[:, 0] = 1.0
    X[:, 1:] = n_hat

    if w is None:
        W = np.eye(N, dtype=float)
    else:
        w = np.asarray(w, dtype=float)
        W = np.diag(w)

    # Solve (X^T W X) beta = X^T W f
    XT_W = X.T @ W
    M    = XT_W @ X
    y    = XT_W @ f
    beta = np.linalg.solve(M, y)
    cov  = np.linalg.inv(M)

    a0   = beta[0]
    b    = beta[1:]
    return a0, b, cov


def dipole_summary_from_b(b):
    """
    From dipole vector b, return amplitude and direction (l,b) in degrees.
    """
    bx, by, bz = b
    amp        = np.linalg.norm(b)

    if amp == 0.0:
        return 0.0, np.nan, np.nan

    # Direction of dipole = direction of +b
    theta = np.arccos(np.clip(bz / amp, -1.0, 1.0))  # colatitude
    phi   = np.arctan2(by, bx)                       # longitude

    if phi < 0:
        phi += 2.0 * np.pi

    # (l,b): l = phi, b = 90° - theta
    l_deg = np.degrees(phi)
    b_deg = 90.0 - np.degrees(theta)
    return amp, l_deg, b_deg

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

def project_dipole_on_cmb_with_cov(b_vec, cov):
    """
    Given dipole vector b_vec (3,) and full covariance matrix cov (4x4) for
    (a0, b_x, b_y, b_z), return:

      amp       = |b_vec|
      amp_par   = component of b_vec along CMB direction
      sigma_par = 1σ uncertainty on amp_par
      frac_par  = amp_par / amp
    """
    b_vec = np.asarray(b_vec, dtype=float)
    amp   = np.linalg.norm(b_vec)

    if amp == 0.0:
        return 0.0, 0.0, 0.0, 0.0

    # 3×3 covariance submatrix for (b_x, b_y, b_z)
    cov_b = np.asarray(cov[1:, 1:], dtype=float)

    # CMB unit vector
    n_cmb = unit_vec_gal(L_CMB_DEG, B_CMB_DEG)

    # Component along CMB direction
    amp_par = float(np.dot(b_vec, n_cmb))

    # Variance of a linear combination: n^T Cov n
    var_par   = float(n_cmb @ cov_b @ n_cmb)
    var_par   = max(var_par, 0.0)   # guard against tiny negative due to rounding
    sigma_par = np.sqrt(var_par)

    frac_par = amp_par / amp if amp > 0.0 else np.nan
    return amp, amp_par, sigma_par, frac_par


def analyze_slice(tag, human_label):
    """
    Load maps for a given tag and fit dipole.
    """
    N_map, zmean_map, NSIDE_loaded, z_min_arr, z_max_arr = load_map_bundle(tag)
    assert int(NSIDE_loaded) == NSIDE, f"NSIDE mismatch in {tag}: {NSIDE_loaded} vs {NSIDE}"

    # Only pixels with enough objects and non-NaN zmean
    good   = (N_map >= N_MIN_PER_PIX) & np.isfinite(zmean_map)
    N_good = int(good.sum())

    # If there are no (or almost no) good pixels, skip this slice
    if N_good < 3:
        print(f"[warn] slice '{human_label}' has N_good={N_good} < 3; skipping dipole fit.")
        return dict(
                    tag      = tag,
                    label    = human_label,
                    N_good   = N_good,
                    amp      = np.nan,
                    l_deg    = np.nan,
                    b_deg    = np.nan,
                    a0       = np.nan,
                    amp_par  = np.nan,
                    frac_par = np.nan,
                    sep_cmb  = np.nan,
                )

    if N_good < 50:
        print(f"[warn] slice '{human_label}' has only {N_good} good pixels; results will be noisy.")

    # Prepare design matrix
    f = zmean_map[good]
    f = f - f.mean()

    # HEALPix pixel indices for good pixels – force to int64 for healpy
    import healpy as hp
    pix        = np.where(good)[0]
    pix        = np.asarray(pix, dtype=np.int64)

    # healpy wants an integer NSIDE as well
    nside_int  = int(NSIDE_loaded)
    theta, phi = hp.pix2ang(nside_int, pix)

    # Convert to Cartesian unit vectors n_hat = (x,y,z)
    x          = np.sin(theta) * np.cos(phi)
    y          = np.sin(theta) * np.sin(phi)
    z          = np.cos(theta)
    n_hat      = np.vstack([x, y, z]).T    # shape (N_good, 3)

    # Weight by N per pixel – more galaxies -> more weight
    w          = N_map[good].astype(float)

    N_good     = np.count_nonzero(good)
    if N_good < 4:
        print(f"[warn] slice '{human_label}' has only {N_good} good pixels; skipping dipole fit.")
        return dict(
                    tag      = tag,
                    label    = human_label,
                    N_good   = N_good,
                    amp      = np.nan,
                    l_deg    = np.nan,
                    b_deg    = np.nan,
                    a0       = np.nan,
                    amp_par  = np.nan,
                    frac_par = np.nan,
                    sep_cmb  = np.nan,
                  )

    # --- Fit dipole and project on CMB direction ---
    a0, b_vec, cov = fit_dipole_linear(f, n_hat, w=w)

    # Basic dipole geometry
    amp, l_deg, b_deg = dipole_summary_from_b(b_vec)
    
    # CMB projection + uncertainty
    amp, amp_par, sigma_par, frac_par = project_dipole_on_cmb_with_cov(b_vec, cov)
    
    sep_cmb = ang_sep_deg(l_deg, b_deg, L_CMB_DEG, B_CMB_DEG) 


    # --- Print full summary (now variables are defined) ---
    print(f"\n[slice] {human_label}")
    print(f"  N_good_pix = {N_good}")
    print(f"  dipole amplitude |b|             = {amp:.4e}")
    print(f"  dipole direction (l,b)           = ({l_deg:.1f}°, {b_deg:.1f}°)")
    print(f"  a0 (offset, after mean-subtract) = {a0:.4e}")
    print(f"  component along CMB dir: amp_par = {amp_par:.4e} ± {sigma_par:.4e}  "
          f"(frac_par = {frac_par:+.3f})")
    print(f"  separation from CMB dipole       ≈ {sep_cmb:.1f}°")


    return dict(
                tag      = tag,
                label    = human_label,
                N_good   = N_good,
                amp      = amp,
                l_deg    = l_deg,
                b_deg    = b_deg,
                a0       = a0,
                amp_par  = amp_par,
                sigma_par= sigma_par,
                frac_par = frac_par,
                sep_cmb  = sep_cmb,
               )

def main():
    print("[quaia] dipole fits from <z> maps")

    results = []

    # 1) Full sample
    results.append(analyze_slice("full", "full sample"))

    # 2) Redshift slices
    for (z_lo, z_hi) in Z_BINS:
        tag = f"z{z_lo:.2f}_{z_hi:.2f}".replace(".", "p")
        label = f"{z_lo:.2f} ≤ z < {z_hi:.2f}"
        results.append(analyze_slice(tag, label))

    # Save summary table
    out_txt = OUT_DIR / "quaia_dipole_summary.txt"
    with open(out_txt, "w") as f:
        f.write("# tag  label  N_good  amp  amp_par  sigma_par  frac_par  "
                "l_deg  b_deg  sep_cmb  a0\n")
        for r in results:
            f.write(
                f"{r['tag']:30s}  {r['label']:20s}  {r['N_good']:6d}  "
                f"{r['amp']:10.4e}  {r['amp_par']:10.4e}  {r['sigma_par']:10.4e}  "
                f"{r['frac_par']:8.3f}  {r['l_deg']:8.3f}  {r['b_deg']:8.3f}  "
                f"{r['sep_cmb']:8.3f}  {r['a0']:10.4e}\n"
            )

    print(f"\n[save] wrote dipole summary -> {out_txt}")
    print("[quaia] done.")


if __name__ == "__main__":
    main()
