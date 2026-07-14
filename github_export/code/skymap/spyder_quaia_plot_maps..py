#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_plot_maps.py  (0V2)

Generate HEALPix skymaps for Quaia:
- <z> maps per slice
- optional residual <z> maps after removing z–N trend
- overlay measured dipole direction and CMB dipole

Uses the existing .npz files built by:
  - quaia_build_zmean_maps.py
  - quaia_build_zmean_maps_bcut.py

Output PNG files go to OUT_DIR.
"""

import numpy as np
import healpy as hp
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for Spyder runs
import matplotlib.pyplot as plt
from pathlib import Path

from quaia_config import OUT_DIR, NSIDE, Z_BINS
from quaia_dipole_fit import fit_dipole_linear, dipole_summary_from_b

# CMB dipole direction (Galactic coords, degrees)
L_CMB_DEG = 264.0
B_CMB_DEG = 48.0

# Minimum N per pixel when plotting / fitting
N_MIN = 3


def load_zmean_bundle(tag: str, bcut: float | None = None):
    """
    Load N and <z> map for a given redshift tag and optional |b| cut.

    tag examples:
      - "full"
      - "z0p10_0p50", "z0p50_1p00", ...

    If bcut is None: use the plain maps:
        quaia_zmean_{tag}_nside{NSIDE}.npz

    If bcut is 10, 20, 30: use the bcut maps:
        quaia_zmean_{tag}_bcutXX_nside{NSIDE}.npz
    """
    if bcut is None:
        fname = f"quaia_zmean_{tag}_nside{NSIDE}.npz"
    else:
        fname = f"quaia_zmean_{tag}_bcut{int(bcut):02d}_nside{NSIDE}.npz"

    path = OUT_DIR / fname
    if not path.exists():
        raise FileNotFoundError(f"Map file not found: {path}")

    data = np.load(path)
    N_map = data["N"]
    zmean = data["zmean"]
    return N_map, zmean


def build_residual_map(N_map: np.ndarray, zmean: np.ndarray):
    """
    Regress zmean on N over good pixels and return the residual map:
        z_resid = zmean - (a + b*N)
    Also return (a, b) for reporting.
    """
    good = (N_map >= N_MIN) & np.isfinite(zmean)
    N_good = N_map[good].astype(float)
    z_good = zmean[good].astype(float)

    A = np.vstack([np.ones_like(N_good), N_good]).T
    beta, *_ = np.linalg.lstsq(A, z_good, rcond=None)
    a, b = beta

    z_model = np.full_like(zmean, np.nan, dtype=float)
    z_model[good] = a + b * N_good
    z_resid = zmean - z_model

    return z_resid, (a, b)


def compute_dipole_from_map(z_map: np.ndarray, mask: np.ndarray):
    """
    Given a z_map and a boolean mask of good pixels, fit a dipole
    using fit_dipole_linear and return (amp, l_deg, b_deg).
    """
    pix = np.where(mask)[0]
    if pix.size < 3:
        return np.nan, np.nan, np.nan

    theta, phi = hp.pix2ang(NSIDE, pix)  # radians

    # Unit vectors for each pixel (shape (3, N))
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    n_hat = np.vstack([x, y, z])  # (3, N) -> fit_dipole_linear expects this

    f = z_map[mask].astype(float)
    f = f - f.mean()  # mean-subtract before dipole fit

    a0, b_vec, _ = fit_dipole_linear(f, n_hat, w=None)
    amp, l_deg, b_deg = dipole_summary_from_b(b_vec)
    return amp, l_deg, b_deg


def mollview_with_dipoles(m, title, unit, out_path,
                          amp=None, l_deg=None, b_deg=None,
                          vmin=None, vmax=None):
    """
    Helper to make a Mollweide plot with optional dipole + CMB markers.
    coord='G' to stay in Galactic coords.
    """
    # Choose color scale if not given
    if vmin is None or vmax is None:
        finite = np.isfinite(m)
        if np.any(finite):
            vmin = np.nanpercentile(m[finite], 5)
            vmax = np.nanpercentile(m[finite], 95)

    hp.mollview(m, coord="G", unit=unit, title=title,
                min=vmin, max=vmax, cmap="viridis")

    # Overlay measured dipole direction (if provided)
    if amp is not None and np.isfinite(amp) and np.isfinite(l_deg) and np.isfinite(b_deg):
        theta_d = np.deg2rad(90.0 - b_deg)
        phi_d = np.deg2rad(l_deg)
        hp.projplot(theta_d, phi_d, marker="o", markersize=7)

    # Overlay CMB dipole direction as a star
    theta_c = np.deg2rad(90.0 - B_CMB_DEG)
    phi_c = np.deg2rad(L_CMB_DEG)
    hp.projplot(theta_c, phi_c, marker="*", markersize=9)

    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"[plot] wrote {out_path}")


def plot_for_tag(tag: str,
                 human_label: str,
                 bcut: float | None = 20.0,
                 make_residual: bool = True):
    """
    For a given tag and |b| cut:
      - plot the raw <z> map
      - optionally plot the residual <z> map
    """
    # ----- raw <z> map -----
    N_map, zmean = load_zmean_bundle(tag, bcut=bcut)
    good = (N_map >= N_MIN) & np.isfinite(zmean)

    m_z = zmean.copy()
    m_z[~good] = np.nan

    amp_z, l_z, b_z = compute_dipole_from_map(m_z, good)

    label_cut = f" (|b|>{int(bcut)}°)" if bcut is not None else ""
    title_z = f"{human_label}{label_cut}\n<z> dipole: amp={amp_z:.2e}, (l,b)=({l_z:.1f}°, {b_z:.1f}°)"

    out_tag = f"{tag}" if bcut is None else f"{tag}_bcut{int(bcut):02d}"
    out_z = OUT_DIR / f"quaia_map_{out_tag}_zmean.png"

    mollview_with_dipoles(m_z, title_z, unit="<z>", out_path=out_z,
                          amp=amp_z, l_deg=l_z, b_deg=b_z)

    # ----- residual <z> map (optional) -----
    if make_residual:
        z_resid, (a, b) = build_residual_map(N_map, zmean)
        m_r = z_resid.copy()
        m_r[~good] = np.nan

        amp_r, l_r, b_r = compute_dipole_from_map(m_r, good)
        title_r = (
            f"{human_label}{label_cut} residual\n"
            f"z_resid = z - (a + bN), a={a:.3f}, b={b:.2e}\n"
            f"dipole: amp={amp_r:.2e}, (l,b)=({l_r:.1f}°, {b_r:.1f}°)"
        )

        out_r = OUT_DIR / f"quaia_map_{out_tag}_zresid.png"
        mollview_with_dipoles(m_r, title_r, unit="z_resid",
                              out_path=out_r,
                              amp=amp_r, l_deg=l_r, b_deg=b_r)


def main():
    # Tags and labels matching your current pipeline
    tags = [("full", "full sample")]
    for (z_lo, z_hi) in Z_BINS:
        tag = f"z{z_lo:.2f}_{z_hi:.2f}".replace(".", "p")
        label = f"{z_lo:.2f} ≤ z < {z_hi:.2f}"
        tags.append((tag, label))

    # Visualise |b|>20° maps (your main FR-friendly mask)
    bcut = 20.0

    for tag, label in tags:
        plot_for_tag(tag, label, bcut=bcut, make_residual=True)


if __name__ == "__main__":
    main()
