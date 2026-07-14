#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SM_SkyMap_FR_KdS.py — Build an anisotropy sky map and overlay FR/KdS vectors

Generates a healpy Mollweide image of the sky, mapping a BAO-derived ratio
    R_BAO = (D_H / r_d) / (D_M / r_d) ≡ D_H / D_M
from a unified PLamb catalog, with optional overlays:
  • Dipole axis inferred from the R_BAO map
  • FR and KdS model axis vectors (provided on CLI)
  • IMF2/IMF3 HHT points (colored by energy) and optional locking filter

Inputs (unified catalog; CSV/Parquet) — suggested columns:
  ra_deg, dec_deg, z, kind,
  DH_over_rd, DM_over_rd,  # or: DH, DM, rd
  imf2_E, imf3_E, hht_lock # optional (0..1 or True/False)

Output: a saved PNG in plamb_runs/SkyMap/figs (configurable).

Notes:
  • Uses HEALPix NSIDE grid (default 64) and inverse-variance weights if
    columns like sigma_DH_over_rd / sigma_DM_over_rd exist.
  • For sparse sky, Gaussian smoothing is applied (configurable FWHM deg).
  • All overlays are plotted in Galactic (l, b) for cross-comparison with CMB.

Author : PLamb / ChatGPT (GPT-5 Thinking)
Version: 0V1 (2025-10-23)
License: MIT
"""

from __future__ import annotations

import argparse
import json
import math
import os
from dataclasses import dataclass
from typing import Iterable, Optional, Tuple

import numpy as np
import pandas as pd
import healpy as hp
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

# ----------------------------
# Utilities & Data Structures
# ----------------------------

@dataclass
class Columns:
    ra_deg           : str = "ra_deg"
    dec_deg          : str = "dec_deg"
    z                : str = "z"
    DH_over_rd       : str = "DH_over_rd"
    DM_over_rd       : str = "DM_over_rd"
    DH               : str = "DH"
    DM               : str = "DM"
    rd               : str = "rd"
    imf2_E           : str = "imf2_E"
    imf3_E           : str = "imf3_E"
    hht_lock         : str = "hht_lock"
    sigma_DH_over_rd : str = "sigma_DH_over_rd"
    sigma_DM_over_rd : str = "sigma_DM_over_rd"


@dataclass
class Args:
    catalog                : str
    out                    : str
    nside                  : int
    coord_in               : str
    zmin                   : float
    zmax                   : float
    gal_mask_abs_b_deg     : float
    fwhm_smooth_deg        : float
    rbao_mode              : str
    rd_const               : Optional[float]
    highlight_quantile     : float
    fr_vec_lb_deg          : Optional[Tuple[float, float]]
    kds_vec_lb_deg         : Optional[Tuple[float, float]]
    lock_min               : Optional[float]
    dipole_search_nside    : int
    cmap                   : str
    vmin                   : Optional[float]
    vmax                   : Optional[float]


# Align assignment columns nicely within logical blocks

def parse_args() -> Args:
    p                            = argparse.ArgumentParser(
        description              = "SkyMap: R_BAO anisotropy with FR/KdS vector overlays",
        formatter_class          = argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument("catalog",                      type=str, help="Unified catalog CSV/Parquet path")
    p.add_argument("--out",                        type=str, default="plamb_runs/SkyMap/figs/SM_SkyMap_FR_KdS.png")
    p.add_argument("--nside",                      type=int, default=64, help="HEALPix NSIDE (power of 2)")
    p.add_argument("--coord-in",                   type=str, default="icrs", choices=["icrs", "galactic"], help="Input coord frame")
    p.add_argument("--zmin",                       type=float, default=0.1)
    p.add_argument("--zmax",                       type=float, default=1.5)
    p.add_argument("--gal-mask-abs-b-deg",         type=float, default=20.0, help="Mask |b|<value (Galactic)")
    p.add_argument("--fwhm-smooth-deg",            type=float, default=8.0, help="Gaussian smoothing FWHM (deg)")

    p.add_argument("--rbao-mode",                  type=str, default="auto", choices=["auto", "DH_over_rd/DM_over_rd", "DH/DM"],
                    help="How to compute R_BAO")
    p.add_argument("--rd-const",                   type=float, default=None, help="Constant r_d if using DH/DM and per-row rd is absent")

    p.add_argument("--highlight-quantile",         type=float, default=0.9, help="Quantile for IMF highlighting")
    p.add_argument("--fr-vec-lb-deg",              type=str, default=None, help="FR vector as 'l_deg,b_deg'")
    p.add_argument("--kds-vec-lb-deg",             type=str, default=None, help="KdS vector as 'l_deg,b_deg'")
    p.add_argument("--lock-min",                   type=float, default=None, help="Minimum HHT locking value to include SN points")

    p.add_argument("--dipole-search-nside",        type=int, default=32, help="NSIDE for coarse dipole search")
    p.add_argument("--cmap",                        type=str, default="coolwarm")
    p.add_argument("--vmin",                        type=float, default=None)
    p.add_argument("--vmax",                        type=float, default=None)

    a                            = p.parse_args()

    fr_vec_lb_deg                = parse_lb(a.fr_vec_lb_deg)
    kds_vec_lb_deg               = parse_lb(a.kds_vec_lb_deg)

    return Args(
        catalog                  = a.catalog,
        out                      = a.out,
        nside                    = a.nside,
        coord_in                 = a.coord_in.lower(),
        zmin                     = a.zmin,
        zmax                     = a.zmax,
        gal_mask_abs_b_deg       = a.gal_mask_abs_b_deg,
        fwhm_smooth_deg          = a.fwhm_smooth_deg,
        rbao_mode                = a.rbao_mode,
        rd_const                 = a.rd_const,
        highlight_quantile       = a.highlight_quantile,
        fr_vec_lb_deg            = fr_vec_lb_deg,
        kds_vec_lb_deg           = kds_vec_lb_deg,
        lock_min                 = a.lock_min,
        dipole_search_nside      = a.dipole_search_nside,
        cmap                     = a.cmap,
        vmin                     = a.vmin,
        vmax                     = a.vmax,
    )


def parse_lb(s: Optional[str]) -> Optional[Tuple[float, float]]:
    if s is None:
        return None
    try:
        l_str, b_str            = s.split(",")
        l_deg                   = float(l_str.strip())
        b_deg                   = float(b_str.strip())
        return (l_deg, b_deg)
    except Exception:
        raise ValueError("Vector must be 'l_deg,b_deg', e.g. '264.0,48.0'")


# ----------------------------
# Coordinate transforms
# ----------------------------

def radec_to_gal_l_b(ra_deg: np.ndarray, dec_deg: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    sc                         = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame="icrs")
    gal                        = sc.galactic
    return gal.l.deg, gal.b.deg


def ensure_galactic(df: pd.DataFrame, coord_in: str, col: Columns) -> pd.DataFrame:
    if coord_in.lower() == "galactic":
        df                      = df.copy()
        if "l_deg" not in df.columns or "b_deg" not in df.columns:
            raise ValueError("When --coord-in=galactic, catalog must have l_deg and b_deg columns.")
        return df

    l_deg, b_deg               = radec_to_gal_l_b(df[col.ra_deg].to_numpy(), df[col.dec_deg].to_numpy())
    out                        = df.copy()
    out["l_deg"]               = l_deg
    out["b_deg"]               = b_deg
    return out


# ----------------------------
# R_BAO construction
# ----------------------------

def compute_rbao(df: pd.DataFrame, col: Columns, mode: str, rd_const: Optional[float]) -> np.ndarray:
    mode                       = mode.lower()

    if mode == "auto":
        has_over               = col.DH_over_rd in df.columns and col.DM_over_rd in df.columns
        has_plain              = col.DH in df.columns and col.DM in df.columns
        if has_over:
            mode               = "dh_over_rd/dm_over_rd"
        elif has_plain:
            mode               = "dh/dm"
        else:
            raise ValueError("Cannot auto-compute R_BAO: missing DH/DM columns.")

    if mode == "dh_over_rd/dm_over_rd":
        num                    = df[col.DH_over_rd].astype(float).to_numpy()
        den                    = df[col.DM_over_rd].astype(float).to_numpy()
        rbao                   = np.divide(num, den, out=np.full_like(num, np.nan), where=den!=0)
        return rbao

    if mode == "dh/dm":
        dh                     = df[col.DH].astype(float).to_numpy()
        dm                     = df[col.DM].astype(float).to_numpy()
        if col.rd in df.columns:
            rd                 = df[col.rd].astype(float).to_numpy()
            with np.errstate(divide='ignore', invalid='ignore'):
                num            = np.divide(dh, rd)
                den            = np.divide(dm, rd)
                rbao           = np.divide(num, den, out=np.full_like(num, np.nan), where=den!=0)
        else:
            if rd_const is None:
                raise ValueError("DH/DM selected but neither per-row rd nor --rd-const provided.")
            # r_d cancels identically; constant not required numerically
            rbao               = np.divide(dh, dm, out=np.full_like(dh, np.nan), where=dm!=0)
        return rbao

    raise ValueError(f"Unknown rbao mode: {mode}")


# ----------------------------
# HEALPix binning & smoothing
# ----------------------------

def make_map_rbao(
    l_deg    : np.ndarray,
    b_deg    : np.ndarray,
    z        : np.ndarray,
    rbao     : np.ndarray,
    nside    : int,
    zmin     : float,
    zmax     : float,
    fwhm_deg : float,
    weights  : Optional[np.ndarray] = None,
) -> np.ndarray:

    sel                        = np.isfinite(rbao) & np.isfinite(l_deg) & np.isfinite(b_deg) & (z >= zmin) & (z < zmax)
    l_rad                      = np.deg2rad(l_deg[sel])
    b_rad                      = np.deg2rad(b_deg[sel])
    theta                      = np.deg2rad(90.0 - np.rad2deg(b_rad))  # colatitude
    phi                        = l_rad

    npix                       = hp.nside2npix(nside)
    m_sum                      = np.full(npix, 0.0)
    w_sum                      = np.full(npix, 0.0)

    pix                        = hp.ang2pix(nside, theta, phi, nest=False)

    rvals                      = rbao[sel].astype(float)
    wvals                      = np.ones_like(rvals) if weights is None else weights[sel].astype(float)

    np.add.at(m_sum, pix, rvals * wvals)
    np.add.at(w_sum, pix, wvals)

    with np.errstate(divide='ignore', invalid='ignore'):
        m_map                  = np.divide(m_sum, w_sum, out=np.full(npix, hp.UNSEEN), where=w_sum>0)

    # Smooth masked map
    m_masked                   = hp.ma(m_map)
    m_masked.mask              = m_map == hp.UNSEEN

    fwhm_rad                   = np.deg2rad(fwhm_deg)
    m_smooth                   = hp.sphtfunc.smoothing(m_masked, fwhm=fwhm_rad, verbose=False)

    # Restore UNSEEN mask after smoothing leakage
    m_out                      = np.array(m_smooth.filled(hp.UNSEEN))
    return m_out


# ----------------------------
# Dipole estimation & vector tools
# ----------------------------

def estimate_dipole_direction(map_rbao: np.ndarray, search_nside: int) -> Tuple[float, float]:
    """Return (l_deg, b_deg) of maximal direction as a proxy for dipole axis."""
    # Coarse up/down-grade to avoid pixel-scale maxima
    if hp.get_nside(map_rbao) != search_nside:
        tmp                    = hp.ud_grade(map_rbao, search_nside, order_in='RING', order_out='RING', power=-2)
    else:
        tmp                    = map_rbao.copy()

    valid                     = tmp != hp.UNSEEN
    if not np.any(valid):
        return (np.nan, np.nan)

    idx_max                   = np.nanargmax(np.where(valid, tmp, np.nan))
    theta, phi                = hp.pix2ang(search_nside, idx_max, nest=False)
    l_deg                     = np.rad2deg(phi)
    b_deg                     = 90.0 - np.rad2deg(theta)
    return (l_deg % 360.0, b_deg)


def lb_to_unit(l_deg: float, b_deg: float) -> np.ndarray:
    l_rad                     = math.radians(l_deg)
    b_rad                     = math.radians(b_deg)
    x                         = math.cos(b_rad) * math.cos(l_rad)
    y                         = math.cos(b_rad) * math.sin(l_rad)
    z                         = math.sin(b_rad)
    return np.array([x, y, z], dtype=float)


def unit_to_lb(v: np.ndarray) -> Tuple[float, float]:
    x, y, z                   = v
    b_rad                     = math.asin(max(-1.0, min(1.0, z)))
    l_rad                     = math.atan2(y, x)
    l_deg                     = math.degrees(l_rad) % 360.0
    b_deg                     = math.degrees(b_rad)
    return (l_deg, b_deg)


def great_circle_between_antipodes(n_hat: np.ndarray, npts: int = 128) -> Tuple[np.ndarray, np.ndarray]:
    """Generate a great-circle path from +n_hat to -n_hat, returned as (l_deg, b_deg)."""
    n_hat                     = n_hat / np.linalg.norm(n_hat)
    # Choose a stable perpendicular direction
    z_axis                    = np.array([0.0, 0.0, 1.0])
    x_axis                    = np.array([1.0, 0.0, 0.0])

    u                         = np.cross(n_hat, z_axis)
    if np.linalg.norm(u) < 1e-8:
        u                     = np.cross(n_hat, x_axis)
    u                         = u / np.linalg.norm(u)

    psi                       = np.linspace(0.0, math.pi, npts)
    pts                       = np.cos(psi)[:, None] * n_hat[None, :] + np.sin(psi)[:, None] * u[None, :]

    l_list, b_list            = [], []
    for p in pts:
        l, b                  = unit_to_lb(p / np.linalg.norm(p))
        l_list.append(l)
        b_list.append(b)
    return np.array(l_list), np.array(b_list)


# ----------------------------
# Plotting
# ----------------------------

def plot_skymap(
    map_rbao       : np.ndarray,
    args           : Args,
    df_gal         : pd.DataFrame,
    col            : Columns,
    dipole_lb      : Optional[Tuple[float, float]] = None,
) -> plt.Figure:

    # Apply Galactic latitude mask to display only (optional)
    if args.gal_mask_abs_b_deg > 0:
        l_mask, b_mask         = hp.pix2ang(hp.get_nside(map_rbao), np.arange(hp.nside2npix(hp.get_nside(map_rbao))), nest=False)
        b_mask_deg             = 90.0 - np.rad2deg(l_mask)
        show_map               = map_rbao.copy()
        pix_to_hide            = np.where(np.abs(b_mask_deg) < args.gal_mask_abs_b_deg)[0]
        show_map[pix_to_hide]  = hp.UNSEEN
    else:
        show_map               = map_rbao

    hp.mollview(
        show_map,
        title       = f"R_BAO (D_H/D_M), z∈[{args.zmin:.2f},{args.zmax:.2f}], NSIDE={args.nside}",
        unit        = "R_BAO",
        norm        = "hist",
        min         = args.vmin,
        max         = args.vmax,
        cmap        = args.cmap,
        cbar        = True,
        coord       = ['G'],
        notext      = False,
        badcolor    = "#202020",
    )

    # Overlay IMF2/IMF3 points (top quantile)
    for name, en_col, marker in [("IMF2", col.imf2_E, "^"), ("IMF3", col.imf3_E, "s")]:
        if en_col in df_gal.columns:
            df_plot            = df_gal.copy()
            if args.lock_min is not None and col.hht_lock in df_plot.columns:
                df_plot        = df_plot[df_plot[col.hht_lock].astype(float) >= float(args.lock_min)]
            ener               = pd.to_numeric(df_plot[en_col], errors='coerce')
            qcut               = ener.quantile(args.highlight_quantile) if ener.notna().any() else np.nan
            sel_pts            = ener >= qcut if np.isfinite(qcut) else np.array([False]*len(df_plot))
            l_pts              = df_plot.loc[sel_pts, "l_deg"].to_numpy()
            b_pts              = df_plot.loc[sel_pts, "b_deg"].to_numpy()
            if len(l_pts) > 0:
                hp.projscatter(l_pts, b_pts, lonlat=True, s=18, marker=marker, alpha=0.8, label=f"{name} top {args.highlight_quantile:.2f}")

    # Overlay dipole axis from R_BAO map
    if dipole_lb is not None and all(np.isfinite(dipole_lb)):
        n_hat                 = lb_to_unit(*dipole_lb)
        l_line, b_line        = great_circle_between_antipodes(n_hat)
        hp.projplot(l_line, b_line, lonlat=True, ls='-', lw=2.0, alpha=0.9, label="R_BAO dipole")
        # Mark the + pole
        hp.projscatter([dipole_lb[0]], [dipole_lb[1]], lonlat=True, marker="o", s=40, alpha=0.9)

    # Overlay FR and KdS axes if provided
    if args.fr_vec_lb_deg is not None:
        n_hat                 = lb_to_unit(*args.fr_vec_lb_deg)
        l_line, b_line        = great_circle_between_antipodes(n_hat)
        hp.projplot(l_line, b_line, lonlat=True, ls='--', lw=1.8, alpha=0.95, label="FR axis")
        hp.projscatter([args.fr_vec_lb_deg[0]], [args.fr_vec_lb_deg[1]], lonlat=True, marker="*", s=60, alpha=0.95)

    if args.kds_vec_lb_deg is not None:
        n_hat                 = lb_to_unit(*args.kds_vec_lb_deg)
        l_line, b_line        = great_circle_between_antipodes(n_hat)
        hp.projplot(l_line, b_line, lonlat=True, ls='-.', lw=1.8, alpha=0.95, label="KdS axis")
        hp.projscatter([args.kds_vec_lb_deg[0]], [args.kds_vec_lb_deg[1]], lonlat=True, marker="D", s=52, alpha=0.95)

    # Legend (matplotlib)
    plt.legend(loc='lower left', fontsize=8, framealpha=0.7)

    fig                      = plt.gcf()
    return fig


# ----------------------------
# Main
# ----------------------------

def main() -> None:
    args                     = parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # Load catalog
    ext                      = os.path.splitext(args.catalog)[1].lower()
    if ext in [".csv", ".txt", ".tsv"]:
        df                   = pd.read_csv(args.catalog)
    elif ext in [".parquet", ".pq"]:
        df                   = pd.read_parquet(args.catalog)
    else:
        raise ValueError("Unsupported catalog format. Use CSV or Parquet.")

    col                      = Columns()

    # Compute R_BAO and ensure Galactic coordinates
    rbao                     = compute_rbao(df, col, args.rbao_mode, args.rd_const)
    df_gal                   = ensure_galactic(df, args.coord_in, col)

    # Basic weights if available (inverse variance of log-ratio approx.)
    if col.sigma_DH_over_rd in df.columns and col.sigma_DM_over_rd in df.columns:
        # crude propagation for ratio a/b → σ ≈ R * sqrt( (σ_a/a)^2 + (σ_b/b)^2 )
        a                    = pd.to_numeric(df[col.DH_over_rd], errors='coerce').to_numpy()
        b                    = pd.to_numeric(df[col.DM_over_rd], errors='coerce').to_numpy()
        sa                   = pd.to_numeric(df[col.sigma_DH_over_rd], errors='coerce').to_numpy()
        sb                   = pd.to_numeric(df[col.sigma_DM_over_rd], errors='coerce').to_numpy()
        with np.errstate(divide='ignore', invalid='ignore'):
            rel2             = np.square(sa / a) + np.square(sb / b)
            sigma_r          = np.abs(rbao) * np.sqrt(rel2)
            weights          = 1.0 / np.square(sigma_r)
            weights[~np.isfinite(weights)] = 0.0
    else:
        weights              = None

    # Build map
    map_rbao                 = make_map_rbao(
        l_deg    = df_gal["l_deg"].to_numpy(),
        b_deg    = df_gal["b_deg"].to_numpy(),
        z        = pd.to_numeric(df_gal[col.z], errors='coerce').to_numpy(),
        rbao     = rbao,
        nside    = args.nside,
        zmin     = args.zmin,
        zmax     = args.zmax,
        fwhm_deg = args.fwhm_smooth_deg,
        weights  = weights,
    )

    # Estimate dipole direction
    dipole_lb                = estimate_dipole_direction(map_rbao, args.dipole_search_nside)

    # Plot and save
    fig                      = plot_skymap(map_rbao, args, df_gal, col, dipole_lb=dipole_lb)
    fig.savefig(args.out, dpi=200, bbox_inches='tight')

    print(f"[ok] saved skymap → {args.out}")


if __name__ == "__main__":
    main()
