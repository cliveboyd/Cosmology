#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SM_SkyMap_FR_KdS.py — Build an anisotropy sky map and overlay FR/KdS vectors

Generates a healpy Mollweide image of the sky, mapping a BAO-derived ratio
    R_BAO = (D_H / r_d) / (D_M / r_d) ≡ D_H / D_M
from a unified catalog, with optional overlays:
  • Dipole axis inferred from the R_BAO map (mask-aware centered fit)
  • FR and KdS model axis vectors (provided on CLI)
  • IMF2/IMF3 HHT points (colored by energy) and optional locking filter

Suggested columns in the catalog (CSV/Parquet):
  ra_deg, dec_deg, z, kind,
  DH_over_rd, DM_over_rd,  # or: DH, DM, rd
  imf2_E, imf3_E, hht_lock # optional (0..1 or True/False)

Example:
  %run '/Users/boyde/.spyder-py3/SM_SkyMap_FR_KdS.py' \
  '/Users/boyde/.spyder-py3/plamb_runs/SkyMap/data/unified_catalog.csv' \
  --out '/Users/boyde/.spyder-py3/plamb_runs/SkyMap/figs/SM_RBAO_z0p1_1p2.png' \
  --nside 64 --zmin 0.1 --zmax 1.2 --fwhm-smooth-deg 8.0 \
  --coord-in icrs --rbao-mode auto --gal-mask-abs-b-deg 20.0 \
  --highlight-quantile 0.90 --lock-min 0.60 \
  --fr-vec-lb-deg "264,48" --kds-vec-lb-deg "300,30" \
  --dipole-search-nside 32 --cmap coolwarm --inset-loc 'lower right'
  

%run '/Users/boyde/.spyder-py3/SM_SkyMap_FR_KdS.py' \
     '/Users/boyde/.spyder-py3/plamb_runs/SkyMap/data/unified_catalog.csv' \
--out '/Users/boyde/.spyder-py3/plamb_runs/SkyMap/figs/SM_RBAO_z0p1_1p2.png' \
--nside 64 \
--zmin 0.1 \
--zmax 1.2 \
--fwhm-smooth-deg 8.0 \
--coord-in icrs \
--rbao-mode auto \
--gal-mask-abs-b-deg 20.0 \
--highlight-quantile 0.90 \
--lock-min 0.60 \
--fr-vec-lb-deg '264,48' \
--kds-vec-lb-deg '300,30' \
--dipole-search-nside 64 \
--cmap coolwarm \
--inset-loc 'lower right' \
--save-stats \
--save-json



Author : PLamb / ChatGPT (GPT-5 Thinking)
Version: 0V3 (2025-10-28) • License: MIT
"""

from __future__ import annotations

import argparse
import math
import sys
import os
from   dataclasses         import dataclass
from   typing              import Optional, Tuple, List
from   astropy.coordinates import SkyCoord
from   math                import pi

import numpy               as np
import pandas              as pd
import healpy              as hp
import matplotlib.pyplot   as plt
import astropy.units       as u


# ----------------------------
# Dataclasses
# ----------------------------

@dataclass
class Columns:
    ra_deg:           str = "ra_deg"
    dec_deg:          str = "dec_deg"
    z:                str = "z"
    DH_over_rd:       str = "DH_over_rd"
    DM_over_rd:       str = "DM_over_rd"
    DH:               str = "DH"
    DM:               str = "DM"
    rd:               str = "rd"
    imf2_E:           str = "imf2_E"
    imf3_E:           str = "imf3_E"
    hht_lock:         str = "hht_lock"
    sigma_DH_over_rd: str = "sigma_DH_over_rd"
    sigma_DM_over_rd: str = "sigma_DM_over_rd"


@dataclass
class Args:
    catalog                 : str
    out                     : str
    nside                   : int
    coord_in                : str
    zmin                    : float
    zmax                    : float
    gal_mask_abs_b_deg      : float
    fwhm_smooth_deg         : float
    rbao_mode               : str
    rd_const                : Optional[float]
    highlight_quantile      : float
    fr_vec_lb_deg           : Optional[Tuple[float, float]]
    kds_vec_lb_deg          : Optional[Tuple[float, float]]
    lock_min                : Optional[float]
    dipole_search_nside     : int
    cmap                    : str
    vmin                    : Optional[float]
    vmax                    : Optional[float]
    title                   : Optional[str]
    inset_loc               : str
    perm_n                  : int
    perm_seed               : int
    save_stats              : bool
    z_bins                  : Optional[List[Tuple[float, float]]] = None
    use_gal_mask_in_analysis: bool = False
    vclip_quantiles         : Optional[str] = None
    use_gal_mask_in_display : bool = False
    display_mask_only       : bool = False
    tomography_min_pix      : int  = 400
    save_json               : bool = False   # <— add this line
    fr_compatible           : bool = False

# ----------------------------
# CLI parsing
# ----------------------------

def parse_args() -> Args:
    p = argparse.ArgumentParser(
        description="SkyMap: R_BAO anisotropy with FR/KdS vector overlays",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument("catalog",               type=str,                     help="Unified catalog CSV/Parquet path")
    p.add_argument("--out",                 type=str,   default="plamb_runs/SkyMap/figs/SM_SkyMap_FR_KdS.png")
    p.add_argument("--nside",               type=int,   default=64,       help="HEALPix NSIDE (power of 2)")
    p.add_argument("--coord-in",            type=str,   default="icrs",
                   choices=["icrs", "galactic"], help="Input coord frame")
    p.add_argument("--zmin",                type=float, default=0.1)
    p.add_argument("--zmax",                type=float, default=1.5)
    p.add_argument("--gal-mask-abs-b-deg",  type=float, default=20.0,     help="Mask |b|<value (Galactic)")
    p.add_argument("--fwhm-smooth-deg",     type=float, default=8.0,      help="Gaussian smoothing FWHM (deg)")
    p.add_argument("--rbao-mode",           type=str,   default="auto",
                   choices=["auto", "DH_over_rd/DM_over_rd", "DH/DM", "none"], help="How to compute R_BAO")
    p.add_argument("--rd-const",            type=float, default=None,     help="Constant r_d if using DH/DM and per-row rd is absent")
    p.add_argument("--highlight-quantile",  type=float, default=0.9,      help="Quantile for IMF highlighting")
    p.add_argument("--fr-vec-lb-deg",       type=str,   default=None,     help="FR  vector as 'l_deg,b_deg'")
    p.add_argument("--kds-vec-lb-deg",      type=str,   default=None,     help="KdS vector as 'l_deg,b_deg'")
    p.add_argument("--lock-min",            type=float, default=None,     help="Minimum HHT locking value to include SN points")

    p.add_argument("--dipole-search-nside", type=int,   default=32,       help="NSIDE for coarse dipole search")
    p.add_argument("--cmap",                type=str,   default="coolwarm")
    p.add_argument("--vmin",                type=float, default=None)
    p.add_argument("--vmax",                type=float, default=None)
    p.add_argument("--inset-loc",           type=str,   default="lower right",
                   choices=["upper left", "upper right", "lower left", "lower right"], help="Inset summary box location")
    p.add_argument("--title",               type=str,   default=None,     help="Custom title for the Mollweide plot")
    p.add_argument("--save-stats",          action="store_true",          help="Save a TXT sidecar with dipole + tomography stats next to --out")

    p.add_argument("--perm-n",              type=int, default=500,        help="Permutation count for null test")
    p.add_argument("--perm-seed",           type=int, default=42,         help="Random seed for permutation test")
    p.add_argument("--z-bins",              type=str, default=None,       help='Custom tomography bins: "0.10,0.40;0.40,0.70;0.70,1.00;1.00,1.20"')

    p.add_argument("--use-gal-mask-in-analysis", action="store_true",     help="Also apply |b| mask to the map before dipole fit, permutations, and tomography.")
    p.add_argument("--vclip-quantiles",     type=str, default=None,       help='Quantile clip for color scale like "0.02,0.98" (over finite map pixels).')
    p.add_argument('--use-gal-mask-in-display',  action='store_true',     help='Also apply the Galactic mask to the rendered map.')
    p.add_argument('--display-mask-only',        action='store_true',     help='Render only the |b|-mask for visual sanity check.')
    p.add_argument("--tomography-min-pix",  type=int, default=400,        help="Skip a z-bin if valid HEALPix pixels < this number (after masking).")
    p.add_argument("--save-json",                action="store_true",     help="Write a JSON summary next to --out")
    p.add_argument("--fr-compatible", action="store_true",                help="Use FR-safe BAO ratio (DH/DM) and annotate outputs accordingly.")

    a = p.parse_args()

    return Args(
        catalog                     = a.catalog,
        out                         = a.out,
        nside                       = a.nside,
        coord_in                    = a.coord_in.lower(),
        zmin                        = a.zmin,
        zmax                        = a.zmax,
        gal_mask_abs_b_deg          = a.gal_mask_abs_b_deg,
        fwhm_smooth_deg             = a.fwhm_smooth_deg,
        rbao_mode                   = a.rbao_mode,
        rd_const                    = a.rd_const,
        highlight_quantile          = a.highlight_quantile,
        fr_vec_lb_deg               = parse_lb(a.fr_vec_lb_deg),
        kds_vec_lb_deg              = parse_lb(a.kds_vec_lb_deg),
        lock_min                    = a.lock_min,
        dipole_search_nside         = a.dipole_search_nside,
        cmap                        = a.cmap,
        vmin                        = a.vmin,
        vmax                        = a.vmax,
        title                       = a.title,
        inset_loc                   = a.inset_loc,
        perm_n                      = a.perm_n,
        perm_seed                   = a.perm_seed,
        z_bins                      = parse_zbins(a.z_bins),
        save_stats                  = a.save_stats,
        use_gal_mask_in_analysis    = a.use_gal_mask_in_analysis,
        vclip_quantiles             = a.vclip_quantiles,
        use_gal_mask_in_display     = a.use_gal_mask_in_display,
        display_mask_only           = a.display_mask_only,
        tomography_min_pix          = a.tomography_min_pix,
        save_json                   = a.save_json,      # <— add this line
        fr_compatible               = a.fr_compatible,
    )


def parse_lb(s: Optional[str]) -> Optional[Tuple[float, float]]:
    if s is None:
        return None
    try:
        l_str, b_str = s.split(",")
        return float(l_str.strip()), float(b_str.strip())
    except Exception as e:
        raise ValueError("Vector must be 'l_deg,b_deg', e.g. '264.0,48.0'") from e


def parse_zbins(s: Optional[str]) -> Optional[List[Tuple[float, float]]]:
    if not s:
        return None
    out = []
    for chunk in s.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        lo_str, hi_str = chunk.split(",")
        lo, hi = float(lo_str), float(hi_str)
        if not (np.isfinite(lo) and np.isfinite(hi)):
            raise ValueError("Non-finite in --z-bins.")
        if hi <= lo:
            raise ValueError("Each z-bin must have hi > lo.")
        out.append((lo, hi))
    return out


# ----------------------------
# Coordinates
# ----------------------------

def radec_to_gal_l_b(ra_deg: np.ndarray, dec_deg: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    sc  = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    gal = sc.galactic
    return gal.l.deg, gal.b.deg


def ensure_galactic(df: pd.DataFrame, coord_in: str, col: Columns) -> pd.DataFrame:
    if coord_in == "galactic":
        if "l_deg" in df.columns and "b_deg" in df.columns:
            return df.copy()
        raise ValueError("When --coord-in=galactic, the catalog must have l_deg,b_deg columns.")
    
    # coord_in == "icrs"
    if col.ra_deg not in df.columns or col.dec_deg not in df.columns:
        raise ValueError("Catalog missing 'ra_deg' and/or 'dec_deg' required for --coord-in=icrs.")
    l_deg, b_deg = radec_to_gal_l_b(df[col.ra_deg].to_numpy(), df[col.dec_deg].to_numpy())
    out          = df.copy(); out["l_deg"], out["b_deg"] = l_deg, b_deg
    return out

# ----------------------------
# R_BAO construction
# ----------------------------

def compute_rbao(df: pd.DataFrame, col: Columns, mode: str, rd_const: Optional[float]) -> np.ndarray:
    """
    Compute R_BAO with a single return.
    - Supports: "none", "auto", "dh_over_rd/dm_over_rd", "dh/dm"
    - Always returns float array; divides safely with masks.
    """
    m   = (mode or "auto").strip().lower()
    n   = len(df)
    out = np.full(n, np.nan, dtype=float)

    if m == "none":
        return out

    has_over  = (col.DH_over_rd in df.columns) and (col.DM_over_rd in df.columns)
    has_plain = (col.DH in df.columns)        and (col.DM in df.columns)

    if m == "auto":
        if   has_over:  m = "dh_over_rd/dm_over_rd"
        elif has_plain: m = "dh/dm"
        else:
            raise ValueError("Cannot auto-compute R_BAO: missing DH/DM columns.")

    if m == "dh_over_rd/dm_over_rd":
        num = pd.to_numeric(df[col.DH_over_rd], errors="coerce").to_numpy(dtype=float)
        den = pd.to_numeric(df[col.DM_over_rd], errors="coerce").to_numpy(dtype=float)
        np.divide(num, den, out=out, where=(den != 0) & np.isfinite(den))

    elif m == "dh/dm":
        dh = pd.to_numeric(df[col.DH], errors="coerce").to_numpy(dtype=float)
        dm = pd.to_numeric(df[col.DM], errors="coerce").to_numpy(dtype=float)
        np.divide(dh, dm, out=out, where=(dm != 0) & np.isfinite(dm))

    else:
        raise ValueError(f"Unknown rbao mode: {mode}")

    # Clean up any infinities that slipped through
    out[~np.isfinite(out)] = np.nan
    return out


def _valid_pix(map_rbao: np.ndarray):
    idx = np.where((map_rbao != hp.UNSEEN) & np.isfinite(map_rbao))[0]
    return idx, map_rbao[idx].astype(float)


def _pix_unitvecs(nside: int, idx: np.ndarray) -> np.ndarray:
    theta, phi = hp.pix2ang(nside, idx, nest=False)
    x          = np.sin(theta) * np.cos(phi)
    y          = np.sin(theta) * np.sin(phi)
    z          = np.cos(theta)
    return np.vstack([x, y, z]).T


# ----------------------------
# Dipole & helpers (mask-aware centered fit)
# ----------------------------

def lb_to_unit(l_deg: float, b_deg: float) -> np.ndarray:
    l_rad = math.radians(l_deg)
    b_rad = math.radians(b_deg)
    return np.array([
        math.cos(b_rad) * math.cos(l_rad),
        math.cos(b_rad) * math.sin(l_rad),
        math.sin(b_rad),
    ], dtype=float)


def unit_to_lb(v: np.ndarray) -> Tuple[float, float]:
    x, y, z = v
    b_rad   = math.asin(max(-1.0, min(1.0, z)))
    l_rad   = math.atan2(y, x)
    return (math.degrees(l_rad) % 360.0, math.degrees(b_rad))


def great_circle_between_antipodes(n_hat: np.ndarray, npts: int = 128) -> Tuple[np.ndarray, np.ndarray]:
    n_hat  = n_hat / np.linalg.norm(n_hat)
    z_axis = np.array([0.0, 0.0, 1.0])
    x_axis = np.array([1.0, 0.0, 0.0])
    u      = np.cross(n_hat, z_axis)
    if np.linalg.norm(u) < 1e-8:
        u = np.cross(n_hat, x_axis)
    u /= np.linalg.norm(u)

    psi = np.linspace(0.0, math.pi, npts)
    pts = np.cos(psi)[:, None] * n_hat[None, :] + np.sin(psi)[:, None] * u[None, :]
    l_list, b_list = [], []
    for p in pts:
        l, b = unit_to_lb(p / np.linalg.norm(p))
        l_list.append(l)
        b_list.append(b)
    return np.array(l_list), np.array(b_list)


def angsep_lb(a: Optional[Tuple[float, float]], b: Optional[Tuple[float, float]]) -> float:
    if a is None or b is None:
        return np.nan
    if not (np.isfinite(a).all() and np.isfinite(b).all()):
        return np.nan
    ua, ub = lb_to_unit(*a), lb_to_unit(*b)
    dot    = float(np.clip(np.dot(ua, ub), -1.0, 1.0))
    return math.degrees(math.acos(dot))


def min_axis_sep(theta_deg: float) -> float:
    return np.nan if not np.isfinite(theta_deg) else min(theta_deg, 180.0 - theta_deg)


def fit_dipole_linear_centered(map_rbao: np.ndarray, ridge: float = 1e-2):
    idx = np.where((map_rbao != hp.UNSEEN) & np.isfinite(map_rbao))[0]
    if idx.size == 0:
        return None

    y          = map_rbao[idx].astype(float)
    y0         = float(np.nanmean(y))
    y_         = y - y0

    nside      = hp.get_nside(map_rbao)
    theta, phi = hp.pix2ang(nside, idx, nest=False)
    V          = np.vstack([np.sin(theta) * np.cos(phi),
                           np.sin(theta) * np.sin(phi),
                           np.cos(theta)]).T

    Vc         = V - np.nanmean(V, axis=0, keepdims=True)
    G          = Vc.T @ Vc
    lam        = ridge * (np.trace(G) / 3.0) if np.isfinite(np.trace(G)) else ridge

    try:
        a = np.linalg.solve(G + lam * np.eye(3), Vc.T @ y_)
    except np.linalg.LinAlgError:
        return None

    a_norm = float(np.linalg.norm(a))
    axis_lb = (np.nan, np.nan)
    if a_norm > 0 and np.isfinite(a_norm):
        axis_lb = unit_to_lb(a / a_norm)

    y_dip = Vc @ a
    amp_eff = 0.5 * (np.nanmax(y_dip) - np.nanmin(y_dip))
    frac_eff = amp_eff / abs(y0) if (np.isfinite(y0) and abs(y0) > 0) else np.nan
    cond = np.linalg.cond(G + lam * np.eye(3))

    return dict(R0=y0, a_vec=a, axis_lb=axis_lb,
                amp_eff=amp_eff, frac_eff=frac_eff, cond=cond)


def hemisphere_contrast_projected(map_rbao: np.ndarray, axis_lb):
    idx, y = _valid_pix(map_rbao)
    if len(idx) == 0 or axis_lb is None or not np.isfinite(axis_lb).all():
        return np.nan, np.nan
    V = _pix_unitvecs(hp.get_nside(map_rbao), idx)
    d = lb_to_unit(*axis_lb)
    s = V @ d
    med = np.median(s)
    pos = y[s >= med]
    neg = y[s < med]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan, np.nan
    delta = float(np.nanmean(pos) - np.nanmean(neg))
    R0 = float(np.nanmean(y))
    return delta, (delta / abs(R0) if np.isfinite(R0) and abs(R0) > 0 else np.nan)


def galactic_b_for_nside(nside: int) -> np.ndarray:
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))  # colatitude
    b = np.degrees(0.5 * pi - theta)
    return b


def build_gal_mask(nside: int, abs_b_deg: float) -> np.ndarray:
    b = galactic_b_for_nside(nside)
    return (np.abs(b) >= abs_b_deg)  # True = keep


def dipole_permutation_pvalue(map_rbao: np.ndarray, n_sims: int = 500, seed: int = 42):
    rng = np.random.default_rng(seed)
    idx = np.where((map_rbao != hp.UNSEEN) & np.isfinite(map_rbao))[0]
    if idx.size == 0:
        return np.nan, np.array([])

    obs = fit_dipole_linear_centered(map_rbao, ridge=1e-2)
    if obs is None:
        return np.nan, np.array([])
    amp_obs    = obs["amp_eff"]

    nside      = hp.get_nside(map_rbao)
    theta, phi = hp.pix2ang(nside, idx, nest=False)
    V          = np.vstack([np.sin(theta) * np.cos(phi),
                 np.sin(theta) * np.sin(phi),
                 np.cos(theta)]).T
    Vc         = V - np.nanmean(V, axis=0, keepdims=True)
    G          = Vc.T @ Vc
    lam        = 1e-2 * (np.trace(G) / 3.0) if np.isfinite(np.trace(G)) else 1e-2

    y = map_rbao[idx].astype(float)
    amps = np.empty(n_sims, dtype=float)
    for i in range(n_sims):
        y_shuf = y[rng.permutation(y.size)]
        y0     = float(np.nanmean(y_shuf))
        y_     = y_shuf - y0
        try:
            a = np.linalg.solve(G + lam * np.eye(3), Vc.T @ y_)
        except np.linalg.LinAlgError:
            amps[i] = np.nan
            continue
        y_dip   = Vc @ a
        amps[i] = 0.5 * (np.nanmax(y_dip) - np.nanmin(y_dip))

    p = (np.sum(amps >= amp_obs) + 1) / (n_sims + 1)

    return p, amps


# ----------------------------
# HEALPix binning & smoothing (mask-aware)
# ----------------------------

def make_map_rbao(
    l_deg           : np.ndarray,
    b_deg           : np.ndarray,
    z               : np.ndarray,
    rbao            : np.ndarray,
    nside           : int,
    zmin            : float,
    zmax            : float,
    fwhm_deg        : float,
    weights         : Optional[np.ndarray] = None,
    analysis_mask   : Optional[np.ndarray] = None,
) -> np.ndarray:
    """Build a smoothed HEALPix map of R_BAO with an optional |b|-analysis mask.
       Uses numerator/denominator smoothing to prevent mask-edge leakage."""
    sel = (
        np.isfinite(rbao)
        & np.isfinite(l_deg)
        & np.isfinite(b_deg)
        & (z >= float(zmin))
        & (z <  float(zmax))
    )

    npix = hp.nside2npix(nside)
    if not np.any(sel):
        return np.full(npix, hp.UNSEEN, dtype=float)

    phi   = np.deg2rad(l_deg[sel])
    theta = np.deg2rad(90.0 - b_deg[sel])

    m_sum = np.zeros(npix, dtype=float)
    w_sum = np.zeros(npix, dtype=float)

    pix   = hp.ang2pix(nside, theta, phi, nest=False)
    rvals = rbao[sel].astype(float)
    wvals = (np.ones_like(rvals) if weights is None else weights[sel].astype(float))

    np.add.at(m_sum, pix, rvals * wvals)
    np.add.at(w_sum, pix, wvals)

    with np.errstate(divide="ignore", invalid="ignore"):
        m_map = np.divide(m_sum, w_sum, out=np.full(npix, hp.UNSEEN), where=w_sum > 0)

    # --- mask-aware smoothing via normalized weights ---
    valid = (m_map != hp.UNSEEN)
    if analysis_mask is not None:
        valid &= analysis_mask.astype(bool)

    m_in     = np.where(valid, m_map, 0.0)
    w_in     = valid.astype(float)

    fwhm_rad = np.deg2rad(fwhm_deg)
    m_num    = hp.sphtfunc.smoothing(m_in, fwhm=fwhm_rad)
    w_den    = hp.sphtfunc.smoothing(w_in, fwhm=fwhm_rad)

    out      = np.full(npix, hp.UNSEEN, dtype=float)
    ok       = w_den > 1e-6
    out[ok]  = m_num[ok] / w_den[ok]

    # Enforce strict mask (no accidental bleed outside keep-region)
    if analysis_mask is not None:
        out[~analysis_mask.astype(bool)] = hp.UNSEEN

    return out


# ----------------------------
# Analysis helpers (Y1m, cross-power, tomography)
# ----------------------------

def y1_dipole_amplitude(map_in: np.ndarray, lmax: int = 1) -> float:
    if not np.any((map_in != hp.UNSEEN) & np.isfinite(map_in)):
        return float("nan")
    
    m         = np.copy(map_in)
    valid     = (m != hp.UNSEEN) & np.isfinite(m)
    m[~valid] = 0.0
    alm       = hp.map2alm(m, lmax=lmax, iter=0)
    a10       = alm[hp.Alm.getidx(lmax, 1, 0)]
    a11       = alm[hp.Alm.getidx(lmax, 1, 1)]
    A         = np.sqrt(np.abs(a10) ** 2 + 2.0 * (a11.real ** 2 + a11.imag ** 2))
    return float(A)


def cross_power_C1(map_a: np.ndarray, map_b: np.ndarray, nest: bool = False) -> float:
    nside_a = hp.get_nside(map_a)
    nside_b = hp.get_nside(map_b)
    nside   = max(nside_a, nside_b)

    def _prep(m, nside_target):
        if hp.get_nside(m) != nside_target:
            m = hp.ud_grade(m, nside_out=nside_target, order_in="RING", order_out="RING", power=-0)
        v = (m != hp.UNSEEN) & np.isfinite(m)
        m = np.where(v, m, 0.0)
        return m

    a  = _prep(map_a, nside)
    b  = _prep(map_b, nside)
    Cl = hp.anafast(a, b, lmax=1, alm=False, iter=0)
    return float(Cl[1])


def make_map_for_bin(df_gal: pd.DataFrame,
                     rbao:          np.ndarray,
                     zlo:           float,
                     zhi:           float,
                     args:          "Args",
                     col:           "Columns",
                     weights:       Optional[np.ndarray],
                     analysis_mask: Optional[np.ndarray]) -> np.ndarray:
    return make_map_rbao(
        l_deg           =df_gal["l_deg"].to_numpy(),
        b_deg           =df_gal["b_deg"].to_numpy(),
        z               =pd.to_numeric(df_gal[col.z], errors="coerce").to_numpy(),
        rbao            =rbao,
        nside           =args.nside,
        zmin            =zlo,
        zmax            =zhi,
        fwhm_deg        =args.fwhm_smooth_deg,
        weights         =weights,
        analysis_mask   =analysis_mask,
    )


# ----------------------------
# Plotting
# ----------------------------

def plot_skymap(
    map_rbao:   np.ndarray,
    args:       Args,
    df_gal:     pd.DataFrame,
    col:        Columns,
    dipole_lb:  Optional[Tuple[float, float]] = None,
) -> plt.Figure:
    # Build display map (mask only if requested)
    show_map = map_rbao.copy()
    if args.use_gal_mask_in_display or args.display_mask_only:
        nside                 = hp.get_nside(map_rbao)
        theta_mask, _phi_mask = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), nest=False)
        b_mask_deg            = 90.0 - np.rad2deg(theta_mask)
        if args.display_mask_only:
            mask_only       = np.full_like(show_map, hp.UNSEEN, dtype=float)
            keep            = np.where(np.abs(b_mask_deg) >= args.gal_mask_abs_b_deg)[0]
            mask_only[keep] = 1.0
            show_map        = mask_only
        else:
            show_map[np.where(np.abs(b_mask_deg) < args.gal_mask_abs_b_deg)[0]] = hp.UNSEEN

    has_field = np.any(show_map != hp.UNSEEN)

    # Optional quantile clipping (without mutating args)
    vmin, vmax   = args.vmin, args.vmax
    if args.vclip_quantiles and (vmin is None and vmax is None):
        qlo, qhi = [float(x) for x in args.vclip_quantiles.split(",")]
        if not (0.0 < qlo < qhi < 1.0):
            raise ValueError("--vclip-quantiles must be two numbers with 0<qlo<qhi<1")
        
        finite   = show_map[(show_map != hp.UNSEEN) & np.isfinite(show_map)]
        if finite.size:
            vmin, vmax = np.quantile(finite, [qlo, qhi])

    title       = args.title if args.title else f"R_BAO (D_H/D_M), z∈[{args.zmin:.2f},{args.zmax:.2f}], NSIDE={args.nside}"
    norm_choice = None if (vmin is not None or vmax is not None) else "hist"

    hp.mollview(
        show_map,
        title       = title   if not args.display_mask_only else "Galactic |b|-mask",
        unit        = "R_BAO" if not args.display_mask_only else "mask",
        norm        = norm_choice,
        min         = vmin, 
        max         = vmax,
        cmap        = args.cmap,
        cbar        = has_field,
        coord       = ["G"], notext=False, badcolor="#202020",
    )

    # Overlay IMF2/IMF3 points
    for name, en_col, marker in [("IMF2", col.imf2_E, "^"), ("IMF3", col.imf3_E, "s")]:
        if en_col in df_gal.columns:
            df_plot = df_gal.copy()
            if args.lock_min is not None and col.hht_lock in df_plot.columns:
                df_plot = df_plot[pd.to_numeric(df_plot[col.hht_lock], errors="coerce") >= float(args.lock_min)]
            ener    = pd.to_numeric(df_plot[en_col], errors="coerce")
            qcut    = ener.quantile(args.highlight_quantile) if ener.notna().any() else np.nan
            sel_pts = ener >= qcut if np.isfinite(qcut) else np.array([False] * len(df_plot))
            l_pts   = df_plot.loc[sel_pts, "l_deg"].to_numpy()
            b_pts   = df_plot.loc[sel_pts, "b_deg"].to_numpy()
            if len(l_pts) > 0:
                hp.projscatter(l_pts, b_pts, lonlat=True, s=18, marker=marker, alpha=0.8,
                               label=f"{name} top {args.highlight_quantile:.2f}")

    # R_BAO dipole line & pole
    if dipole_lb is not None and np.isfinite(dipole_lb).all():
        n_hat          = lb_to_unit(*dipole_lb)
        l_line, b_line = great_circle_between_antipodes(n_hat)
        hp.projplot(l_line, b_line, lonlat=True, ls="-", lw=2.0, alpha=0.9, label="R_BAO dipole")
        hp.projscatter([dipole_lb[0]], [dipole_lb[1]], lonlat=True, marker="o", s=40, alpha=0.9)

    # FR and KdS axes
    if args.fr_vec_lb_deg is not None:
        n_hat          = lb_to_unit(*args.fr_vec_lb_deg)
        l_line, b_line = great_circle_between_antipodes(n_hat)
        hp.projplot(l_line, b_line, lonlat=True, ls="--", lw=1.8, alpha=0.95, label="FR axis")
        hp.projscatter([args.fr_vec_lb_deg[0]], [args.fr_vec_lb_deg[1]], lonlat=True, marker="*", s=60, alpha=0.95)
    if args.kds_vec_lb_deg is not None:
        n_hat          = lb_to_unit(*args.kds_vec_lb_deg)
        l_line, b_line = great_circle_between_antipodes(n_hat)
        hp.projplot(l_line, b_line, lonlat=True, ls="-.", lw=1.8, alpha=0.95, label="KdS axis")
        hp.projscatter([args.kds_vec_lb_deg[0]], [args.kds_vec_lb_deg[1]], lonlat=True, marker="D", s=52, alpha=0.95)

    # Legend (dedupe labels)
    fig = plt.gcf()
    ax  = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    if labels:
        seen, uh, ul = set(), [], []
        for h, l in zip(handles, labels):
            if l not in seen:
                seen.add(l); uh.append(h); ul.append(l)
        ax.legend(uh, ul, loc="lower left", fontsize=8, framealpha=0.7)

    return fig


def render_inset_summary(fig: plt.Figure, lines: List[str], loc: str = "lower right") -> None:
    if not lines:
        return
    ax = fig.axes[0]
    anchors = {
        "upper left":  (0.02, 0.98, "top", "left"),
        "upper right": (0.98, 0.98, "top", "right"),
        "lower left":  (0.02, 0.02, "bottom", "left"),
        "lower right": (0.98, 0.02, "bottom", "right"),
    }
    x, y, va, ha = anchors.get(loc, anchors["lower right"])
    text         = "\n".join(lines)
    ax.text(
        x, y, text,
        transform = ax.transAxes,
        va        = va, 
        ha        = ha,
        fontsize  = 8, color="white",
        bbox      = dict(boxstyle="round,pad=0.35", facecolor="black", alpha=0.42, ec="none"),
        zorder    = 20,
    )


def apply_abs_b_mask_to_map(m: np.ndarray, abs_b_min_deg: float) -> np.ndarray:
    if abs_b_min_deg <= 0:
        return m
    m2                                = m.copy()
    nside                             = hp.get_nside(m2)
    theta, _phi                       = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), nest=False)
    b_deg                             = 90.0 - np.rad2deg(theta)
    m2[np.abs(b_deg) < abs_b_min_deg] = hp.UNSEEN
    return m2


def _coverage(tag: str, m: np.ndarray) -> None:
    v = (m != hp.UNSEEN) & np.isfinite(m)
    print(f"[analysis] {tag}: valid_pix = {v.sum():d}/{v.size:d} ({100*v.mean():.1f}%)")

def cos_template_from_axis(nside, l_deg, b_deg):
    idx = np.arange(hp.nside2npix(nside))
    V   = _pix_unitvecs(nside, idx)
    d   = lb_to_unit(l_deg, b_deg)
    tpl = V @ d
    return tpl

# ----------------------------
# Main
# ----------------------------

def main() -> None:
    args = parse_args()
    
    # --- FR invariants & output path normalisation ---
    # Make DH/DM the recorded mode for FR runs and store an absolute PNG path in JSON
    # --- FR tag & absolute output path (do not force BAO column set) ---
    args.out = os.path.abspath(args.out)
    if args.fr_compatible and "FR" not in (args.title or ""):
        args.title = f"{(args.title or '').strip()} [FR]"

    # Fail fast if NSIDE isn’t a power of two
    if args.nside <= 0 or (args.nside & (args.nside - 1)) != 0:
        raise ValueError("ERROR --nside must be a power of two and > 0.")
        
    dirpath = os.path.dirname(args.out)
    if dirpath:
        os.makedirs(dirpath, exist_ok=True)
        
    # Load catalog
    ext = os.path.splitext(args.catalog)[1].lower()
    if ext in [".csv", ".txt", ".tsv"]:
        df = pd.read_csv(args.catalog)
    elif ext in [".parquet", ".pq"]:
        df = pd.read_parquet(args.catalog)
    else:
        raise ValueError("Unsupported catalog format. Use CSV or Parquet.")

    col    = Columns()

    # Compute R_BAO and ensure Galactic coordinates
    rbao   = compute_rbao(df, col, args.rbao_mode, args.rd_const)
    print(f"[rbao] mode={args.rbao_mode}  finite={np.isfinite(rbao).sum()}/{len(rbao)}")

    df_gal = ensure_galactic(df, args.coord_in, col)

    # Optional inverse-variance weights (approximate)
    if col.sigma_DH_over_rd in df.columns and col.sigma_DM_over_rd in df.columns:
        a  = pd.to_numeric(df[col.DH_over_rd], errors="coerce").to_numpy()
        b  = pd.to_numeric(df[col.DM_over_rd], errors="coerce").to_numpy()
        sa = pd.to_numeric(df[col.sigma_DH_over_rd], errors="coerce").to_numpy()
        sb = pd.to_numeric(df[col.sigma_DM_over_rd], errors="coerce").to_numpy()
        with np.errstate(divide="ignore", invalid="ignore"):
            rel2                           = np.square(sa / a) + np.square(sb / b)
            sigma_r                        = np.abs(rbao) * np.sqrt(rel2)
            weights                        = 1.0 / np.square(sigma_r)
            weights[~np.isfinite(weights)] = 0.0
    else:
        weights = None

    # Build |b|-mask for analysis if requested.
    gal_mask = None
    if args.use_gal_mask_in_analysis and args.gal_mask_abs_b_deg > 0:
        gal_mask = build_gal_mask(args.nside, args.gal_mask_abs_b_deg)
        print(f"[gal-mask] |b|≥{args.gal_mask_abs_b_deg:.1f}°, keep ≈ {100.0*gal_mask.mean():.1f}% of sky")

    # Global map
    map_rbao = make_map_rbao(
        l_deg           = df_gal["l_deg"].to_numpy(),
        b_deg           = df_gal["b_deg"].to_numpy(),
        z               = pd.to_numeric(df_gal[col.z], errors="coerce").to_numpy(),
        rbao            = rbao,
        nside           = args.nside,
        zmin            = args.zmin,
        zmax            = args.zmax,
        fwhm_deg        = args.fwhm_smooth_deg,
        weights         = weights,
        analysis_mask   = gal_mask,
    )
    
    finite = map_rbao[(map_rbao != hp.UNSEEN) & np.isfinite(map_rbao)]
    flat = finite.size and (np.nanstd(finite) < max(1e-8, 1e-6*abs(np.nanmean(finite))))
    if flat:
        print("[guard] Map is numerically flat → skipping dipole/permutation/tomography.")
        fig = plot_skymap(map_rbao, args, df_gal, col, dipole_lb=None)
        render_inset_summary(fig, [], loc=args.inset_loc)
        fig.savefig(args.out, dpi=300, bbox_inches="tight")
        hp.write_map(os.path.splitext(args.out)[0] + ".fits", map_rbao, overwrite=True)
        
        if args.save_json:
            import json
            base = os.path.splitext(args.out)[0]
            payload = {
                "out": args.out,
                "nside": args.nside,
                "zmin": args.zmin, "zmax": args.zmax,
                "fwhm_deg": args.fwhm_smooth_deg,
                "rbao_mode": args.rbao_mode,
                "fr_compatible": bool(args.fr_compatible),
                "dipole": {
                    "l_deg": None, "b_deg": None,
                    "R0": float(np.nanmean(finite)) if finite.size else None,
                    "A_eff": None, "frac": None, "p_perm": None,
                    "R2": None, "rho": None, "dFR_deg": None, "dKdS_deg": None,
                },
                "tomography": [],
                "note": "flat_map_guard"
            }
            with open(base + "_summary.json", "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2)
        return

    
    
    # Force UNSEEN inside masked zone for analysis consistency
    if args.use_gal_mask_in_analysis and args.gal_mask_abs_b_deg > 0:
        map_rbao = apply_abs_b_mask_to_map(map_rbao, args.gal_mask_abs_b_deg)

    _coverage("global", map_rbao)

    # Dipole fit & diagnostics
    dipole_stats: Optional[dict] = None
    dipole_lb = (np.nan, np.nan)
    R2        = np.nan
    rho       = np.nan
    pval      = np.nan

    if args.rbao_mode.lower() != "none":
        dipole_stats = fit_dipole_linear_centered(map_rbao, ridge=1e-2)
        if dipole_stats is not None:
            dipole_lb   = dipole_stats["axis_lb"]

            dH, dH_frac = hemisphere_contrast_projected(map_rbao, dipole_lb)
            pval, _     = dipole_permutation_pvalue(map_rbao, n_sims=args.perm_n, seed=args.perm_seed)

            print(
                f"[dipole-fit ⟂ 1] R0 = {dipole_stats['R0']:.6g} "
                f"A_eff = {dipole_stats['amp_eff']:.6g} "
                f"(frac = {100 * dipole_stats['frac_eff']:.3f}% )  "
                f"cond  ≈ {dipole_stats['cond']:.2e}"
            )
            print(f"[hemi-contrast] Δ = {dH:.6g}  (frac = {100 * dH_frac:.3f}% )")
            print(f"[null|permute] n={args.perm_n}  p(A_null ≥ A_obs) ≈ {pval:.4f}")

            if args.fr_vec_lb_deg is not None:
                raw = angsep_lb(dipole_lb, args.fr_vec_lb_deg)
                print(f"[angle ] dipole–FR   = {raw:7.2f} deg   (axis-min = {min_axis_sep(raw):7.2f} deg)")
            if args.kds_vec_lb_deg is not None:
                raw = angsep_lb(dipole_lb, args.kds_vec_lb_deg)
                print(f"[angle ] dipole–KdS  = {raw:7.2f} deg   (axis-min = {min_axis_sep(raw):7.2f} deg)")
            print(f"[dipole-dir] l = {dipole_lb[0]:7.2f}   b = {dipole_lb[1]:7.2f}")

            # 1D regression along fitted axis
            valid = (map_rbao != hp.UNSEEN) & np.isfinite(map_rbao)
            if np.any(valid):
                pix        = np.where(valid)[0]
                y          = map_rbao[valid].astype(float)
                y0         = float(np.nanmean(y))

                theta, phi = hp.pix2ang(args.nside, pix, nest=False)
                ux, uy, uz = lb_to_unit(*dipole_lb)
                sx         = np.sin(theta) * np.cos(phi)
                sy         = np.sin(theta) * np.sin(phi)
                sz         = np.cos(theta)
                s          = (sx * ux + sy * uy + sz * uz).astype(float)

                s_c        = s - s.mean()
                y_c        = y - y0
                den        = float(np.dot(s_c, s_c))
                beta       = 0.0 if den == 0.0 else float(np.dot(y_c, s_c) / den)
                yhat       = y0 + beta * s_c

                SSE        = float(np.sum((y - yhat) ** 2))
                SST        = float(np.sum((y - y0) ** 2))
                R2         = 0.0 if SST == 0.0 else 1.0 - SSE / SST
                rho        = float(np.corrcoef(y, s)[0, 1]) if np.std(s) > 0 else np.nan

                A_sigma    = abs(beta) * float(np.std(s, ddof=1))
                A_range    = abs(beta) * (float(np.max(s)) - float(np.min(s))) / 2.0
                frac_sigma = (A_sigma / y0) * 100.0 if y0 != 0 else np.nan
                frac_range = (A_range / y0) * 100.0 if y0 != 0 else np.nan

                print(f"[fit-quality] R^2 ≈ {R2:.3f} (≈ corr^2), corr(y,cosθ_axis) ≈ {rho:.3f}")
                print(
                    f"[amplitude]  A_sigma ≈ {A_sigma:.6f}  ({frac_sigma:.3f}% of mean), "
                    f"A_range ≈ {A_range:.6f}  ({frac_range:.3f}% of mean)"
                )
            else:
                print("[fit-quality] skipped (no valid pixels).")
        else:
            print("[info] Dipole fit failed (no valid field).")
    else:
        print("[info] rbao-mode=none → overlays only; dipole stats skipped.")

    # Y1 sanity check
    try:
        A_y1 = y1_dipole_amplitude(map_rbao, lmax=1)
        print(f"[Y1-check] |a_l=1| ≈ {A_y1:.6g}  (mask-biased sanity check)")
    except Exception as e:
        print(f"[Y1-check] skipped: {e}")

    # Cross-power placeholder
    print("[X-power] no template_map provided; skipping cross-power.")

    # Redshift tomography
    z_bins   = args.z_bins if args.z_bins else [(0.10, 0.40), (0.40, 0.70), (0.70, 1.00), (1.00, 1.20)]
    tom_rows = []
    tom_msgs = []   # ← for inset notes about skipped bins
    
    for (zlo, zhi) in z_bins:
        m_bin = make_map_for_bin(df_gal, rbao, zlo, zhi, args, col, weights,
                                 analysis_mask=gal_mask)
    
        # Count finite, unmasked pixels in this bin
        valid = (m_bin != hp.UNSEEN) & np.isfinite(m_bin)
        n_val = int(valid.sum())
        npix  = int(m_bin.size)
        frac  = 100.0 * (n_val / npix if npix else 0.0)
        print(f"[analysis] z∈[{zlo:.2f},{zhi:.2f}]: valid_pix = {n_val}/{npix} ({frac:.1f}%)")
    
        # Guard: skip if too sparse
        if n_val < args.tomography_min_pix:
            msg = (f"  z∈[{zlo:.2f},{zhi:.2f}] → N={n_val} "
                   f"(<{args.tomography_min_pix}), skipped")
            print(f"[tomography] {msg}")
            tom_msgs.append(msg)
            continue
    
        stats = None
        if args.rbao_mode.lower() != "none":
            stats = fit_dipole_linear_centered(m_bin, ridge=1e-2)
    
        if stats is not None and np.all(np.isfinite(stats["axis_lb"])):
            lb = stats["axis_lb"]
            dFR = dKdS = np.nan
            if args.fr_vec_lb_deg is not None:
                dFR = min_axis_sep(angsep_lb(lb, args.fr_vec_lb_deg))
            if args.kds_vec_lb_deg is not None:
                dKdS = min_axis_sep(angsep_lb(lb, args.kds_vec_lb_deg))

            tom_rows.append(dict(
                                zmin     = zlo, zmax=zhi,
                                l_deg    = float(lb[0]), 
                                b_deg    = float(lb[1]),
                                R0       = float(stats["R0"]),
                                A_eff    = float(stats["amp_eff"]),
                                frac     = float(stats["frac_eff"]),
                                cond     = float(stats["cond"]),
                                dFR_deg  = float(dFR),
                                dKdS_deg = float(dKdS),
                            ))
            
            print(f"[tomography] z∈[{zlo:.2f},{zhi:.2f}]: "
                  f"l={lb[0]:.2f} b_deg={lb[1]:.2f}  R0={stats['R0']:.6g}  "
                  f"A_eff={stats['amp_eff']:.6g} (frac={100 * stats['frac_eff']:.2f}%)")
        else:
            msg = f"  z∈[{zlo:.2f},{zhi:.2f}] → no valid dipole"
            print(f"[tomography] {msg}")
            tom_msgs.append(msg)
    
   
    if tom_rows:
        out_csv = os.path.splitext(args.out)[0] + "_tomography.csv"
        try:
            pd.DataFrame(tom_rows).to_csv(out_csv, index=False)
            print(f"[tomography] saved → {out_csv}")
        except Exception as e:
            print(f"[tomography] save skipped: {e}")  
    
    if args.fr_vec_lb_deg is not None:
        m_fr = cos_template_from_axis(args.nside, *args.fr_vec_lb_deg)
        print(f"[X-power] C1(map, FR)  ≈ {cross_power_C1(map_rbao, m_fr):.6g}")

    if args.kds_vec_lb_deg is not None:
        m_kd = cos_template_from_axis(args.nside, *args.kds_vec_lb_deg)
        print(f"[X-power] C1(map, KdS) ≈ {cross_power_C1(map_rbao, m_kd):.6g}")
        
    # Inset summary
    # --- Build summary lines for inset (after tomography loop) ---
    summary_lines = []
    if args.rbao_mode.lower() != "none" and dipole_stats is not None:
        summary_lines += [
            "R_BAO dipole (global)",
            f"  R₀ = {dipole_stats['R0']:.3f}",
            f"  A_eff = {dipole_stats['amp_eff']:.3f}  ({100*dipole_stats['frac_eff']:.1f}%)",
            f"  (l, b) = ({dipole_lb[0]:.1f}°, {dipole_lb[1]:.1f}°)",
        ]
        if args.fr_vec_lb_deg is not None and np.all(np.isfinite(args.fr_vec_lb_deg)):
            dfr = min_axis_sep(angsep_lb(dipole_lb, args.fr_vec_lb_deg))
            summary_lines.append(f"  Δ_FR ≈ {dfr:.1f}°")
        if args.kds_vec_lb_deg is not None and np.all(np.isfinite(args.kds_vec_lb_deg)):
            dkds = min_axis_sep(angsep_lb(dipole_lb, args.kds_vec_lb_deg))
            summary_lines.append(f"  Δ_KdS ≈ {dkds:.1f}°")
        if np.isfinite(pval):
            summary_lines.append(f"  perm p ≈ {pval:.3f}")
        if np.isfinite(R2) or np.isfinite(rho):
            summary_lines.append(f"  R² ≈ {R2:.2f}, ρ ≈ {rho:.2f}")

    summary_lines += [f"tom: {m}" for m in tom_msgs[:3]]  # show a few quick notes



    if args.save_json:
        import json
        base = os.path.splitext(args.out)[0]
        payload = {
            "out"           : args.out,
            "nside"         : args.nside,
            "zmin"          : args.zmin, "zmax": args.zmax,
            "fwhm_deg"      : args.fwhm_smooth_deg,
            "rbao_mode"     : args.rbao_mode,          
            "fr_compatible" : bool(args.fr_compatible),            
            "dipole": {
                "l_deg"     : float(dipole_lb[0]), "b_deg": float(dipole_lb[1]),
                "R0"        : float(dipole_stats["R0"]) if dipole_stats else None,
                "A_eff"     : float(dipole_stats["amp_eff"]) if dipole_stats else None,
                "frac"      : float(dipole_stats["frac_eff"]) if dipole_stats else None,
                "p_perm"    : float(pval) if np.isfinite(pval) else None,
                "R2"        : float(R2) if np.isfinite(R2) else None,
                "rho"       : float(rho) if np.isfinite(rho) else None,
                "dFR_deg"   : (min_axis_sep(angsep_lb(dipole_lb, args.fr_vec_lb_deg))
                            if args.fr_vec_lb_deg else None),
                "dKdS_deg"  : (min_axis_sep(angsep_lb(dipole_lb, args.kds_vec_lb_deg))
                            if args.kds_vec_lb_deg else None),
            },
            "tomography"    : tom_rows,
        }
        with open(base + "_summary.json", "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
    
    if args.fr_vec_lb_deg is None and args.kds_vec_lb_deg is None:
        print("[X-power] no template_map provided; skipping cross-power.")
    
    # --- Render map and save to file ---
    fig = plot_skymap(map_rbao, args, df_gal, col, dipole_lb)
    render_inset_summary(fig, summary_lines, loc=args.inset_loc)
    fig.savefig(args.out, dpi=300, bbox_inches="tight")
    
    if "IPython" in sys.modules:
        plt.show()        # show interactively in Spyder/Jupyter
    else:
        plt.close(fig)    # close only in non-interactive/batch

    hp.write_map(os.path.splitext(args.out)[0] + ".fits", map_rbao, overwrite=True)

    # Optional: save stats
    if args.save_stats:
        txt = os.path.splitext(args.out)[0] + "_stats.txt"
        try:
            with open(txt, "w", encoding="utf-8") as f:
                f.write("SM_SkyMap_FR_KdS stats\n")
                f.write(f"out: {args.out}\n")
                f.write(f"nside: {args.nside}, z∈[{args.zmin},{args.zmax}], fwhm: {args.fwhm_smooth_deg} deg\n")
                if dipole_stats is not None:
                    f.write(f"R0: {dipole_stats['R0']}\n")
                    f.write(f"A_eff: {dipole_stats['amp_eff']}  (frac: {dipole_stats['frac_eff']})\n")
                    f.write(f"axis (l,b): {dipole_lb}\n")
                    if np.isfinite(R2):
                        f.write(f"R2: {R2}\n")
                    if np.isfinite(rho):
                        f.write(f"rho: {rho}\n")
                    if np.isfinite(pval):
                        f.write(f"perm_p: {pval}\n")
                    if args.fr_vec_lb_deg is not None:
                        dfr = min_axis_sep(angsep_lb(dipole_lb, args.fr_vec_lb_deg))
                        f.write(f"Δ_FR (axis-min): {dfr}\n")
                    if args.kds_vec_lb_deg is not None:
                        dkds = min_axis_sep(angsep_lb(dipole_lb, args.kds_vec_lb_deg))
                        f.write(f"Δ_KdS (axis-min): {dkds}\n")
                if tom_rows:
                    f.write("\nTomography (per bin):\n")
                    for r in tom_rows:
                        f.write(
                                f"z∈[{r['zmin']:.2f},{r['zmax']:.2f}], "
                                f"l={r['l_deg']:.2f}, b={r['b_deg']:.2f}, "
                                f"R0={r['R0']:.6g}, A_eff={r['A_eff']:.6g}, frac={r['frac']:.6g}\n"
                               )
            print(f"[ok] saved stats → {txt}")
        except Exception as e:
            print(f"[stats] save skipped: {e}")


if __name__ == "__main__":
    main()
