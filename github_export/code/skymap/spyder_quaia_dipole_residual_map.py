#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 13:12:06 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_dipole_residual_map.py — subtract dipole, inspect higher order (0V1)

- Loads a Quaia CSV: RA [deg], DEC [deg], Z
- Builds a mean-z HEALPix map (Galactic coords)
- Fits a *selection-aware* redshift dipole
- Builds a dipole-model map and residual map
- Computes angular power spectrum of residuals up to lmax

Outputs (optional):
    --out-dipole-map    FITS: z_model (z0 + D·n)
    --out-residual-map  FITS: z_mean - z_model
    
%run /Users/boyde/.spyder-py3/quaia_dipole_residual_map.py \
  --csv              "/Users/boyde/Downloads/quaia_G20p0_basic.csv" \
  --nside            64 \
  --min-per-pixel    10 \
  --sel-map          "/Users/boyde/.spyder-py3/selection_function_NSIDE64_G20.0.fits" \
  --sel-min-frac     0.3 \
  --out-dipole-map   "/Users/boyde/Downloads/quaia_G20p0_nside64_zdipole_model.fits" \
  --out-residual-map "/Users/boyde/Downloads/quaia_G20p0_nside64_zresidual.fits" \
  --lmax             10
"""

import argparse
import os
import sys

import numpy  as np
import pandas as pd
import healpy as hp

from astropy.coordinates import SkyCoord
import astropy.units      as u


def build_parser():
    parser                  = argparse.ArgumentParser(
        description = "Subtract redshift dipole and analyse residuals."
    )

    parser.add_argument(
        "--csv",
        required      = True,
        help          = "Input CSV with columns RA, DEC, Z.",
    )

    parser.add_argument(
        "--nside",
        type          = int,
        default       = 64,
        help          = "HEALPix nside (default: 64).",
    )

    parser.add_argument(
        "--zmin",
        type          = float,
        default       = None,
        help          = "Minimum redshift cut (optional).",
    )

    parser.add_argument(
        "--zmax",
        type          = float,
        default       = None,
        help          = "Maximum redshift cut (optional).",
    )

    parser.add_argument(
        "--min-per-pixel",
        type          = int,
        default       = 10,
        help          = "Minimum #objects per pixel for dipole fit (default: 10).",
    )

    parser.add_argument(
        "--sel-map",
        default       = None,
        help          = "Optional selection-function HEALPix FITS (must match nside).",
    )

    parser.add_argument(
        "--sel-min-frac",
        type          = float,
        default       = 0.3,
        help          = "Minimum selection fraction (sel >= frac * max(sel)) "
                        "for pixels to be used (default: 0.3).",
    )

    parser.add_argument(
        "--out-dipole-map",
        default       = None,
        help          = "Output FITS for dipole model map z_model (optional).",
    )

    parser.add_argument(
        "--out-residual-map",
        default       = None,
        help          = "Output FITS for residual map z_mean - z_model (optional).",
    )

    parser.add_argument(
        "--lmax",
        type          = int,
        default       = 10,
        help          = "Maximum multipole l for residual C_l (default: 10).",
    )

    return parser


def load_catalog(csv_path):
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    print(f"[load] {csv_path}")
    df                      = pd.read_csv(csv_path)

    for col in ("RA", "DEC", "Z"):
        if col not in df.columns:
            raise KeyError(f"Input CSV must contain column '{col}'")

    return df


def apply_redshift_cut(df, zmin = None, zmax = None):
    mask                    = np.ones(len(df), dtype = bool)

    if zmin is not None:
        mask               &= df["Z"].values >= zmin

    if zmax is not None:
        mask               &= df["Z"].values <= zmax

    df_cut                  = df.loc[mask].reset_index(drop = True)

    print(f"[info] total objects after z-cut: {len(df_cut):d}")
    return df_cut


def icrs_to_galactic(ra_deg, dec_deg):
    coords_icrs             = SkyCoord(
        ra    = np.asarray(ra_deg) * u.deg,
        dec   = np.asarray(dec_deg) * u.deg,
        frame = "icrs",
    )
    gal                     = coords_icrs.galactic
    l_rad                   = gal.l.rad
    b_rad                   = gal.b.rad
    return l_rad, b_rad


def build_healpix_mean_z(l_rad, b_rad, z, nside):
    npix                    = hp.nside2npix(nside)
    counts_map              = np.zeros(npix, dtype = int)
    sum_z_map               = np.zeros(npix, dtype = float)

    theta                   = 0.5 * np.pi - b_rad
    phi                     = l_rad

    pix_idx                 = hp.ang2pix(nside, theta, phi, nest = False)

    for p, z_i in zip(pix_idx, z):
        counts_map[p]      += 1
        sum_z_map[p]       += z_i

    mean_z_map              = np.full(npix, hp.UNSEEN, dtype = float)
    good                    = counts_map > 0
    mean_z_map[good]        = sum_z_map[good] / counts_map[good]

    return counts_map, mean_z_map


def fit_redshift_dipole(counts_map,
                        mean_z_map,
                        nside,
                        min_per_pixel = 10,
                        sel_map       = None,
                        sel_min_frac  = 0.3):
    npix                    = len(counts_map)
    assert npix == len(mean_z_map)

    if sel_map is not None and len(sel_map) != npix:
        raise ValueError("Selection map length does not match HEALPix maps.")

    # Base mask: enough counts and finite mean-z
    good_pix                = (counts_map >= min_per_pixel) & (mean_z_map != hp.UNSEEN)

    if sel_map is not None:
        positive            = sel_map > 0.0
        if not np.any(positive):
            raise RuntimeError("Selection map has no positive entries.")
        sel_max             = np.max(sel_map[positive])
        sel_thresh          = sel_min_frac * sel_max

        print(
            f"[sel] applying selection cut: sel >= {sel_thresh:.3e} "
            f"({sel_min_frac:.2f} × max(sel) = {sel_max:.3e})"
        )

        good_pix           &= sel_map >= sel_thresh

    if not np.any(good_pix):
        raise RuntimeError("No pixels pass min-per-pixel + selection cuts.")

    pix_indices             = np.where(good_pix)[0]
    weights                 = counts_map[pix_indices].astype(float)
    z_vals                  = mean_z_map[pix_indices]

    print(f"[fit] using {len(pix_indices):d} pixels in dipole fit")

    # Weighted global mean z
    z_global_mean           = np.average(z_vals, weights = weights)

    theta_c, phi_c          = hp.pix2ang(nside, pix_indices, nest = False)

    n_x                     = np.sin(theta_c) * np.cos(phi_c)
    n_y                     = np.sin(theta_c) * np.sin(phi_c)
    n_z                     = np.cos(theta_c)

    delta_z                 = z_vals - z_global_mean

    w_sqrt                  = np.sqrt(weights)

    A                       = np.vstack((
        n_x * w_sqrt,
        n_y * w_sqrt,
        n_z * w_sqrt,
    )).T
    b                       = delta_z * w_sqrt

    D_vec, _, _, _          = np.linalg.lstsq(A, b, rcond = None)

    D_norm                  = np.linalg.norm(D_vec)
    if D_norm == 0.0:
        raise RuntimeError("Dipole fit returned zero vector.")

    frac_amplitude          = D_norm / z_global_mean

    d_unit                  = D_vec / D_norm
    d_x, d_y, d_z           = d_unit

    b_dip_rad               = np.arcsin(d_z)
    l_dip_rad               = np.arctan2(d_y, d_x)
    if l_dip_rad < 0.0:
        l_dip_rad          += 2.0 * np.pi

    l_dip_deg               = np.degrees(l_dip_rad)
    b_dip_deg               = np.degrees(b_dip_rad)

    return z_global_mean, D_vec, frac_amplitude, (l_dip_deg, b_dip_deg), good_pix


def build_dipole_model_map(z_global_mean, D_vec, nside):
    npix                    = hp.nside2npix(nside)
    theta_all, phi_all      = hp.pix2ang(nside, np.arange(npix), nest = False)

    n_x_all                 = np.sin(theta_all) * np.cos(phi_all)
    n_y_all                 = np.sin(theta_all) * np.sin(phi_all)
    n_z_all                 = np.cos(theta_all)

    n_dot_D                 = (
        n_x_all * D_vec[0]
      + n_y_all * D_vec[1]
      + n_z_all * D_vec[2]
    )

    z_model_map             = z_global_mean + n_dot_D

    return z_model_map


def main(argv = None):
    if argv is None:
        argv                = sys.argv[1:]

    parser                  = build_parser()
    args                    = parser.parse_args(argv)

    df                      = load_catalog(args.csv)
    print(f"[info] total objects before z-cut: {len(df):d}")

    if args.zmin is not None or args.zmax is not None:
        df                  = apply_redshift_cut(df, zmin = args.zmin, zmax = args.zmax)
    else:
        print("[info] no redshift cuts applied")

    if len(df) == 0:
        raise RuntimeError("No objects remain after z-cuts.")

    ra_deg                  = df["RA"].values
    dec_deg                 = df["DEC"].values
    z_vals                  = df["Z"].values

    print("[coord] converting ICRS (RA, DEC) -> Galactic (l, b)")
    l_rad, b_rad            = icrs_to_galactic(ra_deg, dec_deg)

    print(f"[healpix] building mean-z map at nside = {args.nside:d}")
    counts_map, mean_z_map  = build_healpix_mean_z(l_rad, b_rad, z_vals, args.nside)

    # Load selection map if provided
    sel_map                 = None
    if args.sel_map is not None:
        print(f"[sel] loading selection map: {args.sel_map}")
        sel_map             = hp.read_map(args.sel_map)
        nside_sel           = hp.get_nside(sel_map)
        if nside_sel != args.nside:
            raise ValueError(
                f"Selection map nside={nside_sel} does not match requested "
                f"nside={args.nside}."
            )

    print(
        f"[fit] fitting redshift dipole with min-per-pixel = {args.min_per_pixel:d}, "
        f"selection map = {'yes' if sel_map is not None else 'no'}"
    )

    z_mean, D_vec, frac_amp, (l_dip, b_dip), fit_mask = fit_redshift_dipole(
        counts_map      = counts_map,
        mean_z_map      = mean_z_map,
        nside           = args.nside,
        min_per_pixel   = args.min_per_pixel,
        sel_map         = sel_map,
        sel_min_frac    = args.sel_min_frac,
    )

    D_x, D_y, D_z          = D_vec
    D_norm                 = np.linalg.norm(D_vec)

    print("[dipole-fit]")
    print(f"  N_objects         = {len(df):d}")
    print(f"  nside             = {args.nside:d}")
    print(f"  z_global_mean     = {z_mean:.6f}")
    print(f"  |D|               = {D_norm:.6e}")
    print(f"  frac_amplitude    = {frac_amp:.6e}")
    print(f"  l_dip [deg]       = {l_dip:.3f}")
    print(f"  b_dip [deg]       = {b_dip:.3f}")
    print(
        f"  vector (D_x, D_y, D_z) "
        f"= ({D_x:.3e}, {D_y:.3e}, {D_z:.3e})"
    )

    # Build dipole model map and residuals
    print("[model] building dipole model map")
    z_model_map            = build_dipole_model_map(z_mean, D_vec, args.nside)

    residual_map           = np.full_like(mean_z_map, hp.UNSEEN, dtype = float)
    good_pixels            = (mean_z_map != hp.UNSEEN)
    residual_map[good_pixels] = mean_z_map[good_pixels] - z_model_map[good_pixels]

    # Optional outputs
    if args.out_dipole_map is not None:
        print(f"[save] dipole model map -> {args.out_dipole_map}")
        hp.write_map(args.out_dipole_map, z_model_map, overwrite = True)

    if args.out_residual_map is not None:
        print(f"[save] residual map     -> {args.out_residual_map}")
        hp.write_map(args.out_residual_map, residual_map, overwrite = True)

    # --- higher-order: power spectrum of residuals ---

    print(f"[cl] computing C_l of residuals up to lmax = {args.lmax:d}")

    npix                   = hp.nside2npix(args.nside)
    mask                   = np.zeros(npix, dtype = float)
    mask[fit_mask]         = 1.0

    # Fill UNSEEN with 0 for anafast, apply mask
    res_filled             = np.zeros(npix, dtype = float)
    valid_res              = residual_map != hp.UNSEEN
    res_filled[valid_res]  = residual_map[valid_res]

    res_masked             = res_filled * mask

    cl                     = hp.anafast(res_masked, lmax = args.lmax)

    print("[cl] residual power spectrum (l, C_l):")
    for ell in range(min(args.lmax + 1, len(cl))):
        print(f"  l = {ell:2d},  C_l = {cl[ell]:.6e}")

    print("[done]")


if __name__ == "__main__":
    main()
