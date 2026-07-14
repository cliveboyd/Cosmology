#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 13:35:30 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_random_mocks_zdipole_quad.py — Quaia random mocks for dipole + quadrupole (0V1)

Builds many mock realisations using the Quaia random catalog:

- Uses real Quaia G<20 CSV to get the observed Z distribution.
- Uses random_G20.0_10x.fits for angular selection.
- For each mock:
    * Sample N positions from random catalog.
    * Assign Z by resampling from real Z distribution.
    * Build mean-z HEALPix map (Galactic coords).
    * Fit selection-aware redshift dipole (same as main analysis).
    * Subtract dipole, build residual map.
    * Compute C_l of residuals up to lmax, extract C_2 and D_2.

Outputs:
    - prints mean/std of dipole fraction and D_2 across mocks.
    - writes CSV with one row per mock.
"""

import argparse
import os
import sys

import numpy  as np
import pandas as pd
import healpy as hp

from   astropy.table        import Table
from   astropy.coordinates  import SkyCoord
import astropy.units        as u


# ----------------------------------------------------------------------
# Utility functions (mostly copied from quaia_dipole_residual_map.py)
# ----------------------------------------------------------------------

def icrs_to_galactic(ra_deg, dec_deg):
    coords_icrs         = SkyCoord(
        ra    = np.asarray(ra_deg) * u.deg,
        dec   = np.asarray(dec_deg) * u.deg,
        frame = "icrs",
    )
    gal                 = coords_icrs.galactic
    l_rad               = gal.l.rad
    b_rad               = gal.b.rad
    return l_rad, b_rad


def build_healpix_mean_z(l_rad, b_rad, z, nside):
    npix                = hp.nside2npix(nside)
    counts_map          = np.zeros(npix, dtype = int)
    sum_z_map           = np.zeros(npix, dtype = float)

    theta               = 0.5 * np.pi - b_rad
    phi                 = l_rad

    pix_idx             = hp.ang2pix(nside, theta, phi, nest = False)

    for p, z_i in zip(pix_idx, z):
        counts_map[p]  += 1
        sum_z_map[p]   += z_i

    mean_z_map          = np.full(npix, hp.UNSEEN, dtype = float)
    good                = counts_map > 0
    mean_z_map[good]    = sum_z_map[good] / counts_map[good]

    return counts_map, mean_z_map


def fit_redshift_dipole(counts_map,
                        mean_z_map,
                        nside,
                        min_per_pixel = 10,
                        sel_map       = None,
                        sel_min_frac  = 0.3):
    npix                = len(counts_map)
    assert npix == len(mean_z_map)

    if sel_map is not None and len(sel_map) != npix:
        raise ValueError("Selection map length does not match HEALPix maps.")

    good_pix            = (counts_map >= min_per_pixel) & (mean_z_map != hp.UNSEEN)

    if sel_map is not None:
        positive        = sel_map > 0.0
        if not np.any(positive):
            raise RuntimeError("Selection map has no positive entries.")
        sel_max         = np.max(sel_map[positive])
        sel_thresh      = sel_min_frac * sel_max
        good_pix       &= sel_map >= sel_thresh

    if not np.any(good_pix):
        raise RuntimeError("No pixels pass min-per-pixel + selection cuts.")

    pix_indices         = np.where(good_pix)[0]
    weights             = counts_map[pix_indices].astype(float)
    z_vals              = mean_z_map[pix_indices]

    z_global_mean       = np.average(z_vals, weights = weights)

    theta_c, phi_c      = hp.pix2ang(nside, pix_indices, nest = False)

    n_x                 = np.sin(theta_c) * np.cos(phi_c)
    n_y                 = np.sin(theta_c) * np.sin(phi_c)
    n_z                 = np.cos(theta_c)

    delta_z             = z_vals - z_global_mean
    w_sqrt              = np.sqrt(weights)

    A                   = np.vstack((
        n_x * w_sqrt,
        n_y * w_sqrt,
        n_z * w_sqrt,
    )).T
    b                   = delta_z * w_sqrt

    D_vec, _, _, _      = np.linalg.lstsq(A, b, rcond = None)

    D_norm              = np.linalg.norm(D_vec)
    if D_norm == 0.0:
        raise RuntimeError("Dipole fit returned zero vector.")

    frac_amplitude      = D_norm / z_global_mean

    d_unit              = D_vec / D_norm
    d_x, d_y, d_z       = d_unit

    b_dip_rad           = np.arcsin(d_z)
    l_dip_rad           = np.arctan2(d_y, d_x)
    if l_dip_rad < 0.0:
        l_dip_rad      += 2.0 * np.pi

    l_dip_deg           = np.degrees(l_dip_rad)
    b_dip_deg           = np.degrees(b_dip_rad)

    return z_global_mean, D_vec, frac_amplitude, (l_dip_deg, b_dip_deg), good_pix


def build_dipole_model_map(z_global_mean, D_vec, nside):
    npix                = hp.nside2npix(nside)
    theta_all, phi_all  = hp.pix2ang(nside, np.arange(npix), nest = False)

    n_x_all             = np.sin(theta_all) * np.cos(phi_all)
    n_y_all             = np.sin(theta_all) * np.sin(phi_all)
    n_z_all             = np.cos(theta_all)

    n_dot_D             = (
        n_x_all * D_vec[0]
      + n_y_all * D_vec[1]
      + n_z_all * D_vec[2]
    )

    z_model_map         = z_global_mean + n_dot_D
    return z_model_map


# ----------------------------------------------------------------------
# Main mock-generation logic
# ----------------------------------------------------------------------

def build_parser():
    parser              = argparse.ArgumentParser(
        description = "Quaia random mocks: dipole + quadrupole stats."
    )

    parser.add_argument(
        "--data-csv",
        required    = True,
        help        = "Real Quaia G<20 CSV (RA, DEC, Z).",
    )

    parser.add_argument(
        "--random-fits",
        required    = True,
        help        = "Quaia random_G20.0_10x.fits.",
    )

    parser.add_argument(
        "--sel-map",
        required    = True,
        help        = "Selection function HEALPix map (NSIDE64 G20.0).",
    )

    parser.add_argument(
        "--nside",
        type        = int,
        default     = 64,
        help        = "HEALPix nside (default: 64).",
    )

    parser.add_argument(
        "--min-per-pixel",
        type        = int,
        default     = 10,
        help        = "Min #objects per pixel for dipole fit (default: 10).",
    )

    parser.add_argument(
        "--sel-min-frac",
        type        = float,
        default     = 0.3,
        help        = "Selection cut: sel >= frac * max(sel) (default: 0.3).",
    )

    parser.add_argument(
        "--nmocks",
        type        = int,
        default     = 50,
        help        = "Number of mock catalogues (default: 50).",
    )

    parser.add_argument(
        "--lmax",
        type        = int,
        default     = 10,
        help        = "lmax for residual C_l (default: 10).",
    )

    parser.add_argument(
        "--seed",
        type        = int,
        default     = 1234,
        help        = "Random seed (default: 1234).",
    )

    parser.add_argument(
        "--outcsv",
        default     = "quaia_random_mocks_stats.csv",
        help        = "Output CSV with mock statistics.",
    )

    parser.add_argument(
        "--bmin-abs",
        type        = float,
        default     = 0.0,
        help        = "Minimum |b| in degrees for analysis mask (default: 0).",
    )

    return parser


def load_real_catalog(csv_path):
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Real data CSV not found: {csv_path}")

    print(f"[real] load {csv_path}")
    df                  = pd.read_csv(csv_path)

    for col in ("RA", "DEC", "Z"):
        if col not in df.columns:
            raise KeyError(f"Real CSV must contain column '{col}'")

    return df


def load_random_catalog(fits_path):
    if not os.path.exists(fits_path):
        raise FileNotFoundError(f"Random FITS not found: {fits_path}")

    print(f"[rand] load {fits_path}")
    tbl                 = Table.read(fits_path, format = "fits")
    df                  = tbl.to_pandas()

    cols                = {c.lower(): c for c in df.columns}
    for need in ("ra", "dec"):
        if need not in cols:
            raise KeyError(
                f"Random catalog must have RA/DEC column; "
                f"could not find case-insensitive '{need}'. "
                f"Columns are: {list(df.columns)!r}"
            )

    ra_col              = cols["ra"]
    dec_col             = cols["dec"]

    ra_deg              = df[ra_col].values
    dec_deg             = df[dec_col].values

    return ra_deg, dec_deg


def run_mock_suite(args):
    # Load real catalog
    df_real             = load_real_catalog(args.data_csv)
    N_real              = len(df_real)
    z_real              = df_real["Z"].values

    print(f"[real] N = {N_real:d},  z in [{z_real.min():.3f}, {z_real.max():.3f}]")

    # Load random positions
    ra_rand, dec_rand   = load_random_catalog(args.random_fits)
    N_rand              = len(ra_rand)
    print(f"[rand] N_random = {N_rand:d}")

    if N_rand < N_real:
        print(
            "[warn] random catalog has fewer objects than real sample — "
            "mock sampling will use replacement."
        )

    # Selection map
    print(f"[sel] load {args.sel_map}")
    sel_map             = hp.read_map(args.sel_map)
    nside_sel           = hp.get_nside(sel_map)
    if nside_sel != args.nside:
        raise ValueError(
            f"Selection map nside={nside_sel} does not match requested "
            f"nside={args.nside}."
        )

    # Prepare output storage
    rows                = []

    rng                 = np.random.default_rng(args.seed)

    print(f"[mocks] starting {args.nmocks:d} mocks ...")

    for i in range(args.nmocks):
        print(f"[mock {i+1:03d}/{args.nmocks:d}]")

        # --- build mock catalog ---
        if N_rand >= N_real:
            idx_pos     = rng.choice(N_rand, size = N_real, replace = False)
        else:
            idx_pos     = rng.choice(N_rand, size = N_real, replace = True)

        ra_mock         = ra_rand[idx_pos]
        dec_mock        = dec_rand[idx_pos]

        # resample redshifts from real distribution (with replacement)
        z_mock          = rng.choice(z_real, size = N_real, replace = True)

        # mean-z map
        l_rad, b_rad    = icrs_to_galactic(ra_mock, dec_mock)
        counts_map, mean_z_map = build_healpix_mean_z(l_rad, b_rad, z_mock, args.nside)

        # dipole fit
        z_mean, D_vec, frac_amp, (l_dip, b_dip), fit_mask = fit_redshift_dipole(
            counts_map      = counts_map,
            mean_z_map      = mean_z_map,
            nside           = args.nside,
            min_per_pixel   = args.min_per_pixel,
            sel_map         = sel_map,
            sel_min_frac    = args.sel_min_frac,
        )

        D_norm          = np.linalg.norm(D_vec)

        # dipole model + residuals
        z_model_map     = build_dipole_model_map(z_mean, D_vec, args.nside)

        residual_map    = np.full_like(mean_z_map, hp.UNSEEN, dtype = float)
        good_pix        = (mean_z_map != hp.UNSEEN)
        residual_map[good_pix] = mean_z_map[good_pix] - z_model_map[good_pix]

        # residual C_l
        npix            = hp.nside2npix(args.nside)
        mask            = np.zeros(npix, dtype = float)
        mask[fit_mask]  = 1.0
        
        # add |b| cut if requested
        if args.bmin_abs > 0.0:
            theta_all, phi_all = hp.pix2ang(args.nside,
                                            np.arange(npix),
                                            nest = False)
            b_all_rad          = 0.5 * np.pi - theta_all
            b_all_deg          = np.rad2deg(b_all_rad)
            mask[np.abs(b_all_deg) < args.bmin_abs] = 0.0
        
        res_filled      = np.zeros(npix, dtype = float)
        valid_res       = residual_map != hp.UNSEEN
        res_filled[valid_res] = residual_map[valid_res]
        
        res_masked      = res_filled * mask
        
        cl              = hp.anafast(res_masked, lmax = args.lmax)
        ell_arr         = np.arange(len(cl), dtype = float)
        Dl              = ell_arr * (ell_arr + 1.0) * cl / (2.0 * np.pi)


        C2              = cl[2] if len(cl) > 2 else np.nan
        D2              = Dl[2] if len(Dl) > 2 else np.nan

        rows.append({
            "mock_id"        : i,
            "z_global_mean"  : z_mean,
            "dipole_amp"     : D_norm,
            "dipole_frac"    : frac_amp,
            "l_dip_deg"      : l_dip,
            "b_dip_deg"      : b_dip,
            "C2_residual"    : C2,
            "D2_residual"    : D2,
        })

    df_out              = pd.DataFrame(rows)
    print(f"[save] mock stats -> {args.outcsv}")
    df_out.to_csv(args.outcsv, index = False)

    # summary statistics
    print("\n[summary] over mocks:")
    for col in ("dipole_frac", "D2_residual"):
        vals            = df_out[col].values
        mu              = np.mean(vals)
        sig             = np.std(vals)
        print(f"  {col:15s}: mean = {mu:.6e}, std = {sig:.6e}")

    print("[done]")


def main(argv = None):
    if argv is None:
        argv            = sys.argv[1:]

    parser              = build_parser()
    args                = parser.parse_args(argv)

    run_mock_suite(args)


if __name__ == "__main__":
    main()
