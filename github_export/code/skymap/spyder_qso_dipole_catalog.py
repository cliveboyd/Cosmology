#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qso_dipole_catalog.py — Dipole direction from eBOSS QSO catalog (0V1)

Takes a CSV like eBOSS_QSO_quasar_table_z0p1_3p0.csv and:

    * applies optional redshift cuts (zmin, zmax),
    * converts (RA, DEC) [ICRS] -> (l, b) [Galactic],
    * builds a dipole vector from unit position vectors,
    * reports the dipole amplitude and direction in Galactic coords,
    * optionally builds a HEALPix number-count map for skymap2.

Requires:
    numpy
    pandas
    astropy
    healpy     (only if --nside is provided for map output)
"""

import argparse
import sys
import os

import numpy  as np
import pandas as pd

from astropy.coordinates import SkyCoord
import astropy.units      as u

try:
    import healpy as hp
    HAVE_HEALPY  = True
except ImportError:
    HAVE_HEALPY  = False


def load_qso_catalog(csv_path, zmin = None, zmax = None):
    """
    Load quasar CSV and apply optional redshift cuts.
    Expects columns: RA [deg], DEC [deg], Z.
    """
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    df           = pd.read_csv(csv_path)

    if "Z" not in df.columns:
        raise RuntimeError("Input CSV must contain column 'Z' for redshift.")

    mask         = np.ones(len(df), dtype = bool)

    if zmin is not None:
        mask    &= df["Z"].values >= zmin

    if zmax is not None:
        mask    &= df["Z"].values <= zmax

    df_sel       = df.loc[mask].reset_index(drop = True)

    return df_sel


def equatorial_to_galactic(ra_deg, dec_deg):
    """
    Convert (RA, DEC) in degrees (ICRS) to Galactic (l, b) in degrees.
    """
    coords_icrs  = SkyCoord(ra = ra_deg * u.deg, dec = dec_deg * u.deg, frame = "icrs")
    gal          = coords_icrs.galactic

    l_deg        = gal.l.deg
    b_deg        = gal.b.deg

    return l_deg, b_deg


def compute_dipole_direction(l_deg, b_deg):
    """
    Compute dipole vector from Galactic coordinates (l, b) in degrees.

    We approximate the number-density dipole by summing unit vectors:

        D = (1/N) * sum_i r_i

    where r_i is the unit vector pointing toward each source.
    The dipole direction is then the direction of D.
    """
    l_rad        = np.deg2rad(l_deg)
    b_rad        = np.deg2rad(b_deg)

    cos_b        = np.cos(b_rad)
    sin_b        = np.sin(b_rad)

    x            = cos_b * np.cos(l_rad)
    y            = cos_b * np.sin(l_rad)
    z            = sin_b

    # Sum over all sources
    Dx           = np.sum(x)
    Dy           = np.sum(y)
    Dz           = np.sum(z)

    N            = len(l_deg)
    D_norm       = np.sqrt(Dx * Dx + Dy * Dy + Dz * Dz)

    if D_norm == 0.0:
        raise RuntimeError("Dipole vector has zero norm; something is wrong with the input.")

    # Mean unit vector
    Dx_hat       = Dx / D_norm
    Dy_hat       = Dy / D_norm
    Dz_hat       = Dz / D_norm

    # Convert back to Galactic (l, b)
    b_dip_rad    = np.arcsin(Dz_hat)
    l_dip_rad    = np.arctan2(Dy_hat, Dx_hat)

    l_dip_deg    = np.rad2deg(l_dip_rad) % 360.0
    b_dip_deg    = np.rad2deg(b_dip_rad)

    # Define a simple dimensionless amplitude
    dip_amp      = D_norm / float(N)

    result       = {
        "Dx"      : Dx,
        "Dy"      : Dy,
        "Dz"      : Dz,
        "N"       : N,
        "amp"     : dip_amp,
        "l_deg"   : l_dip_deg,
        "b_deg"   : b_dip_deg,
    }

    return result


def build_healpix_map(l_deg, b_deg, nside):
    """
    Build a simple HEALPix number-count map from Galactic (l, b) in degrees.
    """
    if not HAVE_HEALPY:
        raise RuntimeError("healpy is not available; cannot build HEALPix map.")

    theta        = np.deg2rad(90.0 - b_deg)   # colatitude
    phi          = np.deg2rad(l_deg)          # longitude

    npix         = hp.nside2npix(nside)
    pix          = hp.ang2pix(nside, theta, phi)

    counts       = np.bincount(pix, minlength = npix).astype(float)

    return counts


def build_parser():
    parser       = argparse.ArgumentParser(
        description = "Compute quasar dipole direction from eBOSS QSO CSV."
    )

    parser.add_argument(
        "--csv",
        required = True,
        help     = "Path to quasar CSV file (e.g. eBOSS_QSO_quasar_table_z0p1_3p0.csv).",
    )

    parser.add_argument(
        "--zmin",
        type     = float,
        default  = None,
        help     = "Minimum redshift cut (optional, overrides CSV-range).",
    )

    parser.add_argument(
        "--zmax",
        type     = float,
        default  = None,
        help     = "Maximum redshift cut (optional, overrides CSV-range).",
    )

    parser.add_argument(
        "--nside",
        type     = int,
        default  = None,
        help     = "If provided, build a HEALPix number-count map at this nside.",
    )

    parser.add_argument(
        "--map-out",
        default  = None,
        help     = "Output path for HEALPix map FITS (optional). "
                   "If omitted but --nside is given, a default name is used.",
    )

    return parser


def main(argv = None):
    if argv is None:
        argv     = sys.argv[1:]

    parser       = build_parser()
    args         = parser.parse_args(argv)

    print(f"[load] {args.csv}")
    df           = load_qso_catalog(args.csv, zmin = args.zmin, zmax = args.zmax)

    print(f"[info] total objects used for dipole: {len(df):d}")

    # Convert to Galactic
    l_deg, b_deg = equatorial_to_galactic(df["RA"].values, df["DEC"].values)

    # Compute dipole
    dip          = compute_dipole_direction(l_deg, b_deg)

    print("[dipole]")
    print(f"  N           = {dip['N']:d}")
    print(f"  amplitude   = {dip['amp']:.6e}")
    print(f"  l_dip [deg] = {dip['l_deg']:.3f}")
    print(f"  b_dip [deg] = {dip['b_deg']:.3f}")
    print(f"  vector (Dx, Dy, Dz) = ({dip['Dx']:.3e}, {dip['Dy']:.3e}, {dip['Dz']:.3e})")

    # Optional HEALPix map
    if args.nside is not None:
        if not HAVE_HEALPY:
            raise RuntimeError("healpy is required for HEALPix map output (pip install healpy).")

        nside      = int(args.nside)
        print(f"[healpix] building number-count map at nside = {nside:d}")

        hmap       = build_healpix_map(l_deg, b_deg, nside = nside)

        if args.map_out is None:
            base    = os.path.splitext(os.path.basename(args.csv))[0]
            outmap  = f"{base}_nside{nside:d}_counts.fits"
        else:
            outmap  = args.map_out

        print(f"[save] HEALPix map -> {outmap}")
        hp.write_map(outmap, hmap, overwrite = True)

    print("[done]")


if __name__ == "__main__":
    main()
