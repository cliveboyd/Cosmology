#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 17:25:14 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
skymap2_dipole_checks.py  (0V1)

Quick consistency checks for Quaia / skymap2 pipeline:

- Load catalogue (e.g. Quaia)
- Take high-z subset (z >= zmin)
- Compute simple "mean-direction" dipole from:
    (a) RA/Dec (equatorial)
    (b) RA/Dec -> Galactic (via astropy)
    (c) catalogue Galactic columns, if present
    (d) "swap test": treat RA/Dec as if they were (l, b)

For each, print:
- Dipole direction in RA/Dec and (l, b)
- Angle to Galactic Centre (GC)
- Angle to CMB dipole direction

Usage examples:
    python skymap2_dipole_checks.py \
        --fits /path/to/quaia.fits \
        --zmin 1.5 \
        --maxN 200000

You can restrict to a single mode:
    --mode equatorial
    --mode gal_from_ra
    --mode gal_from_cat
    --mode swap
"""

import argparse
import numpy as np

from   pathlib             import Path

from   astropy.io          import fits
from   astropy.coordinates import SkyCoord
import astropy.units as u


def angle_deg(v1, v2):
    v1 = np.asarray(v1, dtype=float)
    v2 = np.asarray(v2, dtype=float)
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    dot = np.clip(np.dot(v1, v2), -1.0, 1.0)
    return np.degrees(np.arccos(dot))

def lonlat_to_vec(lon_deg, lat_deg):
    lon  = np.radians(lon_deg)
    lat  = np.radians(lat_deg)
    cosb = np.cos(lat)
    x    = cosb * np.cos(lon)
    y    = cosb * np.sin(lon)
    z    = np.sin(lat)
    return np.vstack([x, y, z]).T

def vec_to_lonlat(vec):
    v       = np.asarray(vec, dtype=float)
    v      /= np.linalg.norm(v)
    x, y, z = v
    lon     = np.degrees(np.arctan2(y, x)) % 360.0
    lat     = np.degrees(np.arcsin(z))
    return lon, lat

# ----------------------------------------------------------------------
# Bootstrap angular scatter of the dipole direction
# Returns percentiles of angle between bootstrap dipoles and full-sample dipole
# ----------------------------------------------------------------------
def bootstrap_dipole(lon_deg, lat_deg, w=None, nboot=500, seed=12345):
    lon_deg = np.asarray(lon_deg, dtype=float)
    lat_deg = np.asarray(lat_deg, dtype=float)
    N       = lon_deg.size
    if N == 0:
        raise RuntimeError("No objects to bootstrap.")

    # Dipole from full sample (reference direction)
    base_vec = estimate_dipole(lon_deg, lat_deg, w=w)

    rng          = np.random.default_rng(seed)
    boot_angles  = []

    for _ in range(nboot):
        idx   = rng.integers(0, N, size=N)
        lon_b = lon_deg[idx]
        lat_b = lat_deg[idx]
        if w is not None:
            w_b = np.asarray(w)[idx]
        else:
            w_b = None

        vec_b = estimate_dipole(lon_b, lat_b, w=w_b)
        # angle between bootstrap dipole and full-sample dipole
        ang_b = angle_deg(vec_b, base_vec)
        boot_angles.append(ang_b)

    boot_angles = np.asarray(boot_angles)
    p16, p50, p84 = np.percentile(boot_angles, [16, 50, 84])
    return p16, p50, p84


# ----------------------------------------------------------------------
# Simple "dipole" estimator: mean unit vector of positions
# Optional weights w (same length as lon/lat)
# ----------------------------------------------------------------------
def estimate_dipole(lon_deg, lat_deg, w=None):
    xyz = lonlat_to_vec(lon_deg, lat_deg)   # (N, 3)

    if w is None:
        mean_vec = xyz.mean(axis=0)
    else:
        w = np.asarray(w, dtype=float)
        if w.shape[0] != xyz.shape[0]:
            raise ValueError("weights w must have same length as lon/lat.")
        wsum = w.sum()
        if wsum <= 0:
            raise ValueError("sum of weights must be > 0.")
        mean_vec = np.sum(xyz * (w[:, None] / wsum), axis=0)

    norm = np.linalg.norm(mean_vec)
    if norm == 0:
        raise RuntimeError("Mean vector has zero norm; sample is too symmetric or empty.")
    mean_vec /= norm
    return mean_vec


# ----------------------------------------------------------------------
# Load catalogue + columns
# ----------------------------------------------------------------------
def load_catalogue(path, zmin=1.5, maxN=None):
    path    = Path(path)
    with fits.open(path, memmap=True) as hdul:
        # guess table HDU: first BinTable after primary
        hdu = None
        for h in hdul:
            if isinstance(h, fits.BinTableHDU):
                hdu = h
                break
        if hdu is None:
            raise RuntimeError("No BinTableHDU found in FITS file.")

        data        = hdu.data
        colnames    = [c.upper() for c in hdu.columns.names]

    def get_col(candidates, required=False):
        for cand in candidates:
            if cand.upper() in colnames:
                return data[cand]
        if required:
            raise KeyError(f"None of columns {candidates} found in catalogue.")
        return None

    # RA/Dec
    ra          = get_col(["RA", "RA_ICRS", "ALPHA_J2000"], required=True)
    dec         = get_col(["DEC", "DE_ICRS", "DELTA_J2000"], required=True)

    # Redshift
    # Redshift
    # Try a broader set of possibilities to match Quaia-style catalogues:
        # Redshift
    z = get_col(
        [
            "REDSHIFT_QUAIA",  # <- Quaia-specific column
            "Z",
            "Z_QSO",
            "Z_MEAN",
            "ZMEAN",
            "Z_MED",
            "ZMED",
            "REDSHIFT",
            "Z_COSMO",
            "Z_PHOTO",
        ],
        required=False,
    )

    if z is None:
        raise KeyError(
            "No redshift column found. Available columns are:\n"
            f"{colnames}\n"
            "Please add the correct redshift column name to the list in `load_catalogue`."
        )



    # Galactic from catalogue (if any)
    glon_cat    = get_col(["GLON", "GAL_L", "L_GAL", "L"], required=False)
    glat_cat    = get_col(["GLAT", "GAL_B", "B_GAL", "B"], required=False)

    # High-z mask
    mask        = z >= zmin
    ra          = ra[mask]
    dec         = dec[mask]
    z           = z[mask]

    if glon_cat is not None and glat_cat is not None:
        glon_cat = glon_cat[mask]
        glat_cat = glat_cat[mask]

    # Optional thinning for speed
    if maxN is not None and len(ra) > maxN:
        idx     = np.random.default_rng(12345).choice(len(ra), size=maxN, replace=False)
        ra      = ra[idx]
        dec     = dec[idx]
        z       = z[idx]
        if glon_cat is not None and glat_cat is not None:
            glon_cat = glon_cat[idx]
            glat_cat = glat_cat[idx]

    return {
        "ra"        : np.asarray(ra, float),
        "dec"       : np.asarray(dec, float),
        "z"         : np.asarray(z, float),
        "glon_cat"  : None if glon_cat is None else np.asarray(glon_cat, float),
        "glat_cat"  : None if glat_cat is None else np.asarray(glat_cat, float),
    }

def reference_directions():
    # Galactic Centre: l=0, b=0 by definition
    gc_gal      = lonlat_to_vec(0.0, 0.0)

    # GC in equatorial (approx ICRS)
    gc_eq_coord = SkyCoord(l=0*u.deg, b=0*u.deg, frame="galactic").icrs
    gc_eq       = lonlat_to_vec(gc_eq_coord.ra.deg, gc_eq_coord.dec.deg)

    # CMB dipole direction (Planck-ish; close enough for angles)
    cmb_eq      = lonlat_to_vec(167.942, -6.944)
    cmb_gal     = lonlat_to_vec(264.021, 48.253)

    return {
        "gc_eq"   : gc_eq[0],
        "gc_gal"  : gc_gal[0],
        "cmb_eq"  : cmb_eq[0],
        "cmb_gal" : cmb_gal[0],
    }


# ----------------------------------------------------------------------
# Reference directions: GC & CMB dipole
# (numbers are standard; tweak if you prefer other definitions)
# ----------------------------------------------------------------------
def run_mode(mode, cat, zweight="none", nboot=0):
    refs = reference_directions()

    # Choose weights for the dipole estimator
    if zweight == "z":
        w = cat["z"]
    elif zweight == "z2":
        w = cat["z"]**2
    else:
        w = None

    # Build SkyCoord from RA/Dec
    c_eq          = SkyCoord(ra=cat["ra"]*u.deg, dec=cat["dec"]*u.deg, frame="icrs")
    c_gal_from_ra = c_eq.galactic

    results = []

    # --- Equatorial dipole (RA, Dec) ---
    if mode in ("equatorial", "all"):
        vec_eq          = estimate_dipole(cat["ra"], cat["dec"], w=w)
        ra_dip, dec_dip = vec_to_lonlat(vec_eq)
        # convert to galactic for reporting
        dip_coord       = SkyCoord(ra=ra_dip*u.deg, dec=dec_dip*u.deg, frame="icrs").galactic
        l_dip, b_dip    = dip_coord.l.deg, dip_coord.b.deg

        boot = None
        if nboot > 0:
            boot = bootstrap_dipole(cat["ra"], cat["dec"], w=w, nboot=nboot)

        results.append(("equatorial", vec_eq, ra_dip, dec_dip, l_dip, b_dip, boot))

    # --- Galactic dipole via RA/Dec -> Galactic transform ---
    if mode in ("gal_from_ra", "all"):
        l = c_gal_from_ra.l.deg
        b = c_gal_from_ra.b.deg

        vec_gal       = estimate_dipole(l, b, w=w)
        l_dip, b_dip  = vec_to_lonlat(vec_gal)
        dip_coord     = SkyCoord(l=l_dip*u.deg, b=b_dip*u.deg, frame="galactic").icrs
        ra_dip, dec_dip = dip_coord.ra.deg, dip_coord.dec.deg

        boot = None
        if nboot > 0:
            boot = bootstrap_dipole(l, b, w=w, nboot=nboot)

        results.append(("gal_from_ra", vec_gal, ra_dip, dec_dip, l_dip, b_dip, boot))

    # --- Galactic dipole directly from catalogue l,b (if present) ---
    if mode in ("gal_from_cat", "all") and cat["glon_cat"] is not None:
        vec_gal_c      = estimate_dipole(cat["glon_cat"], cat["glat_cat"], w=w)
        l_dip, b_dip   = vec_to_lonlat(vec_gal_c)
        dip_coord      = SkyCoord(l=l_dip*u.deg, b=b_dip*u.deg, frame="galactic").icrs
        ra_dip, dec_dip = dip_coord.ra.deg, dip_coord.dec.deg

        boot = None
        if nboot > 0:
            boot = bootstrap_dipole(cat["glon_cat"], cat["glat_cat"], w=w, nboot=nboot)

        results.append(("gal_from_cat", vec_gal_c, ra_dip, dec_dip, l_dip, b_dip, boot))

    # --- Swap test: treat RA,Dec as if they were (l,b) directly ---
    if mode in ("swap", "all"):
        vec_swap      = estimate_dipole(cat["ra"], cat["dec"], w=w)
        l_dip, b_dip  = vec_to_lonlat(vec_swap)
        dip_coord     = SkyCoord(l=l_dip*u.deg, b=b_dip*u.deg, frame="galactic").icrs
        ra_dip, dec_dip = dip_coord.ra.deg, dip_coord.dec.deg

        boot = None
        if nboot > 0:
            boot = bootstrap_dipole(cat["ra"], cat["dec"], w=w, nboot=nboot)

        results.append(("swap_RA_as_gal", vec_swap, ra_dip, dec_dip, l_dip, b_dip, boot))

    # Print summary
    print("\n[results] high-z dipole estimates")
    print("  (angles to GC and CMB are in degrees)\n")
    for name, vec, ra_dip, dec_dip, l_dip, b_dip, boot in results:
        ang_gc_eq   = angle_deg(vec, refs["gc_eq"])
        ang_gc_gal  = angle_deg(vec, refs["gc_gal"])
        ang_cmb_eq  = angle_deg(vec, refs["cmb_eq"])
        ang_cmb_gal = angle_deg(vec, refs["cmb_gal"])

        print(f"Mode: {name}")
        print(f"  dipole (equatorial): RA={ra_dip:7.3f} deg, Dec={dec_dip:7.3f} deg")
        print(f"  dipole (galactic):   l ={l_dip:7.3f} deg, b ={b_dip:7.3f} deg")
        print(f"  angle to GC   (eq) : {ang_gc_eq:7.3f}")
        print(f"  angle to GC   (gal): {ang_gc_gal:7.3f}")
        print(f"  angle to CMB  (eq) : {ang_cmb_eq:7.3f}")
        print(f"  angle to CMB  (gal): {ang_cmb_gal:7.3f}")

        if boot is not None:
            p16, p50, p84 = boot
            print(f"  bootstrap Δθ (deg) between bootstrap dipoles and full-sample dipole")
            print(f"    median = {p50:6.3f},  68% CI = [{p16:6.3f}, {p84:6.3f}]  (nboot={nboot})")

        print("")


def main():
    parser = argparse.ArgumentParser(description="skymap2 dipole consistency checks")
    parser.add_argument("--fits", type=str, required=True,
                        help="Path to Quaia (or other) FITS catalogue")
    parser.add_argument("--zmin", type=float, default=1.5,
                        help="Minimum redshift for high-z subset (default: 1.5)")
    parser.add_argument("--maxN", type=int, default=200000,
                        help="Max number of objects to use (None for all)")
    parser.add_argument("--mode", type=str, default="all",
                        choices=["all", "equatorial", "gal_from_ra", "gal_from_cat", "swap"],
                        help="Which mode to run (default: all)")
    parser.add_argument("--zweight", type=str, default="none",
                        choices=["none", "z", "z2"],
                        help="Weighting scheme for dipole: none, z, or z^2 (default: none)")
    parser.add_argument("--nboot", type=int, default=0,
                        help="Number of bootstrap resamples for angular error (0 = off)")
    args = parser.parse_args()

    cat = load_catalogue(args.fits, zmin=args.zmin, maxN=args.maxN)
    print(f"[load] selected N={len(cat['ra'])} objects with z>={args.zmin}")

    run_mode(args.mode, cat, zweight=args.zweight, nboot=args.nboot)



if __name__ == "__main__":
    main()
