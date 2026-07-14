#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 08:08:41 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_make_simple_randoms.py  (0V1)

Generate a simple random RA/Dec catalogue matched to the Quaia G20p0
selection map.

Example:
    python quaia_make_simple_randoms.py \
        --sel-map "/Users/boyde/Downloads/quaia_G20p0_nside64_counts.fits" \
        --nside 64 \
        --nrand 2000000 \
        --sel-min-frac 0.5 \
        --bmin-abs 20 \
        --out-fits "/Users/boyde/Downloads/quaia_G20p0_randoms_simple.fits"
"""

import argparse
import numpy      as     np
import healpy     as     hp
from   astropy.io import fits

def build_parser():
    p = argparse.ArgumentParser(description="Build simple random RA/Dec catalogue from a selection map.")
    p.add_argument("--sel-map",      type=str,   required=True,
                   help="HEALPix selection/counts map (FITS).")
    p.add_argument("--nside",        type=int,   required=True,
                   help="nside of the selection map.")
    p.add_argument("--nrand",        type=int,   required=True,
                   help="Number of random points to generate.")
    p.add_argument("--sel-min-frac", type=float, default=0.5,
                   help="Keep pixels with sel >= sel_min_frac * max(sel).")
    p.add_argument("--bmin-abs",     type=float, default=0.0,
                   help="Optional |b| cut in galactic latitude (deg).")
    p.add_argument("--out-fits",     type=str,   required=True,
                   help="Output FITS file for random catalogue.")
    return p

def main():  
    parser           = build_parser()
    args             = parser.parse_args()

    print(f"[load] selection map: {args.sel_map}")
    m                = hp.read_map(args.sel_map, verbose=False)

    # Mask on selection strength
    sel_max         = np.nanmax(m)
    thresh          = args.sel_min_frac * sel_max
    good            = m >= thresh

    print(f"[sel] max(sel) = {sel_max:.3e},  thresh = {thresh:.3e}")
    print(f"[sel] good pixels = {good.sum()} / {len(m)}")

    # Optional |b| cut in galactic coords
    if args.bmin_abs > 0.0:
        ipix_all    = np.arange(len(m))
        theta, phi  = hp.pix2ang(args.nside, ipix_all)
        
        # Galactic lat b = 90 - theta in deg
        b_deg       = np.degrees(0.5*np.pi - theta)
        bmask       = np.abs(b_deg) >= args.bmin_abs
        good       &= bmask
        print(f"[sel] |b| >= {args.bmin_abs:.1f} deg: remaining good pixels = {good.sum()}")

    ipix_good      = np.where(good)[0]
    if len(ipix_good) == 0:
        raise RuntimeError("No good pixels in selection mask with given cuts.")

    print(f"[rand] drawing {args.nrand} random points from {len(ipix_good)} good pixels")
    rng            = np.random.default_rng(12345)

    # Choose pixels uniformly among good ones
    pix_choice     = rng.choice(ipix_good, size=args.nrand, replace=True)

    # Within each pixel, throw a random point uniformly
    # (hp.pix2ang gives the pixel centre; we dither by small random offsets)
    theta_c, phi_c = hp.pix2ang(args.nside, pix_choice)
    # Isotropic dither within the pixel solid angle: approximate by small random shifts
    # of order pixel size
    omega_pix      = 4.0 * np.pi / (12 * args.nside**2)
    dtheta         = np.sqrt(omega_pix) * rng.normal(size=args.nrand) * 0.5
    dphi           = np.sqrt(omega_pix) * rng.normal(size=args.nrand) * 0.5

    theta          = np.clip(theta_c + dtheta, 0.0, np.pi)
    phi            = (phi_c + dphi) % (2.0 * np.pi)

    # Convert to RA/Dec in equatorial: here we assume (for simplicity)
    # that the selection map is already in ICRS. If it's in Galactic,
    # you would convert (theta, phi) via astropy SkyCoord instead.
    # For your current workflow (Quaia basic CSV already has RA/DEC),
    # it's perfectly fine if the randoms are in Galactic and then
    # transformed inside the analysis script. To keep it simple, we treat
    # (theta,phi) as ICRS here.
    ra_deg         = np.degrees(phi) % 360.0
    dec_deg        = 90.0 - np.degrees(theta)

    # Build FITS table
    col_ra         = fits.Column(name="RA",  format="D", array=ra_deg)
    col_dec        = fits.Column(name="DEC", format="D", array=dec_deg)
    cols           = fits.ColDefs([col_ra, col_dec])
    hdu            = fits.BinTableHDU.from_columns(cols)

    hdu.writeto(args.out_fits, overwrite=True)
    print(f"[save] random catalogue -> {args.out_fits}")

if __name__ == "__main__":
    main()
