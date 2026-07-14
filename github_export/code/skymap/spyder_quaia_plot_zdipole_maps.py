#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 23:14:19 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_plot_zdipole_maps.py  (0V1)

Quick helper to make sky maps from the Quaia mean-z, dipole-model,
and residual HEALPix maps.

Example:

    python quaia_plot_zdipole_maps.py \
        --zmean "/Users/boyde/Downloads/quaia_G20p0_nside64_zmean.fits" \
        --zdip  "/Users/boyde/Downloads/quaia_G20p0_nside64_zdipole_model.fits" \
        --zres  "/Users/boyde/Downloads/quaia_G20p0_nside64_zresidual.fits" \
        --outdir "/Users/boyde/Downloads" \
        --prefix "quaia_G20p0_nside64_zdipole"

This will write:
    quaia_G20p0_nside64_zdipole_zmean.png
    quaia_G20p0_nside64_zdipole_zdipole_model.png
    quaia_G20p0_nside64_zdipole_zresidual.png
"""

import argparse
from pathlib import Path

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np


def plot_map(infile, outfile, title, unit=None, vmin=None, vmax=None):
    """Generic Mollweide plot helper."""
    print(f"[plot] {infile} -> {outfile}")
    m = hp.read_map(infile, verbose=False)

    plt.figure(figsize=(8, 5))
    hp.mollview(
        m,
        title=title,
        unit=unit,
        min=vmin,
        max=vmax,
        norm="hist",
    )
    hp.graticule()
    plt.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close()


def main():
    p = argparse.ArgumentParser(description="Plot Quaia z-dipole HEALPix maps")
    p.add_argument("--zmean", type=str, default=None,
                   help="mean-z map FITS (e.g. ..._zmean.fits)")
    p.add_argument("--zdip", type=str, default=None,
                   help="dipole-model map FITS (e.g. ..._zdipole_model.fits)")
    p.add_argument("--zres", type=str, default=None,
                   help="residual map FITS (e.g. ..._zresidual.fits)")
    p.add_argument("--outdir", type=str, default=".",
                   help="output directory for PNGs (default: .)")
    p.add_argument("--prefix", type=str, default="quaia_zdipole",
                   help="prefix for output PNG filenames")
    args = p.parse_args()

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Mean-z map
    if args.zmean is not None:
        zfile = Path(args.zmean)
        m = hp.read_map(zfile, verbose=False)
        vmin = float(np.nanpercentile(m, 5.0))
        vmax = float(np.nanpercentile(m, 95.0))
        outfile = outdir / f"{args.prefix}_zmean.png"
        plot_map(
            zfile,
            outfile,
            title="Quaia mean z (z ≥ 1.5, G ≤ 20)",
            unit="z",
            vmin=vmin,
            vmax=vmax,
        )

    # Dipole model map
    if args.zdip is not None:
        dfile = Path(args.zdip)
        m = hp.read_map(dfile, verbose=False)
        vmax_abs = float(np.nanpercentile(np.abs(m), 95.0))
        outfile = outdir / f"{args.prefix}_zdipole_model.png"
        plot_map(
            dfile,
            outfile,
            title="Quaia redshift dipole model (l=137°, b=21°)",
            unit="Δz (model)",
            vmin=-vmax_abs,
            vmax=vmax_abs,
        )

    # Residual map
    if args.zres is not None:
        rfile = Path(args.zres)
        m = hp.read_map(rfile, verbose=False)
        vmax_abs = float(np.nanpercentile(np.abs(m), 95.0))
        outfile = outdir / f"{args.prefix}_zresidual.png"
        plot_map(
            rfile,
            outfile,
            title="Quaia redshift residuals after dipole subtraction",
            unit="Δz (residual)",
            vmin=-vmax_abs,
            vmax=vmax_abs,
        )

    print("[done] plots written to", outdir)


if __name__ == "__main__":
    main()
