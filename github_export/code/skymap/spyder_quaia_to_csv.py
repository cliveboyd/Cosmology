#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 16 22:25:55 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_to_csv.py — Extract basic RA/DEC/Z table from Quaia FITS (0V1)

Example
-------
    python quaia_to_csv.py \
        --fits  /Users/boyde/.spyder-py3/quaia_G20.0.fits \
        --out   /Users/boyde/Downloads/quaia_G20p0_basic.csv
"""

import argparse
import os
import sys

import numpy  as np
import pandas as pd

from astropy.table import Table


def find_column(df, candidates):
    """
    Given a DataFrame and a list of candidate column names (strings),
    return the first one found (case-insensitive). Raises if none found.
    """
    cols_lower    = {c.lower(): c for c in df.columns}

    for cand in candidates:
        cl       = cand.lower()
        if cl in cols_lower:
            return cols_lower[cl]

    raise KeyError(
        f"Could not find any of the candidate columns {candidates!r} "
        f"in DataFrame with columns: {list(df.columns)!r}"
    )


def load_quaia_basic(fits_path, verbose = True):
    """
    Load Quaia FITS and return a DataFrame with three key columns:
        RA, DEC, Z

    This function tries to be robust to different column name variants.
    """
    if verbose:
        print(f"[load] {fits_path}")

    table          = Table.read(fits_path, format = "fits")
    df             = table.to_pandas()

    # Standardise column names to strings
    df.columns     = [str(c) for c in df.columns]

    # Try to find RA/DEC
    ra_candidates  = ["RA", "ra", "RA_ICRS", "ra_icrs"]
    dec_candidates = ["DEC", "dec", "DEC_ICRS", "dec_icrs"]

    ra_col         = find_column(df, ra_candidates)
    dec_col        = find_column(df, dec_candidates)

    # Try to find a redshift column.
    # Quaia uses an improved redshift estimate; likely candidates:
    #   "Z", "z", "Z_PHOT", "z_phot", "Z_QSO", etc.
    # Try to find a redshift column.
    # For Quaia this is typically "redshift_quaia".
    z_candidates   = [
        "redshift_quaia",
        "redshift",
        "Z", "z",
        "Z_PHOT", "z_phot",
        "Z_QSO", "z_qso",
        "Z_KNN", "z_knn",
    ]


    z_col          = find_column(df, z_candidates)

    if verbose:
        print(f"[info] using RA column : {ra_col}")
        print(f"[info] using DEC column: {dec_col}")
        print(f"[info] using Z column  : {z_col}")

    df_basic       = pd.DataFrame({
        "RA":  df[ra_col].values,
        "DEC": df[dec_col].values,
        "Z":   df[z_col].values,
    })

    return df_basic


def build_parser():
    parser         = argparse.ArgumentParser(
        description = "Convert Quaia FITS (quaia_G20.X.fits) to a basic RA/DEC/Z CSV.",
    )

    parser.add_argument(
        "--fits",
        required  = True,
        help      = "Path to quaia_G20.X.fits",
    )

    parser.add_argument(
        "--out",
        default   = None,
        help      = "Output CSV file (default: same base name with _basic.csv).",
    )

    parser.add_argument(
        "--maxrows",
        type      = int,
        default   = None,
        help      = "Optional: limit to first N rows (for quick tests).",
    )

    return parser


def main(argv = None):
    if argv is None:
        argv      = sys.argv[1:]

    parser        = build_parser()
    args          = parser.parse_args(argv)

    fits_path     = args.fits
    if not os.path.exists(fits_path):
        raise FileNotFoundError(f"FITS file not found: {fits_path}")

    if args.out is None:
        base, ext  = os.path.splitext(fits_path)
        out_path   = base + "_basic.csv"
    else:
        out_path   = args.out

    df_basic      = load_quaia_basic(fits_path, verbose = True)

    if args.maxrows is not None:
        df_basic  = df_basic.iloc[:args.maxrows].copy()
        print(f"[info] limiting to first {len(df_basic):d} rows for output")

    print(f"[info] rows in table: {len(df_basic):d}")
    print(f"[save] {out_path}")
    df_basic.to_csv(out_path, index = False)
    print("[done]")


if __name__ == "__main__":
    main()
