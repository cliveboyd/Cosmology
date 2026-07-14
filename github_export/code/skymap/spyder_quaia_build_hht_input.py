#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  6 18:25:13 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_build_hht_input.py  (0V1)

Builds a 1D series vs redshift from the Quaia <z>-dipole MC summary, suitable
for HHT/CEEMDAN analysis.

It reads:
  OUT_DIR / "quaia_mc_bcut_shuffle_summary.txt"

and for a chosen |b|>bcut it:
  * selects bins in a given z-range,
  * applies simple quality cuts (N_good, p-values),
  * sorts by z_mid,
  * saves:
      - a compact 2-column file: z_mid, signal(z)
      - a richer multi-column file with all diagnostics.

Usage (from Spyder/IPython):
        
        
%run /Users/boyde/.spyder-py3/quaia_build_hht_input.py \
    --bcut 25 \
    --metric Apar \
    --z-min 0.30 \
    --z-max 2.30 \
    --N-min 20 \
    --pA-max 1.0 \
    --pApar-max 1.0
    
"""

import argparse
from   pathlib      import Path

import numpy        as     np

from   quaia_config import OUT_DIR


MC_SUMMARY_NAME = "quaia_mc_bcut_shuffle_summary.txt"


def load_mc_summary(path):
    """
    Load the MC summary table.

    Expected columns (whitespace separated):
      bcut  tag  z_mid  N_good  A_obs  A_par_obs  f_par_obs  p_A  p_|A_par|  p_|f_par|
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"MC summary file not found: {path}")

    data = np.genfromtxt(
                         path,
                         dtype    = None,
                         encoding = None,
                         comments = "#",
                         names=[
                                "bcut",
                                "tag",
                                "z_mid",
                                "N_good",
                                "A_obs",
                                "A_par_obs",
                                "f_par_obs",
                                "p_A",
                                "p_Apar",
                                "p_fpar",
                               ],
                         )
    return data

def build_series(
                 data,
                 bcut,
                 z_min,
                 z_max,
                 N_min,
                 pA_max,
                 pApar_max,
                ):
    """
    Filter the MC summary to a single bcut and redshift range, with quality cuts.

    Returns a dict with arrays and some metadata.
    """
    # Basic selection
    mask = (
            (np.isclose(data["bcut"], bcut))
            & (data["z_mid"] >= z_min)
            & (data["z_mid"] <= z_max)
            & (data["N_good"] >= N_min)
           )

    # Optional p-value cuts (can be liberal initially)
    mask &= (data["p_A"] <= pA_max) & (data["p_Apar"] <= pApar_max)

    sel = data[mask]
    if sel.size == 0:
        raise RuntimeError(
                           f"No bins left after cuts for bcut={bcut}, "
                           f"z∈[{z_min},{z_max}], N_good>={N_min}"
                          )

    # Sort by redshift
    order = np.argsort(sel["z_mid"])
    sel = sel[order]

    return {
            "z_mid": sel["z_mid"],
            "A": sel["A_obs"],
            "A_par": sel["A_par_obs"],
            "f_par": sel["f_par_obs"],
            "N_good": sel["N_good"],
            "p_A": sel["p_A"],
            "p_Apar": sel["p_Apar"],
            "p_fpar": sel["p_fpar"],
            "tag": sel["tag"],
           }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bcut", type=float, default=25.0,
                        help="Galactic latitude cut: use |b| > bcut (deg)")
    parser.add_argument("--z-min", type=float, default=0.30,
                        help="Minimum z_mid to include")
    parser.add_argument("--z-max", type=float, default=2.30,
                        help="Maximum z_mid to include")
    parser.add_argument("--N-min", type=int, default=20,
                        help="Minimum N_good per pixel slice")
    parser.add_argument("--pA-max", type=float, default=1.0,
                        help="Maximum allowed p(A) for inclusion")
    parser.add_argument("--pApar-max", type=float, default=1.0,
                        help="Maximum allowed p(|A_par|) for inclusion")
    parser.add_argument(
                        "--metric",
                        type=str,
                        default="Apar",
                        choices=["A", "Apar", "fpar"],
                        help="Which signal to output for HHT: A, Apar (=A_par), or fpar (=f_par)",
                       )
    args = parser.parse_args()

    mc_path = OUT_DIR / MC_SUMMARY_NAME
    print(f"[hht-build] Reading MC summary: {mc_path}")
    data = load_mc_summary(mc_path)

    series = build_series(
                          data      = data,
                          bcut      = args.bcut,
                          z_min     = args.z_min,
                          z_max     = args.z_max,
                          N_min     = args.N_min,
                          pA_max    = args.pA_max,
                          pApar_max = args.pApar_max,
                         )

    # Choose the metric for the primary 2-column output
    if args.metric.lower() == "a":
        y            = series["A"]
        metric_label = "A_obs"
        
    elif args.metric.lower() == "apar":
        y            = series["A_par"]
        metric_label = "A_par_obs"
        
    else:  # fpar
        y            = series["f_par"]
        metric_label = "f_par_obs"

    z = series["z_mid"]

    # Build filenames
    def ztag(val):
        return str(val).replace(".", "p")

    base = (
            f"quaia_hht_bcut{int(args.bcut)}"
            f"_z{ztag(args.z_min)}_{ztag(args.z_max)}"
            f"_{args.metric.lower()}"
           )

    out_simple = OUT_DIR / f"{base}_simple.txt"
    out_full   = OUT_DIR / f"{base}_full.txt"

    # 2-column: z_mid, signal(z)
    arr_simple = np.column_stack([z, y])
    np.savetxt(
               out_simple,
               arr_simple,
               fmt    = "%.6f",
               header = f"z_mid  {metric_label}  (|b|>{args.bcut} deg)",
              )

    # Full multi-column output
    arr_full = np.column_stack([
                                z,
                                series["A"],
                                series["A_par"],
                                series["f_par"],
                                series["N_good"],
                                series["p_A"],
                                series["p_Apar"],
                                series["p_fpar"],
                              ])
    header_full = (
                   "z_mid  A_obs  A_par_obs  f_par_obs  N_good  "
                   "p_A  p_|A_par|  p_|f_par|  (|b|>"
                   f"{args.bcut} deg)"
                  )
    np.savetxt(
               out_full,
               arr_full,
               fmt = "%.6e",
               header=header_full,
             )

    print(f"[hht-build] Selected {z.size} bins for |b|>{args.bcut}°")
    print(f"[hht-build] Simple series: {out_simple}")
    print(f"[hht-build] Full table  : {out_full}")


if __name__ == "__main__":
    main()
