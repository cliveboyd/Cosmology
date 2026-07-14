#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 09:03:31 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_plot_bcut_grid.py  (0V1)

Read quaia_dipole_bcut_grid.txt (written by quaia_dipole_bcut_grid.py),
reconstruct the z-slice midpoints, compute the CMB-parallel component of
the <z> dipole, and make summary plots:

  - dipole amplitude |b| vs z
  - CMB-parallel amplitude amp_par vs z
  - CMB-parallel fraction frac_par = amp_par / amp vs z

Also print a compact table of (bcut, z_mid, amp, amp_par, frac_par).

Usage (in IPython / Spyder):

    %run /Users/boyde/.spyder-py3/quaia_plot_bcut_grid.py
    
"""

import numpy             as     np
import matplotlib.pyplot as     plt
from   pathlib           import Path

from   quaia_config      import OUT_DIR

# Must match the bcut list used in quaia_dipole_bcut_grid.py
BCUT_LIST = [10.0, 15.0, 20.0, 25.0, 30.0, 35.0]


def parse_grid_file(path: Path):
    """
    Parse quaia_dipole_bcut_grid.txt.

    Returns a list of dicts containing:
      bcut, tag, z_mid, N_good, amp, l_deg, b_deg, sep_cmb
    """
    records = []

    with path.open("r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split()
            if len(cols) < 7:
                continue

            bcut    = float(cols[0])
            tag     = cols[1]
            N_good  = int(cols[2])
            amp     = float(cols[3])
            l_deg   = float(cols[4])
            b_deg   = float(cols[5])
            sep_cmb = float(cols[6])

            # Work out z_mid from tag if it's a z-slice (e.g. 'z0p10_0p30')
            if tag == "full":
                z_mid = None
            else:
                # tag format: z0p10_0p30, z0p30_0p50, ...
                try:
                    z_lo_str, z_hi_str = tag[1:].split("_")
                    z_lo               = float(z_lo_str.replace("p", "."))
                    z_hi               = float(z_hi_str.replace("p", "."))
                    z_mid              = 0.5 * (z_lo + z_hi)
                except Exception as exc:
                    print(f"[plot-grid] WARNING: could not parse tag='{tag}': {exc}")
                    z_mid = None

            records.append(
                            dict(
                                bcut=bcut,
                                tag=tag,
                                N_good=N_good,
                                amp=amp,
                                l_deg=l_deg,
                                b_deg=b_deg,
                                sep_cmb=sep_cmb,
                                z_mid=z_mid,
                               )
                           )

    return records


def build_series(records):
    """
    Group records by bcut and build arrays of
      z_mid, amp, amp_par, frac_par

    Uses sep_cmb to compute:
      frac_par = cos(sep_cmb)
      amp_par  = amp * frac_par
    """
    series = {b: [] for b in BCUT_LIST}

    for r in records:
        if r["z_mid"] is None:
            # Skip the 'full sample' entries for z-based plots
            continue

        bcut = r["bcut"]
        if bcut not in series:
            continue

        z_mid    = r["z_mid"]
        amp      = r["amp"]
        theta    = np.deg2rad(r["sep_cmb"])
        frac_par = np.cos(theta)          # signed, between -1 and +1
        amp_par  = amp * frac_par

        series[bcut].append((z_mid, amp, amp_par, frac_par))

    # Convert lists to sorted arrays
    out = {}
    for bcut, rows in series.items():
        if not rows:
            out[bcut] = dict(
                            z        = np.array([]),
                            amp      = np.array([]),
                            amp_par  = np.array([]),
                            frac_par = np.array([]),
                            )
            continue

        arr = np.array(rows)  # shape (nslice, 4)
        idx = np.argsort(arr[:, 0])

        out[bcut] = dict(
                         z        = arr[idx, 0],
                         amp      = arr[idx, 1],
                         amp_par  = arr[idx, 2],
                         frac_par = arr[idx, 3],
                        )
    return out


def print_table(series):
    """
    Print a compact table for notes / LaTeX.
    """
    print("# bcut  z_mid   amp        amp_par    frac_par")
    for bcut in BCUT_LIST:
        s = series[bcut]
        for z, amp, amp_par, frac in zip(
            s["z"], s["amp"], s["amp_par"], s["frac_par"]
        ):
            print(
                  f"{bcut:4.1f}  {z:5.2f}  "
                  f"{amp:9.3e}  {amp_par:9.3e}  {frac:+7.3f}"
                 )

def plot_amp(series):
    """
    Plot dipole amplitude |b| vs z for each bcut.
    """
    plt.figure(figsize=(7, 5))
    markers = ["o", "s", "^", "D", "v", "x"]
    
    for i, bcut in enumerate(BCUT_LIST):
        s = series[bcut]
        if s["z"].size == 0:
            continue
        plt.plot(
                 s["z"],
                 s["amp"],
                 marker    = markers[i % len(markers)],
                 linestyle = "-",
                 label     = f"|b|>{int(bcut)}°",
                )

    plt.xlabel("z (slice midpoint)")
    plt.ylabel("dipole amplitude |b|")
    plt.title( "<z> dipole amplitude vs z")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_DIR / "quaia_bcut_amp_vs_z.png", dpi=200)


def plot_amp_par(series):
    """
    Plot CMB-parallel component amp_par vs z.
    """
    plt.figure(figsize=(7, 5))
    markers = ["o", "s", "^", "D", "v", "x"]

    for i, bcut in enumerate(BCUT_LIST):
        s = series[bcut]
        if s["z"].size == 0:
            continue
        plt.plot(
                 s["z"], s["amp_par"],
                 marker    = markers[i % len(markers)],
                 linestyle = "-",
                 label     = f"|b|>{int(bcut)}°",
                )

    plt.axhline(0.0, ls="--", color="k", alpha=0.4)
    plt.xlabel("z (slice midpoint)")
    plt.ylabel("CMB-parallel dipole amp_par")
    plt.title( "CMB-parallel component of <z> dipole vs z")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_DIR / "quaia_bcut_amp_par_vs_z.png", dpi=200)


def plot_frac_par(series):
    """
    Plot frac_par = amp_par / amp vs z (i.e. cos(sepCMB)).
    """
    plt.figure(figsize=(7, 5))
    markers = ["o", "s", "^", "D", "v", "x"]

    for i, bcut in enumerate(BCUT_LIST):
        s = series[bcut]
        if s["z"].size == 0:
            continue
        plt.plot(
                 s["z"], s["frac_par"],
                 marker    = markers[i % len(markers)],
                 linestyle = "-",
                 label     = f"|b|>{int(bcut)}°",
                )

    plt.axhline(0.0, ls="--", color="k", alpha=0.4)
    plt.ylim(-1.05, 1.05)
    plt.xlabel("z (slice midpoint)")
    plt.ylabel("frac_par = amp_par / amp")
    plt.title( "CMB-parallel fraction of <z> dipole vs z")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_DIR / "quaia_bcut_frac_par_vs_z.png", dpi=200)


def main():
    grid_path = OUT_DIR / "quaia_dipole_bcut_grid.txt"
    print(f"[plot-grid] reading {grid_path}")

    records = parse_grid_file(grid_path)
    series  = build_series(records)

    # Print table for notes
    print_table(series)

    # Make plots
    plot_amp(series)
    plot_amp_par(series)
    plot_frac_par(series)

    # Let the current matplotlib backend handle showing / saving
    plt.show()


if __name__ == "__main__":
    main()
