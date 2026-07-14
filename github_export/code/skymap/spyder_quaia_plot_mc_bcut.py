#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 18:37:56 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_plot_mc_bcut.py  (0V1)

Simple visualisation of the Monte Carlo distributions produced by
quaia_mc_bcut_shuffle.py.

Usage (from Spyder/IPython):

    %run /Users/boyde/.spyder-py3/quaia_plot_mc_bcut.py --bcut 20 --zlo 0.10 --zhi 0.30

This will look for:
    OUT_DIR / f"quaia_mc_bcut_shuffle_z0p10_0p30_bcut20.npz"

and plot histograms for A, A_par, and f_par with vertical lines at the observed
values.
"""

import argparse
import numpy             as     np
import matplotlib.pyplot as     plt

from   quaia_config      import OUT_DIR


def build_tag(z_lo, z_hi):
    return f"z{z_lo:.2f}_{z_hi:.2f}".replace(".", "p")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bcut", type=float, required=True)
    parser.add_argument("--zlo",  type=float, required=True)
    parser.add_argument("--zhi",  type=float, required=True)
    args      = parser.parse_args()

    tag       = build_tag(args.zlo, args.zhi)
    path      = OUT_DIR / f"quaia_mc_bcut_shuffle_{tag}_bcut{int(args.bcut)}.npz"

    data      = np.load(path, allow_pickle=True)
    amps      = data["amps"]
    amp_pars  = data["amp_pars"]
    frac_pars = data["frac_pars"]
    A_obs     = float(data["obs_amp"])
    A_par_obs = float(data["obs_amp_par"])
    f_obs     = float(data["obs_frac_par"])
    z_mid     = float(data["z_mid"])
    bcut      = float(data["bcut"])

    print(f"[plot-mc] file: {path}")
    print(f"[plot-mc] |b|>{bcut:.1f}, tag={tag}, z_mid={z_mid:.2f}")
    print(f"[plot-mc] A_obs={A_obs:.3e}, A_par_obs={A_par_obs:.3e}, f_par_obs={f_obs:+.3f}")

    # --- A ---
    plt.figure(figsize=(7, 5))
    plt.hist(amps, bins=40, alpha=0.7)
    plt.axvline(A_obs, color="r", linestyle="--", label="observed")
    plt.xlabel("dipole amplitude A")
    plt.ylabel("count")
    plt.title(f"A MC distribution |b|>{bcut:.0f}, z~{z_mid:.2f}")
    plt.legend()
    plt.tight_layout()

    # --- A_par ---
    plt.figure(figsize=(7, 5))
    plt.hist(amp_pars, bins=40, alpha=0.7)
    plt.axvline(A_par_obs, color="r", linestyle="--", label="observed")
    plt.xlabel("CMB-parallel amplitude A_par")
    plt.ylabel("count")
    plt.title(f"A_par MC distribution |b|>{bcut:.0f}, z~{z_mid:.2f}")
    plt.legend()
    plt.tight_layout()

    # --- f_par ---
    plt.figure(figsize=(7, 5))
    plt.hist(frac_pars, bins=40, alpha=0.7)
    plt.axvline(f_obs, color="r", linestyle="--", label="observed")
    plt.xlabel("fraction f_par = A_par / A")
    plt.ylabel("count")
    plt.title(f"f_par MC distribution |b|>{bcut:.0f}, z~{z_mid:.2f}")
    plt.legend()
    plt.tight_layout()

    print("[plot-mc] done (figures are in the Plots pane / backend).")


if __name__ == "__main__":
    main()
