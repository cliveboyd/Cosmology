#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 20:27:42 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_plot_mc_bcut_rotate.py  (0V1)

Plot rotation-based Monte Carlo results for the <z>–dipole.

This script expects NPZ files created by:
    quaia_mc_bcut_rotate.py

Filename pattern:
    quaia_mc_bcut_rotate_{tag}_bcut{BCUT}.npz

Each file contains:
    amps        : MC distribution of |b_vec|
    amp_pars    : MC distribution of A_parallel (projection on CMB direction)
    frac_pars   : MC distribution of f_parallel = A_parallel / A

    obs_amp     : observed |b_vec|
    obs_amp_par : observed A_parallel
    obs_frac_par: observed f_parallel

    z_mid, label, bcut, N_good_obs
    p_amp, p_par_abs, p_frac_abs : precomputed p-values

Usage examples (Spyder/IPython):

    %run /Users/boyde/.spyder-py3/quaia_plot_mc_bcut_rotate.py \
         --bcut 20 --tag z0p10_0p30

    %run /Users/boyde/.spyder-py3/quaia_plot_mc_bcut_rotate.py \
         --bcut 30 --tag z1p00_1p50 --save

"""

import argparse
from   pathlib           import Path

import numpy             as     np
import matplotlib.pyplot as     plt

from   quaia_config      import OUT_DIR


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def load_mc_file(bcut, tag):
    """
    Load the NPZ file for a given (bcut, tag) combination.
    """
    path = OUT_DIR / f"quaia_mc_bcut_rotate_{tag}_bcut{int(bcut)}.npz"
    if not path.exists():
        raise FileNotFoundError(f"MC file not found: {path}")

    data = np.load(path, allow_pickle=True)
    return data, path

def pdf_hist(x, bins=40, range=None):
    n, edges = np.histogram(x, bins=bins, range=range)
    
    if n.sum() == 0:
        return n, edges  # empty; caller should handle
    db       = np.diff(edges)
    pdf      = n / db / n.sum()
    return pdf, edges


def add_hist_panel(ax, samples, obs, xlabel, title, pval, zero_line=False):
    """
    Draw a 1D histogram + vertical line at observed value, with a p-value label.
    """
    samples = np.asarray(samples, dtype=float)

    # Basic histogram
    ax.hist(samples, bins=40,   histtype="stepfilled", alpha=0.6, density=True)
    ax.axvline(obs,  color="k", linestyle="--", linewidth=2, label="observed")

    if zero_line:
        ax.axvline(0.0, color="k", linestyle=":", linewidth=1, label="0")

    ax.set_xlabel(xlabel)
    ax.set_ylabel("PDF")
    ax.set_title(title)

    # Text box with observed value and p-value
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    x_text       = x_min + 0.05 * (x_max - x_min)
    y_text       = y_max * 0.80

    ax.text(
            x_text,
            y_text,
            f"obs = {obs:.3e}\n"
            f"p = {pval:.3f}",
            fontsize = 10,
            bbox     = dict(boxstyle="round", facecolor="white", alpha=0.8),
           )

    ax.legend(loc="best")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
                        "--bcut",
                        type=float,
                        required=True,
                        help="Galactic latitude cut |b|>bcut (must match MC file)",
                       )
    parser.add_argument(
                        "--tag",
                        type=str,
                        required=True,
                        help="redshift-slice tag, e.g. z0p10_0p30, z0p75_1p00, z1p00_1p50",
                       )
    parser.add_argument(
                        "--save",
                        action="store_true",
                        help="save PNG to OUT_DIR instead of (or as well as) showing on screen",
                       )

    args = parser.parse_args()

    data, path = load_mc_file(args.bcut, args.tag)

    amps         = data["amps"]
    amp_pars     = data["amp_pars"]
    frac_pars    = data["frac_pars"]

    obs_amp      = float(data["obs_amp"])
    obs_amp_par  = float(data["obs_amp_par"])
    obs_frac_par = float(data["obs_frac_par"])

    z_mid        = float(data["z_mid"])
    label        = str(data["label"])
    bcut         = float(data["bcut"])
    N_good_obs   = int(data["N_good_obs"])

    p_amp        = float(data["p_amp"])
    p_par_abs    = float(data["p_par_abs"])
    p_frac_abs   = float(data["p_frac_abs"])

    print(f"[plot-rot] loaded {path}")
    print(
          f"[plot-rot] |b|>{bcut:.1f}, {label}, z_mid={z_mid:.2f}, "
          f"N_good_obs={N_good_obs}"
         )
    print(
          f"[plot-rot] obs A={obs_amp:.3e}, A_par={obs_amp_par:.3e}, "
          f"f_par={obs_frac_par:+.3f}"
         )
    print(
          f"[plot-rot] p(A)={p_amp:.3f}, "
          f"p(|A_par|>=|obs|)={p_par_abs:.3f}, "
          f"p(|f_par|>=|obs|)={p_frac_abs:.3f}"
         )

    # 3-panel figure
    fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=False)
    fig.suptitle(
                 rf"Rotation MC for $\langle z\rangle$ dipole "
                 rf"(|b|>{bcut:.0f}^\circ, {label}, $z_\mathrm{{mid}}={z_mid:.2f}$)",
                 fontsize=13,
                )

    add_hist_panel(
                   axes[0],
                   amps,
                   obs_amp,
                   xlabel    = r"$A = |{\bf b}|$",
                   title     = r"Dipole amplitude $A$",
                   pval      = p_amp,
                   zero_line = False,
                  )

    add_hist_panel(
                   axes[1],
                   amp_pars,
                   obs_amp_par,
                   xlabel    = r"$A_\parallel$ (CMB–parallel component)",
                   title     = r"CMB–parallel component $A_\parallel$",
                   pval      = p_par_abs,
                   zero_line = True,
                  )

    add_hist_panel(
                   axes[2],
                   frac_pars,
                   obs_frac_par,
                   xlabel    = r"$f_\parallel = A_\parallel / A$",
                   title     = r"CMB–parallel fraction $f_\parallel$",
                   pval      = p_frac_abs,
                   zero_line = True,
                  )
    
    fig.tight_layout(rect=[0, 0.0, 1, 0.95])

    if args.save:
        out_png = OUT_DIR / f"quaia_mc_bcut_rotate_{args.tag}_bcut{int(bcut)}.png"
        fig.savefig(out_png, dpi=150)
        print(f"[plot-rot] saved figure -> {out_png}")

    # Show on screen (Spyder/interactive)
    plt.show()


if __name__ == "__main__":
    main()
