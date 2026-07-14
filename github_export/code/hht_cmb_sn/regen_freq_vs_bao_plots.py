#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
regen_freq_vs_bao_plots.py
-----------------------------------
Regenerate the "IMF{n} freq vs BAO proxy — colored by <meta>" plots with
clear coloring and titles. Supports:
  • Single-proxy plotting (DM_over_rd OR DH_over_rd), color-coded by a meta field
    (or by L if no meta is provided)
  • Dual-proxy mode: side-by-side panels for DM_over_rd and DH_over_rd using
    the same color scale and a single shared colorbar.
  • Optional robust smoothing overlay to emphasize the trend.

Inputs expected from a prior run:
  - <outdir>/imf_inst_freq.csv   # has columns: L, IMF0, IMF1, ...
  - Your BAO CSV (with columns z, DM_over_rd, DH_over_rd)
  - Optional meta CSV (with at least columns: z and the chosen meta column)

Usage (examples)
----------------
Single-proxy, color by meta:
    python regen_freq_vs_bao_plots.py \
      --freq-csv "/Users/boyde/.spyder-py3/plamb_runs/hht_sn_BAO_ALL_rd_DR2LyA/imf_inst_freq.csv" \
      --bao "/Users/boyde/.spyder-py3/plamb_runs/bao_eBOSS/bao_scales_ALL_rd_DR2LyA.csv" \
      --bao-cols z DM_over_rd DH_over_rd \
      --meta "/path/to/model_meta.csv" --meta-cols z noether_strength \
      --proxy DM_over_rd --imfs 1 2 3 4 \
      --outdir "/Users/boyde/.spyder-py3/plamb_runs/hht_sn_BAO_ALL_rd_DR2LyA"

Dual-proxy (two panels), color by L (if no meta):
    python regen_freq_vs_bao_plots.py \
      --freq-csv "/Users/boyde/.spyder-py3/plamb_runs/hht_sn_BAO_ALL_rd_DR2LyA/imf_inst_freq.csv" \
      --bao "/Users/boyde/.spyder-py3/plamb_runs/bao_eBOSS/bao_scales_ALL_rd_DR2LyA.csv" \
      --bao-cols z DM_over_rd DH_over_rd \
      --dual --imfs 1 2 3 4 \
      --outdir "/Users/boyde/.spyder-py3/plamb_runs/hht_sn_BAO_ALL_rd_DR2LyA"

Notes
-----
- We intentionally keep the matplotlib look simple (no seaborn).
- If meta is provided, the scatter is colored by that field and a labeled colorbar is drawn.
- If meta is not provided, we color by L (log(1+z)) so the "colored by ..." label is unambiguous.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

def _load_freq_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "L" not in df.columns:
        raise ValueError("imf_inst_freq.csv must contain column 'L'")
    return df

def _load_bao(path: Path, cols: list[str]) -> pd.DataFrame:
    df = pd.read_csv(path)
    for c in cols:
        if c not in df.columns:
            raise ValueError(f"Missing BAO column: {c}")
    return df

def _load_meta(meta_path: Path | None, meta_cols: list[str] | None) -> pd.DataFrame | None:
    if meta_path is None:
        return None
    if meta_cols is None or len(meta_cols) < 2:
        raise ValueError("--meta-cols must be provided as: z <meta_key>")
    mdf = pd.read_csv(meta_path)
    for c in meta_cols:
        if c not in mdf.columns:
            raise ValueError(f"Missing meta column: {c}")
    return mdf

def _build_proxy_on_L(Lgrid: np.ndarray, bao_df: pd.DataFrame, proxy_col: str) -> np.ndarray:
    Lb = np.log1p(bao_df["z"].astype(float).to_numpy())
    vb = bao_df[proxy_col].astype(float).to_numpy()
    f  = interp1d(Lb, vb, kind="linear", bounds_error=False, fill_value=np.nan, assume_sorted=False)
    return f(Lgrid)

def _build_meta_on_L(Lgrid: np.ndarray, meta_df: pd.DataFrame, z_col: str, meta_key: str) -> np.ndarray:
    Lm = np.log1p(meta_df[z_col].astype(float).to_numpy())
    vm = meta_df[meta_key].astype(float).to_numpy()
    f  = interp1d(Lm, vm, kind="nearest", bounds_error=False, fill_value=np.nan, assume_sorted=False)
    return f(Lgrid)

def _smooth_xy(x: np.ndarray, y: np.ndarray, frac: float = 0.06) -> tuple[np.ndarray, np.ndarray]:
    """Return a simple smoothed curve for trend viz: sort by x, then SG filter."""
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 15:
        return np.array([]), np.array([])
    xs, ys = x[m], y[m]
    o = np.argsort(xs)
    xs, ys = xs[o], ys[o]
    n = len(xs)
    win = max(7, int(frac * n) | 1)  # odd window
    win = min(win, n - (1 - n % 2))  # ensure <= n and odd
    if win < 7:
        return xs, ys
    try:
        ys_s = savgol_filter(ys, window_length=win, polyorder=2, mode="interp")
        return xs, ys_s
    except Exception:
        return xs, ys

def make_single_panel(imf_idx: int, L: np.ndarray, f_imf: np.ndarray,
                      bao_on_L: np.ndarray, color_vals: np.ndarray | None,
                      color_label: str, outdir: Path, proxy_name: str):
    plt.figure(figsize=(8, 5))
    m = np.isfinite(bao_on_L) & np.isfinite(f_imf)
    x, y = bao_on_L[m], f_imf[m]
    if color_vals is not None:
        c = color_vals[m]
        norm = Normalize(vmin=np.nanmin(c), vmax=np.nanmax(c))
        sc = plt.scatter(x, y, s=12, c=c, cmap="viridis", norm=norm)
        cbar = plt.colorbar(sc)
        cbar.set_label(color_label)
        subtitle = f"colored by {color_label}"
    else:
        # Color by L for clarity
        c = L[m]
        norm = Normalize(vmin=np.nanmin(c), vmax=np.nanmax(c))
        sc = plt.scatter(x, y, s=12, c=c, cmap="plasma", norm=norm)
        cbar = plt.colorbar(sc)
        cbar.set_label("log(1+z)")
        subtitle = "colored by log(1+z)"

    # Smooth overlay
    xs, ys = _smooth_xy(x, y, frac=0.06)
    if len(xs):
        plt.plot(xs, ys, lw=1.5, alpha=0.9)

    plt.xlabel("BAO proxy (arb.)")
    plt.ylabel("IMF instantaneous freq (cycles/sample)")
    plt.title(f"IMF{imf_idx} freq vs {proxy_name} — {subtitle}")
    out = outdir / f"imf{imf_idx}_freq_vs_{proxy_name}_color.png"
    plt.tight_layout()
    plt.savefig(out, dpi=160, bbox_inches="tight")
    plt.close()

def make_dual_panel(imf_idx: int, L: np.ndarray, f_imf: np.ndarray,
                    bao_dm: np.ndarray, bao_dh: np.ndarray,
                    color_vals: np.ndarray | None, color_label: str, outdir: Path):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    panels = [("DM_over_rd", bao_dm), ("DH_over_rd", bao_dh)]
    # Color setup
    if color_vals is None:
        c_full = L
        label  = "log(1+z)"
        cmap   = "plasma"
    else:
        c_full = color_vals
        label  = color_label
        cmap   = "viridis"
    vmin, vmax = np.nanmin(c_full), np.nanmax(c_full)

    for ax, (name, bx) in zip(axes, panels):
        m = np.isfinite(bx) & np.isfinite(f_imf)
        x, y = bx[m], f_imf[m]
        c    = c_full[m]
        sc = ax.scatter(x, y, s=12, c=c, cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax))
        xs, ys = _smooth_xy(x, y, frac=0.06)
        if len(xs):
            ax.plot(xs, ys, lw=1.5, alpha=0.9)
        ax.set_xlabel(f"{name} (arb.)")
        ax.set_title(f"IMF{imf_idx} vs {name}")

    axes[0].set_ylabel("IMF instantaneous freq (cycles/sample)")
    fig.suptitle(f"IMF{imf_idx} freq vs BAO proxies — colored by {label}", y=0.98)

    # One shared colorbar
    cax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cb  = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap=cmap),
                       cax=cax)
    cb.set_label(label)

    out = outdir / f"imf{imf_idx}_freq_vs_bao_color_dual.png"
    plt.tight_layout(rect=[0, 0, 0.90, 0.95])
    plt.savefig(out, dpi=160, bbox_inches="tight")
    plt.close()

def main():
    ap = argparse.ArgumentParser(description="Regenerate IMF freq vs BAO proxy plots with clear coloring.")
    ap.add_argument("--freq-csv", required=True, type=Path, help="Path to imf_inst_freq.csv produced by your run.")
    ap.add_argument("--bao", required=True, type=Path, help="Path to BAO CSV.")
    ap.add_argument("--bao-cols", nargs="+", default=["z","DM_over_rd","DH_over_rd"],
                    help="Columns in BAO CSV: z DM_over_rd DH_over_rd (at least z + 1 proxy).")
    ap.add_argument("--meta", type=Path, default=None, help="Optional path to meta CSV (for coloring).")
    ap.add_argument("--meta-cols", nargs="+", default=None, help="Columns in meta CSV: z <meta_key>.")
    ap.add_argument("--proxy", type=str, default="DM_over_rd", choices=["DM_over_rd","DH_over_rd"],
                    help="Which BAO proxy to use for single-panel mode.")
    ap.add_argument("--dual", action="store_true", help="If set, make side-by-side panels for both proxies.")
    ap.add_argument("--imfs", nargs="+", type=int, default=[1,2,3,4], help="IMF indices to plot.")
    ap.add_argument("--outdir", type=Path, default=None, help="Where to save figures (default: same as freq-csv).")
    args = ap.parse_args()

    outdir = args.outdir if args.outdir is not None else args.freq_csv.parent
    outdir.mkdir(parents=True, exist_ok=True)

    # Load inputs
    freq_df = _load_freq_csv(args.freq_csv)
    bao_df  = _load_bao(args.bao, args.bao_cols)
    meta_df = _load_meta(args.meta, args.meta_cols)

    L = freq_df["L"].to_numpy(float)

    # Build color field
    color_vals = None
    color_label = ""
    if meta_df is not None:
        z_col, meta_key = args.meta_cols[0], args.meta_cols[1]
        color_vals  = _build_meta_on_L(L, meta_df, z_col, meta_key)
        color_label = meta_key

    # Build BAO proxies on L
    proxies_present = set(args.bao_cols[1:]) & {"DM_over_rd","DH_over_rd"}
    have_dm = "DM_over_rd" in proxies_present
    have_dh = "DH_over_rd" in proxies_present

    dm_on_L = _build_proxy_on_L(L, bao_df, "DM_over_rd") if have_dm else None
    dh_on_L = _build_proxy_on_L(L, bao_df, "DH_over_rd") if have_dh else None

    for imf_idx in args.imfs:
        col = f"IMF{imf_idx}"
        if col not in freq_df.columns:
            print(f"[warn] Column '{col}' not in freq file; skipping.")
            continue
        f_imf = freq_df[col].to_numpy(float)

        if args.dual and (have_dm and have_dh):
            make_dual_panel(imf_idx, L, f_imf, dm_on_L, dh_on_L, color_vals, color_label, outdir)
        else:
            proxy_name = args.proxy
            bx = dm_on_L if proxy_name == "DM_over_rd" else dh_on_L
            if bx is None:
                print(f"[warn] Proxy '{proxy_name}' not present in BAO file; trying the other one.")
                if proxy_name == "DM_over_rd" and have_dh:
                    bx = dh_on_L; proxy_name = "DH_over_rd"
                elif proxy_name == "DH_over_rd" and have_dm:
                    bx = dm_on_L; proxy_name = "DM_over_rd"
                else:
                    print("[error] No BAO proxy available to plot; skipping.")
                    continue
            make_single_panel(imf_idx, L, f_imf, bx, color_vals, color_label or "log(1+z)", outdir, proxy_name)

    print(f"[done] Plots written to: {outdir}")

if __name__ == "__main__":
    main()
