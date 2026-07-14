from __future__ import annotations
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 31 15:26:22 2026

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_imf_diagnostics.py

Load Quaia FITS quasar catalogue, build a 1D signal vs angle bin (galactic longitude by default),
then decompose with CEEMDAN/EMD and plot IMFs.

Typical use:
  python quaia_imf_diagnostics.py --fits /path/quaia_G20.0.fits --out outdir --signal counts

Outputs:
  - signal.csv         (x, signal)
  - imfs.npy           (imfs array)
  - residue.npy        (residue array)
  - imf_stack.png      (stacked plot: signal + IMFs + residue)
  - imf_overlay.png    (overlay plot)

Notes:
  - IMFs are 1D. We reduce the sky catalogue to a 1D angular signal by binning in longitude (or RA).
  - For dipole / anisotropy work, galactic longitude is usually the right axis.
"""


import argparse
from   pathlib import Path
import numpy   as     np

def _clean_path(s: str) -> str:
    # protects against Spyder runfile passing quotes through to argparse
    return s.strip().strip('"').strip("'")

def _moving_average_circular(y: np.ndarray, win: int) -> np.ndarray:
    if win <= 1:
        return y
    win = int(win)
    if win % 2 == 0:
        win += 1
    pad = win // 2
    w               = np.ones(win, dtype=float) / float(win)
    y2              = np.r_[y[-pad:], y, y[:pad]]
    ys              = np.convolve(y2, w, mode="same")[pad:-pad]
    return ys

def _load_fits_cols(fits_path: Path, cols: list[str]) -> dict[str, np.ndarray]:
    # Try fitsio first (fast). Fall back to astropy.
    try:
        import fitsio  # type: ignore
        arr         = fitsio.read(str(fits_path), columns=cols, ext=1)
        out         = {c: np.asarray(arr[c]) for c in cols if c in arr.dtype.names}
        return out
    except Exception:
        pass

    try:
        from astropy.table import Table  # type: ignore
        tab         = Table.read(str(fits_path), hdu=1, memmap=True)
        out         = {}
        for c in cols:
            if c in tab.colnames:
                out[c] = np.asarray(tab[c])
        return out
    except Exception as e:
        raise RuntimeError(
            "Could not read FITS. Install fitsio or astropy.\n"
            "Try: conda install -c conda-forge fitsio astropy"
        ) from e

def _pick_first_existing(candidates: list[str], available: set[str]) -> str | None:
    for c in candidates:
        if c in available:
            return c
    return None

def _icrs_to_galactic_lonlat_deg(ra_deg: np.ndarray, dec_deg: np.ndarray, chunk: int = 200_000):
    try:
        from astropy.coordinates import SkyCoord  # type: ignore
        import astropy.units as u  # type: ignore
    except Exception as e:
        raise RuntimeError("Need astropy for ICRS -> Galactic conversion. Install astropy.") from e

    n            = ra_deg.size
    l               = np.empty(n, dtype=float)
    b               = np.empty(n, dtype=float)
    for i in range(0, n, chunk):
        j           = min(i + chunk, n)
        c           = SkyCoord(ra=ra_deg[i:j] * u.deg, dec=dec_deg[i:j] * u.deg, frame="icrs")
        g        = c.galactic
        l[i:j]      = g.l.deg
        b[i:j]   = g.b.deg
    return l, b

def build_signal(
                 x_deg       : np.ndarray,
                 z           : np.ndarray | None,
                 nbins       : int,
                 signal_mode : str,
                ) -> tuple[np.ndarray, np.ndarray]:
    # x bins [0, 360)
    x               = np.mod(x_deg, 360.0)
    edges           = np.linspace(0.0, 360.0, nbins + 1)
    centers         = 0.5 * (edges[:-1] + edges[1:])

    if signal_mode == "counts":
        counts, _   = np.histogram(x, bins=edges)
        y           = counts.astype(float)
        return centers, y

    if signal_mode == "zmean":
        if z is None:
            raise ValueError("signal_mode=zmean requires a redshift column.")
        idx         = np.digitize(x, edges) - 1
        good        = (idx >= 0) & (idx < nbins) & np.isfinite(z)
        idx         = idx[good]
        zv          = z[good].astype(float)

        sumz        = np.bincount(idx, weights=zv, minlength=nbins).astype(float)
        cnt         = np.bincount(idx, minlength=nbins).astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            y       = sumz / cnt
        # fill empty bins with 0 (or you can np.nan_to_num with nan=0)
        y           = np.nan_to_num(y, nan=0.0, posinf=0.0, neginf=0.0)
        return centers, y

    raise ValueError(f"Unknown signal_mode: {signal_mode}")

def emd_decompose(y: np.ndarray, ce_trials: int, noise: float, seed: int):
    # Prefer CEEMDAN if available, else EMD fallback
    try:
        from PyEMD import CEEMDAN  # type: ignore
        ce              = CEEMDAN(trials=int(ce_trials), noise_width=float(noise))
        ce.random_seed  = int(seed)
        imfs            = ce.ceemdan(y)
        residue         = y - np.sum(imfs, axis=0)
        return imfs, residue, "CEEMDAN"
    except Exception:
        try:
            from PyEMD import EMD  # type: ignore
            emd         = EMD()
            imfs        = emd.emd(y)
            residue     = y - np.sum(imfs, axis=0)
            return imfs, residue, "EMD"
        
        except Exception as e:
            raise RuntimeError(
                                "No PyEMD available. Install with:\n"
                                "  pip install EMD-signal\n"
                                "(import is `from PyEMD import CEEMDAN`)\n"
                              ) from e

def plot_imfs(x, y, imfs, residue, out_dir: Path, dpi: int, title: str):
    import matplotlib.pyplot as plt

    n_imf = imfs.shape[0]
    # stacked plot
    rows            = 1 + n_imf + 1
    fig, axes = plt.subplots(rows, 1, figsize=(12, 2.0 * rows), sharex=True)
    axes[0].plot(x, y, linewidth=1.0)
    axes[0].set_ylabel("signal")
    axes[0].set_title(title)

    for k in range(n_imf):
        axes[1 + k].plot(x, imfs[k], linewidth=0.9)
        axes[1 + k].set_ylabel(f"IMF{k+1}")

    axes[-1].plot(x, residue, linewidth=0.9)
    axes[-1].set_ylabel("resid")
    axes[-1].set_xlabel("angle (deg)")

    fig.tight_layout()
    fig.savefig(out_dir / "imf_stack.png", dpi=int(dpi))
    plt.close(fig)

    # overlay plot (first few imfs)
    fig2            = plt.figure(figsize=(12, 6))
    plt.plot(x, y, linewidth=1.2, label="signal")
    top             = min(6, n_imf)
    for k in range(top):
        plt.plot(x, imfs[k], linewidth=0.9, label=f"IMF{k+1}")
    plt.xlabel("angle (deg)")
    plt.title(title + f" (overlay first {top} IMFs)")
    plt.legend(loc="best", fontsize=9)
    plt.tight_layout()
    fig2.savefig(out_dir / "imf_overlay.png", dpi=int(dpi))
    plt.close(fig2)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fits", required=True, type=str, help="Quaia FITS file (e.g., quaia_G20.0.fits)")
    ap.add_argument("--out", required=True, type=str, help="Output directory")

    ap.add_argument("--xcoord", default="galactic_l", choices=["galactic_l", "ra"],
                    help="Axis used to build 1D signal.")
    ap.add_argument("--bcut", default=0.0, type=float,
                    help="Apply |b| > bcut in degrees (only when xcoord=galactic_l).")

    ap.add_argument("--nbins", default=720, type=int, help="Number of bins in 0..360")
    ap.add_argument("--signal", default="counts", choices=["counts", "zmean"],
                    help="What to bin into the 1D signal.")

    ap.add_argument("--zcol", default="", type=str, help="Optional override: redshift column name")
    ap.add_argument("--racol", default="", type=str, help="Optional override: RA column name")
    ap.add_argument("--deccol", default="", type=str, help="Optional override: DEC column name")
    ap.add_argument("--lcol", default="", type=str, help="Optional override: galactic l column name")
    ap.add_argument("--bcol", default="", type=str, help="Optional override: galactic b column name")

    ap.add_argument("--smooth", default=1, type=int, help="Circular moving-average window (odd).")
    ap.add_argument("--ce-trials", default=200, type=int, help="CEEMDAN trials (if available).")
    ap.add_argument("--noise", default=0.2, type=float, help="CEEMDAN noise_width.")
    ap.add_argument("--seed", default=7, type=int, help="Random seed.")
    ap.add_argument("--max-imfs", default=14, type=int, help="Limit number of IMFs saved/plotted.")
    ap.add_argument("--dpi", default=180, type=int)

    args = ap.parse_args()

    fits_path       = Path(_clean_path(args.fits)).expanduser()
    out_dir         = Path(_clean_path(args.out)).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not fits_path.exists():
        raise FileNotFoundError(f"FITS not found: {fits_path}")

    # Common Quaia-ish candidates (robust guessing)
    ra_candidates   = [args.racol] if args.racol else ["RA", "ra", "ALPHA_J2000", "alpha", "RA_ICRS"]
    dec_candidates  = [args.deccol] if args.deccol else ["DEC", "dec", "DELTA_J2000", "delta", "DEC_ICRS"]

    # Many catalogues use some z phot column names; this is just best-effort
    z_candidates    = [args.zcol] if args.zcol else ["z", "Z", "Z_PHOT", "z_phot", "redshift", "ZBEST", "Z_BEST", "Z_QSO"]

    l_candidates    = [args.lcol] if args.lcol else ["l", "L", "GLON", "gal_l", "l_gal", "GAL_L"]
    b_candidates    = [args.bcol] if args.bcol else ["b", "B", "GLAT", "gal_b", "b_gal", "GAL_B"]

    # Load a superset; we'll pick what exists
    want_cols       = sorted(set(ra_candidates + dec_candidates + z_candidates + l_candidates + b_candidates))
    want_cols       = [c for c in want_cols if c]  # drop blanks

    cols            = _load_fits_cols(fits_path, want_cols)
    avail           = set(cols.keys())
    if not avail:
        raise RuntimeError("No requested columns were found in FITS. Check ext=1 and column names.")

    ra_col          = _pick_first_existing(ra_candidates, avail)
    dec_col         = _pick_first_existing(dec_candidates, avail)
    z_col           = _pick_first_existing(z_candidates, avail)
    l_col           = _pick_first_existing(l_candidates, avail)
    b_col           = _pick_first_existing(b_candidates, avail)

    print("\n[quaia-imf] FITS:", fits_path)
    print("[quaia-imf] Available (loaded) cols:", sorted(avail))
    print("[quaia-imf] Picked ra/dec:", ra_col, dec_col)
    print("[quaia-imf] Picked z:", z_col)
    print("[quaia-imf] Picked l/b:", l_col, b_col)
    print("[quaia-imf] xcoord:", args.xcoord, "signal:", args.signal, "nbins:", args.nbins)

    if args.xcoord == "ra":
        if ra_col is None:
            raise RuntimeError("xcoord=ra requested but RA column not found.")
        x_deg       = np.asarray(cols[ra_col], dtype=float)
        b_deg       = None
    else:
        # galactic_l
        if l_col is not None and b_col is not None:
            x_deg           = np.asarray(cols[l_col], dtype=float)
            b_deg           = np.asarray(cols[b_col], dtype=float)
        else:
            if ra_col is None or dec_col is None:
                raise RuntimeError("Need RA/DEC to compute galactic coordinates, but columns not found.")
            ra              = np.asarray(cols[ra_col], dtype=float)
            dec             = np.asarray(cols[dec_col], dtype=float)
            x_deg, b_deg    = _icrs_to_galactic_lonlat_deg(ra, dec)

    z                       = np.asarray(cols[z_col], dtype=float) if (z_col is not None and args.signal == "zmean") else None

    # Apply bcut if relevant
    if args.xcoord == "galactic_l" and float(args.bcut) > 0.0:
        if b_deg is None:
            raise RuntimeError("bcut requested but no galactic latitude available.")
            
        m = np.isfinite(x_deg) & np.isfinite(b_deg)
        
        if z is not None:
            m      &= np.isfinite(z)
            
        m          &= (np.abs(b_deg) > float(args.bcut))
        x_deg       = x_deg[m]
        
        if z is not None:
            z       = z[m]
        print(f"[quaia-imf] bcut |b|>{args.bcut} kept {x_deg.size} rows")

    # Build 1D signal
    x, y            = build_signal(x_deg, z, int(args.nbins), args.signal)
    y               = _moving_average_circular(y, int(args.smooth))

    # Save signal
    import pandas as pd
    pd.DataFrame({"x_deg": x, "signal": y}).to_csv(out_dir / "signal.csv", index=False)

    # Decompose
    imfs, residue, method = emd_decompose(y, args.ce_trials, args.noise, args.seed)

    # Limit IMFs for plotting/saving
    max_imfs        = int(args.max_imfs)
    if imfs.shape[0] > max_imfs:
        imfs_use    = imfs[:max_imfs]
    else:
        imfs_use    = imfs

    np.save(out_dir / "imfs.npy", imfs_use)
    np.save(out_dir / "residue.npy", residue)

    title = f"Quaia IMF decomposition ({method}) | x={args.xcoord} | signal={args.signal} | nbins={args.nbins} | bcut={args.bcut}"
    plot_imfs(x, y, imfs_use, residue, out_dir, args.dpi, title)

    print("\n[quaia-imf] Done. Wrote:")
    print(" ", out_dir / "signal.csv")
    print(" ", out_dir / "imfs.npy")
    print(" ", out_dir / "residue.npy")
    print(" ", out_dir / "imf_stack.png")
    print(" ", out_dir / "imf_overlay.png")

if __name__ == "__main__":
    main()
