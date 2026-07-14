#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  6 19:49:17 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hht_quaia_residuals.py  (0V1)

HHT on Quaia dipole/residual series y(z):

Input:
    2-column ASCII file, typically:
        col 0: z (or z_mid)
        col 1: signal (e.g. A_par(z) aligned with CMB, or f_par)

Pipeline:
    * Sort by z (or log(1+z))
    * Optional polynomial de-trend in L = log(1+z)
    * Resample to a uniform grid in L
    * Apply cosine taper at edges
    * EMD/EEMD decomposition
    * Hilbert transform -> instantaneous frequency & amplitude
    * Surrogate significance (phase-random + circular block bootstrap)
    * Save detail NPZ + stats CSV/NPZ + summary plot

Typical usage (Spyder/IPython):

    %run '/Users/boyde/.spyder-py3/hht_quaia_residuals.py' \
         '/Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_hht_bcut25_z0p30_2p30_apar_simple.txt' \
         'plamb_runs/hht_quaia/bcut25_z0p30_2p30'


OUTPUT...

script = "/Users/boyde/.spyder-py3/hht_quaia_residuals.py"
bcuts  = [10, 15, 20, 25, 30, 35]

for bcut_target in bcuts:
    series = f"/Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_hht_bcut{int(bcut_target)}_allz_apar_simple.txt"
    outdir = f"plamb_runs/hht_quaia/bcut{int(bcut_target)}_allz_apar"

    print("Running:", script, series, outdir)
    %run $script $series $outdir
    
Running: /Users/boyde/.spyder-py3/hht_quaia_residuals.py /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_hht_bcut10_allz_apar_simple.txt plamb_runs/hht_quaia/bcut10_allz_apar
[hht] series file: /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_hht_bcut10_allz_apar_simple.txt
[hht] outdir     : /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut10_allz_apar
[clip] kept 1.00 of points.
[grid] uniform log(1+z) bins = 512, L ∈ [0.4055,3.8395]
[save] HHT detail -> /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut10_allz_apar/hht_quaia_detail.npz
[save] IMF stats -> /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut10_allz_apar/hht_quaia_imf_stats.csv, /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut10_allz_apar/hht_quaia_imf_stats.npz
[files] /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut10_allz_apar/hht_quaia.png
Running: /Users/boyde/.spyder-py3/hht_quaia_residuals.py /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_hht_bcut15_allz_apar_simple.txt plamb_runs/hht_quaia/bcut15_allz_apar
[hht] series file: /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_hht_bcut15_allz_apar_simple.txt
[hht] outdir     : /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut15_allz_apar
[clip] kept 1.00 of points.
[grid] uniform log(1+z) bins = 512, L ∈ [0.4055,3.8395]
[save] HHT detail -> /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut15_allz_apar/hht_quaia_detail.npz
[save] IMF stats -> /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut15_allz_apar/hht_quaia_imf_stats.csv, /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut15_allz_apar/hht_quaia_imf_stats.npz
[files] /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut15_allz_apar/hht_quaia.png
<Figure size 864x576 with 0 Axes>
Running: /Users/boyde/.spyder-py3/hht_quaia_residuals.py /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_hht_bcut20_allz_apar_simple.txt plamb_runs/hht_quaia/bcut20_allz_apar
[hht] series file: /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_hht_bcut20_allz_apar_simple.txt
[hht] outdir     : /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut20_allz_apar
[clip] kept 1.00 of points.
[grid] uniform log(1+z) bins = 512, L ∈ [0.4055,3.8395]
[save] HHT detail -> /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut20_allz_apar/hht_quaia_detail.npz
[save] IMF stats -> /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut20_allz_apar/hht_quaia_imf_stats.csv, /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut20_allz_apar/hht_quaia_imf_stats.npz
[files] /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut20_allz_apar/hht_quaia.png
<Figure size 864x576 with 0 Axes>
Running: /Users/boyde/.spyder-py3/hht_quaia_residuals.py /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_hht_bcut25_allz_apar_simple.txt plamb_runs/hht_quaia/bcut25_allz_apar
[hht] series file: /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_hht_bcut25_allz_apar_simple.txt
[hht] outdir     : /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut25_allz_apar
[clip] kept 1.00 of points.
[grid] uniform log(1+z) bins = 512, L ∈ [0.4055,3.8395]
[save] HHT detail -> /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut25_allz_apar/hht_quaia_detail.npz
[save] IMF stats -> /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut25_allz_apar/hht_quaia_imf_stats.csv, /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut25_allz_apar/hht_quaia_imf_stats.npz
[files] /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut25_allz_apar/hht_quaia.png
<Figure size 864x576 with 0 Axes>
Running: /Users/boyde/.spyder-py3/hht_quaia_residuals.py /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_hht_bcut30_allz_apar_simple.txt plamb_runs/hht_quaia/bcut30_allz_apar
[hht] series file: /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_hht_bcut30_allz_apar_simple.txt
[hht] outdir     : /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut30_allz_apar
[clip] kept 1.00 of points.
[grid] uniform log(1+z) bins = 512, L ∈ [0.4055,3.8395]
[save] HHT detail -> /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut30_allz_apar/hht_quaia_detail.npz
[save] IMF stats -> /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut30_allz_apar/hht_quaia_imf_stats.csv, /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut30_allz_apar/hht_quaia_imf_stats.npz
[files] /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut30_allz_apar/hht_quaia.png
<Figure size 864x576 with 0 Axes>
Running: /Users/boyde/.spyder-py3/hht_quaia_residuals.py /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_hht_bcut35_allz_apar_simple.txt plamb_runs/hht_quaia/bcut35_allz_apar
[hht] series file: /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_hht_bcut35_allz_apar_simple.txt
[hht] outdir     : /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut35_allz_apar
[clip] kept 1.00 of points.
[grid] uniform log(1+z) bins = 512, L ∈ [0.4055,3.8395]
[save] HHT detail -> /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut35_allz_apar/hht_quaia_detail.npz
[save] IMF stats -> /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut35_allz_apar/hht_quaia_imf_stats.csv, /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut35_allz_apar/hht_quaia_imf_stats.npz
[files] /Users/boyde/.spyder-py3/plamb_runs/hht_quaia/bcut35_allz_apar/hht_quaia.png
<Figure size 864x576 with 0 Axes>
<Figure size 864x576 with 0 Axes>

"""

import os
import sys
import json
import numpy             as     np
import matplotlib.pyplot as     plt
from   datetime          import datetime

plt.ioff()

# ---------- optional deps ----------
try:
    from PyEMD import EMD as _EMD
    HAS_EMD = True
except Exception:
    HAS_EMD = False

try:
    from PyEMD import EEMD as _EEMD
    HAS_EEMD = True
except Exception:
    HAS_EEMD = False

try:
    from scipy.signal import hilbert as _hilbert
    HAS_HILBERT = True
    
except Exception:
    HAS_HILBERT = False


# ---------- helpers ----------
def _instantaneous_freq_phase(signal, x):
    """
    Instantaneous frequency over x (monotone, typically x = log(1+z)).
    Prefers analytic-signal Hilbert; falls back to finite differences.
    """
    y  = np.asarray(signal, float)
    x  = np.asarray(x,      float)
    dx = np.gradient(x)
    dx = np.where(np.abs(dx) < 1e-9, 1e-9, dx)

    if HAS_HILBERT:
        a         = _hilbert(y)
        phase     = np.unwrap(np.angle(a))
        inst_freq = np.gradient(phase) / np.maximum(dx, 1e-12)
        amp_env   = np.abs(a)
        return inst_freq, phase, amp_env

    # fallback: pseudo-phase from cumulative sum
    y_norm    = (y - y.mean()) / (y.std() + 1e-12)
    phase     = np.unwrap(np.cumsum(y_norm) * dx / np.max(dx))
    inst_freq = np.gradient(phase) / np.maximum(dx, 1e-12)
    amp_env   = np.abs(y_norm)
    return inst_freq, phase, amp_env


def _phase_randomize(y):
    """Fourier phase-randomization surrogate (preserves power spectrum)."""
    y   = np.asarray(y, float)
    Y   = np.fft.rfft(y)
    mag = np.abs(Y)
    ph  = np.angle(Y)
    rnd = np.random.uniform(-np.pi, np.pi, size=ph.shape)
    Ys  = mag * np.exp(1j * rnd)
    ys  = np.fft.irfft(Ys, n=y.size)
    return ys


def _cbb_surrogate(y, block=48):
    """
    Circular Block Bootstrap surrogate: preserves local correlation structure.
    """
    y      = np.asarray(y, float)
    n      = y.size
    k      = int(np.ceil(n / block))
    starts = np.random.randint(0, n, size=k)
    segs   = []
    for s in starts:
        e = s + block
        if e <= n:
            segs.append(y[s:e])
        else:
            wrap = (s + block) % n
            segs.append(np.r_[y[s:], y[:wrap]])
    ys    = np.concatenate(segs)[:n]
    return ys


def main(series_path,
         outdir="plamb_runs/hht_quaia",
         n_surr=200,
         seed=123,
         poly_degree=5,
         use_log1pz=True):
    """
    series_path : path to 2-col ASCII file (z, signal)
    outdir      : output directory
    n_surr      : number of surrogates
    seed        : RNG seed
    poly_degree : polynomial degree for de-trend in L
    use_log1pz  : if True, x -> L=log(1+z); else L=z
    """
    np.random.seed(int(seed))
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    print(f"[hht] series file: {series_path}")
    print(f"[hht] outdir     : {outdir}")

    # ---- load series ----
    arr = np.loadtxt(series_path)
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError("Input must be at least 2 columns: z, signal")

    z_raw = arr[:, 0]
    y_raw = arr[:, 1]

    # ---- choose axis: L = log(1+z) or z ----
    if use_log1pz:
        L = np.log1p(z_raw)
        xlabel = "log(1+z)"
    else:
        L = z_raw
        xlabel = "z"

    # sort by L
    order = np.argsort(L)
    L     = L[order]
    y     = y_raw[order]

        # ---- polynomial de-trend in L (NumPy polyfit, no sklearn) ----
    DO_POLY_DETREND = True

    if DO_POLY_DETREND:
        # limit degree so we don't overfit tiny samples
        n_pts     = len(L)
        deg       = int(poly_degree)
        deg       = max(1, min(deg, n_pts - 2))   # at least 1, at most N-2

        coef      = np.polyfit(L, y, deg)
        y_trend   = np.polyval(coef, L)
        resid     = y - y_trend

        detrend_info = {
                        "method"    : "np_polyfit_in_L",
                        "degree"    : deg,
                        "coef"      : coef.tolist(),
                        "intercept" : float(0.0),  # included in coef[-1]
                       }
    else:
        y_trend      = np.zeros_like(y)
        resid        = y.copy()
        detrend_info = {"method": "none"}


    # ---- robust clip on residuals ----
    m_resid   = np.nanmedian(resid)
    s_resid   = np.nanmedian(np.abs(resid - m_resid)) * 1.4826
    good      = np.abs(resid - m_resid) <= 3.5 * max(s_resid, 1e-6)
    keep_frac = float(np.count_nonzero(good)) / float(resid.size)
    
    if keep_frac < 0.6:
        print(f"[clip] skipped (keep_frac={keep_frac:.2f} < 0.60).")
    
    else:
        L       = L[good]
        resid   = resid[good]
        y_trend = y_trend[good]
        y       = y[good]
        print(f"[clip] kept {keep_frac:.2f} of points.")

    # ---- resample to a uniform grid in L ----
    N_uni      = max(512, int(L.size // 1.5))
    Lmin, Lmax = L.min(), L.max()
    L_uni      = np.linspace(Lmin, Lmax, N_uni)

    try:
        from scipy.interpolate import CubicSpline as _CS
        resid_spline = _CS(L, resid)(L_uni)
    
    except Exception:
        resid_spline = np.interp(L_uni, L, resid)

    # ---- cosine taper ----
    M_tap          = max(12, N_uni // 50)
    taper          = np.ones_like(resid_spline)
    w              = 0.5 * (1 - np.cos(np.linspace(0, np.pi, M_tap)))
    taper[:M_tap]  = w
    taper[-M_tap:] = w[::-1]
    rs             = resid_spline * taper

    Ls             = L_uni
    print(f"[grid] uniform {xlabel} bins = {Ls.size}, "
          f"L ∈ [{Ls.min():.4f},{Ls.max():.4f}]")

    # ---- EMD/EEMD ----
    if (not HAS_EMD) and (not HAS_EEMD):
        print("[warn] PyEMD not found; using single-IMF proxy.")
        IMFs = rs[None, :]
    
    else:
        if HAS_EEMD:
            try:
                eemd = _EEMD(trials=100, noise_width=0.2)
                IMFs = eemd.eemd(rs, Ls)
            
            except Exception:
                emd  = _EMD()
                IMFs = emd.emd(rs)
        else:
            emd = _EMD()
            IMFs = emd.emd(rs)

        if IMFs is None or IMFs.size == 0:
            IMFs = rs[None, :]

    # ---- Hilbert: inst. frequency & envelope ----
    inst_freqs, envelopes = [], []
    for k in range(IMFs.shape[0]):
        f, ph, A = _instantaneous_freq_phase(IMFs[k], Ls)
        inst_freqs.append(f)
        envelopes.append(A)
        
    inst_freqs = np.array(inst_freqs, float)
    envelopes  = np.array(envelopes,  float)

    # ---- IMF energy spectrum ----
    energies   = np.sum(IMFs**2, axis=1)
    ranks      = np.argsort(energies)[::-1]

    # ---- surrogate significance on max IMF energy ----
    surr_peaks = []
    for s in range(n_surr):
        ys = _phase_randomize(rs) if (s % 2 == 0) else _cbb_surrogate(rs, block=48)

        if HAS_EEMD:
            try:
                eemd_s = _EEMD(trials=50, noise_width=0.2)
                imfs_s = eemd_s.eemd(ys, Ls)
            except Exception:
                imfs_s = _EMD().emd(ys) if HAS_EMD else ys[None, :]
        
        elif HAS_EMD:
            imfs_s = _EMD().emd(ys)
        
        else:
            imfs_s = ys[None, :]

        if imfs_s is None or imfs_s.size == 0:
            imfs_s = ys[None, :]

        e_s = np.sum(imfs_s**2, axis=1).max()
        surr_peaks.append(e_s)

    surr_peaks = np.array(surr_peaks)
    e_max      = energies.max()
    p_val      = float((np.sum(surr_peaks >= e_max) + 1) / (n_surr + 1))

    # ---- IMF stats ----
    rank_of = np.empty(IMFs.shape[0], dtype=int)
    for r, k in enumerate(ranks, 1):
        rank_of[k] = r

    imf_stats = []
    for k in range(IMFs.shape[0]):
        f     = inst_freqs[k]
        f_med = np.nanmedian(f[np.isfinite(f)]) if np.any(np.isfinite(f)) else np.nan
        imf_stats.append((k + 1, energies[k], f_med, int(rank_of[k])))

    imf_stats = np.array(
                         imf_stats,
                         dtype=[
                                ("IMF", "i4"),
                                ("energy", "f8"),
                                ("f_med", "f8"),
                                ("rank", "i4"),
                               ],
                         )

    # ---- save detailed arrays ----
    detail_npz_path = os.path.join(outdir, "hht_quaia_detail.npz")
    np.savez_compressed(
                        detail_npz_path,
                        Ls=Ls,
                        IMFs=IMFs,
                        inst_freqs=inst_freqs,
                        envelopes=envelopes,
                        energies=energies,
                        ranks=ranks,
                        detrend_info=json.dumps(detrend_info),
                        z_raw=z_raw,
                        y_raw=y_raw,
                       )
    
    print(f"[save] HHT detail -> {detail_npz_path}")

    # ---- save IMF stats ----
    import csv
    csv_path = os.path.join(outdir, "hht_quaia_imf_stats.csv")
    npz_path = os.path.join(outdir, "hht_quaia_imf_stats.npz")

    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
                  [ "IMF",
                    "energy",
                    "median_freq",
                    "rank",
                    "p_val",
                    "series_path",
                    "poly_degree",
                    "use_log1pz"]
                  )
        
        for row in imf_stats:
            w.writerow(
                       [row["IMF"],
                        f"{row['energy']:.10g}",
                        f"{row['f_med']:.10g}",
                        row["rank"],
                        f"{p_val:.6g}",
                        os.path.basename(series_path),
                        poly_degree,
                        int(use_log1pz)]
                      )

    run_id = datetime.now().strftime("%Y%m%d-%H%M%S")
    np.savez_compressed(
                        npz_path,
                        imf_stats   = imf_stats,
                        p_val       = p_val,
                        series_path = series_path,
                        poly_degree=poly_degree,
                        use_log1pz=use_log1pz,
                        run_id=run_id,
                       )
    print(f"[save] IMF stats -> {csv_path}, {npz_path}")

    # ---- plots ----
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 7.8))

    # raw series + trend (on clipped set)
    axs[0, 0].plot(L, y, "o", ms=4, alpha=0.6, label="raw")
    
    if DO_POLY_DETREND:
        axs[0, 0].plot(L, y_trend, "r-", lw=1.0, label="trend")
    
    axs[0, 0].set_xlabel(xlabel)
    axs[0, 0].set_ylabel("signal")
    axs[0, 0].set_title("Quaia series (clipped) + trend")
    axs[0, 0].legend(fontsize=8)

    # detrended / HHT input
    axs[0, 1].plot(Ls, rs, "k-", lw=0.8)
    axs[0, 1].axhline(0, ls=":", lw=1)
    axs[0, 1].set_xlabel(xlabel)
    axs[0, 1].set_ylabel("detrended signal")
    axs[0, 1].set_title("HHT input (uniform grid)")

    # IMF decomposition (offset)
    off = 0.0
    for k in ranks:
        axs[1, 0].plot(Ls, IMFs[k] + off, lw=1.0, label=f"IMF{k+1}")
        off += 0.07
    
    axs[1, 0].set_xlabel(xlabel)
    axs[1, 0].set_ylabel("IMFs (offset)")
    axs[1, 0].set_title("EMD decomposition")
    axs[1, 0].legend(fontsize=7, ncol=2)

    # energy spectrum + surrogate median
    axs[1, 1].bar(np.arange(1, IMFs.shape[0] + 1), energies, width=0.7)
    axs[1, 1].axhline(np.median(surr_peaks), ls="--", lw=1.0, label="surrogate median")
    axs[1, 1].set_xlabel("IMF index")
    axs[1, 1].set_ylabel("energy")
    axs[1, 1].set_title(f"IMF energy spectrum (p≈{p_val:.3f})")
    axs[1, 1].legend(fontsize=8)

    fig.tight_layout()
    png_path = os.path.join(outdir, "hht_quaia.png")
    fig.savefig(png_path, dpi=160, bbox_inches="tight")
    print(f"[files] {png_path}")
    plt.close(fig)

    # quick summary
    np.savez_compressed(
                        os.path.join(outdir, "hht_quaia_summary.npz"),
                        Ls=Ls,
                        resid=rs,
                        energies=energies,
                        ranks=ranks,
                        p_val=p_val,
                       )


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: hht_quaia_residuals.py <series_2col.txt> [outdir]")
        sys.exit(1)

    series_path = sys.argv[1]
    outdir      = sys.argv[2] if len(sys.argv) > 2 else "plamb_runs/hht_quaia"
    main(series_path, outdir=outdir)
