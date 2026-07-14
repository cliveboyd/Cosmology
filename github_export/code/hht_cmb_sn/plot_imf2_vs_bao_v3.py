
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
IMF2 & BAO plotting utility (v3)
--------------------------------
Adds an Energy vs z figure:
  - Computes IMF2 instantaneous amplitude via Hilbert transform and defines
    energy E(x) = |A(x)|^2 on x = log(1+z).
  - Plots raw normalized energy and a binned-median curve.
  - Optionally overlays BAO z positions as vertical dashed lines.

Also retains:
  - IMF2 oscillations vs log(1+z)
  - IMF2 instantaneous frequency vs BAO 1/k proxy (normalized)
"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt

# ========== User settings ==========
# Update these paths before running:
NPZ_PATH = r"/Users/boyde/.spyder-py3/plamb_runs/hht_sn_FRPBH/hht_sn_detail.npz"
BAO_CSV  = r"/Users/boyde/Downloads/bao_long.csv"

# Optional features
PLOT_FREQ_OVERLAY = True      # make IMF2 freq vs BAO 1/k proxy
PLOT_ENERGY       = True      # make IMF2 energy vs log(1+z)
SHOW_BAO_MARKERS  = True      # draw BAO vertical lines on energy plot
NBINS_ENERGY      = 40        # bins for median energy curve

# Output folder for figures
OUTDIR = os.path.join(os.path.dirname(__file__) if '__file__' in globals() else os.getcwd(),
                      "imf2_figs_out")
os.makedirs(OUTDIR, exist_ok=True)

# ========== Helpers ==========
def load_hht_npz(path):
    npz = np.load(path, allow_pickle=True)
    # IMFs
    imfs = None
    for k in ["imfs", "IMFs", "imf"]:
        if k in npz:
            imfs = np.array(npz[k])
            break
    if imfs is None:
        raise KeyError("Could not find IMFs in NPZ (tried keys: imfs, IMFs, imf)." )

    # x-axis (log1pz expected)
    x = None
    for k in ["log1pz", "x", "t", "log_1pz"]:
        if k in npz:
            x = np.array(npz[k])
            break
    if x is None:
        N = imfs.shape[-1]
        x = np.linspace(0.0, 1.0, N)

    # instantaneous frequency, if present
    instfreq = None
    for k in ["inst_freq", "omega", "instantaneous_frequency", "w"]:
        if k in npz:
            instfreq = np.array(npz[k])
            break
    return x, imfs, instfreq

def get_imf_index_by_energy(imfs, prefer_index=1):
    n = imfs.shape[0]
    idx = min(max(prefer_index, 0), n-1)
    return idx

def compute_instfreq_hilbert(x, imf, smooth=True):
    """Compute instantaneous frequency from IMF via Hilbert transform."""
    try:
        from scipy.signal import hilbert, savgol_filter
        analytic = hilbert(imf)
        phase = np.unwrap(np.angle(analytic))
        omega = np.gradient(phase, x)
        if smooth and len(omega) >= 21:
            win = 21 if len(omega) >= 21 else (len(omega)//2)*2 + 1
            omega = savgol_filter(omega, window_length=win, polyorder=3, mode="interp")
        return omega, True, None
    except Exception as e:
        return None, False, str(e)

def norm_series(v):
    v = np.array(v, dtype=float)
    m = np.nanmean(v)
    s = np.nanstd(v)
    out = (v - m) / s if s and np.isfinite(s) and s > 0 else v - m
    return out

def load_bao_csv(path):
    """Expect columns: z, kind, value, sigma. Creates inv_k_eff_proxy = 1/value."""
    import pandas as pd
    df = pd.read_csv(path)
    df["inv_k_eff_proxy"] = 1.0 / df["value"]
    return df

def hilbert_amplitude(imf):
    """Return instantaneous amplitude |analytic_signal| with light smoothing if available."""
    try:
        from scipy.signal import hilbert, savgol_filter
        analytic = hilbert(imf)
        amp = np.abs(analytic)
        if len(amp) >= 51:
            win = 51 if len(amp) >= 51 else (len(amp)//2)*2 + 1
            amp = savgol_filter(amp, window_length=win, polyorder=3, mode="interp")
        return amp, None
    except Exception as e:
        # Fallback RMS window if SciPy unavailable
        w = max(11, int(0.02*len(imf))//2*2 + 1)
        pad = w//2
        xpad = np.pad(imf, (pad, pad), mode="reflect")
        amp = np.sqrt(np.convolve(xpad**2, np.ones(w)/w, mode="valid"))
        return amp, f"[warn] SciPy hilbert unavailable: used RMS window {w} as amplitude proxy."

def binned_stat(x, y, nbins=40, stat="median"):
    x = np.asarray(x); y = np.asarray(y)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    bins = np.linspace(x.min(), x.max(), nbins+1)
    xc = 0.5*(bins[:-1] + bins[1:])
    out = np.full_like(xc, np.nan, dtype=float)
    for i in range(nbins):
        sel = (x >= bins[i]) & (x < bins[i+1])
        if np.any(sel):
            out[i] = np.median(y[sel]) if stat == "median" else np.mean(y[sel])
    return xc, out

# ========== Main ==========
def main():
    try:
        x, imfs, instfreq = load_hht_npz(NPZ_PATH)
    except Exception as e:
        print(f"[error] Failed to load NPZ: {e}")
        sys.exit(1)

    # Choose IMF2 (index 1, zero-based)
    imf_idx = get_imf_index_by_energy(imfs, prefer_index=1)
    imf2 = imfs[imf_idx]

    # --- Figure 1: IMF2 oscillations ---
    plt.figure(figsize=(7.5, 4.5), dpi=300)
    plt.plot(x, imf2, label="IMF2")
    plt.xlabel("log(1 + z)")
    plt.ylabel("IMF2 amplitude (arb. units)")
    plt.title("IMF2 Oscillations vs log(1+z)")
    plt.legend()
    plt.tight_layout()
    f1_png = os.path.join(OUTDIR, "imf2_oscillations.png")
    f1_svg = os.path.join(OUTDIR, "imf2_oscillations.svg")
    plt.savefig(f1_png); plt.savefig(f1_svg); plt.close()
    print(f"[save] {f1_png}"); print(f"[save] {f1_svg}")

    # --- Figure 2: IMF2 instantaneous frequency vs BAO proxy ---
    if PLOT_FREQ_OVERLAY:
        have_freq = False
        freq = None
        if instfreq is not None and instfreq.ndim >= 2 and instfreq.shape[0] > imf_idx:
            freq = instfreq[imf_idx]; have_freq = True
            print("[info] Using instantaneous frequency from NPZ.")
        else:
            freq, ok, err = compute_instfreq_hilbert(x, imf2, smooth=True)
            if ok:
                have_freq = True
                print("[info] Computed instantaneous frequency via Hilbert transform.")
            else:
                print(f"[warn] Could not compute instantaneous frequency: {err}")
        if have_freq:
            try:
                import pandas as pd
                df = load_bao_csv(BAO_CSV)
                if (df["z"] < -1e-6).any():
                    raise ValueError("BAO z values appear invalid (negative)." )
                df["x_from_z"] = np.log(1.0 + df["z"])
                good = np.isfinite(df["x_from_z"]) & np.isfinite(df["inv_k_eff_proxy"])
                zx = df.loc[good, "x_from_z"].values
                yk = df.loc[good, "inv_k_eff_proxy"].values
                if len(zx) > 3:
                    order = np.argsort(zx); zx = zx[order]; yk = yk[order]
                    yk_interp = np.interp(x, zx, yk, left=np.nan, right=np.nan)
                else:
                    yk_interp = np.full_like(x, np.nan)
                y1 = norm_series(freq); y2 = norm_series(yk_interp)
                plt.figure(figsize=(7.5, 4.5), dpi=300)
                plt.plot(x, y1, label="IMF2 instantaneous frequency (normalized)")
                plt.plot(x, y2, label="BAO 1/k_eff proxy (normalized)")
                plt.xlabel("log(1 + z)\"); plt.ylabel("Normalized units")
                plt.title("IMF2 Instantaneous Frequency vs BAO 1/k_eff (proxy)")
                plt.legend(); plt.tight_layout()
                f2_png = os.path.join(OUTDIR, "imf2_freq_vs_bao.png")
                f2_svg = os.path.join(OUTDIR, "imf2_freq_vs_bao.svg")
                plt.savefig(f2_png); plt.savefig(f2_svg); plt.close()
                print(f"[save] {f2_png}"); print(f"[save] {f2_svg}")
            except Exception as e:
                print(f"[warn] Skipping frequency vs BAO overlay: {e}")

    # --- Figure 3: IMF2 energy vs log(1+z) ---
    if PLOT_ENERGY:
        amp, note = hilbert_amplitude(imf2)
        if note: print(note)
        energy = amp**2
        # Normalize for plotting
        eplot = (energy - np.nanmin(energy)) / (np.nanmax(energy) - np.nanmin(energy) + 1e-12)
        xb, eb = binned_stat(x, eplot, nbins=NBINS_ENERGY, stat="median")

        # Optional BAO markers
        x_bao = None
        if SHOW_BAO_MARKERS:
            try:
                import pandas as pd
                dfb = load_bao_csv(BAO_CSV)
                x_bao = np.log1p(dfb["z"].values)
                x_bao = x_bao[np.isfinite(x_bao)]
            except Exception:
                x_bao = None

        plt.figure(figsize=(7.5, 4.5), dpi=300)
        plt.plot(x, eplot, alpha=0.5, label="IMF2 energy (normalized, instantaneous)")
        plt.plot(xb, eb, linewidth=2.0, label=f"IMF2 energy (binned median, n={NBINS_ENERGY})")
        if x_bao is not None:
            for xv in x_bao:
                plt.axvline(x=xv, linestyle="--", linewidth=0.8)
            plt.text(0.01, 0.95, "BAO z positions shown as dashed lines", transform=plt.gca().transAxes)

        plt.xlabel("log(1 + z)\"); plt.ylabel("Normalized energy")
        plt.title("IMF2 Energy vs log(1+z)")
        plt.legend(); plt.tight_layout()
        f3_png = os.path.join(OUTDIR, "imf2_energy_vs_log1pz.png")
        f3_svg = os.path.join(OUTDIR, "imf2_energy_vs_log1pz.svg")
        plt.savefig(f3_png); plt.savefig(f3_svg); plt.close()
        print(f"[save] {f3_png}"); print(f"[save] {f3_svg}")

if __name__ == "__main__":
    main()
