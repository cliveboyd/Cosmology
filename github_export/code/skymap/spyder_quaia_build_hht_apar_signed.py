#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_build_hht_apar_signed.py  (0V1)

Run a 1D HHT (CEEMDAN + Hilbert) on the *signed* A_par(z) series
for all |b| cuts, using the files produced by
`quaia_build_hht_apar_signed.py`.

Outputs (for each bcut, into plamb_runs/hht_quaia/bcutXX_allz_apar_signed):

  - hht_quaia_detail.npz   (Ls, IMFs, inst_freqs, envelopes, energies, ranks, z_raw, y_raw)
  - hht_quaia_imf_stats.npz
  - hht_quaia_imf_stats.csv
"""

from   pathlib      import Path
import numpy        as     np

from   PyEMD        import CEEMDAN
from   scipy.signal import hilbert
from   math         import pi

# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------

BASE_IN                = Path("/Users/boyde/.spyder-py3/skymap2/quaia_outputs")
BASE_OUT               = Path("/Users/boyde/.spyder-py3/plamb_runs/hht_quaia")

BCUTS                  = [10, 15, 20, 25, 30, 35]
N_GRID                 = 512   # uniform bins in L = log(1+z)

# ----------------------------------------------------------------------
# Core HHT runner
# ----------------------------------------------------------------------

def run_hht_series(series_path, outdir):
    series_path        = Path(series_path)
    outdir             = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("[hht] series file:", series_path)
    print("[hht] outdir     :", outdir)

    # Load (z_mid, A_par_signed)
    z, y               = np.loadtxt(series_path, unpack=True)

    # Sort by z
    order              = np.argsort(z)
    z                  = z[order]
    y                  = y[order]

    # Map to uniform grid in L = log(1+z)
    L                  = np.log1p(z)
    L_min              = L.min()
    L_max              = L.max()
    L_grid             = np.linspace(L_min, L_max, N_GRID)
    y_grid             = np.interp(L_grid, L, y)

    print(f"[grid] uniform log(1+z) bins = {N_GRID}, "
          f"L ∈ [{L_min:.4f},{L_max:.4f}]")

    # CEEMDAN decomposition
    ce                 = CEEMDAN()
    IMFs               = ce(y_grid)          # shape (n_imf, N_GRID)
    n_imf, N           = IMFs.shape

    # Hilbert-based instantaneous frequency + envelopes
    inst_freqs         = np.zeros_like(IMFs)
    envelopes          = np.zeros_like(IMFs)
    dL                 = L_grid[1] - L_grid[0]

    for k in range(n_imf):
        analytic       = hilbert(IMFs[k])
        envelopes[k]   = np.abs(analytic)
        phase          = np.unwrap(np.angle(analytic))
        dphase         = np.gradient(phase, dL)
        inst_freqs[k]  = dphase / (2.0 * pi)

    # Energies and median frequency per IMF
    energies           = np.sum(IMFs**2, axis=1)
    f_med              = np.median(np.abs(inst_freqs), axis=1)

    # Rank: 1 = lowest-frequency (trend), larger rank = higher freq
    order_f            = np.argsort(f_med)        # small f_med -> low freq
    ranks              = np.empty_like(order_f)
    ranks[order_f]     = np.arange(1, n_imf + 1)

    # Structured IMF stats
    imf_dtype          = np.dtype([
                                   ("IMF",    np.int32),
                                   ("energy", np.float64),
                                   ("f_med",  np.float64),
                                   ("rank",   np.int32) ]    
                                 )
    imf_stats           = np.zeros(n_imf, dtype=imf_dtype)
    imf_stats["IMF"]    = np.arange(1, n_imf + 1)
    imf_stats["energy"] = energies
    imf_stats["f_med"]  = f_med
    imf_stats["rank"]   = ranks

    # Save detail NPZ (matches what your diagnostic overlay expects)
    detail_path = outdir / "hht_quaia_detail.npz"
    np.savez(
            detail_path,
            Ls              = L_grid,
            IMFs            = IMFs,
            inst_freqs      = inst_freqs,
            envelopes       = envelopes,
            energies        = energies,
            ranks           = ranks,
            detrend_info    = "none",
            z_raw           = z,
            y_raw           = y,
          )
    print("[save] HHT detail ->", detail_path)

    # Save stats NPZ + CSV
    stats_npz = outdir / "hht_quaia_imf_stats.npz"
    np.savez(
        stats_npz,
        imf_stats   = imf_stats,
        p_val       = np.nan,
        series_path = str(series_path),
        poly_degree = 0,
        use_log1pz  = True,
        run_id      = "apar_signed",
    )

    stats_csv = outdir / "hht_quaia_imf_stats.csv"
    header    = "IMF,energy,f_med,rank"
    np.savetxt(
        stats_csv,
        np.column_stack([
            imf_stats["IMF"],
            imf_stats["energy"],
            imf_stats["f_med"],
            imf_stats["rank"],
        ]),
        delimiter = ",",
        header    = header,
        comments  = "",
    )
    print("[save] IMF stats ->", stats_csv, ",", stats_npz)


# ----------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------

def main():
    for bcut in BCUTS:
        series = BASE_IN / f"quaia_hht_bcut{bcut}_allz_apar_signed.txt"
        outdir = BASE_OUT / f"bcut{bcut}_allz_apar_signed"
        print(f"\nRunning A_par(signed) for bcut={bcut}°")
        run_hht_series(series, outdir)


if __name__ == "__main__":
    main()
