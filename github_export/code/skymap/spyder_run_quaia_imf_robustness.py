#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_quaia_imf_robustness.py

Grid-run quaia_imf_diagnostics.py, select IMF by period band, compute chirp metrics,
compute quick phase-randomisation null p (screen only), compute harmonic phase diagnostics,
and compute ENSEMBLE-level (stacked across seeds) p_shift (aligned) and p_align (raw).
"""

import sys
import shutil
import itertools
import subprocess
from   pathlib           import  Path

import numpy             as      np
import pandas            as      pd
import matplotlib.pyplot as      plt
from   scipy.signal      import  hilbert, savgol_filter


# =========================
# USER SETTINGS
# =========================
SCRIPT                      = Path("/Users/boyde/.spyder-py3/quaia_imf_diagnostics.py")
FITS                        = Path("/Users/boyde/Downloads/quaia_G20.0.fits")
BASE                        = Path("/Users/boyde/.spyder-py3/quaia_imf_grid")

BASE.mkdir(parents = True, exist_ok = True)

SEEDS                       = [  7,   77,  123]
CE_TRIALS                   = [500, 1000, 2000]
NBINS                       = [360,  720, 1080]
BCUTS                       = [10,    20,   30, 40]

SMOOTH                      = 7
NOISE                       = 0.2
MAX_IMFS                    = 14
XCOORD                      = "galactic_l"
SIGNAL                      = "counts"

# --- IMF selection by period band ---
TARGET_PERIOD_DEG           = (140.0, 240.0)

# --- nulls ---
PHASE_NULL_NSURR            = 300      # quick screen only (floors at 1/(N+1))
ENSEMBLE_SHIFT_NULL         = 10000    # ensemble-level circular-shift null

# --- direction/phase settings ---
K_HARM                      = 2
TARGET_DIR_DEG              = 264.0         # compare in k=2 space => mod 180 internally

EDGE_TRIM_DEG               = 20.0


# Overide Options........
BCUTS         = [20]
NBINS         = [360]      
CE_TRIALS     = [2000]
SEEDS         = [7, 77, 123]

CHIRP_KW                    = dict(
                                   edge_trim_deg = EDGE_TRIM_DEG,
                                   amp_floor_q   = 0.15,
                                   smooth_ifreq  = True,
                                   sg_window     = 21,
                                   sg_poly       = 2,
                                  )

CHIRP_R2_MIN                = 0.15   # optional; tune (0.0 to disable)


REQUIRED_FILES = ["signal.csv", "imfs.npy", "residue.npy"]

def write_run_logs(outdir: Path, stdout: str, stderr: str) -> None:
    (outdir / "run_stdout.txt").write_text(stdout or "", encoding = "utf-8")
    (outdir / "run_stderr.txt").write_text(stderr or "", encoding = "utf-8")
    

def assert_run_artifacts(outdir: Path) -> None:
    present                 = sorted([p.name for p in outdir.glob("*")])

    # common “wrong name” possibilities we should detect explicitly
    alt_imfs                = [p.name for p in outdir.glob("imf*.*")]
    alt_npys                = [p.name for p in outdir.glob("*.npy")]

    missing                 = [f for f in REQUIRED_FILES if not (outdir / f).exists()]
    if missing:
        raise FileNotFoundError(
                                f"{outdir} missing {missing}.\n"
                                f"Present: {present}\n"
                                f"Alt imf-ish files: {alt_imfs}\n"
                                f"All npy: {alt_npys}\n"
                                f"Check logs:\n  {outdir/'run_stdout.txt'}\n  {outdir/'run_stderr.txt'}"
                               )

    imfs                    = np.load(outdir / "imfs.npy")
    
    if imfs.ndim != 2 or imfs.shape[0] < 1 or imfs.shape[1] < 3:
        raise ValueError(f"{outdir}/imfs.npy looks wrong: shape={imfs.shape}")

# =========================
# HELPERS
# =========================
def run_one(outdir: Path, bcut: int, nbins: int, ce_trials: int, seed: int) -> None:
    outdir.mkdir(parents = True, exist_ok = True)

    cmd = [
            sys.executable, str(SCRIPT),
            "--fits",      str(FITS),
            "--out",       str(outdir),
            "--xcoord",    XCOORD,
            "--bcut",      str(bcut),
            "--nbins",     str(nbins),
            "--signal",    SIGNAL,
            "--smooth",    str(SMOOTH),
            "--ce-trials", str(ce_trials),
            "--noise",     str(NOISE),
            "--seed",      str(seed),
            "--max-imfs",  str(MAX_IMFS),
            "--dpi",       "180",
           ]

    p                       = subprocess.run(cmd, capture_output=True, text=True)
    write_run_logs(outdir, p.stdout, p.stderr)

    if p.returncode != 0:
        raise RuntimeError(
                            f"Run failed (bcut={bcut}, nbins={nbins}, ce={ce_trials}, seed={seed})\n"
                            f"See logs:\n  {outdir/'run_stdout.txt'}\n  {outdir/'run_stderr.txt'}"
                           )

    # CRITICAL: prevent “success” without imfs.npy
    assert_run_artifacts(outdir)

def select_imf_best_in_band(imfs: np.ndarray, x: np.ndarray, band: tuple[float, float],
                            chirp_kw: dict, r2_min: float = 0.0):
    
    per                         = period_per_imf(imfs, x)
    lo, hi                      = band
    cand = np.where((per >= lo) & (per <= hi))[0]

    # fallback: closest to band midpoint
    mid                         = 0.5*(lo + hi)
    if cand.size == 0:
        k0                      = int(np.nanargmin(np.abs(per - mid)))
        return k0, per

    best_k                      = None
    best_score                  = -np.inf

    for k in cand:
        met                     = imf_chirp_metrics(x, imfs[k], **chirp_kw)
        sc                      = met["score"]
        r2                      = met["r2"]

        if not np.isfinite(sc):
            continue
        
        if np.isfinite(r2) and r2 < r2_min:
            continue

        if sc > best_score:
            best_score          = sc
            best_k              = int(k)

    if best_k is None:
        # nothing passed r2/finite filters → revert to closest-to-mid
        best_k                  = int(cand[np.argmin(np.abs(per[cand] - mid))])

    return best_k, per


def load_outputs(outdir: Path):
    sig                         = pd.read_csv(outdir / "signal.csv")
    imfs                        = np.load(outdir     / "imfs.npy")
    res                         = np.load(outdir     / "residue.npy")

    if "x_deg" not in sig.columns:
        raise ValueError(f"{outdir}/signal.csv missing x_deg column: {sig.columns.tolist()}")

    if "signal" in sig.columns:
        ycol                    = "signal"
    else:
        cands                   = [c for c in sig.columns if c != "x_deg"]
        
        if not cands:
            raise ValueError(f"{outdir}/signal.csv has no signal column candidates")
        ycol                    = cands[0]

    x                           = sig["x_deg"].to_numpy()
    y                           = sig[ycol].to_numpy()

    if imfs.ndim != 2:
        raise ValueError(f"imfs.npy expected 2D, got {imfs.shape}")

    if res.shape[0] != y.shape[0] or imfs.shape[1] != y.shape[0]:
        raise ValueError(
                         f"Length mismatch in {outdir}: signal={y.shape}, imfs={imfs.shape}, res={res.shape}"
                        )
    return x, y, imfs, res


def est_period_zero_cross(y: np.ndarray, x: np.ndarray):
    zc                          = np.where(np.signbit(y[:-1]) != np.signbit(y[1:]))[0]  
    
    if len(zc) < 4:
        return np.nan
    dx                          = np.diff(x[zc])
    return 2.0 * np.nanmean(dx)


def period_per_imf(imfs: np.ndarray, x: np.ndarray):
    return np.array([est_period_zero_cross(imfs[k], x) for k in range(imfs.shape[0])], dtype = float)


def select_imf_by_period(imfs: np.ndarray, x: np.ndarray, band: tuple[float, float]):
    per                         = period_per_imf(imfs, x)
    lo, hi                      = band
    ok                          = np.where((per >= lo) & (per <= hi))[0]
    mid                         = 0.5 * (lo + hi)
    
    if len(ok) == 0:
        k0                      = int(np.nanargmin(np.abs(per - mid)))
        return k0, per
    
    k0                          = int(ok[np.argmin(np.abs(per[ok] - mid))])
    return k0, per


def circ_dist(x, y, period):
    d                           = (x - y + period/2) % period - period/2
    return float(abs(d))


def fit_harmonic_phase(l_deg, y, k = 2, include_const = True):
    theta                       = np.deg2rad(l_deg)
    cosk                        = np.cos(k * theta)
    sink                        = np.sin(k * theta)

    if include_const:
        X                       = np.column_stack([np.ones_like(theta), cosk, sink])
    else:
        X                       = np.column_stack([cosk, sink])

    beta, *_                    = np.linalg.lstsq(X, y, rcond = None)
    
    if include_const:
        c, a, b                 = beta
        
    else:
        c, a, b                 = 0.0, beta[0], beta[1]

    amp                         = float(np.hypot(a, b))
    phi                         = float(np.arctan2(b, a))
    l0                          = (np.rad2deg(phi) / k) % (360.0 / k)
    
    return dict(c = float(c), a = float(a), b = float(b), amp = amp, phi_rad = phi, l0_deg = float(l0))


def imf_chirp_metrics(x_deg             : np.ndarray,
                      imf               : np.ndarray,
                      edge_trim_deg     : float = 20.0,
                      amp_floor_q       : float = 0.15,
                      smooth_ifreq      : bool  = True,
                      sg_window         : int   = 21,
                      sg_poly           : int   = 2):
    """
    Chirp metrics using Hilbert instantaneous frequency, but robustified by:
      - edge trim
      - amplitude masking (drop low-amplitude points where inst-freq is unstable)
      - optional Savitzky–Golay smoothing of inst-freq
    """

    dx                          = float(np.median(np.diff(x_deg)))
    z                           = hilbert(imf)
    A                           = np.abs(z)
    ph                          = np.unwrap(np.angle(z))

    inst_f                      = np.diff(ph) / (2*np.pi*dx)             # cycles/deg
    x_mid                       = 0.5 * (x_deg[:-1] + x_deg[1:])
    A_mid                       = 0.5 * (A[:-1] + A[1:])                 # amplitude aligned to inst_f samples

    # edge trim mask
    m_edge                      = (x_mid >= edge_trim_deg) & (x_mid <= (360.0 - edge_trim_deg))

    # amplitude mask (drop bottom amp_floor_q quantile inside edge window)
    A_use                       = A_mid[m_edge]
    
    if A_use.size < 10:
        return dict(m=np.nan, r2=np.nan, f_med=np.nan, f_iqr=np.nan, score=np.nan)

    thr                         = np.quantile(A_use, amp_floor_q)
    m_amp                       = A_mid >= thr

    msk                         = m_edge & m_amp
    
    if msk.sum() < 20:
        return dict(m=np.nan, r2=np.nan, f_med=np.nan, f_iqr=np.nan, score=np.nan)

    xm                          = x_mid[msk]
    fm                          = inst_f[msk]

    # optional smoothing (helps reduce spiky inst-freq without destroying broad trend)
    if smooth_ifreq and fm.size >= sg_window and sg_window >= 5 and (sg_window % 2 == 1):
        fm = savgol_filter(fm, window_length=sg_window, polyorder=sg_poly, mode="interp")

    # fit slope
    m, c                        = np.polyfit(xm, fm, 1)
    fhat                        = m*xm + c

    ss_res                      = np.sum((fm - fhat)**2)
    ss_tot                      = np.sum((fm - np.mean(fm))**2)
    r2                          = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else np.nan

    f_med                       = float(np.median(fm))
    f_iqr                       = float(np.percentile(fm, 75) - np.percentile(fm, 25))
    score                       = float(np.abs(m) / (f_iqr + 1e-12))

    return dict(m=float(m), r2=float(r2), f_med=f_med, f_iqr=f_iqr, score=score)






def phase_randomised_surrogate(y: np.ndarray, rng: np.random.Generator):
    Y                           = np.fft.rfft(y)
    mag                         = np.abs(Y)

    rand_ph                     = rng.uniform(-np.pi, np.pi, size = len(Y))
    rand_ph[0]                  = 0.0
    
    if (len(y) % 2) == 0:
        rand_ph[-1]             = 0.0

    Y2                          = mag * np.exp(1j * rand_ph)
    
    return np.fft.irfft(Y2, n = len(y))


def chirp_pvalue_phase_null(x_deg, imf, nsurr=300, seed=0):
    rng                         = np.random.default_rng(seed)
    obs                         = imf_chirp_metrics(x_deg, imf, **CHIRP_KW)["score"]
    
    if not np.isfinite(obs):
        return np.nan

    scores = []
    for _ in range(nsurr):
        ys                      = phase_randomised_surrogate(imf, rng)
        sc                      = imf_chirp_metrics(x_deg, ys, **CHIRP_KW)["score"]
        
        if np.isfinite(sc):
            scores.append(sc)

    scores                      = np.array(scores, dtype=float)
    if scores.size < 20:
        return np.nan

    return float((np.sum(scores >= obs) + 1.0) / (scores.size + 1.0))


def best_circular_alignment(ref: np.ndarray, y: np.ndarray):
    a                           = ref - ref.mean()
    b                           = y   - y.mean()
    fa                          = np.fft.rfft(a)
    fb                          = np.fft.rfft(b)
    cc                          = np.fft.irfft(fa * np.conj(fb), n = len(a))
    k                           = int(np.argmax(cc))
    return k, float(cc[k])


def pvalue_circular_shift_null(y, x, n_shifts, rng):
    obs                         = imf_chirp_metrics(x, y, **CHIRP_KW)["score"]
    
    if not np.isfinite(obs):
        return np.nan, obs

    n                           = len(y)
    ge                          = 0
    
    for _ in range(n_shifts):
        k                       = int(rng.integers(0, n))
        ys                      = np.roll(y, k)
        sc                      = imf_chirp_metrics(x, ys, **CHIRP_KW)["score"]
        
        if np.isfinite(sc) and sc >= obs:
            ge += 1

    p                           = (ge + 1.0) / (n_shifts + 1.0)
    return float(p), float(obs)


def pvalue_phase_alignment_shift_null(y: np.ndarray, x_deg: np.ndarray, k: int, target_deg: float, n_shifts: int, rng: np.random.Generator):
    period                      = 360.0 / k
    obs_fit                     = fit_harmonic_phase(x_deg, y, k=k)
    obs_l0                      = obs_fit["l0_deg"] % period
    obs_d                       = circ_dist(obs_l0, target_deg % period, period)

    n                           = len(y)
    ge                          = 0
    for _ in range(n_shifts):
        shift                   = int(rng.integers(0, n))
        ys                      = np.roll(y, shift)
        l0s                     = fit_harmonic_phase(x_deg, ys, k=k)["l0_deg"] % period
        ds                      = circ_dist(l0s, target_deg % period, period)
        if ds <= obs_d:
            ge += 1

    p                           = (ge + 1.0) / (n_shifts + 1.0)
    
    return float(p), float(obs_l0), float(obs_d)


def quick_plot_ensemble(x, y_mean, imf_raw_mean, imf_al_mean, outpng: Path, title: str):
    plt.figure(figsize=(10, 5))
    plt.plot(x, y_mean,       label="mean signal",        linewidth = 2)
    plt.plot(x, imf_raw_mean, label="mean IMF (raw)",     alpha     = 0.55)
    plt.plot(x, imf_al_mean,  label="mean IMF (aligned)", linewidth = 2)
    plt.axhline(0, linestyle = "--", alpha = 0.4)
    plt.title(title)
    plt.xlabel("x_deg")
    plt.ylabel("counts / IMF units")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpng, dpi = 180)
    plt.close()


# =========================
# MAIN GRID RUN
# =========================
def main():
    if not SCRIPT.exists():
        raise FileNotFoundError(f"Script not found: {SCRIPT}")
        
    if not FITS.exists():
        raise FileNotFoundError(f"FITS not found: {FITS}")

    results                     = []
    runs                        = list(itertools.product(BCUTS, 
                                                         NBINS, 
                                                         CE_TRIALS, 
                                                         SEEDS))
    
    print(f"[grid] total runs: {len(runs)} -> output base: {BASE}")

    for i, (bcut, nbins, ce, seed) in enumerate(runs, start=1):
        tag                     = f"b{bcut}_n{nbins}_ce{ce}_s{seed}"
        outdir                  = BASE / tag

        print(f"[{i}/{len(runs)}] running {tag}")
        try:
            run_one(outdir, bcut=bcut, nbins=nbins, ce_trials=ce, seed=seed)
            x, y, imfs, res     = load_outputs(outdir)

            k0, per             = select_imf_best_in_band(imfs, x, TARGET_PERIOD_DEG, CHIRP_KW, r2_min=CHIRP_R2_MIN)
            imf                 = imfs[k0]
            met                 = imf_chirp_metrics(x, imf, **CHIRP_KW)

            picked_period_deg   = float(per[k0]) if np.isfinite(per[k0]) else np.nan

            

            p_null               = chirp_pvalue_phase_null(
                                                           x, 
                                                           imf,
                                                           edge_trim_deg = EDGE_TRIM_DEG,
                                                           nsurr         = PHASE_NULL_NSURR,
                                                           seed          = seed
                                                          )

            period_h            = 360.0 / K_HARM
            target_h            = TARGET_DIR_DEG % period_h
            hf                  = fit_harmonic_phase(x, imf, k = K_HARM)
            l0                  = hf["l0_deg"] % period_h
            d_to_target         = circ_dist(l0, target_h, period_h)
            amp_h               = hf["amp"]

            results.append({
                            "tag"                   : tag,
                            "bcut"                  : bcut,
                            "nbins"                 : nbins,
                            "ce_trials"             : ce,
                            "seed"                  : seed,
                            "n_imfs"                : imfs.shape[0],
                            "picked_imf_idx_1based" : int(k0 + 1),
                            "picked_period_deg"     : picked_period_deg,
            
                            "chirp_slope_m"         : met["m"],
                            "chirp_r2"              : met["r2"],
                            "f_med"                 : met["f_med"],
                            "f_iqr"                 : met["f_iqr"],
                            "chirp_score"           : met["score"],
                            "p_phase_null"          : p_null,
            
                            "harm_k"                : int(K_HARM),
                            "harm_period_deg"       : float(period_h),
                            "harm_l0_deg_modP"      : float(l0),
                            "harm_d_to_target_deg"  : float(d_to_target),
                            "harm_amp"              : float(amp_h),
            
                            "outdir"                : str(outdir),
                          })

        except Exception as e:
            print(f"[WARN] {tag} failed: {e}")
        
            # if artifacts missing, remove the partial folder so reruns are clean
            if "missing" in str(e) and outdir.exists():
                shutil.rmtree(outdir,  ignore_errors = True)
        
            results.append({
                            "tag"       : tag, 
                            "bcut"      : bcut, 
                            "nbins"     : nbins, 
                            "ce_trials" : ce, 
                            "seed"      : seed,
                            "error"     : str(e), 
                            "outdir"    : str(outdir),
                        })


    df                          = pd.DataFrame(results)
    out_csv                     = BASE / "grid_metrics.csv"
    df.to_csv(out_csv, index=  False)
    print(f"[grid] wrote metrics: {out_csv}")

    # =========================
    # ENSEMBLES (same bcut/nbins/ce, average over seeds)
    # =========================
    ok                          = df[df["error"].isna()].copy() if "error" in df.columns else df.copy()
    
    if ok.empty:
        print("[ensemble] no successful runs, skipping.")
        return

    ensemble_rows = []
    group_cols                  = ["bcut", "nbins", "ce_trials"]

    for (bcut, nbins, ce), g in ok.groupby(group_cols):
        xs, ys, picked          = [], [], []

        for _, row in g.iterrows():
            outdir              = Path(row["outdir"])
            x, y, imfs, _       = load_outputs(outdir)

            k0                  = int(row["picked_imf_idx_1based"]) - 1
            k0                  = max(0, min(k0, imfs.shape[0] - 1))

            xs.append(x)
            ys.append(y)
            picked.append(imfs[k0])

        # Require identical x grids
        x0                      = xs[0]
        
        if not all(np.allclose(x0, xi) for xi in xs[1:]):
            print(f"[ensemble] grid mismatch for b{bcut}_n{nbins}_ce{ce}, skipping.")
            continue

        y_mean                  = np.mean(np.vstack(ys), axis=0)

        Yraw                    = np.vstack(picked)
        imf_raw_mean            = Yraw.mean(axis = 0)

        # aligned mean (shape/chirp)
        ref                     = Yraw[0]
        aligned                 = [ref]
        shifts                  = [0]
        
        for j in range(1, Yraw.shape[0]):
            sh, _               = best_circular_alignment(ref, Yraw[j])
            
            aligned.append(np.roll(Yraw[j], sh))
            shifts.append(sh)

        Yal                     = np.vstack(aligned)
        imf_al_mean             = Yal.mean(axis = 0)

        ens_dir                 = BASE / f"ENSEMBLE_b{bcut}_n{nbins}_ce{ce}"
        ens_dir.mkdir(parents   = True, exist_ok = True)

        # Save ensemble products
        pd.DataFrame({"x_deg": x0, "signal_mean": y_mean}).to_csv(ens_dir / "ensemble_signal.csv", index=False)
        np.save(ens_dir / "ensemble_imf_mean_raw.npy",     imf_raw_mean)
        np.save(ens_dir / "ensemble_imf_mean_aligned.npy", imf_al_mean)

        quick_plot_ensemble(
                            x0, y_mean, 
                            imf_raw_mean, 
                            imf_al_mean,
                            outpng = ens_dir / "ensemble_signal_imf_raw_vs_aligned.png",
                            title  = f"Ensemble (bcut={bcut}, nbins={nbins}, ce={ce})"
                           )

        # Ensemble metrics / p-values
        met_raw                 = fit_harmonic_phase(x0, imf_raw_mean, k=K_HARM)
        period_h                = 360.0 / K_HARM
        l0_raw                  = met_raw["l0_deg"] % period_h
        d_raw                   = circ_dist(l0_raw, TARGET_DIR_DEG % period_h, period_h)

        rngE                    = np.random.default_rng(12345 + bcut*1000 + nbins + ce)

        p_shift, chirp_obs      = pvalue_circular_shift_null(imf_al_mean, x0, ENSEMBLE_SHIFT_NULL, rngE)
        p_align, _, _           = pvalue_phase_alignment_shift_null(imf_raw_mean, x0, K_HARM, TARGET_DIR_DEG, ENSEMBLE_SHIFT_NULL, rngE)

        met_chirp               = imf_chirp_metrics(x0, imf_al_mean, edge_trim_deg=EDGE_TRIM_DEG)

        # write a readable txt
        with open(ens_dir / "ensemble_metrics.txt", "w") as f:
            f.write(f"chirp_score_aligned : {met_chirp['score']}\n")
            f.write(f"chirp_slope_m       : {met_chirp['m']}\n")
            f.write(f"chirp_r2            : {met_chirp['r2']}\n")
            f.write(f"p_shift_aligned     : {p_shift}\n")
            f.write("\n")
            f.write(f"harm_k              : {K_HARM}\n")
            f.write(f"l0_raw_modP         : {l0_raw}\n")
            f.write(f"d_to_target_raw     : {d_raw}\n")
            f.write(f"harm_amp_raw        : {met_raw['amp']}\n")
            f.write(f"p_align_raw         : {p_align}\n")
            f.write("\n")
            f.write(f"aligned_shifts_bins : {shifts}\n")

        ensemble_rows.append({
                              "bcut"                : bcut, "nbins": nbins, "ce_trials": ce,
                              "n_runs"              : len(g),
                              "chirp_score_aligned" : met_chirp["score"],
                              "p_shift_aligned"     : p_shift,
                              "l0_raw_modP"         : l0_raw,
                              "d_to_target_raw"     : d_raw,
                              "harm_amp_raw"        : met_raw["amp"],
                              "p_align_raw"         : p_align,
                              "ens_dir"             : str(ens_dir),
                            })

        print(f"[ensemble] wrote: {ens_dir}")

    if ensemble_rows:
        ens_csv             = BASE / "ensemble_summary.csv"
        pd.DataFrame(ensemble_rows).to_csv(ens_csv, index = False)
        print(f"[ensemble] wrote summary: {ens_csv}")

    # =========================
    # QUICK SUMMARY PLOT (run-level; will still show p_phase_null floors)
    # =========================
    if "chirp_score" in ok.columns and "p_phase_null" in ok.columns:
        tmp                 = ok.dropna(subset = ["chirp_score", "p_phase_null"])
        
        if not tmp.empty:
            plt.figure(figsize = (9, 5))
            plt.scatter(tmp["chirp_score"], tmp["p_phase_null"], s = 15)
            plt.xlabel("chirp_score (|slope| / IQR(freq))")
            plt.ylabel("p_phase_null (phase randomisation)")
            plt.title("Chirp score vs phase-null p (all runs)")
            plt.tight_layout()
            plt.savefig(BASE / "chirp_score_vs_p_phase_null.png", dpi = 180)
            plt.close()
            print(f"[summary] wrote: {BASE/'chirp_score_vs_p_phase_null.png'}")


if __name__ == "__main__":
    main()
