#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  1 08:30:33 2026

@author: boyde


"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stack_quaia_imfs.py

Stack / average IMFs across robustness runs, selecting IMFs by period band
(not IMF index), phase-aligning them, and re-testing chirp significance
with a higher-resolution circular-shift null.

Assumes run folders live under:
  /Users/boyde/.spyder-py3/quaia_imf_grid/<tag>/
and each contains:
  signal.csv, imfs.npy, residue.npy
and a metrics CSV in:
  /Users/boyde/.spyder-py3/quaia_imf_grid/grid_metrics_with_shiftnull.csv
(or grid_metrics.csv)
"""

import numpy                    as     np
import pandas                   as     pd
import matplotlib.pyplot        as     plt
from   pathlib                  import Path
from   scipy.signal             import hilbert
import run_quaia_imf_robustness as     rq


BASE                    = Path("/Users/boyde/.spyder-py3/quaia_imf_grid")
METRICS1                = BASE / "grid_metrics_with_shiftnull.csv"
METRICS0                = BASE / "grid_metrics.csv"

OUT                     = BASE / "STACKED"

OUT.mkdir(parents = True, exist_ok = True)

MIN_RUNS_PER_GROUP      = 2   # set to 1 if you want a “stack” row even for a single run


# ---------- knobs you can change ----------
TARGET_PERIOD_DEG       = (140.0, 240.0)        # choose scale band you care about (deg)
N_SHIFTS_NULL           = 10000                 # increase to avoid p-value floor
ALIGN_METHOD            = "corr"                # "corr" (circular correlation)
MAX_PLOT_GROUPS         = 30                    # cap plots if you have many groups
# ------------------------------------------

AMP_FLOOR_Q             = 0.15
SMOOTH_IFREQ            = True
SG_WINDOW               = 21
SG_POLY                 = 2
EDGE_TRIM_DEG           = 20.0

# ---- chirp metric knobs (must match run_quaia_imf_robustness) ----
CHIRP_KW                = dict(
                               edge_trim_deg = EDGE_TRIM_DEG,
                               amp_floor_q   = AMP_FLOOR_Q,
                               smooth_ifreq  = SMOOTH_IFREQ,
                               sg_window     = SG_WINDOW,
                               sg_poly       = SG_POLY,
                              )

CHIRP_R2_MIN            = 0.15   # set 0.0 to disable; helps reject junk inst-freq fits

N_PHASE_NULL            = 200000  # or 500000 if you want it very tight

# quick sanity
if SMOOTH_IFREQ:
    if SG_WINDOW % 2 != 1:
        raise ValueError("SG_WINDOW must be odd when SMOOTH_IFREQ=True")
        
    if SG_WINDOW < 5:
        raise ValueError("SG_WINDOW too small; use >= 5")
        
    if SG_POLY >= SG_WINDOW:
        raise ValueError("SG_POLY must be < SG_WINDOW")


def circ_dist(x, y, period):
    """Circular distance on [0, period)."""
    d                   = (x - y + period/2) % period - period/2
    
    return abs(d)


def fit_harmonic_phase(l_deg, y, k = 2, include_const = True):
    """
    Fit y(l) ≈ c + a cos(kθ) + b sin(kθ), θ=l*pi/180.
    Returns the principal maximum longitude l0 in [0, 360/k).
    """
    theta               = np.deg2rad(l_deg)
    cosk                = np.cos(k * theta)
    sink                = np.sin(k * theta)

    if include_const:
        X               = np.column_stack([np.ones_like(theta), cosk, sink])
    else:
        X               = np.column_stack([cosk, sink])

    beta, *_            = np.linalg.lstsq(X, y, rcond = None)
    
    if include_const:
        c, a, b         = beta
        
    else:
        c, a, b         = 0.0, beta[0], beta[1]

    amp                 = np.hypot(a, b)
    phi                 = np.arctan2(b, a)  # phase of cos term

    # Max occurs when kθ - phi = 0 => θ0 = phi/k
    l0                  = (np.rad2deg(phi) / k) % (360.0 / k)

    return dict(c       = c, 
                a       = a, 
                b       = b, 
                amp     = amp, 
                phi_rad = phi, 
                l0_deg  = l0)


def load_metrics():
    if METRICS1.exists():
        return pd.read_csv(METRICS1)
    
    return pd.read_csv(METRICS0)


def pick_xy_cols(sig_df: pd.DataFrame):
    # try common names
    x_candidates        = ["x_deg", "x", "l", "lon", "galactic_l", "bin_center", "theta"]
    y_candidates        = ["signal", "counts", "y",  "value", "n", "N"]
    xcol                = next((c for c in x_candidates if c in sig_df.columns), sig_df.columns[0])
    ycol                = next((c for c in y_candidates if c in sig_df.columns), sig_df.columns[1])
    return xcol, ycol

def est_period_zero_cross(y: np.ndarray, x: np.ndarray):
    zc                  = np.where(np.signbit(y[:-1]) != np.signbit(y[1:]))[0]
    
    if len(zc) < 4:
        return np.nan
    
    dx                  = np.diff(x[zc])
    
    return 2.0 * np.nanmean(dx)


def period_per_imf(imfs: np.ndarray, x: np.ndarray):
    return np.array([est_period_zero_cross(imfs[k], x) for k in range(imfs.shape[0])], dtype = float)



def circular_shift(y: np.ndarray, k: int):
    return np.roll(y, k)


def best_circular_alignment(ref: np.ndarray, y: np.ndarray):
    """
    Find shift k that maximizes circular correlation corr(ref, roll(y,k)).
    Uses FFT for speed.
    """
    # demean
    a                   = ref - ref.mean()
    b                   = y   - y.mean()
    fa                  = np.fft.rfft(a)
    fb                  = np.fft.rfft(b)
    
    # circular cross-correlation via FFT
    cc                  = np.fft.irfft(fa * np.conj(fb), n = len(a))
    k                   = int(np.argmax(cc))
    
    return k, cc[k]


def pvalue_circular_shift_null(y: np.ndarray, x: np.ndarray, n_shifts: int, rng: np.random.Generator):
    obs                 = rq.imf_chirp_metrics(x, y, **CHIRP_KW)["score"]

    if not np.isfinite(obs):
        return np.nan, obs

    n                   = len(y)
    ge                  = 0
    
    for _ in range(n_shifts):
        k               = int(rng.integers(0, n))
        ys              = np.roll(y, k)
        sc              = rq.imf_chirp_metrics(x, ys, **CHIRP_KW)["score"]
        
        if np.isfinite(sc) and sc >= obs:
            ge += 1
    
    return float((ge + 1.0) / (n_shifts + 1.0)), float(obs)


def pvalue_phase_null(y: np.ndarray, x: np.ndarray, nsurr: int, rng: np.random.Generator):
    obs                 = rq.imf_chirp_metrics(x, y, **CHIRP_KW)["score"]
    
    if not np.isfinite(obs):
        return np.nan, obs

    ge                  = 0
    n_good              = 0
    
    for _ in range(nsurr):
        ys              = rq.phase_randomised_surrogate(y, rng)
        sc              = rq.imf_chirp_metrics(x, ys, **CHIRP_KW)["score"]
        
        if np.isfinite(sc):
            n_good += 1
            if sc >= obs:
                ge += 1

    p                   = (ge + 1.0) / (n_good + 1.0)
    
    return float(p), float(obs)


def load_run(tag: str):
    rdir                = BASE / tag
    sig_path            = rdir / "signal.csv"
    imf_path            = rdir / "imfs.npy"
    
    if not (sig_path.exists() and imf_path.exists()):
        return None

    sig                 = pd.read_csv(sig_path)
    xcol, ycol          = pick_xy_cols(sig)
    x                   = sig[xcol].to_numpy()
    y                   = sig[ycol].to_numpy()
    imfs                = np.load(imf_path)
    
    return x, y, imfs, xcol, ycol, rdir


def select_imf_best_in_band(imfs     : np.ndarray, 
                            x        : np.ndarray, 
                            band     : tuple[float, float],
                            chirp_kw : dict, 
                            r2_min   : float = 0.0):
    """
    Among IMFs whose estimated period lies in band, choose the one with the
    highest chirp score (optionally requiring r2 >= r2_min). Falls back to
    closest-to-midpoint if nothing passes filters.
    """
    per                 = period_per_imf(imfs, x)
    
    if not np.any(np.isfinite(per)):
        
        # total failure of period estimation; pick IMF 0 as a safe default
        met0            = rq.imf_chirp_metrics(x, imfs[0], **chirp_kw)
        return 0, per, met0

    lo, hi              = band
    cand                = np.where((per >= lo) & (per <= hi))[0]

    mid                 = 0.5*(lo + hi)

    # fallback if no IMFs in band
    if cand.size == 0:
        k0              = int(np.nanargmin(np.abs(per - mid)))
        return k0, per, None  # None => no per-imf chirp diag

    best_k              = None
    best_sc             = -np.inf
    best_met            = None

    for k in cand:
        met             = rq.imf_chirp_metrics(x, imfs[k], **chirp_kw)
        sc              = met.get("score", np.nan)
        r2              = met.get("r2", np.nan)

        if not np.isfinite(sc):
            continue
        
        if np.isfinite(r2) and (r2 < r2_min):
            continue

        if sc > best_sc:
            best_sc     = float(sc)
            best_k      = int(k)
            best_met    = met

    # if everything filtered out, revert to midpoint choice within cand
    if best_k is None:
        best_k          = int(cand[np.argmin(np.abs(per[cand] - mid))])
        
        best_met        = rq.imf_chirp_metrics(x, 
                                               imfs[best_k], 
                                               **chirp_kw)

    return best_k, per, best_met


def pvalue_phase_alignment_shift_null(y, x_deg, k, target_deg, n_shifts, rng):
    period              = 360.0 / k
    obs_fit             = fit_harmonic_phase(x_deg, y, k = k)
    obs_l0              = obs_fit["l0_deg"] % period
    obs_d               = circ_dist(obs_l0, target_deg % period, period)

    n                   = len(y)
    ge                  = 0
    
    for _ in range(n_shifts):
        shift           = int(rng.integers(0, n))
        y_s             = np.roll(y, shift)
        l0_s            = fit_harmonic_phase(x_deg, 
                                             y_s, 
                                             k = k)["l0_deg"] % period
        d_s             = circ_dist(l0_s, 
                                    target_deg % period, 
                                    period)
        
        if d_s <= obs_d:
            ge += 1

    p                               = (ge + 1.0) / (n_shifts + 1.0)
    
    return p, obs_l0, obs_d


def main():
    df                              = load_metrics()
    if "error" in df.columns:
        df                          = df[df["error"].isna()].copy()

    if "tag" not in df.columns:
        raise RuntimeError("metrics file must include a 'tag' column matching run folder names.")

    # group knobs
    keys                            = [k for k in ["bcut", "nbins", "ce_trials"] if k in df.columns]
    
    if len(keys) == 0:
        raise RuntimeError("metrics file needs columns like bcut/nbins/ce_trials to group runs.")

    groups                          = df.groupby(keys)
    rng                             = np.random.default_rng(12345)

    rows                            = []
    plot_count                      = 0

    for gkey, gdf in groups:
        tags                        = gdf["tag"].astype(str).tolist()

        X                           = None

        # 1) select IMF per run (raw), record its chirp diagnostics
        Y_runs                      = []          # list of selected IMFs (raw)
        tag_runs                    = []
        meta_runs                   = []       # per-run meta for debugging: k, period, score, r2
        
        for tag in tags:
            loaded                  = load_run(tag)
            
            if loaded is None:
                continue
        
            x, y, imfs, xcol, ycol, rdir = loaded
        
            if X is None:
                X                   = x.copy()
        
            if len(x) != len(X) or not np.allclose(x, X):
                continue
        
            k, per, met_k           = select_imf_best_in_band(imfs, 
                                                              X, 
                                                              TARGET_PERIOD_DEG, 
                                                              CHIRP_KW, 
                                                              r2_min = CHIRP_R2_MIN)
            
            # QC gate: if the best-in-band IMF fails r2_min, skip this run entirely
            r2                      = met_k.get("r2", np.nan) if met_k is not None else np.nan
            
            if (not np.isfinite(r2)) or (r2 < CHIRP_R2_MIN):
                continue

            yk                      = imfs[k].copy()
        
            Y_runs.append(yk)
            
            tag_runs.append(tag)
        
            sc                      = np.nan if met_k is None else met_k.get("score", np.nan)
            r2                      = np.nan if met_k is None else met_k.get("r2",    np.nan)
            
            meta_runs.append((k, 
                              float(per[k]) if np.isfinite(per[k]) else np.nan, 
                              float(sc), float(r2)))
        
        if len(Y_runs) < MIN_RUNS_PER_GROUP:
            continue
        
        # 2) choose reference IMF as the one with highest chirp_score (finite), else max RMS
        scores                      = np.array([m[2] for m in meta_runs], dtype=float)
        
        if np.any(np.isfinite(scores)):
            ref_idx                 = int(np.nanargmax(scores))
            
        else:
            rms                     = np.array([np.sqrt(np.mean(v*v)) for v in Y_runs], dtype=float)
            ref_idx                 = int(np.argmax(rms))
        
        ref                         = Y_runs[ref_idx]
        
        # 3) build raw + aligned stacks
        Ys_raw                      = [v.copy() for v in Y_runs]
        
        Ys_aligned                  = []
        shifts                      = []
        for v in Y_runs:
            if ALIGN_METHOD == "corr":
                sh, _               = best_circular_alignment(ref, v)
                Ys_aligned.append(np.roll(v, sh))
                shifts.append(int(sh))
                
            else:
                Ys_aligned.append(v.copy())
                shifts.append(0)
        
        if len(Ys_raw) < MIN_RUNS_PER_GROUP:
            continue

        
        Yraw                        = np.vstack(Ys_raw)
        Yal                         = np.vstack(Ys_aligned)
        
        y_mean_raw                  = Yraw.mean(axis = 0)
        y_std_raw                   = Yraw.std( axis = 0)
        
        y_mean_al                   = Yal.mean( axis = 0)
        y_std_al                    = Yal.std(  axis = 0)
        
        # -----------------------------
        # (A) Harmonic / direction test
        # Use RAW mean (absolute longitude phase preserved)
        # -----------------------------
        K_HARM                      = 2  # k=2 => 180° periodic component
        
        fit_raw                     = fit_harmonic_phase(X, 
                                                         y_mean_raw, 
                                                         k = K_HARM)
        amp_h                       = fit_raw["amp"]
        
        p_align, l0_raw, d_to_target = pvalue_phase_alignment_shift_null(
                                                                         y_mean_raw,
                                                                         X,
                                                                         k          = K_HARM,
                                                                         target_deg = 264.0,          # function mods by 180 internally for k=2
                                                                         n_shifts   = N_SHIFTS_NULL,
                                                                         rng=rng
                                                                        )
        
        # phi0 for RAW stack (optional diagnostic)
        ph_raw                      = np.angle(hilbert(y_mean_raw))
        phi0_raw                    = float(np.degrees(ph_raw[0]) % 360.0)
        
        # -----------------------------
        # (B) Chirp significance test
        # Use ALIGNED mean (better SNR; phase is not absolute)
        # -----------------------------
        p_shift, chirp_obs          = pvalue_circular_shift_null(y_mean_al, 
                                                                 X, 
                                                                 N_SHIFTS_NULL, 
                                                                 rng)
        
        rng_phase                   = np.random.default_rng(54321)   # fixed seed for reproducibility
        
        p_phase_stacked, _          = pvalue_phase_null(y_mean_al, 
                                                        X, 
                                                        N_PHASE_NULL, 
                                                        rng_phase)

        ph_al                       = np.angle(hilbert(y_mean_al))
        phi0_al                     = float(np.degrees(ph_al[0]) % 360.0)





        rec                         = dict(zip(keys, gkey if isinstance(gkey, tuple) else (gkey,)))
        rec.update({
                    "n_runs_loaded"       : len(Ys_raw),

                    "chirp_score_stacked" : chirp_obs,     # from aligned stack
                    "p_shift_stacked"     : p_shift,       # from aligned stack
                    
                    "phi0_deg_raw"        : phi0_raw,
                    "phi0_deg_aligned"    : phi0_al,
                    
                    "l0_deg_raw"          : l0_raw,
                    "d_to_target_deg_raw" : d_to_target,
                    "harm_amp_raw"        : amp_h,
                    "p_align_raw"         : p_align,  
                    "picked_k_list"       : str([m[0] for m in meta_runs]),
                    "picked_period_list"  : str([m[1] for m in meta_runs]),
                    "picked_score_list"   : str([m[2] for m in meta_runs]),
                    "picked_r2_list"      : str([m[3] for m in meta_runs]),
                    "align_shifts_bins"   : str(shifts),
                    "ref_tag"             : tag_runs[ref_idx],
                    "p_phase_stacked"     : p_phase_stacked,

                  })
        
        rows.append(rec)

        # plots
        if plot_count < MAX_PLOT_GROUPS:
            label                         = "_".join([f"{k}{v}" for k, v in zip(keys, (gkey if isinstance(gkey, tuple) else (gkey,)))])
            fig                           = plt.figure(figsize = (10,4))
            plt.plot(X, y_mean_raw, label = "stack mean (raw)", alpha = 0.5)
            plt.fill_between(X, 
                             y_mean_raw - y_std_raw, 
                             y_mean_raw + y_std_raw, 
                             alpha = 0.10)
            
            plt.plot(X, 
                     y_mean_al, 
                     label      = "stack mean (aligned)", 
                     linewidth  = 2)
            
            plt.fill_between(X, 
                             y_mean_al - y_std_al, 
                             y_mean_al + y_std_al, 
                             alpha = 0.25)
            
            plt.title(
                      f"Stacked IMF (period band {TARGET_PERIOD_DEG[0]}-{TARGET_PERIOD_DEG[1]} deg)\n"
                      f"{label} | n={len(Ys_raw)} | chirp={chirp_obs:.4g} | p_shift={p_shift:.3g} | "
                      f"l0_raw={l0_raw:.1f}° | p_align={p_align:.3g}"
                     )

            plt.xlabel("x_deg (galactic longitude)")
            plt.ylabel("IMF amplitude")
            plt.legend()
            plt.tight_layout()
            figpath                     = OUT / f"stacked_{label}.png"
            plt.savefig(figpath, dpi    = 180)
            plt.close(fig)
            plot_count += 1

    # ---- output (never write a blank/no-header CSV) ----
    out_csv                 = OUT / "stacked_summary.csv"
    
    OUT_COLS                = keys + [
                                      "n_runs_loaded",
                                      "chirp_score_stacked",
                                      "p_shift_stacked",
                                      "phi0_deg_raw",
                                      "phi0_deg_aligned",
                                      "l0_deg_raw",
                                      "d_to_target_deg_raw",
                                      "harm_amp_raw",
                                      "p_align_raw",
                                      
                                      # debug / provenance
                                      "picked_k_list",
                                      "picked_period_list",
                                      "picked_score_list",
                                      "picked_r2_list",
                                      "align_shifts_bins",
                                      "ref_tag",
                                     ]

    
    if len(rows) == 0:
        # header-only file so pandas can read it
        out_df              = pd.DataFrame(columns = OUT_COLS)
        
    else:
        out_df              = pd.DataFrame(rows)
    
    out_df.to_csv(out_csv, index = False)
    
    print(f"[stack] wrote: {out_csv}")
    print(f"[stack] plots in: {OUT}")


if __name__ == "__main__":
    main()
