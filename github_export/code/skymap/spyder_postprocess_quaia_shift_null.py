#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  1 08:03:48 2026

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy             as     np
import pandas            as     pd
from   pathlib           import Path
import matplotlib.pyplot as     plt
from   scipy.signal      import hilbert

GRID                    = Path("/Users/boyde/.spyder-py3/quaia_imf_grid")
METRICS                 = GRID / "grid_metrics.csv"

# --- knobs ---
N_SHIFT                 = 10000          # shift-null surrogates per run
IMF_TO_SCORE            = 6              # "IMF6" as you’ve been looking at (1-based)
EPS                     = 1e-12

def chirp_score_from_series(x_deg, y):
    """
    chirp_score = |slope(freq vs x)| / IQR(freq)
    freq estimated from instantaneous phase derivative (Hilbert).
    """
    x                   = np.asarray(x_deg)
    y                   = np.asarray(y)

    # Ensure strictly increasing x (should already be)
    # Compute dx robustly
    dx                  = float(np.median(np.diff(x)))
    if dx <= 0:
        return np.nan

    # analytic signal + phase
    analytic            = hilbert(y)
    phase               = np.unwrap(np.angle(analytic))

    # inst freq in cycles per degree (length n-1)
    inst                = np.diff(phase) / (2*np.pi*dx)
    x_mid               = 0.5*(x[1:] + x[:-1])

    # robust-ish slope via ordinary least squares on midpoints
    # (you can swap in robust regression later if needed)
    slope               = np.polyfit(x_mid, inst, 1)[0]

    q75, q25            = np.percentile(inst, [75, 25])
    iqr                 = (q75 - q25)
    
    return float(np.abs(slope) / (iqr + EPS))

def dipole_phase_deg(x_deg, y):
    """
    Fit y(l) ~ a*cos(l) + b*sin(l); phase phi = atan2(b, a)
    Returns phi in degrees in [0, 360).
    """
    l                   = np.deg2rad(np.asarray(x_deg))
    y                   = np.asarray(y)
    X                   = np.column_stack([np.cos(l), np.sin(l)])
    
    # least squares
    coef, *_            = np.linalg.lstsq(X, y, rcond = None)
    a, b                = coef
    phi                 = np.arctan2(b, a)                       # radians
    phi_deg             = (np.rad2deg(phi) + 360.0) % 360.0
    
    return float(phi_deg)

def load_run_outputs(run_dir: Path):
    sig                 = pd.read_csv(run_dir   / "signal.csv")
    imfs                = np.load(run_dir       / "imfs.npy")     # (n_imf, n)
    res                 = np.load(run_dir       / "residue.npy")  # (n,)

    # expected columns from your pipeline
    if "x_deg" in sig.columns:
        x               = sig["x_deg"].to_numpy()
    else:
        # fallback: first column
        x               = sig.iloc[:, 0].to_numpy()

    if "signal" in sig.columns:
        y               = sig["signal"].to_numpy()
    else:
        # fallback: second column
        y               = sig.iloc[:, 1].to_numpy()

    return x, y, imfs, res


def inst_freq_from_series(x_deg, y):
    x                   = np.asarray(x_deg)
    y                   = np.asarray(y)

    dx                  = float(np.median(np.diff(x)))
    
    if not np.isfinite(dx) or dx <= 0:
        return None, None

    analytic            = hilbert(y)
    phase               = np.unwrap(np.angle(analytic))
    inst                = np.diff(phase) / (2*np.pi*dx)
    x_mid               = 0.5*(x[1:] + x[:-1])

    # sanity check
    if inst.size < 5 or not np.all(np.isfinite(inst)):
        return None, None

    return x_mid, inst



def shift_null_pvalue(x, series, score_obs, n_shift = 2000, rng = None):
    if rng is None:
        rng                     = np.random.default_rng(0)
        
    n                           = len(series)
    
    # draw random shifts excluding 0
    shifts                      = rng.integers(1,   n, size = n_shift)
    null                        = np.empty(n_shift, dtype = float)
    
    for i, s in enumerate(shifts):
        y_shift                 = np.roll(series, s)
        null[i]                 = chirp_score_from_series(x, y_shift)
    
    # empirical p (>=)
    p                           = (1.0 + np.sum(null >= score_obs)) / (1.0 + n_shift)
    return float(p), null


def shift_null_pvalue_from_inst(x_mid, inst, n_shift, rng, eps = 1e-12):
    x_mid                       = np.asarray(x_mid)
    inst                        = np.asarray(inst)

    # OLS slope quickly: slope = dot(xc, y)/dot(xc, xc)
    xc                          = x_mid - x_mid.mean()
    den                         = float(np.dot(xc, xc)) + eps

    slope_obs                   = float(np.dot(xc, inst) / den)
    iqr                         = float(np.subtract(*np.percentile(inst, [75, 25])))
    score_obs                   = float(np.abs(slope_obs) / (iqr + eps))

    L                           = len(inst)
    max_unique                  = L - 1  # shifts 1..L-1
    
    if max_unique <= 0:
        return np.nan, np.nan

    # sample shifts WITHOUT replacement; if n_shift >= max_unique, do exhaustive
    if (n_shift is None) or (int(n_shift) >= max_unique):
        shifts                  = np.arange(1, L, dtype=int)  # exhaustive
        denom                   = max_unique + 1             # +1 for observed arrangement
    
    else:
        shifts                  = rng.choice(np.arange(1, 
                                                       L, 
                                                       dtype   = int), 
                                                       size    = int(n_shift), 
                                                       replace = False)
        denom                   = len(shifts) + 1

    ge                          = 0
    for s in shifts:
        inst_s                  = np.roll(inst, s)
        slope_s                 = float(np.dot(xc, inst_s) / den)
        score_s                 = float(np.abs(slope_s)    / (iqr + eps))
        ge                     += (score_s >= score_obs)

    p = (ge + 1.0) / float(denom)
    
    return float(p), float(score_obs)


def main():
    df                          = pd.read_csv(METRICS)
    
    if "tag" not in df.columns:
        raise RuntimeError("grid_metrics.csv must have a 'tag' column")

    rng                         = np.random.default_rng(12345)
    out_rows                    = []

    for i, row in df.iterrows():
        tag                     = row["tag"]
        run_dir                 = GRID / tag
        
        if not run_dir.exists():
            continue

        try:
            x, y_sig, imfs, res = load_run_outputs(run_dir)
            
        except Exception as e:
            print("skip", tag, "load failed:", e)
            continue

        # pick IMF index per run (period-band selection from grid_metrics.csv)
        k0_1based               = int(row.get("picked_imf_idx_1based", IMF_TO_SCORE))
        k                       = max(0, min(k0_1based - 1, imfs.shape[0] - 1))
        imf                     = imfs[k]

        # fast shift-null using instantaneous frequency series (no hilbert inside loop)
        x_mid, inst             = inst_freq_from_series(x, imf)
        
        if x_mid is None or inst is None:
            continue
        
        p_shift, score_obs      = shift_null_pvalue_from_inst(x_mid, 
                                                              inst, 
                                                              N_SHIFT, 
                                                              rng)

        phi                     = dipole_phase_deg(x, imf)

        out                     = dict(row)
        
        out.update({
                    "picked_imf_idx_1based_used"    : k + 1,
                    "chirp_score_recalc"            : score_obs,
                    "p_shift_null"                  : p_shift,
                    "phi_deg"                       : phi,
                    "n_shift_null"                  : N_SHIFT,
                  })
        
        out_rows.append(out)

    out_df                      = pd.DataFrame(out_rows)
    out_path                    = GRID / "grid_metrics_with_shiftnull.csv"
    
    out_df.to_csv(out_path, index = False)
    print("wrote:", out_path)


    plt.figure(figsize = (9, 5))
    plt.scatter(out_df["chirp_score_recalc"], out_df["p_shift_null"], s = 15)
    plt.yscale("log")
    plt.xlabel("chirp_score_recalc")
    plt.ylabel("p_shift_null (circular shift, log)")
    plt.title("Chirp score vs circular-shift null p-value")
    plt.tight_layout()
    plt.savefig(GRID / "chirp_vs_p_shiftnull.png", dpi = 180)
    plt.close()

    plt.figure(figsize = (9, 4))
    plt.hist(out_df["phi_deg"], bins = 36)
    plt.xlabel("phi_deg (phase angle, degrees)")
    plt.ylabel("count")
    plt.title("Phase angle distribution across runs")
    plt.tight_layout()
    plt.savefig(GRID / "phi_imf_hist.png", dpi = 180)
    plt.close()


if __name__ == "__main__":
    main()
