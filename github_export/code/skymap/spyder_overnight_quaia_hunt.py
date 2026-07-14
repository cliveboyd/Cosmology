#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 23:17:55 2026

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
overnight_quaia_hunt.py

Overnight scan for interesting stacked-IMF chirp signals.

What it does:
  - Runs missing robustness runs (quaia_imf_diagnostics.py) for a chosen grid
  - Builds a dedicated metrics CSV for *this scan* (does NOT overwrite your main metrics)
  - Runs your stack_quaia_imfs.py (configured to use that metrics CSV + separate output dir)
  - Performs subset search per group:
      * full stack
      * leave-one-out (LOO)
      * best pairs (screen with small nsurr, refine the most promising with large nsurr)
  
    - Writes:
      <BASE>/OVERNIGHT_<stamp>/metrics_scan.csv
      <BASE>/OVERNIGHT_<stamp>/STACKED/stacked_summary.csv
      <BASE>/OVERNIGHT_<stamp>/subset_hits.csv
      <BASE>/OVERNIGHT_<stamp>/subset_hits_top.csv
      <BASE>/OVERNIGHT_<stamp>/STACKED/plots_*/...

This is designed around the tension you observed:
  - p_shift (circular longitude roll null) vs
  - p_phase (phase-randomisation null)

So we compute BOTH and flag “phase-only” hits explicitly.

"""

import os
import time
import itertools
from   datetime import datetime
from   pathlib  import Path

import numpy    as np
import pandas   as pd

# Ensure no interactive matplotlib windows if run from Terminal
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

import run_quaia_imf_robustness as rq
import stack_quaia_imfs         as sq


# =========================
# CONFIG (edit these)
# =========================

         
# --- selection + chirp knobs (must match your stack script intent) ---
TARGET_PERIOD_DEG   = (140.0, 240.0)

EDGE_TRIM_DEG       = 20.0
AMP_FLOOR_Q         = 0.15
SMOOTH_IFREQ        = True
SG_WINDOW           = 21
SG_POLY             = 2
CHIRP_R2_MIN        = 0.15       # QC gate

# --- null sizes ---
SHIFT_NULL_N        =  50000     # for pvalue_circular_shift_null (stack-level)
PHASE_NULL_SCREEN   =   5000     # cheap screen for many candidate subsets
PHASE_NULL_BIG      = 200000     # refine only best candidates

# --- subset search ---
MAX_PAIRS_PER_GROUP = 30        # only refine top-N pairs by chirp_score
TOP_PLOTS_PER_GROUP =  8        # how many subset plots to save per group

# Flag thresholds (used for ranking “interesting”)
P_PHASE_STRONG      = 1e-4
P_SHIFT_STRONG      = 1e-3
PHASE_BIG_TRIGGER   = 1e-3


# --- scan scope (allow %run -i overrides) ---
BCUTS     = globals().get("BCUTS",     [20])
NBINS     = globals().get("NBINS",     [360])
CE_TRIALS = globals().get("CE_TRIALS", [2000])
SEEDS     = globals().get("SEEDS",     list(range(1, 201)))

#OVERIDES
BCUTS     = [10, 20, 30]
NBINS     = [360, 720]
CE_TRIALS = [2000]
SEEDS     = list(range(1, 151))     # 150 per group

SHIFT_NULL_N      = 50000
PHASE_NULL_SCREEN = 3000
PHASE_NULL_BIG    = 200000

MAX_PAIRS_PER_GROUP = 40
CHIRP_R2_MIN        = 0.15
TARGET_PERIOD_DEG   = (140.0, 240.0)




# =======================================
#BCUTS               = [ 10,   20,    30]
#NBINS               = [360,  720,  2000]
#SEEDS               = list(range(1,  81))


def stamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def ensure_chirp_kw():
    if SMOOTH_IFREQ:
        if SG_WINDOW % 2 != 1:
            raise ValueError("SG_WINDOW must be odd when SMOOTH_IFREQ=True")
            
        if SG_WINDOW < 5:
            raise ValueError("SG_WINDOW too small; use >= 5")
            
        if SG_POLY >= SG_WINDOW:
            raise ValueError("SG_POLY must be < SG_WINDOW")


def chirp_kw_dict():
    return dict(
                edge_trim_deg = EDGE_TRIM_DEG,
                amp_floor_q   = AMP_FLOOR_Q,
                smooth_ifreq  = SMOOTH_IFREQ,
                sg_window     = SG_WINDOW,
                sg_poly       = SG_POLY,
               )


def run_exists(outdir: Path) -> bool:
    try:
        rq.assert_run_artifacts(outdir)
        return True
    
    except Exception:
        return False


def compute_run_selection(tag: str, band, chirp_kw, r2_min: float):
    loaded              = sq.load_run(tag)
    
    if loaded is None:
        return   None

    x, y, imfs, *_      = loaded
    k, per, met         = sq.select_imf_best_in_band(imfs, x, band, chirp_kw, r2_min=r2_min)
   
    if met is  None:
        return None

    sc                  = float(met.get("score", np.nan))
    r2                  = float(met.get("r2", np.nan))

    return dict(
                tag    = tag,
                x      = x,
                imf    = imfs[k].copy(),
                k      = int(k),
                period = float(per[k]) if np.isfinite(per[k]) else np.nan,
                score  = sc,
                r2     = r2,
               )


def phase_null_p(x, y, nsurr: int, rng: np.random.Generator, chirp_kw):
    obs                 = rq.imf_chirp_metrics(x, 
                                               y, 
                                               **chirp_kw).get("score", np.nan)
    
    if not np.isfinite(obs):
        return np.nan, obs

    ge                  = 0
    n_good              = 0
    
    for _ in range(nsurr):
        ys              = rq.phase_randomised_surrogate(y, rng)
        sc              = rq.imf_chirp_metrics(x, ys, **chirp_kw).get("score", np.nan)
        
        if np.isfinite(sc):
            n_good += 1
            
            if sc >= obs:
                ge += 1

    return float((ge + 1.0) / (n_good + 1.0)), float(obs)


def align_to_ref(ref, v):
    sh, _               = sq.best_circular_alignment(ref, v)
    return np.roll(v, sh), int(sh)


def build_aligned_mean(imfs_list):
    # ref = max score handled outside; this just aligns to the provided ref
    raise NotImplementedError


def plot_subset(outpng: Path, x, y_mean_raw, y_mean_al, title: str):
    plt.figure(figsize = (10, 4))
    plt.plot(x, y_mean_raw, label = "mean (raw)", alpha         = 0.55)
    plt.plot(x, y_mean_al,  label = "mean (aligned)", linewidth = 2)
    plt.axhline(0, linestyle = "--", alpha = 0.4)
    plt.title(title)
    plt.xlabel("x_deg")
    plt.ylabel("IMF amplitude")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpng, dpi = 180)
    plt.close()


def scan_group_subsets(group_key, tags, outdir: Path, band, chirp_kw, r2_min: float):
    """
    For a group (bcut, nbins, ce_trials): load selected IMFs, QC-gate,
    then evaluate:
      - FULL stack
      - LOO
      - best pairs (screen + refine)
    """
    rows                = []
    
    if len(tags) == 0:
        raise RuntimeError("No runs matched. Check bcut/nbins/ce_trials and status values.")

    # Load + QC gate
    selected = []
    for tag in tags:
        rec             = compute_run_selection(tag, band, chirp_kw, r2_min=r2_min)
        if rec is None:
            continue
        
        if (not np.isfinite(rec["score"])) or (not np.isfinite(rec["r2"])) or (rec["r2"] < r2_min):
            continue
        
        selected.append(rec)

    if len(selected) < 2:
        return rows  # nothing to do

    # Choose reference = max score (your script-style)
    scores              = np.array([r["score"] for r in selected], dtype = float)
    ref_idx             = int(np.nanargmax(scores))
    ref_tag             = selected[ref_idx]["tag"]
    ref_imf             = selected[ref_idx]["imf"]
    
    x0                  = selected[0]["x"]

    def eval_subset(name, subset):
        
        
        
        # choose subset-specific reference (script-style): max score (else max RMS)
        sub_scores      = np.array([r["score"] for r in subset], dtype = float)
        if np.any(np.isfinite(sub_scores)):
            sub_ref_idx = int(np.nanargmax(sub_scores))
            
        else:
            rms         = np.array([np.sqrt(np.mean(r["imf"] * r["imf"])) for r in subset], dtype=float)
            sub_ref_idx = int(np.argmax(rms))

        ref_tag_local   = subset[sub_ref_idx]["tag"]
        ref_imf_local   = subset[sub_ref_idx]["imf"]

        # Raw mean (absolute phase)
        Yraw            = np.vstack([r["imf"] for r in subset])
        y_mean_raw      = Yraw.mean(axis=0)

        # Aligned mean (shape/chirp)
        aligned         = []
        shifts          = []
        for r in subset:
            ya, sh      = align_to_ref(ref_imf_local, r["imf"])
            aligned.append(ya)
            shifts.append(sh)

        Yal             = np.vstack(aligned)
        y_mean_al       = Yal.mean(axis = 0)

        # save FULL series for later reproducibility (NOW we have y_mean_al + shifts)
        if name == "FULL":
            npz_dir     = outdir / "STACK_SERIES"
            
            npz_dir.mkdir(parents = True, exist_ok = True)
            
            npz_path    = npz_dir / f"stack_series_{group_key}.npz"
            np.savez(
                     npz_path,
                     X            = x0,
                     y_mean_raw   = y_mean_raw,
                     y_mean_al    = y_mean_al,
                     tags         = np.array([r["tag"]   for r in subset], dtype=object),
                     picked_score = np.array([r["score"] for r in subset], dtype=float),
                     picked_r2    = np.array([r["r2"]    for r in subset], dtype=float),
                     shifts       = np.array(shifts, dtype = int),
                     ref_tag      = ref_tag_local,
                    )

        # p_shift (your key null)
        rng_shift           = np.random.default_rng(123)
        p_shift, chirp_obs  = sq.pvalue_circular_shift_null(y_mean_al, x0, SHIFT_NULL_N, rng_shift)

        # p_phase (screen, then possibly big)
        rng_phase           = np.random.default_rng(321)
        p_phase_s, _        = phase_null_p(x0, y_mean_al, PHASE_NULL_SCREEN, rng_phase, chirp_kw)

        p_phase_b           = np.nan
        
        if np.isfinite(p_phase_s) and (p_phase_s <= 5e-3):
            
            rng_phase2      = np.random.default_rng(777)
            p_phase_b, _    = phase_null_p(x0, y_mean_al, PHASE_NULL_BIG, rng_phase2, chirp_kw)

        # harmonic / directional test on RAW mean
        K_HARM                          =   2
        TARGET_DIR_DEG                  = 264.0
        rng_align                       = np.random.default_rng(999)
        
        p_align, l0_raw, d_to_target    = sq.pvalue_phase_alignment_shift_null(
                                                                               y_mean_raw, 
                                                                               x0, 
                                                                               K_HARM, 
                                                                               TARGET_DIR_DEG, 
                                                                               
                                                                               20000, 
                                                                               rng_align
                                                                              )
        amp_h               = sq.fit_harmonic_phase(x0, y_mean_raw, k=K_HARM)["amp"]

        # plot
        label               = f"{group_key}_{name}".replace(" ", "_").replace("/", "_")
        plot_dir            = outdir / "plots_subsets"
        plot_dir.mkdir(parents = True, exist_ok = True)
        plot_subset(
                    plot_dir / f"{label}.png",
                    x0, y_mean_raw, y_mean_al,
                    title=(f"{group_key} | {name} | n={len(subset)} | "
                           f"chirp={chirp_obs:.6g} | p_shift={p_shift:.3g} | "
                           f"p_phase_s={p_phase_s:.3g} | p_phase_b={p_phase_b:.3g} | "
                           f"p_align={p_align:.3g}")
                    )

        return dict(
                    group          = str(group_key),
                    subset         = name,
                    n              = len(subset),
                    ref_tag        = ref_tag_local,
                    tags           = [r["tag"]   for r in subset],
                    picked_scores  = [r["score"] for r in subset],
                    picked_r2      = [r["r2"]    for r in subset],
                    chirp          = float(chirp_obs),
                    p_shift        = float(p_shift),
                    p_phase_screen = float(p_phase_s) if np.isfinite(p_phase_s) else np.nan,
                    p_phase_big    = float(p_phase_b) if np.isfinite(p_phase_b) else np.nan,
                    p_align_raw    = float(p_align),
                    l0_raw         = float(l0_raw),
                    d_to_target    = float(d_to_target),
                    harm_amp       = float(amp_h),
                   )

    # FULL
    rows.append(eval_subset("FULL", selected))

    # LOO
    for drop in selected:
        subset          = [r for r in selected if r["tag"] != drop["tag"]]
        
        if len(subset) >= 2:
            rows.append(eval_subset(f"LOO_drop_{drop['tag']}", subset))

    # Pair screen: rank by chirp score only (fast), then refine a few
    # Build pair candidates
    pair_rows           = []
    for i in range(len(selected)):
        for j in range(i+1, len(selected)):
            a, b        = selected[i], selected[j]
            
            # quick chirp on aligned mean
            ya, _       = align_to_ref(ref_imf, a["imf"])
            yb, _       = align_to_ref(ref_imf, b["imf"])
            y_mean      = 0.5*(ya + yb)
            sc          = rq.imf_chirp_metrics(x0, y_mean, **chirp_kw).get("score", np.nan)
            
            pair_rows.append((float(sc) if np.isfinite(sc) else -np.inf, a["tag"],  b["tag"]))

    pair_rows.sort(reverse = True, key = lambda t: t[0])
    top_pairs           = pair_rows[:MAX_PAIRS_PER_GROUP]

    # Evaluate top pairs with nulls
    for rank, (sc, t1, t2) in enumerate(top_pairs, start = 1):
        subset          = [r for r in selected if r["tag"] in (t1, t2)]
        rows.append(eval_subset(f"PAIR_rank{rank}_{t1}__{t2}", subset))

    return rows


def main():
    ensure_chirp_kw()
    chirp_kw            = chirp_kw_dict()

    base                = rq.BASE  # /Users/boyde/.spyder-py3/quaia_imf_grid
    out_root            = base / f"OVERNIGHT_{stamp()}"
    out_root.mkdir(parents = True, exist_ok = True)

    metrics_csv = out_root / "metrics_scan.csv"

    runs                = list(itertools.product(BCUTS, NBINS, CE_TRIALS, SEEDS))
    print(f"[overnight] scan out: {out_root}")
    print(f"[overnight] runs to consider: {len(runs)}")

    rows                = []
    t0                  = time.time()

    for idx, (bcut, nbins, ce, seed) in enumerate(runs, start=1):
        tag             = f"b{bcut}_n{nbins}_ce{ce}_s{seed}"
        outdir          = base / tag

        if run_exists(outdir):
            status      = "skip"
        
        else:
            status      = "run"   # means: we attempted to run it now
            try:
                rq.run_one(outdir, bcut=bcut, nbins=nbins, ce_trials=ce, seed=seed)
            
            except Exception as e:
                rows.append(dict(
                                 tag       = tag, 
                                 bcut      = bcut, 
                                 nbins     = nbins, 
                                 ce_trials = ce, 
                                 seed      = seed,
                                 status    = "fail", 
                                 error     = str(e), 
                                 outdir    = str(outdir)
                               ))
                
                pd.DataFrame(rows).to_csv(metrics_csv, index = False)
                
                print(f"[{idx}/{len(runs)}] {tag} FAIL: {e}")
                continue


        # record minimal row (stack script only needs grouping + tag)
        rows.append(dict(tag=tag, bcut=bcut, nbins=nbins, ce_trials=ce, seed=seed, status=status, outdir=str(outdir)))

        if (idx % 10) == 0:
            pd.DataFrame(rows).to_csv(metrics_csv, index = False)
            dt              = time.time() - t0
            print(f"[{idx}/{len(runs)}] checkpoint -> {metrics_csv} (elapsed {dt:.1f}s)")

    pd.DataFrame(rows).to_csv(metrics_csv, index = False)
    print(f"[overnight] wrote metrics: {metrics_csv}")

    # Configure stack module to use this scan's metrics and its own output folder
    sq.METRICS1             = metrics_csv
    sq.METRICS0             = metrics_csv

    sq.OUT                  = out_root / "STACKED"
    
    sq.OUT.mkdir(parents = True, exist_ok = True)

    sq.TARGET_PERIOD_DEG    = TARGET_PERIOD_DEG
    sq.N_SHIFTS_NULL        = SHIFT_NULL_N

    sq.EDGE_TRIM_DEG        = EDGE_TRIM_DEG
    sq.AMP_FLOOR_Q          = AMP_FLOOR_Q
    sq.SMOOTH_IFREQ         = SMOOTH_IFREQ
    sq.SG_WINDOW            = SG_WINDOW
    sq.SG_POLY              = SG_POLY
    sq.CHIRP_R2_MIN         = CHIRP_R2_MIN
    sq.CHIRP_KW             = chirp_kw_dict()

    print("[overnight] running stack_quaia_imfs.main()...")
    sq.main()

    # Subset scan
    dfm                     = pd.read_csv(metrics_csv)
    dfm                     = dfm[dfm["status"].isin(["run", "skip"])].copy()

    group_cols              = [c for c in ["bcut", "nbins", "ce_trials"] if c in dfm.columns]
    hits                    = []

    for gkey, gdf in dfm.groupby(group_cols):
        tags                = gdf["tag"].astype(str).tolist()
        key_str             = "bcut{}_nbins{}_ce{}".format(*gkey) if isinstance(gkey, tuple) else str(gkey)

        group_hits          = scan_group_subsets(
                                                 group_key=key_str,
                                                 tags     = tags,
                                                 outdir   = out_root,
                                                 band     = TARGET_PERIOD_DEG,
                                                 chirp_kw = chirp_kw,
                                                 r2_min   = CHIRP_R2_MIN,


                                                )
        hits.extend(group_hits)

    hits_df                 = pd.DataFrame(hits)
    hits_csv                = out_root / "subset_hits.csv"
    hits_df.to_csv(hits_csv, index = False)
    print(f"[overnight] wrote subset hits: {hits_csv}")

    if not hits_df.empty:
        # Rank: strongest phase evidence first, then (optionally) shift evidence
        hits_df["p_phase_rank"]     = hits_df["p_phase_big"].fillna(hits_df["p_phase_screen"])
        hits_df["phase_only_flag"]  = (hits_df["p_phase_rank"] <= P_PHASE_STRONG) & (hits_df["p_shift"] > 0.05)
        hits_df["shift_flag"]       = (hits_df["p_shift"] <= P_SHIFT_STRONG)

        top                         = hits_df.sort_values(["p_phase_rank", "p_shift", "chirp"], ascending=[True, True, False]).head(80)
        top_csv                     = out_root / "subset_hits_top.csv"
        
        top.to_csv(top_csv, index = False)
        print(f"[overnight] wrote top list: {top_csv}")

    print(f"[overnight] DONE. Review folder: {out_root}")


if __name__ == "__main__":
    main()
