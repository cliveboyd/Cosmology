#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 01:20:50 2026

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
overnight_quaia_pipeline.py

Runs the Quaia IMF pipeline end-to-end overnight:
  1) run_quaia_imf_robustness.py
  2) postprocess_quaia_shift_null.py
  3) stack_quaia_imfs.py

Archives outputs into a timestamped folder and prints quick summaries.
"""

import sys
import shutil
import subprocess
from   pathlib  import Path
from   datetime import datetime

import pandas   as     pd


ROOT                    = Path("/Users/boyde/.spyder-py3")
BASE                    = ROOT / "quaia_imf_grid"

ROBUST                  = ROOT / "run_quaia_imf_robustness.py"
SHIFT                   = ROOT / "postprocess_quaia_shift_null.py"
STACK                   = ROOT / "stack_quaia_imfs.py"

# What to run
DO_ROBUST               = True
DO_SHIFT                = True
DO_STACK                = True

# Files/folders we care about archiving after the run
ARTIFACTS = [
             BASE / "grid_metrics.csv",
             BASE / "grid_metrics_with_shiftnull.csv",
             BASE / "chirp_score_vs_p.png",
             BASE / "STACKED",
            ]

# Also archive all ENSEMBLE_* folders (if present)
ARCHIVE_ENSEMBLES       = True


def run_script(path: Path, args=None, cwd=None) -> None:
    args                = [] if args is None else list(args)
    
    if not path.exists():
        raise FileNotFoundError(f"Missing script: {path}")

    cmd = [sys.executable, str(path)] + args
    print("\n" + "=" * 90)
    print("RUN:", " ".join(cmd))
    print("CWD:", str(cwd or ROOT))
    print("=" * 90)

    p = subprocess.run(cmd, cwd=str(cwd or ROOT), text = True)
    if p.returncode != 0:
        raise RuntimeError(f"Script failed (returncode={p.returncode}): {path.name}")


def safe_copy(src: Path, dst: Path) -> None:
    if not src.exists():
        return
    dst.parent.mkdir(parents = True, exist_ok = True)
    
    if src.is_dir():
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(src, dst)
        
    else:
        shutil.copy2(src, dst)


def snapshot_existing(archive_dir: Path) -> None:
    """
    If key files already exist, copy them to archive_dir/_preexisting
    so you don't lose the previous run outputs.
    """
    pre                     = archive_dir / "_preexisting"
    pre.mkdir(parents = True, exist_ok = True)

    # copy any existing key files/folders before overwriting
    for a in ARTIFACTS:
        safe_copy(a, pre / a.name)

    if ARCHIVE_ENSEMBLES and BASE.exists():
        for d in BASE.glob("ENSEMBLE_b*_n*_ce*"):
            safe_copy(d, pre / d.name)


def archive_outputs(archive_dir: Path) -> None:
    archive_dir.mkdir(parents = True, exist_ok = True)

    for a in ARTIFACTS:
        safe_copy(a, archive_dir / a.name)

    if ARCHIVE_ENSEMBLES and BASE.exists():
        ens_dir             = archive_dir / "ENSEMBLES"
        
        ens_dir.mkdir(parents = True, exist_ok = True)
        for d in BASE.glob("ENSEMBLE_b*_n*_ce*"):
            safe_copy(d, ens_dir / d.name)


def summarize_grid_metrics():
    p                           = BASE / "grid_metrics.csv"
    
    if not p.exists():
        print("\n[summary] grid_metrics.csv not found.")
        return

    df                          = pd.read_csv(p)
    df_ok                       = df[df.get("error").isna()] if "error" in df.columns else df
    print("\n[summary] grid_metrics.csv rows:", len(df), "ok:", len(df_ok))

    if "p_phase_null" in df_ok.columns:
        print("[summary] p_phase_null min:", df_ok["p_phase_null"].min())
        print(df_ok["p_phase_null"].describe())

    cols = [c for c in ["tag","bcut","nbins","ce_trials","seed","chirp_score","p_phase_null","picked_period_deg","harm_l0_deg_modP","harm_d_to_target_deg"] if c in df_ok.columns]
    if "chirp_score" in df_ok.columns and cols:
        top                     = df_ok.sort_values("chirp_score", ascending=False)[cols].head(15)
        
        print("\n[summary] top 15 runs by chirp_score:")
        print(top.to_string(index = False))


def summarize_shift_metrics():
    p                           = BASE / "grid_metrics_with_shiftnull.csv"
    if not p.exists():
        print("\n[summary] grid_metrics_with_shiftnull.csv not found (shift postprocess may not have run).")
        return

    df = pd.read_csv(p)
    print("\n[summary] grid_metrics_with_shiftnull.csv rows:", len(df))

    # Try to find a shift-null p column
    pcols                       = [c for c in df.columns if c.lower().startswith("p_shift")]
    if pcols:
        for c in pcols:
            s                   = df[c].dropna()
            
            if len(s):
                print(f"[summary] {c} min:", s.min(), "median:", s.median(), "max:", s.max())


def summarize_stacked():
    p                           = BASE / "STACKED" / "stacked_summary.csv"
    
    if not p.exists():
        print("\n[summary] stacked_summary.csv not found (stack step may not have run).")
        return

    df = pd.read_csv(p)
    print("\n[summary] stacked_summary.csv rows:", len(df))

    sort_cols                   = []
    if "p_shift_stacked" in df.columns:
        sort_cols.append("p_shift_stacked")
        
    if "chirp_score_stacked" in df.columns:
        sort_cols.append("chirp_score_stacked")

    if sort_cols:
        df2                     = df.sort_values(sort_cols, ascending = [True] + [False]*(len(sort_cols)-1))
        cols                    = [c for c in ["bcut","nbins","ce_trials","n_runs_loaded","chirp_score_stacked","p_shift_stacked","phi0_deg_stacked","l0_deg_stacked","d_to_target_deg","p_align_stacked"] if c in df2.columns]
        print("\n[summary] top 20 stacked groups:")
        print(df2[cols].head(20).to_string(index = False))


def main():
    BASE.mkdir(parents = True, exist_ok = True)

    stamp                       = datetime.now().strftime("%Y%m%d_%H%M%S")
    archive_root                = ROOT / "quaia_overnight_archives"
    archive_dir                 = archive_root / stamp

    print("[overnight] archive dir:", archive_dir)

    # Snapshot anything that already exists before we overwrite it
    snapshot_existing(archive_dir)

    # 1) Robustness grid (writes BASE/grid_metrics.csv and ENSEMBLE_* folders)
    if DO_ROBUST:
        run_script(ROBUST, cwd = ROOT)

    # 2) Shift-null postprocess (writes BASE/grid_metrics_with_shiftnull.csv)
    if DO_SHIFT:
        run_script(SHIFT, cwd = ROOT)

    # 3) Stacking (writes BASE/STACKED/stacked_summary.csv + plots)
    if DO_STACK:
        run_script(STACK, cwd = ROOT)

    # Archive outputs
    archive_outputs(archive_dir)

    # Summaries
    summarize_grid_metrics()
    summarize_shift_metrics()
    summarize_stacked()

    print("\n[overnight] DONE.")
    print("[overnight] archived to:", archive_dir)


if __name__ == "__main__":
    main()
