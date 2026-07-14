#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
run_hht_seedgrid_v3.py

Seed-grid runner that is:
  - environment-driven (defaults, but ENV overrides win)
  - resumable (skips seeds already logged for same param-set)
  - robust to CSV schema drift (parses from stdout as fallback)
  - incremental (appends to aggregate after each seed)
  - prints one metric per line per seed
  - writes a per-run summary across seeds

Designed for Spyder (%run -i) or terminal python.


import os
os.environ["TRIALS"]     = "200"
os.environ["N_SEEDS"]    =  "40"
os.environ["SEED_START"] =   "3"
os.environ["SEED_STEP"]  =   "2"
os.environ["USE_SEED_GENERATOR"] = "1"
%run -i /Users/boyde/.spyder-py3/run_hht_seedgrid_v3.py
"""

import os
import re
import csv
import time
import runpy
import hashlib
import numpy    as     np
import pandas   as     pd

from pathlib    import Path
from io         import StringIO
from datetime   import datetime
from contextlib import redirect_stdout, redirect_stderr


# ----------------------- helpers: env parsing ----------------------------

def env_str(key: str, default: str) -> str:
    v = os.environ.get(key, "")
    return v if v.strip() != "" else default

def env_int(key: str, default: int) -> int:
    try:
        return int(env_str(key, str(default)))
    
    except Exception:
        return default

def env_float(key: str, default: float) -> float:
    try:
        return float(env_str(key, str(default)))
    
    except Exception:
        return default

def env_list_int(key: str, default_list):
    s = os.environ.get(key, "").strip()
    if not s:
        return list(default_list)
    
    out = []
    for part in s.split(","):
        part = part.strip()
        
        if part:
            out.append(int(part))
    return out

def stamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")

def safe_float(x):
    try:
        if x is None:
            return np.nan
        
        if isinstance(x, (int, float, np.floating)):
            return float(x)
        
        s = str(x).strip().lower()
        if s in ("", "na", "nan", "none"):
            return np.nan
        
        return float(s)
    
    except Exception:
        return np.nan

def is_missing(x) -> bool:
    if x is None:
        return True
    
    if isinstance(x, float) and np.isnan(x):
        return True
    
    s = str(x).strip().lower()
    return s in ("", "na", "nan", "none")


# ----------------------- user config (defaults) --------------------------

HHT_SCRIPT = Path(env_str("HHT_SCRIPT", "/Users/boyde/.spyder-py3/HHT_Surrogates_and_Locking_Compare_V4.py"))

BASE_ENV_DEFAULTS = {
                     "HHT_SANITY":       "0",
                     "HHT_CALENDAR": "solar",
                     "TRIALS":         "100",
                     "N_SPIN":        "2000",
                     "N_PERM":       "20000",
                     "N_SHIFT":      "10000",
                     "N_IAAFT":        "800",
                     "BLOCK_LEN":       "32",
                     # optionally add: "NOISE": "0.2" if your engine reads it
                    }

# Seeds: either provide SEEDS="3,7,11,..." or N_SEEDS/SEED_START/SEED_STEP
SEEDS_DEFAULT = [3, 7, 11, 19, 23, 29, 31, 37, 41, 47]
SEEDS         = env_list_int( "SEEDS", SEEDS_DEFAULT)

if not os.environ.get("SEEDS", "").strip():
    n_seeds   = env_int("N_SEEDS", len(SEEDS_DEFAULT))
    seed0     = env_int("SEED_START", 3)
    step      = env_int("SEED_STEP",  4)
    
    # only use generator mode if user explicitly asks for it
    if os.environ.get("USE_SEED_GENERATOR", "").strip() == "1":
        SEEDS = [seed0 + i * step for i in range(n_seeds)]

OUT_DIR       = Path(env_str("OUT_DIR", "/Users/boyde/.spyder-py3/plamb_runs/hht_bridge/FR_log1pz_hht_export_v4/seedgrid_runs"))
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ----------------------- robust parsing from stdout ----------------------

RE_PLV_PERM   = re.compile(r"PLV perm p \(raw\)\s*=\s*([0-9.eE+-]+)\s+PLV perm p \(det\)\s*=\s*([0-9.eE+-]+)")
RE_RAYLEIGH   = re.compile(r"\[Rayleigh test\]\s*[\r\n]+raw p\s*=\s*([0-9.eE+-]+)\s+det p\s*=\s*([0-9.eE+-]+)", re.MULTILINE)
RE_WALD       = re.compile(r"Wald p\s*=\s*([0-9.eE+-]+)")
RE_RUNL_BEST  = re.compile(r"Running-L best L_center\s*=\s*([0-9.eE+-]+)\s+amp\s*=\s*([0-9.eE+-]+)\s+p\s*=\s*([0-9.eE+-]+)")
RE_RUNL_SPIN  = re.compile(r"\[Running-[λL] spin null\]\s*global p\s*=\s*([0-9.eE+-]+)")
WRITE_RE      = re.compile(r"^\[write\]\s+(.+?)\s+\(schema=", re.MULTILINE)

def parse_metrics_from_text(text: str) -> dict:
    out       = {}
    m         = RE_PLV_PERM.search(text)
    if m:
        out["PLV_p_raw"] = m.group(1); out["PLV_p_det"] = m.group(2)
    m         = RE_RAYLEIGH.search(text)
    
    if m:
        out["Rayleigh_p_raw"] = m.group(1); out["Rayleigh_p_det"] = m.group(2)
    m         = RE_WALD.search(text)
    
    if m:
        out["wald_p_det"] = m.group(1)
    m         = RE_RUNL_BEST.search(text)
    
    if m:
        out["runL_best_center"] = m.group(1); out["runL_best_amp"] = m.group(2); out["runL_best_p"] = m.group(3)
    m        = RE_RUNL_SPIN.search(text)
    
    if m:
        out["runningL_spin_global_p_det"] = m.group(1)
    return out


# ----------------------- summary CSV tail read ---------------------------

def read_header_and_last_row(csv_path: Path, *, max_tail_bytes=256_000) -> dict:
    csv_path             = Path(csv_path)
    
    with open(csv_path, "r", encoding="utf-8", errors="replace", newline="") as f:
        rdr              = csv.reader(f)
        header           = next(rdr, None)
    if not header:
        return {}

    size                 = csv_path.stat().st_size
    nbytes               = min(max_tail_bytes, size)
    
    with open(csv_path, "rb") as f:
        f.seek(-nbytes, 2)
        tail             = f.read().decode("utf-8", errors="replace")

    lines = [ln for ln in tail.splitlines() if ln.strip()]
    if not lines:
        return {}

    last_line = lines[-1]
    try:
        row               = next(csv.reader([last_line]))
    except Exception:
        return {}

    if len(row) < len(header):
        row += [""] * (len(header) - len(row))
        
    elif len(row) > len(header):
        header            = header + [f"extra_{i+1}" for i in range(len(row) - len(header))]

    return dict(zip(header, row))


# ----------------------- run identity + resume ---------------------------

def param_fingerprint(base_env: dict) -> str:
    # stable per-parameter-set fingerprint, used for resume/skip
    items                 = sorted((k, str(v)) for k, v in base_env.items())
    h                     = hashlib.sha1(repr(items).encode("utf-8")).hexdigest()[:10]
    return h

def seed_already_done(seed: int, fp: str) -> bool:
    # look for an existing log matching this seed+fingerprint
    pat                    = f"hht_run_seed{seed}_FP{fp}_"
    return any(p.name.startswith(pat) and p.suffix == ".log" for p in OUT_DIR.glob(f"{pat}*.log"))


# ----------------------- core runner -------------------------------------

def run_once(seed: int, base_env: dict, fp: str) -> dict:
    # apply env (ENV overrides already baked into base_env below)
    for k, v in base_env.items():
        os.environ[k]      = str(v)
    os.environ["SEED"]     = str(seed)
    os.environ["HHT_SEED"] = str(seed)

    np.random.seed(seed)

    buf                    = StringIO()
    t0                     = time.time()

    init_globals           = {"SEED": seed}
    
    with redirect_stdout(buf), redirect_stderr(buf):
        runpy.run_path(str(HHT_SCRIPT), run_name="__main__", init_globals=init_globals)

    dt                     = time.time() - t0
    text                   = buf.getvalue()

    log_path               = OUT_DIR / f"hht_run_seed{seed}_FP{fp}_{stamp()}.log"
    log_path.write_text(text, encoding="utf-8")

    m                      = WRITE_RE.search(text)
    summary_path           = Path(m.group(1).strip()) if m else None
    
    if summary_path is None or not summary_path.exists():
        fallback           = Path("/Users/boyde/.spyder-py3/plamb_runs/hht_bridge/FR_log1pz_hht_export_v4/hht_seasonal_robustness_summary_v4.csv")
        summary_path       = fallback if fallback.exists() else None

    row = {
            "run_fp"       : fp,
            "seed"         : seed,
            "runtime_s"    : float(dt),
            "log_path"     : str(log_path),
            "summary_path" : str(summary_path) if summary_path else "",
            **{f"env_{k}"  : str(v) for k, v in base_env.items()},
          }

    parsed                 = parse_metrics_from_text(text)

    if summary_path and summary_path.exists():
        last               = read_header_and_last_row(summary_path)
        row.update(last)

    for k, v in parsed.items():
        if is_missing(row.get(k)):
            row[k]        = v

    return row


def append_row_csv(path: Path, row: dict):
    path                  = Path(path)
    exists                = path.exists()

    # fixed column ordering: preserve prior header if file already exists
    if exists:
        with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
            header        = next(csv.reader(f), [])
            
        cols              = header if header else list(row.keys())
    else:
        cols              = list(row.keys())

    with open(path, "a", encoding="utf-8", newline="") as f:
        w                 = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        
        if not exists:
            w.writeheader()
        
        # ensure all fields exist
        out               = {k: row.get(k, "") for k in cols}
        w.writerow(out)


def fisher_combined_p(ps):
    ps                    = [p for p in ps if np.isfinite(p) and p > 0.0 and p <= 1.0]
    if len(ps) == 0:
        return np.nan
    
    # Fisher: X = -2 Σ ln(p) ~ χ²(2k)
    x                     = -2.0 * np.sum(np.log(ps))
    try:
        # use scipy if available; otherwise approximate via survival from mpmath
        import scipy.stats as st
        return float(st.chi2.sf(x, 2 * len(ps)))
    
    except Exception:
        try:
            import mpmath as mp
            k             = len(ps)
            # survival function for chi-square
            return float(mp.gammainc(k, x/2, mp.inf) / mp.gamma(k))
        
        except Exception:
            return np.nan


def main():
    # bake ENV overrides into base_env (override defaults)
    base_env              = {k: env_str(k, v) for k, v in BASE_ENV_DEFAULTS.items()}
    fp                    = param_fingerprint(base_env)

    agg_path              = OUT_DIR / f"seedgrid_aggregate_FP{fp}_{stamp()}.csv"

    print("\n[seedgrid] starting")
    print(f"script   = {HHT_SCRIPT}")
    print(f"out_dir  = {OUT_DIR}")
    print(f"run_fp   = {fp}")
    print(f"seeds    = {SEEDS}")
    print("env      = " + "  ".join([f"{k}={v}" for k, v in base_env.items()]))

    rows                  = []
    for seed in SEEDS:
        print(f"\nseed = {seed}")

        if seed_already_done(seed, fp):
            print("  [skip] existing log found for this seed+param-set")
            continue

        r                 = run_once(seed, base_env, fp)
        rows.append(r)

        # one metric per line
        metrics = [
                   ("PLV_p_raw",                  safe_float(r.get("PLV_p_raw"))),
                   ("PLV_p_det",                  safe_float(r.get("PLV_p_det"))),
                   ("Rayleigh_p_raw",             safe_float(r.get("Rayleigh_p_raw"))),
                   ("Rayleigh_p_det",             safe_float(r.get("Rayleigh_p_det"))),
                   ("wald_p_det",                 safe_float(r.get("wald_p_det"))),
                   ("runL_best_p",                safe_float(r.get("runL_best_p"))),
                   ("runningL_spin_global_p_det", safe_float(r.get("runningL_spin_global_p_det"))),
                  ]
        
        for k, v in metrics:
            sv = f"{v:.6g}" if np.isfinite(v) else "NA"
# %%
            print(f"{k:<28} = {sv}")


        append_row_csv(agg_path, r)

    # final summary (read back the aggregate for consistency)
    if agg_path.exists():
        df                    = pd.read_csv(agg_path)
        # coerce p-cols
        for c in ["PLV_p_raw","PLV_p_det","Rayleigh_p_raw","Rayleigh_p_det","runL_best_p"]:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")

        def frac_sig(col):
            if col not in df.columns:
                return np.nan
            x                 = df[col].to_numpy()
            x                 = x[np.isfinite(x)]
            
            if len(x) == 0:
                return np.nan
            return float(np.mean(x < 0.05))

        print("\n[seedgrid] done")
        print(f"aggregate_csv = {agg_path}")

        # headline summary
        for col in ["PLV_p_raw","PLV_p_det","Rayleigh_p_raw","Rayleigh_p_det","runL_best_p"]:
            if col in df.columns:
                vals          = df[col].to_numpy()
                vals          = vals[np.isfinite(vals)]
                med           = np.median(vals) if len(vals) else np.nan
                print(f"{col:<16}  frac<0.05={frac_sig(col):.3f}  median={med:.6g}")

        # combined p across seeds (detrended metrics typically most meaningful)
        if "Rayleigh_p_det" in df.columns:
            comb              = fisher_combined_p(df["Rayleigh_p_det"].to_numpy())
            print(f"Fisher combined p (Rayleigh_p_det) = {comb:.6g}")

if __name__ == "__main__":
    main()
