# =============================================================================
# USAGE / RUN NOTES (drop-in remark)
#
# This script runs the HHT/CEEMDAN surrogate + phase-locking pipeline on two CSVs:
#   - RAW residuals (sn_residuals_enriched.csv)
#   - DETRENDED residuals (sn_residuals_detrended.csv)
#
# INPUT CSV REQUIREMENTS
#   RAW must contain:  RA, DEC, z, resid        (w_mu optional; defaults to 1)
#   DET must contain:  RA, DEC, z, resid_detrended   (or resid_detrended_whitened)
#   Optional for lunar modes: an MJD column (default name 'mjd', configurable)
#
# RUN PARAMETER PRECEDENCE (highest → lowest)
#   1) CLI args (e.g. --trials 10000)
#   2) Environment vars (e.g. TRIALS=10000)
#   3) Globals already defined in the runtime (e.g. TRIALS=10000 in a driver loop)
#   4) Hard defaults in this script
#
# CORE RUN KNOBS (CLI + env + globals)
#   --trials / TRIALS                 : global depth knob (default 1000)
#   --ce-trials / CE_TRIALS           : CEEMDAN ensemble size (defaults to TRIALS)
#   --n-spin / N_SPIN                 : sky-rotation spin null draws (default 1000)
#   --n-spin-runl / N_SPIN_RUNL        : running-λ global spin null draws (default N_SPIN)
#   --n-perm / N_PERM                 : PLV permutation null draws (default 20000)
#   --n-shift / N_SHIFT               : circular-shift / block-λ null draws (default 10000)
#   --n-sur-e2 / N_SUR_E2             : E2 surrogate count (defaults to TRIALS)
#   --n-sur-e23 / N_SUR_E23           : E23 surrogate count (defaults to TRIALS)
#   --seed / SEED                     : main RNG seed (default 7)
#   --spin-seed / SPIN_SEED           : seed for running-λ spin null (default SEED)
#   --block-len / BLOCK_LEN           : seasonal block length for block nulls
#                                      (set BLOCK_LEN=None via env/globals to disable)
#
# CALENDAR / LUNAR CONTROLS (env only)
#   HHT_CALENDAR     : 'solar' | 'lunar-phase' | 'lunar-lambda' | 'none'   (default 'solar')
#   HHT_MJD_COL      : name of MJD column in CSVs (default 'mjd')
#   HHT_WRITE_LUNAR  : '1' to write enriched-lunar CSVs alongside originals (default '0')
#
# SURROGATE DIAGNOSTICS (env)
#   HHT_VERBOSE      : '1' prints failure summaries / scale checks
#   HHT_COLLECT_VALID: '1' try to collect exactly ns valid surrogates
#   HHT_MAX_ATTEMPTS : cap attempts when collecting valid (default 5*ns)
#   IAAFT_ITERS      : IAAFT iterations (default 60)
#
# OPTIONAL DEBUG / OUTPUTS (env)
#   HHT_SANITY       : '1' runs synthetic locked-dataset sanity_check() then exits
#   HHT_DECOMP_LOG   : '1' prints IMF decomposition ladder (PLV/energy/period)
#   HHT_DECOMP_TOPK  : ladder rows (default 10)
#   HHT_IMF_PLOTS    : '1' write IMF power heatmaps + energy bars (RAW + DET)
#   HHT_IMF_PLOT_TOPK, HHT_IMF_PLOT_SMOOTH, HHT_IMF_PLOT_HILBERT, HHT_IMF_PLOT_DPI
#   N_BOOT           : bootstrap draws for dipole direction CI (default 1000)
#
# TYPICAL RUN EXAMPLES
#   (A) Command line:
#       python HHT_Surrogates_and_Locking_Compare_V4.py \
#         --seed 7 --trials 10000 --ce-trials 200 \
#         --n-spin 10000 --n-perm 20000 --n-shift 10000 \
#         --n-sur-e2 10000 --n-sur-e23 10000 --block-len 32
#
#   (B) Spyder/IPython driver loop (globals take effect via parse_known_args()):
#       %env HHT_CALENDAR=solar
#       for seed in (3,7,11):
#           for B in (None, 24, 32, 48):
#               SEED=seed; BLOCK_LEN=B
#               TRIALS=10000; CE_TRIALS=200
#               N_SPIN=10000; N_SHIFT=10000; N_PERM=20000
#               %run -i "/path/to/HHT_Surrogates_and_Locking_Compare_V4.py"
#
# OUTPUTS
#   - Appends a row to: hht_seasonal_robustness_summary_V5.csv (in p_det.parent)
#   - Optional diagnostic PNGs if HHT_IMF_PLOTS=1
# ============================================================================= 



import os
import hashlib
import argparse
import tempfile
import numpy               as     np
import pandas              as     pd
import astropy.units       as     u
import importlib.util

from   datetime            import datetime
from   collections         import Counter
from   astropy.coordinates import SkyCoord
from   pathlib             import Path
from   scipy.interpolate   import CubicSpline
from   scipy.signal        import hilbert
from   PyEMD               import CEEMDAN
from   numpy.linalg        import lstsq
from   astropy.coordinates import BarycentricTrueEcliptic


# --- V5 Summary Ledger ----------------------------------------------------

SCHEMA_VERSION  = 5
SUMMARY_COLS_V5 = [
                    # E2
                    "E2_raw","E2_det","p_raw","p_det",
                    
                    # PLV & locks
                    "PLV_raw_lambda",         "PLV_det_lambda",   "PLV_phase_stab",
                    "PLV_p_raw",              "PLV_p_det",        "PLV_p_det_spin",   "PLV_p_det_shift",
                    "PLV_p_det_block_lambda", "PLV_spin_med_det", "PLV_spin_q95_det",
                    
                    # Rayleigh
                    "Rayleigh_p_raw", "Rayleigh_p_det", "Rayleigh_p_det_shift",
                    
                    # Regression RAW
                    "amp_ecl","wald_F", "wald_p", "phi_ecl_raw_deg", "se_phi_raw_deg",
                    
                    # IAAFT/BLOCK/AR1
                    "E23_raw_iaaft", "p23_raw_iaaft", "E23_det_iaaft", "p23_det_iaaft",
                    "E23_raw_block", "p23_raw_block", "E23_det_block", "p23_det_block",
                    "E23_raw_ar1",   "p23_raw_ar1",   "E23_det_ar1",   "p23_det_ar1",
                    
                    # Regression DET
                    "amp_ecl_det","wald_F_det","wald_p_det","phi_ecl_det_deg","se_phi_det_deg",
                    
                    # Coords + diagnostics
                    "phi_det_gal_l_deg", "phi_det_gal_b_deg",
                    "runL_best_center",  "runL_best_amp",     "runL_best_p",
                    "phi_det_boot_deg",  "phi_det_ci68_lo",   "phi_det_ci68_hi",
                    
                    # Null sizes + controls + provenance
                    "n_perm",    "n_spin",   "n_shift",        "n_iaaft",
                    "BLOCK_LEN", "seed",     "trials",         "noise",
                    "data_raw",  "data_det", "provenance_hash",
                    "calendar",
                    "schema_version",        "run_timestamp",
                    
                    # --- Running-λ (ecliptic sector scan) spin-global (DET) ---
                    "runningL_obs_amp_det",
                    "runningL_obs_L_deg_det",
                    "runningL_spin_global_p_det",
                    "runningL_spin_q50_det",
                    "runningL_spin_q95_det",
                   ]


def _coerce_int(x, default):
    try: return int(x)
    except: return default

   
parser = argparse.ArgumentParser(add_help = False)

parser.add_argument("--trials",                          type = int)
parser.add_argument("--n-spin",      dest="n_spin",      type = int)
parser.add_argument("--seed",                            type = int)
parser.add_argument("--n-perm",      dest="n_perm",      type = int)
parser.add_argument("--n-shift",     dest="n_shift",     type = int)
parser.add_argument("--spin-seed",   dest="spin_seed",   type = int)
parser.add_argument("--n-spin-runl", dest="n_spin_runl", type = int)
parser.add_argument("--n-sur-e2",    dest="n_sur_e2",    type = int)
parser.add_argument("--n-sur-e23",   dest="n_sur_e23",   type = int)
parser.add_argument("--ce-trials",   dest="ce_trials",   type = int)
parser.add_argument("--block-len",   dest="block_len",   type = str)

args, _ = parser.parse_known_args()


def _pick_int(name, arg_val, env_name, default):
    g                         = globals()
    if arg_val is not None:
        return int(arg_val), f"arg:{name}"
    
    if env_name in os.environ:
        try:
            return int(os.environ[env_name]), f"env:{env_name}"
        
        except Exception:
            return int(default), f"env:{env_name}(bad)->default"

    
    if g.get(env_name) is not None:
        return int(g.get(env_name)),      f"global:{env_name}"
    
    return int(default), "default"

# --- resolve run knobs (order matters) ------------------------------------

# pick TRIALS first (this is your global default depth knob)
TRIALS, TRIALS_SRC            = _pick_int("trials", getattr(args, "trials", None), "TRIALS", 1000)


# CE_TRIALS drives CEEMDAN ensemble size; default it to TRIALS unless explicitly set
CE_TRIALS, CE_TRIALS_SRC      = _pick_int(
                                         "ce_trials",
                                         getattr(args, "ce_trials", None),
                                         "CE_TRIALS",
                                         TRIALS
                                        )

# now resolve the rest
N_SPIN,  N_SPIN_SRC            = _pick_int("n_spin",  args.n_spin,  "N_SPIN",  1000)
SEED,    SEED_SRC              = _pick_int("seed",    args.seed,    "SEED",       7)

N_PERM,  N_PERM_SRC            = _pick_int("n_perm",  getattr(args, "n_perm",  None), "N_PERM",  20000)
N_SHIFT, N_SHIFT_SRC           = _pick_int("n_shift", getattr(args, "n_shift", None), "N_SHIFT", 10000)

SPIN_SEED,   SPIN_SEED_SRC     = _pick_int("spin_seed",   args.spin_seed,   "SPIN_SEED", SEED)
N_SPIN_RUNL, N_SPIN_RUNL_SRC   = _pick_int("n_spin_runl", args.n_spin_runl, "N_SPIN_RUNL", N_SPIN)

# these depend on TRIALS (so must come after TRIALS is known)
N_SUR_E2,  N_SUR_E2_SRC        = _pick_int("n_sur_e2",  getattr(args, "n_sur_e2",  None), "N_SUR_E2",  TRIALS)
N_SUR_E23, N_SUR_E23_SRC       = _pick_int("n_sur_e23", getattr(args, "n_sur_e23", None), "N_SUR_E23", TRIALS)

_B = args.block_len if args.block_len is not None else (os.environ.get("BLOCK_LEN") or globals().get("BLOCK_LEN"))
if _B is None:
    BLOCK_LEN = None
else:
    s = str(_B).strip().lower()
    BLOCK_LEN = None if s in ("none", "") else int(float(_B))
 

print(f"seed      = {SEED} ({SEED_SRC})  "
      f"TRIALS    = {TRIALS} ({TRIALS_SRC})  "
      f"CE_TRIALS = {CE_TRIALS} ({CE_TRIALS_SRC})  "
      f"N_SUR_E2  = {N_SUR_E2} ({N_SUR_E2_SRC})  "
      f"N_SUR_E23 = {N_SUR_E23} ({N_SUR_E23_SRC})  "
      f"N_SPIN    = {N_SPIN} ({N_SPIN_SRC})   BLOCK_LEN = {BLOCK_LEN}", flush = True)

    
def ecliptic_lambda_rad(ra_deg, dec_deg):
    c                   = SkyCoord(ra    = ra_deg*u.deg, 
                                   dec   = dec_deg*u.deg, 
                                   frame = "icrs")
    
    return c.transform_to(BarycentricTrueEcliptic()).lon.to_value(u.rad)





def load_module_from_path(py_path: str, module_name="_hht_target"):
    spec = importlib.util.spec_from_file_location(module_name, py_path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod

def detect_csv_paths(mod):
    """Return dict of detected CSV paths keyed by filename."""
    wanted = {
        "sn_residuals_enriched.csv": None,
        "sn_residuals_detrended.csv": None,
    }

    # 1) Prefer explicit variables if present
    explicit_keys = ("p_raw", "p_det", "RAW_CSV", "DET_CSV", "raw_csv", "det_csv", "P_RAW", "P_DET")
    for k in explicit_keys:
        if hasattr(mod, k):
            v = getattr(mod, k)
            if isinstance(v, (str, Path)):
                s = str(v)
                for fname in wanted:
                    if s.endswith(fname):
                        wanted[fname] = s

    # 2) Fallback: scan all module globals for Path/str ending in those filenames
    for _, v in mod.__dict__.items():
        if isinstance(v, (str, Path)):
            s = str(v)
            for fname in wanted:
                if wanted[fname] is None and s.endswith(fname):
                    wanted[fname] = s

    return {k: v for k, v in wanted.items() if v is not None}




def maybe_detect_csv_paths_from_target():
    """
    If HHT_TARGET_SCRIPT is set, load it and detect CSV paths.
    Returns a dict: {filename: fullpath}.
    Returns {} if not enabled or nothing detected.
    """
    target = os.environ.get("HHT_TARGET_SCRIPT", "").strip()
    if not target:
        return {}

    try:
        mod = load_module_from_path(target)
        found = detect_csv_paths(mod)

        if not found:
            print("[WARN] Could not detect RAW/DET CSV paths from target module.", flush=True)
            return {}

        for fname, fullpath in found.items():
            print(f"Detected: {fname} -> {fullpath}", flush=True)
            if not Path(fullpath).exists():
                print(f"  [WARN] {fullpath} does not exist on disk.", flush=True)

        return found

    except Exception as e:
        print(f"[WARN] Target-script detect failed: {type(e).__name__}: {e}", flush=True)
        return {}



def runningL_max_amp(resid, lam, centers = None, half_width = np.deg2rad(20), base_weights = None):
    if centers is None:
        centers         = np.linspace(0, 2*np.pi, 180, endpoint = False)  # every 2°
        
    amp_max             = -np.inf
    best_center         = np.nan
    
    for Lc in centers:
        dphi            = np.angle(np.exp(1j*(lam - Lc)))               # wrap to [-π, π]
        
        w               = np.zeros_like(dphi)
        m               = np.abs(dphi) < half_width
        w[m]            = np.cos((dphi[m] / half_width) * (np.pi/2))    # cosine taper
        
        if base_weights is not None:
            w           = w * np.clip(base_weights, 0, None)
            
        sel             = w > 0
        
        if sel.sum() < 10:
            continue
        
        X               = np.column_stack([np.ones(sel.sum()), 
                                           np.sin(lam[sel]), 
                                           np.cos(lam[sel])])
        
        y               = resid[sel]
        sw              = np.sqrt(w[sel])
        beta, *_        = np.linalg.lstsq(X * sw[:, None], 
                                          y * sw, 
                                          rcond = None)
        
        amp             = np.hypot(beta[1], beta[2])
        
        if amp > amp_max:
            amp_max, best_center = amp, Lc
            
    return amp_max, best_center
    

def runningL_spin_global_p(resid,
                           lam,
                           base_weights = None,
                           N            = 2000,
                           seed         = 7,
                           centers      = None,
                           half_width   = np.deg2rad(20),
                           rng          = None,
                           return_null  = False,
                          ):
    """
    Global spin (random rotation) test for the running-L statistic.

    Idea:
      - Compute observed max amplitude over running-L windows using (resid, lam).
      
      - Create a null distribution by adding a random global rotation to lam:
            lam_rot = (lam + U[0, 2π)) mod 2π
        and recompute the max amplitude each time.
        
      - One-sided p-value: P(null >= observed), with +1 smoothing.


    Parameters
    ----------
    resid : array-like
        Residuals associated with each object (same length as lam).
    
    lam : array-like
        Longitudes (radians), same length as resid. Will be wrapped to [0, 2π).
    
    base_weights : array-like or None
        Optional weights passed through to runningL_max_amp.
   
    N : int
        Number of random spins (null draws). Use your driver-loop N_SPIN/TRIALS here.
    
    seed : int
        Seed used if rng is not provided.
    
    centers : array-like or None
        Centers passed through to runningL_max_amp.
    
    half_width : float
        Half-width (radians) passed through to runningL_max_amp.
    
    rng : np.random.Generator or None
        If provided, used directly (recommended for reproducibility across calls).
    
    return_null : bool
        If True, also return the full null distribution (max_stats).

    Returns
    -------
    obs_amp : float
        Observed maximum amplitude.
    
    obs_L : float
        Observed L (or center) where the maximum occurs (as returned by runningL_max_amp).
    
    p : float
        One-sided Monte Carlo p-value with +1 smoothing.
    
    q50 : float
        Median of null distribution.
    
    q95 : float
        95th percentile of null distribution.
  
    (optional) max_stats : np.ndarray
        Null distribution of max amplitudes (length N), if return_null=True.
    """
    
    if rng is None:
        rng = np.random.default_rng(seed)

    lam = np.asarray(lam, dtype=float) % (2 * np.pi)
    if N is None:
        raise ValueError("N must be an integer > 0 (pass N_SPIN/TRIALS explicitly).")
    
    N = int(N)
    if N <= 0:
        raise ValueError("N must be an integer > 0.")

    # Observed statistic
    obs_amp, obs_L = runningL_max_amp(resid,
                                      lam,
                                      centers      = centers,
                                      half_width   = half_width,
                                      base_weights = base_weights,
                                     )

    # Null distribution
    max_stats = np.empty(N, dtype = float)

    # Pre-draw shifts (slightly faster and makes RNG use explicit)
    shifts = rng.uniform(0.0, 2 * np.pi, size = N)
    for i, sh in enumerate(shifts):
        lam_rot         = (lam + sh) % (2 * np.pi)
        max_stats[i], _ = runningL_max_amp(resid,
                                           lam_rot,
                                           centers      = centers,
                                           half_width   = half_width,
                                           base_weights = base_weights,
                                          )

    # One-sided p-value with +1 smoothing (avoids p=0)
    p                   = (np.count_nonzero(max_stats >= obs_amp) + 1.0) / (N + 1.0)

    q50, q95            = np.quantile(max_stats, [0.50, 0.95])

    if return_null:
        return obs_amp, obs_L, p, q50, q95, max_stats
    
    return obs_amp, obs_L, p, q50, q95


def _sha256_of_paths(paths):
    h                    = hashlib.sha256()
    for p in map(str, paths):
        
        try:
            with open(p, "rb") as f:
                for chunk in iter(lambda: f.read(1<<20), b""):
                    h.update(chunk)
                    
        except FileNotFoundError:
            h.update(f"[missing]{p}".encode())
            
    return h.hexdigest()


def write_summary_V5(row: dict, out_dir, filename="hht_seasonal_robustness_summary_V5.csv"):
        out_dir               = Path(out_dir)
        out_dir.mkdir(parents = True, exist_ok = True)
        summary_path          = out_dir / filename
    
        row["schema_version"] = SCHEMA_VERSION
        row["run_timestamp"]  = datetime.utcnow().isoformat(timespec="seconds") + "Z"
    
        # enforce column order & fill missing with NaN
        row_ordered           = {c: row.get(c, np.nan) for c in SUMMARY_COLS_V5}
    
        append                = summary_path.exists()
        fd, tmp_path          = tempfile.mkstemp(prefix = "hht_sum_", 
                                                 suffix = ".csv", 
                                                 dir    = str(out_dir))
        os.close(fd)
        
        pd.DataFrame([row_ordered]).to_csv(tmp_path, 
                                           index  = False, 
                                           header = True)

        if append:
            with open(summary_path, "ab") as dst, open(tmp_path, "rb") as src:
                _             = src.readline()  # now this really is the header
                dst.write(src.read())
            os.remove(tmp_path)
        
        else:
            os.replace(tmp_path, summary_path)

        print(f"[write] {summary_path} (schema={SCHEMA_VERSION})", flush = True)
        return str(summary_path)
    

__all__             = ["extended_regression", 
                       "plv_perm", 
                       "plv",
                       "interpolate_angle_over_L", 
                       "to_uniform_L",
                       "ecliptic_angles_from_radec",
                       "e23_fraction",
                      ]


# === calendar controls (env-driven; no argparse) =========================
# HHT_CALENDAR: 'solar' | 'lunar-phase' | 'lunar-lambda' | 'none'
CALENDAR_MODE       = os.environ.get("HHT_CALENDAR", "solar").strip().lower()

# HHT_MJD_COL: name of MJD column if present (default 'mjd')
MJD_COL             = os.environ.get("HHT_MJD_COL", "mjd")

# HHT_WRITE_LUNAR: '1' to persist enriched CSVs alongside originals
WRITE_LUNAR_COLS    = os.environ.get("HHT_WRITE_LUNAR", "0") == "1"
    
# ---- DEBUG header (safe to keep) -----------------------------------------
def _ping(stage: str) -> None:
    print(f"[HHT_Surrogates_and_Locking_Compare] {stage}", flush = True)


def safe_target_angle_on_L(df, L_grid, mode):
    try:
        return target_angle_on_L(df, L_grid, mode=mode)
    except Exception as e:
        # fallback to solar
        print(f"[calendar] target angle fallback -> solar (requested={mode}, reason={e})", flush=True)
        return target_angle_on_L(df, L_grid, mode="solar")


# add near other helpers
def target_angle_on_L(df, L_grid, *, mode = "solar"):
    """
    Returns angle series on the L_grid to use as the PLV/Rayleigh target,
    depending on calendar mode: 'solar' (ecliptic λ), 'lunar-phase' (φ_moon),
    or 'lunar-lambda' (λ_moon).
    """
    L_src           = np.log1p(pd.to_numeric(df["z"], errors = "coerce").to_numpy(float))
    order           = np.argsort(L_src)

    if mode == "solar":
        lam, _      = ecliptic_angles_from_radec(df["RA"], df["DEC"])
        return interpolate_angle_over_L(L_src[order], lam[order], L_grid)

    if mode == "lunar-phase":
        if "phi_moon" not in df.columns:
            raise KeyError("phi_moon missing; ensure MJD present or precomputed.")
            
        phi         = pd.to_numeric(df["phi_moon"], errors = "coerce").to_numpy(float)
        return interpolate_angle_over_L(L_src[order], phi[order], L_grid)

    if mode == "lunar-lambda":
        if "lam_moon" not in df.columns:
            raise KeyError("lam_moon missing; enable astropy or lam_sun fallback.")
            
        lm          = pd.to_numeric(df["lam_moon"], errors = "coerce").to_numpy(float)
        return interpolate_angle_over_L(L_src[order], lm[order], L_grid)

    raise ValueError(f"Unknown calendar mode: {mode}")

    
# ---------------- Uniform L-grid with duplicate collapse ---------------- #
def to_uniform_L(z, y, pad = 0.15, ngrid = 384):
    L               = np.log1p(np.asarray(z, dtype = float))
    y               = np.asarray(y,          dtype = float)

    m               = np.isfinite(L) & np.isfinite(y)
    L, y            = L[m], y[m]

    # if nothing usable, return a stable NaN grid (or raise)
    if L.size == 0:
        Lu_grid     = np.linspace(0.0, 1.0, int(ngrid))
        return Lu_grid, np.full_like(Lu_grid, np.nan, dtype = float)

    o               = np.argsort(L, kind="mergesort")
    L, y            = L[o], y[o]

    df              = pd.DataFrame({"L": L, "y": y})
    g               = df.groupby("L", as_index = False).mean()
    Lu, yu_in       = g["L"].to_numpy(), g["y"].to_numpy()

    if Lu.size == 0:
        Lu_grid     = np.linspace(0.0, 1.0, int(ngrid))
        return Lu_grid, np.full_like(Lu_grid, np.nan, dtype = float)

    span            = max(Lu.max() - Lu.min(), 1e-12)
    Lu_grid         = np.linspace(Lu.min() - pad*span, 
                                  Lu.max() + pad*span, 
                                  int(ngrid))

    if Lu.size < 4:
        yu_grid     = np.interp(Lu_grid, Lu, yu_in)
        return Lu_grid, yu_grid

    try:
        cs          = CubicSpline(Lu, yu_in, bc_type = "natural")
        yu_grid     = cs(Lu_grid)
        
    except ValueError:
        yu_grid     = np.interp(Lu_grid, Lu, yu_in)

    return Lu_grid, yu_grid

# --- add near other helpers ---
def pick_det_ycol(df: pd.DataFrame) -> str:
    return "resid_detrended_whitened" if "resid_detrended_whitened" in df.columns else "resid_detrended"


def e23_fraction(y, imfs):
    recon           = np.sum(imfs, axis=0)
    resid           = y - recon
    den             = np.sum(imfs**2)     + np.sum(resid**2) + 1e-30
    num             = (np.sum(imfs[1]**2) + np.sum(imfs[2]**2)) if imfs.shape[0] >= 3 else np.nan
    return float(num / den)


def _energy_fraction(y, imfs, band):
    if band == "imf23":
        return e23_fraction(y, imfs)
    
    elif band == "imf2":
        return imf2_energy_ratio(y, imfs)
    
    return e23_fraction(y, imfs)


# ---------------- CEEMDAN + IMF2 energy and helpers --------------------- #
def ceemdan_imfs(Lu, yu, trials = 100, noise = 0.2, seed = 7):
    ce              = CEEMDAN()
    ce.trials       = int(trials)
    ce.noise_width  = float(noise)
    
    # Ensure determinism if PyEMD exposes RNG control
    for attr in ("random_state", "seed"):
        if hasattr(ce, attr):
            try:
                setattr(ce, attr, int(seed))
            
            except Exception:
                pass
    
   # print("[decomp]", 
   #       type(ce).__name__,
   #                "trials=",       getattr(ce,  "trials",       None),
   #                 "noise_width=", getattr(ce,  "noise_width", None),
   #                 "epsilon=",     getattr(ce,  "epsilon",     None))

    return ce.ceemdan(yu, Lu)

def imf2_energy_ratio(y, imfs):
    denom           = float(np.sum(y**2) + 1e-30)
    if imfs is None or imfs.shape[0] < 2:
        return np.nan
    
    e2              = float(np.sum(imfs[1]**2))
    r               = e2 / denom
    return min(max(r, 0.0), 1.0)

def imf23_energy_ratio(y, imfs):
    denom           = float(np.sum(y**2) + 1e-30)
    
    if imfs is None or imfs.shape[0] < 3:
        return np.nan
    
    e23             = float(np.sum(imfs[1]**2) + np.sum(imfs[2]**2))
    r               = e23 / denom
    return min(max(r, 0.0), 1.0)

def block_iaaft_surrogate(y, block = 32, rng = None, iters = 80):
    """
    Apply IAAFT within contiguous blocks of the series y.
    
    Preserves local distribution + local spectrum approximately, while breaking
    longer-range/global coherence across blocks.
    """
    rng             = rng or np.random.default_rng()
    y               = np.asarray(y, float)
    n               = y.size
    B               = int(block) if block is not None else acf_block_length(y)
    B               = max(8, min(B, n))

    out             = y.copy()
    
    for s in range(0, n, B):
        e           = min(s + B, n)
        seg         = y[s:e]
        
        if seg.size < 16:
            # too short: just shuffle amplitudes (minimal disruption)
            out[s:e] = rng.permutation(seg)
            
        else:
            out[s:e] = iaaft_surrogate(seg, iters = iters, rng = rng)
    return out


def keyed_block_iaaft_surrogate(y, key, block = 32, rng = None, iters = 80):
    """
    Keyed version: sort by key, apply block_iaaft_surrogate, then unsort.
    
    Use when cal_key is supplied.
    """
    rng             = rng or np.random.default_rng()
    y               = np.asarray(y, float)
    n               = y.size

    if key is None:
        return block_iaaft_surrogate(y, block = block, rng = rng, iters = iters)

    key = np.asarray(key, float)
    
    if key.size != n:
        return block_iaaft_surrogate(y, block = block, rng = rng, iters = iters)

    idx             = np.argsort(key, kind = "mergesort")
    y_sort          = y[idx]

    y_sort_s        = block_iaaft_surrogate(y_sort, 
                                            block = block, 
                                            rng   = rng, 
                                            iters = iters)

    y_s             = np.empty_like(y_sort_s)
    y_s[idx]        = y_sort_s
    return y_s

def iaaft_surrogate(y, iters = 100, rng = None, x0 = None):
    """
    Proper IAAFT surrogate:
      - Preserves the sample amplitude distribution of y
      - Approximately  preserves the power spectrum of y
    """
    rng             = rng or np.random.default_rng()
    y               = np.asarray(y, float)
    n               = y.size

    y0              = y - np.mean(y)

    # Target Fourier magnitudes (DO NOT SORT)
    Y               = np.fft.rfft(y0)
    amp_target      = np.abs(Y)

    # Target amplitude distribution
    y_sorted        = np.sort(y0)

    # Initial guess
    if x0 is None:
        x           = rng.permutation(y0)
    
    else:
        x           = np.asarray(x0, float).copy()
        x           = x - np.mean(x)

    for _ in range(int(iters)):
        # Enforce target spectrum
        X           = np.fft.rfft(x)
        ph          = np.exp(1j * np.angle(X))
        x_spec      = np.fft.irfft(amp_target * ph, n = n)

        # Enforce target distribution via rank-order mapping
        ranks       = np.argsort(np.argsort(x_spec))
        x           = y_sorted[ranks]

    # Return with mean restored (optional; your pipeline re-centers anyway)
    return x + np.mean(y)


def acf_block_length(y):
    y0              = np.asarray(y, dtype=float) - np.mean(y)
    n               = len(y0)
    
    if n < 32:
        return max(8, n // 4)
    
    ac              = np.correlate(y0, y0, mode="full")[n-1:] / np.dot(y0, y0)
    thr             = np.exp(-1.0)
    i1              = np.where(ac < thr)[0]
    i2              = np.where(np.diff(np.signbit(ac)))[0]
    cand            = []
    
    if i1.size > 0: cand.append(int(i1[0]))
    if i2.size > 0: cand.append(int(i2[0]))
    
    return max(8, min(64, (min(cand) if cand else 24)))


# === lunar utilities ======================================================
_TWOPI              = 2.0 * np.pi

def _frac(x):
    return x - np.floor(x)

def lunar_phase_from_mjd(mjd):
    """
    Synodic lunar phase φ_moon in radians, 0=new moon, 2π next synodic month.
    Accuracy is sufficient for calendar-keyed nulls/regressors.
    """
    synodic         = 29.530588853
    return _TWOPI * _frac((np.asarray(mjd, float) - 51544.5) / synodic)

def add_lunar_columns(df, *, mjd_col="mjd", lam_sun_col = None):
    """
    Adds columns: phi_moon, sin_phi_moon, cos_phi_moon.
    If astropy is available (or lam_sun_col provided), also adds lam_moon (+sin/cos).

    Notes:
      - Safe with NaN MJD values.
      - Uses astropy.coordinates.get_body("moon", ...) as the primary ephemeris path.
    """
    if mjd_col not in df.columns:
        return df  # graceful no-op

    # --- numeric MJD (may include NaNs) ---
    mjd                = pd.to_numeric(df[mjd_col], errors="coerce").to_numpy(dtype=np.float64)

    # --- synodic phase (NaNs propagate) ---
    phi                = lunar_phase_from_mjd(mjd)
    df["phi_moon"]     = phi
    df["sin_phi_moon"] = np.sin(phi)
    df["cos_phi_moon"] = np.cos(phi)

    lam_moon = None

    # --- astropy path (optional) ---
    try:
        from   astropy.time        import Time
        import astropy.units       as     u
        from   astropy.coordinates import get_body, GeocentricTrueEcliptic

        good = np.isfinite(mjd)
        if np.any(good):
            t               = Time(mjd[good], format="mjd", scale="utc")

            moon_icrs       = get_body("moon", t)             # SkyCoord
            
            # If equinox=t ever causes issues, replace with equinox=Time("J2000")
            moon_ecl        = moon_icrs.transform_to(GeocentricTrueEcliptic(equinox=t))

            lam_tmp         = np.full(mjd.shape, np.nan, dtype=np.float64)
            lam_tmp[good]   = np.mod(moon_ecl.lon.to_value(u.rad), 2.0 * np.pi)
            lam_moon        = lam_tmp

    except Exception:
        lam_moon = None

    # --- fallback: approximate lam_moon from lam_sun + synodic phase ---
    if lam_moon is None and lam_sun_col is not None and lam_sun_col in df.columns:
        lam_sun             = pd.to_numeric(df[lam_sun_col], errors="coerce").to_numpy(dtype=np.float64)
        lam_moon            = np.mod(lam_sun + phi, 2.0 * np.pi)

    if lam_moon is not None:
        df["lam_moon"]      = lam_moon
        df["sin_lam_moon"]  = np.sin(lam_moon)
        df["cos_lam_moon"]  = np.cos(lam_moon)

    return df


def calendar_key_from_df(df, *, mode = "solar"):

    if mode == "none":
        return None

    if mode == "solar":

        for cand in ("lam","lambda","lambda_ecl","lambda_sun"):
            if cand in df.columns:
                x     = pd.to_numeric(df[cand], errors="coerce").to_numpy(np.float64)
                
                if np.nanmax(x) > 2*np.pi + 0.5:
                    x = np.deg2rad(x)
                    
                return np.mod(x, _TWOPI)
    
        # fallback: compute from RA/DEC
        lam, _        = ecliptic_angles_from_radec(df["RA"], df["DEC"])
        return np.mod(lam, _TWOPI)


    if mode == "lunar-phase":
        if "phi_moon" not in df.columns:
            raise KeyError("phi_moon not found; call add_lunar_columns first.")
        return np.mod(pd.to_numeric(df["phi_moon"], errors="coerce").to_numpy(np.float64), _TWOPI)

    if mode == "lunar-lambda":
        if "lam_moon" not in df.columns:
            raise KeyError("lam_moon not found; call add_lunar_columns first.")
        return np.mod(pd.to_numeric(df["lam_moon"], errors="coerce").to_numpy(np.float64), _TWOPI)

    raise ValueError(f"Unknown calendar mode: {mode}")


def keyed_block_shuffle(y, key, block, rng):
    y                 = np.asarray(y, float); n = y.size
    
    if key is None:
        return block_shuffle(y, block=block, rng=rng)
    
    key               = np.asarray(key, float)
    if key.size != n:
        print(f"[keyed_block_shuffle] key length {key.size} != series length {n}; "
              f"falling back to unkeyed block shuffle.", flush=True)
        return block_shuffle(y, block=block, rng=rng)

    key               = np.asarray(key, float)
    

    idx               = np.argsort(key, kind="mergesort")
    y_sort            = y[idx]
    B                 = int(block) if (block is not None and int(block) > 0) else max(1, int(round(n/12)))
    starts            = np.arange(0, n, B, dtype=int)
    blocks            = [y_sort[s:e] for s, e in zip(starts, np.r_[starts[1:], n])]
    order             = np.arange(len(blocks))
    
    rng.shuffle(order)
    
    y_sort_s          = np.concatenate([blocks[i] for i in order], axis=0)
    y_s               = np.empty_like(y_sort_s)
    y_s[idx]          = y_sort_s
    return y_s



def keyed_ar1_block_surrogate(y, key, block, rng):
    return keyed_block_shuffle(ar1_boot_surrogate(y, block=block, rng=rng), key, block, rng)

def check_scale(tag, y_obs, sur_vals):
    v_obs             = float(np.var(y_obs))
    v_med             = float(np.median([np.var(s) for s in sur_vals])) if len(sur_vals) else np.nan
    print(f"[scale] {tag}: Var(obs)={v_obs:.6g}  median Var(sur)={v_med:.6g}", flush=True)


def block_shuffle(y, block=None, rng=None, varmatch=True):
    """
    Block bootstrap (row-shuffle blocks) with optional variance matching.

    Parameters
    ----------
    y : array-like
    block : int or None
        If None, choose via acf_block_length(y). Otherwise, use given size.
    rng : np.random.Generator or None
    varmatch : bool
        If True, rescale surrogate to match original std.

    Returns
    -------
    yb : np.ndarray
        Block-shuffled surrogate, same length as y.
    """
     
    rng               = rng or np.random.default_rng()
    y                 = np.asarray(y, dtype=float)
    n                 = len(y)

    # choose block size
    B                 = int(block) if block is not None else acf_block_length(y)
    B                 = max(2, min(B, n))  # keep sane and ≤ n

    k                 = n // B
    
    if block is None and os.environ.get("HHT_VERBOSE","0") == "1":
        print(f"[block_shuffle] auto block B = {B} (n={n})")
    
    if k == 0:
        # not enough data for one full block — return a copy
        return y.copy()

    core              = y[:k*B].copy().reshape(k, B)
    rng.shuffle(core, axis = 0)

    yb                = core.ravel()
    if k*B < n:
        yb            = np.concatenate([yb, y[k*B:]])

    if varmatch:
        s0            = float(np.std(y,  ddof=1)) or 1.0
        sb            = float(np.std(yb, ddof=1)) or 1.0
        
        if sb > 0:
            yb *= (s0 / sb)

    return yb

def random_rotation_matrix(rng=None):
    rng              = rng or np.random.default_rng()
    u1,u2,u3          = rng.random(3)
    q1                = np.sqrt(1-u1)*np.sin(2*np.pi*u2)
    q2                = np.sqrt(1-u1)*np.cos(2*np.pi*u2)
    q3                = np.sqrt(u1)*np.sin(2*np.pi*u3)
    q4                = np.sqrt(u1)*np.cos(2*np.pi*u3)
    
    # quaternion to rotation matrix
    a,b,c,d           = q4,q1,q2,q3
    
    return np.array([
                     [a*a+b*b-c*c-d*d, 2*(b*c-a*d),     2*(b*d+a*c)],
                     [2*(b*c+a*d),     a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c),     2*(c*d+a*b),     a*a-b*b-c*c+d*d]
                   ])


def rotate_radec(ra_deg, dec_deg, R):
    ra                = np.deg2rad(np.asarray(ra_deg, float))
    dec               = np.deg2rad(np.asarray(dec_deg,float))
    x                 = np.cos(dec)*np.cos(ra)
    y                 = np.cos(dec)*np.sin(ra)
    z                 = np.sin(dec)
    v                 = np.vstack([x,y,z])
    vr                = R @ v
    xr, yr, zr        = vr
    ra_r              = np.arctan2(yr, xr); ra_r[ra_r<0]+=2*np.pi
    dec_r             = np.arcsin(np.clip(zr,-1,1))
    return np.rad2deg(ra_r), np.rad2deg(dec_r)


def _circ_shift_p(phi, lam, n=1000, rng=None, seed=None):
    if rng is None:
        rng           = np.random.default_rng(seed if seed is not None else 0)
    obs               = plv(phi, lam)
    ge                = 0
    
    for _ in range(int(n)):
        k             = rng.integers(0, len(phi))
        
        if plv(np.roll(phi, k), lam) >= obs:
            ge += 1
    return (ge + 1.0) / (int(n) + 1.0)

def rayleigh_shift_p(phi, ang, n=10000, seed=0, rng=None):
    """
    Shift-calibrated Rayleigh test:
      statistic z = N * R^2, where R = |mean(exp(i*(phi-ang)))|.
    Null via circularly shifting phi (preserves autocorrelation structure).
    """
    phi               = np.asarray(phi, float)
    ang               = np.asarray(ang, float)
    if rng is None:
        rng           = np.random.default_rng(seed)

    delta             = np.mod(phi - ang + np.pi, 2*np.pi) - np.pi
    R                 = np.abs(np.mean(np.exp(1j*delta)))
    z_obs             = len(delta) * (R*R)

    ge                = 0
    for _ in range(int(n)):
        k             = rng.integers(0, len(phi))
        d             = np.mod(np.roll(phi, k) - ang + np.pi, 2*np.pi) - np.pi
        Rk            = np.abs(np.mean(np.exp(1j*d)))
        zk            = len(d) * (Rk*Rk)
        if zk >= z_obs:
            ge += 1

    return (ge + 1.0) / (int(n) + 1.0)


def spin_test_plv(df, phi_on_L, L_grid, n=1000, seed=0, preserve_L_rank=True):

    rng               = np.random.default_rng(seed)
    plv_vals          = []

    RA                = pd.to_numeric(df["RA"],  errors="coerce")
    DEC               = pd.to_numeric(df["DEC"], errors="coerce")
    Z                 = pd.to_numeric(df["z"],   errors="coerce")
    m                 = np.isfinite(RA) & np.isfinite(DEC) & np.isfinite(Z)
    RA, DEC, Z        = RA[m], DEC[m], Z[m]

    for _ in range(int(n)):
        R             = random_rotation_matrix(rng)
        ra_r, dec_r   = rotate_radec(RA, DEC, R)
        lam_r, _      = ecliptic_angles_from_radec(ra_r, dec_r)

        L_src         = np.log1p(Z.to_numpy(float))
        order         = np.argsort(L_src)
        
        if preserve_L_rank:
            lam_on_L  = interpolate_angle_over_L(L_src[order], lam_r[order], L_grid)
        else:
            # break strict rank preservation, but keep local smoothness via block shuffle
            B         = max(16, len(order)//16)
            lam_shuf  = block_shuffle(lam_r[order], block=B, rng=rng, varmatch=False)
            lam_on_L  = interpolate_angle_over_L(L_src[order], lam_shuf, L_grid)
            

        plv_vals.append(plv(phi_on_L, lam_on_L))

    return np.asarray(plv_vals, dtype = float)


def circular_ci_deg(samples_deg, center_deg=None, lo=16, hi=84):
    a                 = np.deg2rad(np.asarray(samples_deg, float))
    if center_deg is None:
        center       = np.angle(np.mean(np.exp(1j*a)))
    else:
        center       = np.deg2rad(center_deg)

    # shortest signed angular distance from center
    d                = np.angle(np.exp(1j*(a - center)))
    lo_d, hi_d       = np.percentile(d, [lo, hi])

    lo_deg           = (np.rad2deg(center + lo_d)) % 360.0
    hi_deg           = (np.rad2deg(center + hi_d)) % 360.0
    ctr_deg          = (np.rad2deg(center))        % 360.0
    return ctr_deg, lo_deg, hi_deg


def bootstrap_phi(df, n = 1000, seed = 1):
    rng               = np.random.default_rng(seed)
    phis              = []
    for _ in range(n):
        ix            = rng.integers(0, len(df), len(df))
        reg           = extended_regression(df.iloc[ix].reset_index(drop=True))
        phis.append(reg["phi_ecl_deg"])
    
    phis              = np.array(phis)
    
    ctr_deg, lo_deg, hi_deg = circular_ci_deg(phis, center_deg = None, lo = 16, hi = 84)
    
    return ctr_deg, (lo_deg, hi_deg)

def ecliptic_to_galactic(lam_deg, beta_deg):
    ra, dec           = ecliptic_to_icrs(lam_deg, beta_deg)
    
    cg                = SkyCoord(ra    = np.asarray(ra)*u.deg, 
                                 dec   = np.asarray(dec)*u.deg, 
                                 frame = "icrs").galactic
    
    return float(cg.l.deg), float(cg.b.deg)


def blockshuffle_iaaft_surrogate(y, block_len, rng, cal_key=None, iters=60, nbins=12):
    y                 = np.asarray(y, float)
    n                 = len(y)

    if block_len is None or block_len < 8 or block_len >= n:
        return iaaft_surrogate(y, rng = rng, iters = iters)

    m                 = n // block_len
    cut               = m * block_len
    y0                = y[:cut]

    blocks            = [y0[i*block_len:(i+1)*block_len] for i in range(m)]
    order             = np.arange(m)

    if cal_key is None:
        rng.shuffle(order)
    else:
        ck            = np.asarray(cal_key[:cut], float).reshape(m, block_len)
        block_key     = np.nanmedian(ck, axis = 1)
        bins          = np.floor((np.mod(block_key, 2*np.pi) / (2*np.pi)) * nbins).astype(int)

        for b in range(nbins):
            idx       = np.where(bins == b)[0]
            if idx.size > 1:
                tmp        = idx.copy()
                rng.shuffle(tmp)
                order[idx] = tmp   # keep your FIXED line

    out = np.empty_like(y0)
    for j, bi in enumerate(order):
        blk                              = blocks[bi]
        out[j*block_len:(j+1)*block_len] = iaaft_surrogate(blk, 
                                                           rng   = rng, 
                                                           iters = iters)

    if cut < n:
        tail          = y[cut:]
        tail_s        = iaaft_surrogate(tail, rng=rng, iters=iters) if len(tail) >= 8 else tail.copy()
        return np.concatenate([out, tail_s])

    return out

# --- parameter intake (must be BEFORE any use of these values) ---
def _get_int(name, default):
    try:
        v = globals().get(name, default)
        return int(v)
    
    except Exception:
        return int(default)

def e_band_surrogate_test(
                          Lu, yu, *,
                          band         = "imf23",
                          mode         = "iaaft",
                          ns           = 800,
                          trials       = 100,
                          noise        = 0.2,
                          seed         =   7,
                          block_len    = None,
                          standardize  = True,
                          use_fraction = True,
                          cal_key      = None,
                          tag=""
                        ):
    """
    Energy-band surrogate test with sane normalization and OPTIONAL calendar key.

    Key points:
      - p-value is computed using *valid* surrogate draws, not requested ns.
      - If HHT_COLLECT_VALID=1, we attempt to collect exactly 'ns' valid draws
        (up to a cap HHT_MAX_ATTEMPTS or 5*ns).
      - When cal_key is provided, 'block' and 'ar1' nulls use keyed block shuffles.

    Env controls:
      - HHT_COLLECT_VALID = '1' : collect exactly ns valid surrogates (recommended)
      - HHT_MAX_ATTEMPTS        : cap on attempts (default 5*ns)
      - HHT_VERBOSE      = '1' : print failure summary line
      - IAAFT_ITERS            : iterations for IAAFT (default 60)
    """
    rng                = np.random.default_rng(seed)
    fail               = Counter()
    n_attempt          = 0

    target_valid     = int(ns)
    if target_valid <= 0:
        return dict(
            E_obs      = np.nan, p = np.nan, s_q05 = np.nan, s_q50 = np.nan, s_q95 = np.nan,
            n_valid    = 0, n_attempt = 0, n_fail = 0,
            metric    = ("fraction_std" if (standardize and use_fraction) else "legacy_ratio")
        )

    iaaft_iters        = _get_int("IAAFT_ITERS", 60)                             # 40–80 is usually plenty
    block_len_eff      = block_len if block_len is not None else acf_block_length(yu)

    collect_valid      = (os.environ.get("HHT_COLLECT_VALID", "0") == "1")
    max_attempts       = int(os.environ.get("HHT_MAX_ATTEMPTS", str(5 * target_valid)))
    max_attempts       = max(max_attempts, target_valid)                         # never below target

    def _z(x):
        x              = np.asarray(x, float)
        x              = x - np.mean(x)
        s              = np.std(x, ddof = 1) or 1.0
        return x / s

    # observed
    y_obs              = _z(yu) if standardize else np.asarray(yu, float)
    imfs_obs           = ceemdan_imfs(Lu, 
                                      y_obs, 
                                      trials = trials,  
                                      noise  = noise,  
                                      seed   = seed)
    
    e_obs              = _energy_fraction(y_obs, imfs_obs, band if band in ("imf23", "imf2") else "imf23")

    ge, vals, sur_pool = 0, [], []

    def _make_surrogate():
        """Generate a surrogate series ys (unstandardized), according to mode."""
        if mode == "iaaft":
            return iaaft_surrogate(yu, rng=rng, iters=iaaft_iters)

        elif mode == "block":
            if cal_key is not None:
                return keyed_block_shuffle(yu, cal_key, block_len, rng)
            return block_shuffle(yu, block=block_len, rng=rng)

        elif mode == "ar1":
            if cal_key is not None:
                return keyed_ar1_block_surrogate(yu, cal_key, block_len, rng)
            return ar1_boot_surrogate(yu, block=block_len, rng=rng)

        elif mode == "blockiaaft":
            return blockshuffle_iaaft_surrogate(
                                                yu,
                                                block_len=block_len_eff,
                                                rng=rng,
                                                cal_key=cal_key,
                                                iters=iaaft_iters
                                               )

        else:
            return phase_randomize(yu, rng)

    def _attempt_once():
        """One surrogate attempt; returns True if a valid energy value was recorded."""
        nonlocal n_attempt, ge
        n_attempt += 1

        try:
            ys         = _make_surrogate()
            ys0        = _z(ys) if standardize else np.asarray(ys, float)

            imfs_s     = ceemdan_imfs(
                                      Lu, 
                                      ys0,
                                      trials = trials,
                                      noise  = noise,
                                      seed   = int(rng.integers(1, 1 << 31))
                                     )

            es         = _energy_fraction(ys0, imfs_s, band if band in ("imf23", "imf2") else "imf23")

            if not np.isfinite(es):
                fail["nonfinite_es"] += 1
                return False

            vals.append(float(es))
            if es >= e_obs:
                ge += 1

            if len(sur_pool) < 50:
                sur_pool.append(ys0)

            return True

        except Exception as e:
            fail[type(e).__name__] += 1
            return False

    # main loop: either fixed attempts (=target_valid) or collect_valid until vals == target_valid
    if collect_valid:
        while (len(vals) < target_valid) and (n_attempt < max_attempts):
            _attempt_once()
    
    else:
        for _ in range(target_valid):
            _attempt_once()

    vals               = np.asarray(vals, float)
    n_val              = int(len(vals))
    n_fail             = int(n_attempt - n_val)

    valid              = max(n_val, 1)
    pval               = (ge + 1.0) / (valid + 1.0)

    if valid > 3:
        q05            = float(np.nanpercentile(vals, 5))
        q50            = float(np.nanpercentile(vals, 50))
        q95            = float(np.nanpercentile(vals, 95))
    else:
        q05            = q50 = q95 = np.nan

    # diagnostics
    check_scale(f"{mode}:{tag}", y_obs, sur_pool)

    # failure summary (optional)
    if os.environ.get("HHT_VERBOSE", "0") == "1":
        top             = ", ".join([f"{k}:{v}" for k, v in fail.most_common(5)]) or "none"
        print(
            f"[E23:{mode}:{tag}] requested={target_valid}  valid={n_val}  "
            f"attempts={n_attempt}  fail={n_fail}  top_failures={top}",
            flush=True
        )
        if collect_valid and n_val < target_valid:
            print(
                f"[E23:{mode}:{tag}] WARNING: could not collect target_valid={target_valid} "
                f"within max_attempts={max_attempts}. Using valid={n_val} for p-value depth.",
                flush=True
            )

    # conservative-null note
    if np.isfinite(q50) and np.isfinite(e_obs) and e_obs > 0:
        ratio            = q50 / e_obs
        if ratio > 5.0 or ratio < 0.2:
            print(
                  f"[note] surrogate median differs markedly from observed "
                  f"(mode={mode}, med={q50:.4g}, obs={e_obs:.4g}); "
                  f"this indicates a conservative null rather than a scaling issue.",
                  flush = True
                 )

    return dict(
        E_obs            = float(e_obs),
        p                = float(pval),
        s_q05            = float(q05),
        s_q50            = float(q50),
        s_q95            = float(q95),
        n_valid          = int(valid),
        n_attempt        = int(n_attempt),
        n_fail           = int(n_fail),
        metric           = ("fraction_std" if (standardize and use_fraction) else "legacy_ratio"),
    )


    
def p_and_effect(obs, null):
    null                 = np.asarray(null, float)
    p                    = ((null >= obs).sum() + 1.0) / (len(null) + 1.0)
    med                  = float(np.median(null))
    mad                  = float(np.median(np.abs(null - med))) or float(np.std(null, ddof = 1))
    
    # robust z via MAD→σ (1.4826)
    zrob                 = (obs - med) / ((mad * 1.4826) if mad > 0 else 1.0)
    return p, med, mad, zrob


# ---------------- Surrogates (phase randomization) ---------------------- #
def phase_randomize(y, rng):
    n                     = len(y)
    y0                    = y - np.mean(y)
    F                     = np.fft.rfft(y0)
    amps                  = np.abs(F)
    phs                   = np.angle(F)
    rand_ph               = rng.uniform(0.0, 2.0*np.pi, size = len(F))
    rand_ph[0]            = phs[0]
    
    if (n % 2) == 0:
        rand_ph[-1]       = phs[-1]
    
    F_new                 = amps * np.exp(1j * rand_ph)
    y_new                 = np.fft.irfft(F_new, n = n)
    return y_new

def e2_surrogate_test(Lu, yu, trials = 100, noise = 0.2, ns = 1000, seed = 7):
    rng                   = np.random.default_rng(seed)
    
    imfs_obs              = ceemdan_imfs(Lu, 
                                         yu, 
                                         trials = trials, 
                                         noise  = noise, 
                                         seed   = seed)
    
    e_obs                 = imf2_energy_ratio(yu, imfs_obs)

    ge                    = 0
    vals                  = []

    # NEW: failure diagnostics
    fail                = Counter()
    n_attempt             = 0

    for _ in range(int(ns)):
        n_attempt      += 1
        ys                = phase_randomize(yu, rng)

        try:
            imfs_s        = ceemdan_imfs(
                                         Lu, 
                                         ys, 
                                         trials = trials, 
                                         noise  = noise, 
                                        seed    = rng.integers(1, 1 << 31)
                                        )
            
            es            = imf2_energy_ratio(ys, imfs_s)

            if not np.isfinite(es):
                fail["nonfinite_es"] += 1
                continue

            vals.append(es)
            if es >= e_obs:
                ge += 1

        except Exception as e:
            fail[type(e).__name__] += 1
            continue

    vals                  = np.asarray(vals, dtype = float)
    n_valid               = int(len(vals))
    n_fail                = int(n_attempt - n_valid)

    # NEW: print once (gate with env if you prefer)
    if os.environ.get("HHT_VERBOSE", "0") == "1":
        top = ", ".join([f"{k}:{v}" for k, v in fail.most_common(5)]) or "none"
        print(f"[E2] ns={int(ns)}  valid={n_valid}  fail={n_fail}  top_failures={top}", flush = True)

    # IMPORTANT: p-value uses VALID, not requested ns
    valid                 = max(n_valid, 1)
    pval                  = (ge + 1.0) / (valid + 1.0)
    q05, q50, q95         = (
                             (np.nanpercentile(vals,  5), 
                              np.nanpercentile(vals, 50), 
                              np.nanpercentile(vals, 95))
                              if valid > 3 else (np.nan, np.nan, np.nan)
                            )

    return dict(
                E_obs     = float(e_obs),
                p         = float(pval),
                s_q05     = float(q05),
                s_q50     = float(q50),
                s_q95     = float(q95),
                n_valid   = int(valid),
                
                # Optional extra fields (handy in your CSV later)
                n_attempt = int(n_attempt),
                n_fail    = int(n_fail),
               )


def rayleigh_p(delta):
    n                   = len(delta)
    R                   = np.abs(np.mean(np.exp(1j*delta)))
    z                   = n * R * R
    p                   = np.exp(-z) * (1 + (2*z - z*z)/(4*n))
    p                   = float(np.clip(p, 0.0, 1.0))
    
    return max(p, np.finfo(float).tiny)                                        # avoid printing 0 due to underflow




 
def ceemdan_grid(L, y, combos = [(50, 0.15), (100, 0.20), (200, 0.20), (100, 0.30)]):
    out                 = []
    
    for tr, nw in combos:
        imfs            = ceemdan_imfs(L, 
                                       y, 
                                       trials = tr, 
                                       noise  = nw, 
                                       seed   = 7)
        
        phi             = inst_phase_from_imf2(imfs)
        out.append({"trials"  : tr,
                    "noise"   : nw,
                    "phi_std" : np.std(phi)})
         
    return pd.DataFrame(out)

def fdr_bh(p):
    p                   = np.asarray(p, float); n = len(p)
    o                   = np.argsort(p)
    q                   = np.empty_like(p)
    minv                = 1.0
    
    for k,i in enumerate(o[::-1], start=1):
        v               = p[i]*n/k
        minv            = min(minv, v)
        q[i]            = minv
    
    return q

# -------------------- PLV utilities ----------------------------------------- #
def inst_phase_from_imf2(imfs):
    sig                 = imfs[1]                                              # IMF2
    an                  = hilbert(sig)
    phi                 = np.unwrap(np.angle(an))
    
    return phi

def ecliptic_to_icrs(lam_deg, beta_deg):
    lam                 = np.asarray(lam_deg,  float)
    beta                = np.asarray(beta_deg, float)
    
    try:
        ecl             = SkyCoord(lon   = lam*u.deg, 
                                   lat   = beta*u.deg, 
                                   frame = BarycentricTrueEcliptic())
        
        icrs            = ecl.transform_to("icrs")
        
        return icrs.ra.deg, icrs.dec.deg
    
    except Exception:
        # J2000 fallback
        eps             = np.deg2rad(23.43928)
        L               = np.deg2rad(lam)
        B               = np.deg2rad(beta)
        x_eq            = np.cos(B)*np.cos(L)
        y_eq            = np.cos(B)*np.sin(L)*np.cos(eps) - np.sin(B)*np.sin(eps)
        z_eq            = np.cos(B)*np.sin(L)*np.sin(eps) + np.sin(B)*np.cos(eps)
        ra              = np.degrees(np.arctan2(y_eq, x_eq)); ra[ra < 0] += 360.0
        dec             = np.degrees(np.arcsin(np.clip(z_eq, -1.0, 1.0)))
        
        return ra, dec



def ecliptic_angles_from_radec(ra_deg, dec_deg):
    ra_arr              = np.asarray(ra_deg,  dtype = float)
    dec_arr             = np.asarray(dec_deg, dtype = float)
    
    # Try Astropy transform first (no obstime to stay compatible with older builds)
    try:
        sc              = SkyCoord(ra=ra_arr * u.deg, dec = dec_arr * u.deg, frame = "icrs")
        try:
            ecl         = sc.transform_to(BarycentricTrueEcliptic())           # instance call
        
        except TypeError:
            ecl         = sc.transform_to(BarycentricTrueEcliptic)             # class call (older API)
        
        lam             = ecl.lon.to(u.rad).value
        lam             = np.mod(lam, 2.0 * np.pi)
        beta            = ecl.lat.to(u.rad).value
        return lam, beta
    
    except Exception:
        # J2000 obliquity fallback (manual rotation)
        ra              = np.deg2rad(ra_arr)
        dec             = np.deg2rad(dec_arr)
        eps             = np.deg2rad(23.43928)  # J2000 mean obliquity
        x               = np.cos(dec) * np.cos(ra)
        y               = np.cos(dec) * np.sin(ra)
        z               = np.sin(dec)
        x_e             = x
        y_e             =  y * np.cos(eps) + z * np.sin(eps)
        z_e             = -y * np.sin(eps) + z * np.cos(eps)
        lam             = np.arctan2(y_e, x_e)
        lam[lam < 0.0] += 2.0 * np.pi
        beta            = np.arcsin(np.clip(z_e, -1.0, 1.0))
        return lam, beta



def plv(phi_a, phi_b):
    delta               = phi_a - phi_b
    v                   = np.exp(1j * delta)
    return float(np.abs(np.mean(v)))

def interpolate_angle_over_L(L_src, ang_src, L_tgt):
    import numpy             as     np
    import pandas            as     pd
    from   scipy.interpolate import CubicSpline

    L_src               = np.asarray(L_src,  dtype  = float)
    ang_src             = np.asarray(ang_src, dtype = float)

    o                   = np.argsort(L_src, kind = "mergesort")
    Ls                  = L_src[o]
    ang                 = ang_src[o]

    u                   = np.exp(1j * ang)
    df                  = pd.DataFrame({"L": Ls, "re": np.real(u), "im": np.imag(u)})
    g                   = df.groupby("L", as_index = False).mean()
    Lu                  = g["L"].to_numpy()
    re                  = g["re"].to_numpy()
    im                  = g["im"].to_numpy()

    if Lu.size == 0:
        return np.full_like(L_tgt, np.nan, dtype = float)

    if Lu.size < 4:
        u_tgt           = np.interp(L_tgt, Lu, re) + 1j * np.interp(L_tgt, Lu, im)
        return np.angle(u_tgt)

    try:
        re_spl          = CubicSpline(Lu, re, bc_type = "natural")
        im_spl          = CubicSpline(Lu, im, bc_type = "natural")
        u_tgt           = re_spl(L_tgt) + 1j * im_spl(L_tgt)
    
    except ValueError:
        u_tgt           = np.interp(L_tgt, Lu, re) + 1j * np.interp(L_tgt, Lu, im)

    return np.angle(u_tgt)

def plv_lambda_block_shuffle(lam_on_L, block_len = None, n = None, seed = None):
    """
    Shuffle ecliptic λ along L in contiguous blocks (seasonal-block λ null).
    Returns an (n, nL) array of block-permuted λ series.
    """
    lam_on_L = np.asarray(lam_on_L, dtype=float)
    nL                  = lam_on_L.size

    # pick defaults from your run configuration
    if n is None:
        # For this null, N_SHIFT is typically the right knob; fall back to TRIALS, then 1000.
        n              = int(globals().get("N_SHIFT", globals().get("TRIALS", 1000)))
        
    else:
        n              = int(n)

    if seed is None:
        seed           = int(globals().get("SEED", 7))
    else:
        seed           = int(seed)

    if n <= 0:
        return np.empty((0, nL), dtype = float)

    if nL < 2:
        return np.full((n, nL), np.nan, dtype = float)

    rng                = np.random.default_rng(seed)

    # choose block length
    if block_len is None:
        B              = max(8, nL // 16)
    else:
        B              = int(block_len)

    B                  = max(2, min(B, nL))

    # build blocks (include remainder)
    starts             = range(0, nL, B)
    blocks             = [lam_on_L[s:s + B] for s in starts]
    nb                 = len(blocks)

    if nb < 2:
        return np.tile(lam_on_L, (n, 1))

    lam_null = np.empty((n, nL), dtype = float)
    
    for i in range(n):
        order          = rng.permutation(nb)
        lam_null[i]    = np.concatenate([blocks[j] for j in order])

    return lam_null

def plv_scan(phi_by_imf, theta, imf_idx):
    # returns per-IMF PLVs and max-stat
    plvs               = []
    
    for i in imf_idx:
        plvs.append(plv(phi_by_imf[i], theta))           # <-- FIX: plv not PLV
        
    plvs               = np.array(plvs, dtype = float)
    k                  = int(np.argmax(plvs))
    
    return plvs, float(plvs[k]), int(imf_idx[k])



def phases_by_imf(imfs):
    """Return wrapped phase ([-pi,pi]) per IMF (list length n_imf)."""
    imfs              = np.asarray(imfs, float)
    return [_imf_phase(imfs[i]) for i in range(imfs.shape[0])]

def pick_scan_imfs(imfs, *, start_imf = 2, efrac_min = 0.0):
    """
    Choose IMFs to include in scan.
    start_imf is 1-based IMF index to start from (default 2 => exclude IMF1).
    efrac_min is minimum energy fraction (relative to total IMF energy).
    Returns (imf_idx0based, efrac_per_imf).
    """
    imfs              = np.asarray(imfs, float)
    e                 = np.sum(imfs**2, axis = 1)
    den               = float(np.sum(e) + 1e-30)
    efrac             = e / den

    start0            = max(0, int(start_imf) - 1)
    idx               = [i for i in range(start0, imfs.shape[0]) if efrac[i] >= float(efrac_min)]

    if len(idx) == 0:
        idx           = list(range(start0, imfs.shape[0]))

    return idx, efrac

def plv_scan_perm_p(phi_by_imf, theta, imf_idx, *, n = 20000, seed = 0):
    rng               = np.random.default_rng(seed)
    _, obs_max, best  = plv_scan(phi_by_imf, theta, imf_idx)

    null              = np.empty(int(n), dtype = float)
    
    for r in range(int(n)):
        ix            = rng.permutation(len(theta))
        _, tmax, _    = plv_scan(phi_by_imf, theta[ix], imf_idx)
        null[r]       = tmax

    p                 = pvalue_maxstat(null, obs_max)
    return obs_max, best, p, null

def plv_scan_shift_p(phi_by_imf, theta, imf_idx, *, n = 10000, seed = 0):
    """
    Circular shift null: shift all IMFs by the same k each replicate
    (preserves cross-IMF structure while breaking alignment to theta).
    """
    rng               = np.random.default_rng(seed)
    _, obs_max, best  = plv_scan(phi_by_imf, theta, imf_idx)

    null              = np.empty(int(n), dtype = float)
    N                 = len(theta)
    
    for r in range(int(n)):
        k             = int(rng.integers(0, N))
        phi_shifted   = [np.roll(phi, k) for phi in phi_by_imf]
        _, tmax, _    = plv_scan(phi_shifted, theta, imf_idx)
        null[r]       = tmax

    p                 = pvalue_maxstat(null, obs_max)
    return obs_max, best, p, null

def plv_scan_blocklam_p(phi_by_imf, lam_block_shufs, imf_idx, *, obs_theta,):
    """
    Block-λ null: lam_block_shufs is (n_rep, nL) array.
    obs_theta is the observed theta (e.g. lam_d_on_Ld).
    """
    _, obs_max, best = plv_scan(phi_by_imf, obs_theta, imf_idx)

    lam_block_shufs  = np.asarray(lam_block_shufs, float)
    null             = np.empty(lam_block_shufs.shape[0], dtype = float)
    
    for r in range(lam_block_shufs.shape[0]):
        _, tmax, _   = plv_scan(phi_by_imf, lam_block_shufs[r], imf_idx)
        null[r]      = tmax

    p                = pvalue_maxstat(null, obs_max)
    return obs_max, best, p, null

def spin_test_plv_scan(df, phi_by_imf, L_grid, imf_idx, *, n = 1000, seed = 0, preserve_L_rank = True):
    """
    Spin null for scan: for each random sky rotation, compute lam_on_L and then
    return max PLV across IMFs in imf_idx.
    """
    rng              = np.random.default_rng(seed)
    phi_by_imf       = list(phi_by_imf)

    RA               = pd.to_numeric(df["RA"],  errors="coerce")
    DEC              = pd.to_numeric(df["DEC"], errors="coerce")
    Z                = pd.to_numeric(df["z"],   errors="coerce")
    m                = np.isfinite(RA) & np.isfinite(DEC) & np.isfinite(Z)
    RA, DEC, Z       = RA[m], DEC[m], Z[m]

    L_src            = np.log1p(Z.to_numpy(float))
    order            = np.argsort(L_src)

    out              = np.empty(int(n), dtype = float)
    
    for r in range(int(n)):
        Rm           = random_rotation_matrix(rng)
        ra_r, dec_r  = rotate_radec(RA, DEC, Rm)
        lam_r, _     = ecliptic_angles_from_radec(ra_r, dec_r)

        if preserve_L_rank:
            lam_on_L = interpolate_angle_over_L(L_src[order], lam_r[order], L_grid)
        else:
            B        = max(16, len(order)//16)
            lam_shuf = block_shuffle(lam_r[order], block=B, rng=rng, varmatch=False)
            lam_on_L = interpolate_angle_over_L(L_src[order], lam_shuf, L_grid)

        _, tmax, _ = plv_scan(phi_by_imf, lam_on_L, imf_idx)
        out[r] = tmax

    return out


def pvalue_maxstat(null_maxstats, obs_maxstat):
    null_maxstats     = np.asarray(null_maxstats)
    return (1.0 + np.sum(null_maxstats >= obs_maxstat)) / (1.0 + null_maxstats.size)


# ---------------- Extended regression (ecliptic + equatorial) ----------- #
def extended_regression(df: pd.DataFrame, *, calendar_mode: str = None) -> dict:
    import statsmodels.api as sm

    if calendar_mode is None:
        calendar_mode = globals().get("CALENDAR_MODE", "solar")

    df = df.copy()

    # Ensure weight column exists
    if "w_mu" not in df.columns:
        df["w_mu"]   = 1.0

    # Parse required columns
    RA               = pd.to_numeric(df.get("RA"),    errors = "coerce")
    DEC              = pd.to_numeric(df.get("DEC"),   errors = "coerce")
    y                = pd.to_numeric(df.get("resid"), errors = "coerce")
    z                = pd.to_numeric(df.get("z"),     errors = "coerce")       # no fillna()
    w                = pd.to_numeric(df.get("w_mu"),  errors = "coerce")

    # Keep only valid rows (z must be finite because it is a regressor)
    m0  = (
           np.isfinite(RA) & np.isfinite(DEC) &
           np.isfinite(y)  & np.isfinite(z)   &
           np.isfinite(w)  & (w > 0)
          )
    df2 = df.loc[m0].copy().reset_index(drop=True)

    # If nothing left, fail cleanly
    if len(df2) < 20:
        return {
                "amp_ecl"     : np.nan, 
                "phi_ecl_deg" : np.nan, 
                "se_phi_deg"  : np.nan,
                "p_sin_lam"   : np.nan, 
                "p_cos_lam"   : np.nan,
                "wald_F"      : np.nan, 
                "wald_p"      : np.nan,
                "n_used"      : int(len(df2)),
               }

    RA2                  = pd.to_numeric(df2["RA"],    errors = "coerce").to_numpy(float)
    DEC2                 = pd.to_numeric(df2["DEC"],   errors = "coerce").to_numpy(float)
    y2                   = pd.to_numeric(df2["resid"], errors = "coerce").to_numpy(float)
    z2                   = pd.to_numeric(df2["z"],     errors = "coerce").to_numpy(float)

    w2_raw               = pd.to_numeric(df2.get("w_mu", 1.0), errors="coerce").to_numpy(float)
    w2                   = np.where(np.isfinite(w2_raw) & (w2_raw > 0.0), w2_raw, 1.0)

    # Ecliptic angles
    lam, beta = ecliptic_angles_from_radec(RA2, DEC2)

    # Design matrix (ASCII names)
    X                    = pd.DataFrame({
                                         "const"    : 1.0,
                                         "sin_lam"  : np.sin(lam),
                                         "cos_lam"  : np.cos(lam),
                                         "sin_beta" : np.sin(beta),
                                         "z"        : z2,
                                         "sin_ra"   : np.sin(np.deg2rad(RA2)),
                                         "cos_ra"   : np.cos(np.deg2rad(RA2)),
                                       })

    # Optional lunar regressors
    if {"sin_phi_moon", "cos_phi_moon"}.issubset(df2.columns):
        X["sin_phi_moon"] = pd.to_numeric(df2["sin_phi_moon"], errors = "coerce").to_numpy(float)
        X["cos_phi_moon"] = pd.to_numeric(df2["cos_phi_moon"], errors = "coerce").to_numpy(float)

    if calendar_mode == "lunar-lambda" and {"sin_lam_moon", "cos_lam_moon"}.issubset(df2.columns):
        X["sin_lam_moon"] = pd.to_numeric(df2["sin_lam_moon"], errors = "coerce").to_numpy(float)
        X["cos_lam_moon"] = pd.to_numeric(df2["cos_lam_moon"], errors = "coerce").to_numpy(float)

    # Final finite-mask (covers any optional columns)
    mask                  = np.isfinite(y2) & np.isfinite(w2) & np.isfinite(X.to_numpy(dtype = float)).all(axis=1)
    y3, w3, X3            = y2[mask], w2[mask], X.loc[mask].reset_index(drop = True)

    # Guard again
    if len(y3) <= (X3.shape[1] + 5):
        return {
                "amp_ecl"       : np.nan, 
                "phi_ecl_deg"   : np.nan, 
                "se_phi_deg"    : np.nan,
                "p_sin_lam"     : np.nan, 
                "p_cos_lam"     : np.nan,
                "wald_F"        : np.nan, 
                "wald_p"        : np.nan,
                "n_used"        : int(len(y3)),
               }

    # Fit
    res                 = sm.WLS(y3, X3, weights=w3).fit()
    res                 = res.get_robustcov_results(cov_type="HC3")

    names               = list(X3.columns)
    ix                  = {n: i for i, n in enumerate(names)}
    params              = np.asarray(res.params).ravel()
    pvals               = np.asarray(res.pvalues).ravel()

    b_sin               = params[ix["sin_lam"]]
    b_cos               = params[ix["cos_lam"]]

    amp_ecl             = float(np.hypot(b_sin, b_cos))
    phi_ecl_deg         = float(np.degrees(np.arctan2(b_sin, b_cos)))
    
    if phi_ecl_deg < 0:
        phi_ecl_deg    += 360.0

    V                   = np.asarray(res.cov_params(), float)
    vss                 = float(V[ix["sin_lam"], ix["sin_lam"]])
    vcc                 = float(V[ix["cos_lam"], ix["cos_lam"]])
    vsc                 = float(V[ix["sin_lam"], ix["cos_lam"]])

    denom               = (b_sin*b_sin + b_cos*b_cos)
    se_phi_deg          = (
                           float(np.degrees(np.sqrt(max((b_cos*b_cos*vss + b_sin*b_sin*vcc - 2*b_sin*b_cos*vsc) / (denom*denom), 0.0))))
                           if denom > 0 else float("nan")
                          )

    R                   = np.zeros((2, len(names)), float)
    R[0, ix["sin_lam"]] = 1.0
    R[1, ix["cos_lam"]] = 1.0
    wald                = res.wald_test(R, scalar=True)

    p_sin_lam           = float(pvals[ix["sin_lam"]])
    p_cos_lam           = float(pvals[ix["cos_lam"]])

    out = {
            "amp_ecl"     : amp_ecl,
            "phi_ecl_deg" : phi_ecl_deg,
            "se_phi_deg"  : se_phi_deg,
    
            # canonical ASCII keys
            "p_sin_lam"   : p_sin_lam,
            "p_cos_lam"   : p_cos_lam,
    
            "wald_F"      : float(np.asarray(wald.statistic)),
            "wald_p"      : float(np.asarray(wald.pvalue)),
            "n_used"      : int(len(y3)),
          }

    # backward-compat keys used elsewhere in the script
    out["p_sinλ"] = p_sin_lam
    out["p_cosλ"] = p_cos_lam

    return out




def plv_perm(phi, lam_on_L, n = 2000, seed = 11):
    rng                    = np.random.default_rng(seed)
    obs                    = plv(phi, lam_on_L)
    ge                     = 0
    
    for _ in range(n):
        ix                 = rng.permutation(len(lam_on_L))
        
        if plv(phi, lam_on_L[ix]) >= obs:
            ge            += 1
    
    return (ge + 1.0) / (n + 1.0)


def ar1_fit(y):
    y                      = np.asarray(y, float)
    Y                      = y[1:]
    X                      = y[:-1][:, None]
    phi                    = lstsq(X, Y, rcond = None)[0][0]
    eps                    = Y - phi*X.ravel()
    return float(phi), eps


def ar1_boot_surrogate(y, block=None, rng = None):
    rng                    = rng or np.random.default_rng()
    phi, eps               = ar1_fit(y)
    eps_b                  = block_shuffle(eps, block = block, rng = rng)
    yb                     = np.empty_like(y, float)
    yb[0]                  = y[0]
    
    for t in range(1, len(y)):
        yb[t]              = phi*yb[t-1] + eps_b[t-1]
        
    # variance match to original
    yb                    *= (np.std(y, ddof=1) or 1.0)/(np.std(yb, ddof=1) or 1.0)
    return yb

def build_on_grid_key(df, L_grid, *, mode="solar", lam_on_L=None):
    """
    Returns a calendar key defined on the L-grid (same length as y on that grid),
    or None if not available. For 'solar', lam_on_L is already the desired key.
    """
    if mode == "none":
        return None
    
    if mode == "solar":
        return np.asarray(lam_on_L, float) if lam_on_L is not None else None
    
    if mode == "lunar-phase":
        if "phi_moon" not in df.columns:
            return None
        
        Lp                  = np.log1p(pd.to_numeric(df["z"], errors="coerce").to_numpy(float))
        o                   = np.argsort(Lp)
        phi                 = pd.to_numeric(df["phi_moon"],   errors="coerce").to_numpy(float)
        return interpolate_angle_over_L(Lp[o], phi[o], L_grid)
    
    if mode == "lunar-lambda":
        if "lam_moon" not in df.columns:
            return None
        
        Lp                  = np.log1p(pd.to_numeric(df["z"], errors="coerce").to_numpy(float))
        o                   = np.argsort(Lp)
        lm                  = pd.to_numeric(df["lam_moon"],   errors="coerce").to_numpy(float)
        return interpolate_angle_over_L(Lp[o], lm[o], L_grid)
    raise ValueError(f"Unknown calendar mode: {mode}")


def running_ecl_amp(df, Lwin=0.35, Lstep=0.05):
    L                      = np.log1p(pd.to_numeric(df["z"], errors="coerce").to_numpy(float))
    ok                     = np.isfinite(L)
    df                     = df.loc[ok, ["resid","RA","DEC","z","w_mu"]].reset_index(drop = True)
    L                      = np.log1p(pd.to_numeric(df["z"], errors="coerce").to_numpy(float))
    Lmin, Lmax             = float(np.min(L)), float(np.max(L))
    centers                = np.arange(Lmin+Lwin/2, Lmax-Lwin/2, Lstep)
    rows                   = []
    
    for c in centers:
        sel                = (L>=c-Lwin/2)&(L<=c+Lwin/2)
        
        if sel.sum() < 150: 
            continue
        
        reg                = extended_regression(df.loc[sel])
        rows.append({"L_center":c, "amp_ecl":reg["amp_ecl"], "wald_p":reg["wald_p"], "phi":reg["phi_ecl_deg"]})
    
    return pd.DataFrame(rows)

    
# --- Quick synthetic sanity check ---------------------------------------
def sanity_check():
    """
    Build a synthetic dataset that is deliberately 'locked' to ecliptic λ,
    then run the key steps and assert they behave sensibly.
    """
    print("\n[SANITY] building synthetic locked dataset...")
    
    # inside sanity_check() only
    print("[SANITY] spin null: preserve_L_rank=False")

    rng             = np.random.default_rng(42)

    # monotonic L grid and corresponding z
    n               = 2000
    L               = np.linspace(np.log1p(0.01), np.log1p(1.0), n)
    z               = np.expm1(L)

    # target ecliptic coordinates deliberately locked to L (monotonic λ)
    lam_true        = np.linspace(0, 2*np.pi, n, endpoint = False)
    lam_true       += 0.35*np.sin(3*lam_true)
    beta_true       = rng.uniform(-np.deg2rad(35), np.deg2rad(35), size=n)
    
    # RA/DEC from true ecliptic coords (robust; no obstime usage)
    RA, DEC         = ecliptic_to_icrs(np.degrees(lam_true), np.degrees(beta_true))
    lam_back, _     = ecliptic_angles_from_radec(RA, DEC)
    err             = np.max(np.abs(np.angle(np.exp(1j*(lam_back - lam_true)))))
    print(f"[SANITY] max |λ round-trip error| = {err:.2e} rad")

    # should be clearly > 0 if spins randomize λ relative to original
    R0              = random_rotation_matrix(np.random.default_rng(0))  # SANITY ONLY (intentional fixed seed)
    ra_s, dec_s     = rotate_radec(RA, DEC, R0)
    lam_s, _        = ecliptic_angles_from_radec(ra_s, dec_s)
    dlam_std        = float(np.std(np.angle(np.exp(1j*(lam_s - lam_back)))))
    print(f"[CHECK]  std Δλ after one random spin = {dlam_std:.3f} rad")

    # residuals with λ signal + mild z trend + noise
    y_true          = 0.9*np.sin(lam_back) + 0.6*np.cos(lam_back) + 0.3*(z - z.mean())
    resid           = y_true + rng.normal(0, 0.25, n)

    df_raw          = pd.DataFrame({"RA":RA, "DEC":DEC, "z":z, "resid":resid, "w_mu":1.0})

    # simple linear detrend on z to create "detrended"
    X               = np.c_[np.ones(n), z]
    B               = lstsq(X, resid, rcond = None)[0]
    resid_det       = resid - X @ B
    df_det          = pd.DataFrame({"RA"              : RA, 
                                    "DEC"             : DEC, 
                                    "z"               : z,
                                    "resid_detrended" : resid_det, 
                                    "w_mu"            : 1.0})
    
    
    if CALENDAR_MODE in ("lunar-phase", "lunar-lambda"):
        if MJD_COL not in df_raw.columns or MJD_COL not in df_det.columns:
            raise KeyError(
                           f"HHT_CALENDAR={CALENDAR_MODE} requires an MJD column. "
                           f"Expected '{MJD_COL}' (override with HHT_MJD_COL)."
                          )
        df_raw      = add_lunar_columns(df_raw, mjd_col=MJD_COL)
        df_det      = add_lunar_columns(df_det, mjd_col=MJD_COL)


    # build uniform L grids
    Lr, yr          = to_uniform_L(df_raw["z"], 
                                   df_raw["resid"], 
                                   pad   = 0.1,  
                                   ngrid = 512)
    
    Ld, yd          = to_uniform_L(df_det["z"], 
                                   df_det["resid_detrended"], 
                                   pad   = 0.1, 
                                   ngrid = 512)

    # CEEMDAN phases (lighter settings for sanity)
    imfs_r          = ceemdan_imfs(Lr, yr, 
                                   trials = 50, 
                                   noise  = 0.2, 
                                   seed   = 7)
    
    imfs_d          = ceemdan_imfs(Ld, yd, 
                                   trials = 50, 
                                   noise  = 0.2, 
                                   seed   = 7)
    
    phi_r           = inst_phase_from_imf2(imfs_r)
    phi_d           = inst_phase_from_imf2(imfs_d)

    # λ on L-grids
    lam_r_on_Lr     = interpolate_angle_over_L(np.log1p(z), lam_back, Lr)
    lam_d_on_Ld     = interpolate_angle_over_L(np.log1p(z), lam_back, Ld)


    # --- Calendar keys on the L-grids (so key length matches y length) ---
    CAL_KEY_RAW_L   = CAL_KEY_DET_L = None
    
    CAL_MODE_EFF = CALENDAR_MODE

    if CAL_MODE_EFF in ("lunar-phase", "lunar-lambda"):
        need = "phi_moon" if CAL_MODE_EFF == "lunar-phase" else "lam_moon"
        if need not in df_raw.columns or need not in df_det.columns:
            print(f"[calendar] requested {CALENDAR_MODE} but '{need}' missing; falling back to solar.", flush=True)
            CAL_MODE_EFF = "solar"    

    
    if CAL_MODE_EFF == "solar":
        CAL_KEY_RAW_L       = lam_r_on_Lr
        CAL_KEY_DET_L       = lam_d_on_Ld
    
    elif CAL_MODE_EFF == "lunar-phase":
        
        
        if "phi_moon" in df_raw.columns and "phi_moon" in df_det.columns:
            phi_raw         = pd.to_numeric(df_raw["phi_moon"],          errors = "coerce").to_numpy(float)
            phi_det         = pd.to_numeric(df_det["phi_moon"],          errors = "coerce").to_numpy(float)
            L_raw_pts       = np.log1p(df_raw["z"].to_numpy(float));     o_raw  = np.argsort(L_raw_pts)
            L_det_pts       = np.log1p(df_det["z"].to_numpy(float));     o_det  = np.argsort(L_det_pts)
            CAL_KEY_RAW_L   = interpolate_angle_over_L(L_raw_pts[o_raw], phi_raw[o_raw], Lr)
            CAL_KEY_DET_L   = interpolate_angle_over_L(L_det_pts[o_det], phi_det[o_det], Ld)

    elif CAL_MODE_EFF == "lunar-lambda":
        if "lam_moon" in df_raw.columns and "lam_moon" in df_det.columns:
            lm_raw          = pd.to_numeric(df_raw["lam_moon"],          errors = "coerce").to_numpy(float)
            lm_det          = pd.to_numeric(df_det["lam_moon"],          errors = "coerce").to_numpy(float)
            L_raw_pts       = np.log1p(df_raw["z"].to_numpy(float));     o_raw  = np.argsort(L_raw_pts)
            L_det_pts       = np.log1p(df_det["z"].to_numpy(float));     o_det  = np.argsort(L_det_pts)
            
            CAL_KEY_RAW_L   = interpolate_angle_over_L(L_raw_pts[o_raw], lm_raw[o_raw], Lr)
            CAL_KEY_DET_L   = interpolate_angle_over_L(L_det_pts[o_det], lm_det[o_det], Ld)

    # core metrics
    plv_r                   = plv(phi_r, lam_r_on_Lr)
    plv_d                   = plv(phi_d, lam_d_on_Ld)

    # --- Use the keys: PLV seasonal block-λ keyed null + keyed BLOCK energy test ---
    try:
        shufs               = plv_lambda_block_shuffle(
                                                       lam_d_on_Ld,
                                                       block_len =  32,
                                                       n         = 512,
                                                       seed      =   7
                                                      )

        obs                 = plv_d
        ge                  = sum(plv(phi_d, s) >= obs for s in shufs)
        p_block_keyed       = (ge + 1.0) / (len(shufs) + 1.0)
        print(f"[SANITY] keyed block-λ PLV p (det) = {p_block_keyed:.3g}  [n={len(shufs)}]")
    
    except Exception as e:
        print(f"[SANITY] keyed block-λ PLV skipped: {e}")
    
    try:
        # light settings to keep sanity_check() fast
        block_chk           = e_band_surrogate_test(
                                                    Ld, yd, 
                                                    band         = "imf23", 
                                                    mode         = "block",
                                                    ns           = 80, 
                                                    trials       = 30, 
                                                    noise        = 0.2, 
                                                    seed         = 7,
                                                    block_len    = 32, 
                                                    standardize  = True, 
                                                    use_fraction = True,
                                                    cal_key      = CAL_KEY_DET_L, 
                                                    tag          = "det[sanity]"
                                                   )
        print(f"[SANITY] keyed BLOCK E(IMF2+3) p = {block_chk['p']:.3g}  metric={block_chk['metric']}")
    
    except Exception as e:
        print(f"[SANITY] keyed BLOCK surrogate skipped: {e}")

    # --- spin null: break L–rank to ensure a hard negative in sanity ---
    print("[SANITY] spin null: preserve_L_rank=False")
    sp                      = spin_test_plv(
                                            df_det.rename(columns = {"resid_detrended":"resid"}),
                                                                      phi_on_L        = phi_d,
                                                                      L_grid          = Ld,
                                                                      n               = 4000,      # higher N to reduce flakiness
                                                                      seed            = 7,
                                                                      preserve_L_rank = False
                                                                     )
    
    spin_p                  = ((sp >= plv_d).sum() + 1.0) / (len(sp) + 1.0)
    ray_p                   = rayleigh_p(np.mod(phi_d - lam_d_on_Ld + np.pi, 2*np.pi) - np.pi)

    reg_det                 = extended_regression(pd.DataFrame({
                                                                "resid" : resid_det, 
                                                                "RA"    : RA, 
                                                                "DEC"   : DEC, 
                                                                "z"     : z, 
                                                                "w_mu"  : 1.0
                                                              }))

    print("\n[SANITY results]")
    print(f" PLV raw        = {plv_r:.3f}")
    print(f" PLV det        = {plv_d:.3f}")
    print(f" spin p (det)   = {spin_p:.3g}  [n={len(sp)}]")
    print(f" Rayleigh p     = {ray_p:.3g}")
    print(f" Wald p (det)   = {reg_det['wald_p']:.3g}, amp_ecl={reg_det['amp_ecl']:.3f}, φ_ecl={reg_det['phi_ecl_deg']:.1f}°")

    # stricter, but still robust to tiny numerical drift
    assert 0.0 <= plv_d <= 1.0
    assert spin_p < 0.01,             "Spin p should be < 0.01 on this locked synthetic set"
    assert ray_p  < 1e-6,             "Rayleigh p should be essentially zero on this set"
    assert reg_det["wald_p"] < 1e-10, "Directional regression must be ultra-significant"

    # (optional) quick uniformity canary for spun λ
    # from scipy.stats import kstest
    # lam_u = (np.mod(lam_d_on_Ld, 2*np.pi) / (2*np.pi))
    # print("[CHECK] KS p (λ uniform under spin) =", kstest(lam_u, "uniform").pvalue)

    # Log a row into the V4 ledger so we can spot regressions over time
    try:
        write_summary_V5({
                          "E2_raw"           : np.nan, 
                          "E2_det"           : np.nan, 
                          "p_raw"            : np.nan, 
                          "p_det"            : np.nan,
                          "PLV_raw_lambda"   : plv_r, 
                          "PLV_det_lambda"   : plv_d,
                          "PLV_p_det_spin"   : spin_p,
                          "PLV_spin_med_det" : float(np.median(sp)),
                          "PLV_spin_q95_det" : float(np.quantile(sp, 0.95)),
                          "Rayleigh_p_det"   : ray_p,
                          "amp_ecl_det"      : reg_det["amp_ecl"], 
                          "wald_p_det"       : reg_det["wald_p"],
                          "phi_ecl_det_deg"  : reg_det["phi_ecl_deg"], 
                          "se_phi_det_deg"   : reg_det["se_phi_deg"],
                          "n_spin"           : len(sp), 
                          "seed"             : 7, 
                          "trials"           : 50, 
                          "noise"            : 0.2,
                          "data_raw"         : "[sanity]", 
                          "data_det"         : "[sanity]"
                         }, out_dir=Path.cwd())
    
    except Exception as _e:
        print(f"[SANITY] summary write skipped: {_e}")

    print("[SANITY] OK")
    

def _get_int_or_none(name,  default = None):
    v               = globals().get(name, default)

    if v is None:
        return None

    if isinstance(v, str) and v.strip().lower() == "none":
        return None

    try:
        return int(v)
    
    except Exception:
        return None if default is None else int(default)


def _dominant_period_in_L(imf, L):
    """
    Rough dominant period of an IMF in units of L (log1p(z)).
    Returns np.inf if no clear frequency.
    """
    x                   = np.asarray(imf, float)
    L                   = np.asarray(L, float)
    n                   = x.size
    if n < 8:
        return np.inf

    x                   = x - np.mean(x)
    dL                  = np.median(np.diff(L))
    if not np.isfinite(dL) or dL <= 0:
        return np.inf

    # FFT on uniform L-grid
    F                   = np.fft.rfft(x)
    A                   = np.abs(F)
    
    if A.size <= 1:
        return np.inf

    A[0]                = 0.0  # remove DC
    freqs               = np.fft.rfftfreq(n, d=dL)  # cycles per unit L

    k                   = int(np.argmax(A))
    f                   = float(freqs[k])
    
    if not np.isfinite(f) or f <= 0:
        return np.inf
    return 1.0 / f

def _imf_phase(imf):
    """Wrapped phase from analytic signal; unwrap not required for PLV."""
    an                  = hilbert(np.asarray(imf, float))
    return np.angle(an)

def decomposition_report(tag, L_grid, y_grid, imfs, target_ang = None, topk = 8):
    """
    Prints an IMF “ladder”:
      - energy fraction per IMF (and residual)
      - PLV(IMFi phase vs target_ang) if provided
      - dominant period in L units
    """
    imfs                = np.asarray(imfs,   float)
    y                   = np.asarray(y_grid, float)
    L                   = np.asarray(L_grid, float)

    if imfs.ndim != 2 or imfs.shape[1] != y.size:
        print(f"[DECOMP] {tag}: shape mismatch: imfs={imfs.shape}, y={y.shape}")
        return

    n_imf               = imfs.shape[0]
    recon               = np.sum(imfs, axis=0)
    resid               = y - recon

    # Use decomposition-consistent energy accounting (IMFs + resid)
    e_imf               = np.sum(imfs**2, axis=1)
    e_res               = float(np.sum(resid**2))
    den                 = float(np.sum(e_imf) + e_res + 1e-30)

    frac_imf            = e_imf / den
    frac_res            = e_res / den

    rows                = []
    for i in range(n_imf):
        periodL         = _dominant_period_in_L(imfs[i], L)

        plv_i = np.nan
        if target_ang is not None:
            try:
                phi_i   = _imf_phase(imfs[i])
                plv_i   = plv(phi_i, np.asarray(target_ang, float))
            except Exception:
                plv_i   = np.nan

        rows.append((i+1, frac_imf[i], plv_i, periodL))

    # Sort primarily by PLV if available, else by energy
    if target_ang is not None and np.any(np.isfinite([r[2] for r in rows])):
        rows_sorted     = sorted(rows, key=lambda r: (-(r[2] if np.isfinite(r[2]) else -np.inf), -r[1]))
        best            = rows_sorted[0]
        best_note       = f"best_by_PLV=IMF{best[0]} (PLV={best[2]:.4f}, Efrac={best[1]*100:.2f}%)"
    else:
        rows_sorted     = sorted(rows, key=lambda r: -r[1])
        best            = rows_sorted[0]
        best_note       = f"best_by_energy=IMF{best[0]} (Efrac={best[1]*100:.2f}%)"

    print(f"\n[DECOMP] {tag}: n_imf={n_imf}  resid_Efrac={frac_res*100:.2f}%  {best_note}")
    hdr                 = " idx   Efrac(%)    PLV(vs target)     dom_period_L"
    print(hdr)
    print("-"*len(hdr))

    for (idx, ef, plv_i, pL) in rows_sorted[:min(topk, len(rows_sorted))]:
        ef_pct          = 100.0 * ef
        plv_txt         = f"{plv_i:>10.4f}" if np.isfinite(plv_i)            else "      (na)"
        pL_txt = f"{pL:>12.4f}"             if np.isfinite(pL) and pL < 1e12 else f"{'(inf)':>12}"
        print(f" IMF{idx:<2d}  {ef_pct:>8.2f}     {plv_txt}       {pL_txt}")
    
def _movavg(x, w):
    x                   = np.asarray(x, float)
    w                   = int(w)
    if w <= 1 or x.size < 3:
        return x
    
    w                   = min(w, x.size)
    k                   = np.ones(w, dtype=float) / w
    return np.convolve(x, k, mode="same")


def plot_imf_energy_diagnostics(tag, L_grid, y_grid, imfs, out_dir,
                                *, smooth = 9, topk = 12, dpi = 160, use_hilbert = False):
    """
    Saves:
      (1) heatmap of instantaneous IMF power vs L (great for chirps / energy packets)
      (2) bar chart of total energy fraction per IMF (+ residual)

    use_hilbert:
      - False: power ~ imf^2
      - True : power ~ |hilbert(imf)|^2 (envelope power)
    """
    from   pathlib           import Path
    import matplotlib
    
    matplotlib.use("Agg", force = True)
    import matplotlib.pyplot as plt

    L                   = np.asarray(L_grid, float)
    y                   = np.asarray(y_grid, float)
    imfs                = np.asarray(imfs, float)

    if imfs.ndim != 2 or imfs.shape[1] != y.size:
        print(f"[IMF_PLOT] {tag}: shape mismatch: imfs = {imfs.shape}, y = {y.shape}", flush=True)
        return None

    n_imf               = imfs.shape[0]
    recon               = np.sum(imfs, axis = 0)
    resid               = y - recon

    # instantaneous power per IMF (optionally via Hilbert envelope)
    if use_hilbert:
        # envelope power = |analytic|^2
        pow_imf         = np.empty_like(imfs)
        
        for i in range(n_imf):
            an          = hilbert(imfs[i])
            pow_imf[i]  = (np.abs(an) ** 2)
            
    else:
        pow_imf         = imfs  ** 2

    pow_res             = resid ** 2

    # optional smoothing along L
    if smooth and smooth > 1:
        for i in range(n_imf):
            pow_imf[i]  = _movavg(pow_imf[i], smooth)
        pow_res         = _movavg(pow_res, smooth)

    # total energy fractions (decomp-consistent: IMFs + residual)
    e_imf               = np.sum(pow_imf, axis=1)
    e_res               = float(np.sum(pow_res))
    den                 = float(np.sum(e_imf) + e_res + 1e-30)
    frac_imf            = e_imf / den
    frac_res            = e_res / den

    # limit displayed IMFs (top by total energy)
    topk                = int(min(max(1, topk), n_imf))
    order               = np.argsort(-frac_imf)[:topk]  # indices into imfs (0-based)
    pow_show            = pow_imf[order]

    # heatmap uses log power for dynamic range
    eps                 = 1e-30
    Z                   = np.log10(pow_show + eps)

    out_dir             = Path(out_dir)
    out_dir.mkdir(parents = True, exist_ok = True)

    base = f"hht_imf_energy_{tag.lower()}_seed{SEED}_B{BLOCK_LEN}_cal{CALENDAR_MODE}"
    if use_hilbert:
        base += "_hilb"

    # ---------- (1) HEATMAP ----------
    fig                 = plt.figure(figsize=(11, 4.8))
    ax                  = fig.add_subplot(111)

    # y-axis is "rank within topk" but label as IMF number (original index)
    extent              = [float(L.min()), float(L.max()), 0.5, topk + 0.5]
    im = ax.imshow(Z, aspect = "auto", origin = "lower", extent = extent)

    ax.set_xlabel("L = log(1+z)")
    ax.set_ylabel("IMF (ranked by total energy)")
    ax.set_title(f"IMF instantaneous power heatmap ({tag})  smooth={smooth}  "
                 f"{'Hilbert env^2' if use_hilbert else 'imf^2'}")

    # tick labels show actual IMF numbers
    ax.set_yticks(np.arange(1, topk + 1))
    ax.set_yticklabels([f"IMF{int(i+1)}" for i in order])

    cbar              = fig.colorbar(im, ax = ax)
    cbar.set_label("log10(power)")

    heat_path         = out_dir / f"{base}_heat.png"
    fig.tight_layout()
    fig.savefig(heat_path, dpi = int(dpi))
    plt.close(fig)

    # ---------- (2) BAR (TOTAL ENERGY FRACTIONS) ----------
    fig                        = plt.figure(figsize=(11, 4.2))
    ax                         = fig.add_subplot(111)

    x                          = np.arange(1, n_imf + 1)
    ax.bar(x,  100.0 * frac_imf)
    ax.axhline(100.0 * frac_res, linestyle="--")
    ax.set_xlabel("IMF index")
    ax.set_ylabel("Total energy fraction (%)")
    ax.set_title(f"Total IMF energy fractions ({tag})  resid = {100.0*frac_res:.2f}%")

    bar_path = out_dir / f"{base}_bar.png"
    fig.tight_layout()
    fig.savefig(bar_path, dpi=int(dpi))
    plt.close(fig)

    print(f"[IMF_PLOT] wrote {heat_path}", flush=True)
    print(f"[IMF_PLOT] wrote {bar_path}",  flush=True)

    return str(heat_path), str(bar_path)

    
# --------------------------------- MAIN ---------------------------------- #
if __name__ == "__main__":

    _ping("\n\nentered __main__")
    print("...", flush = True)

    found     = maybe_detect_csv_paths_from_target()

    # — quick synthetic smoke test toggle —
    if os.environ.get("HHT_SANITY", "0") == "1":
        sanity_check()
        raise SystemExit
    
    NOISE     = float(globals().get("NOISE" , 0.2))
    PAD       = float(globals().get("PAD"   , 0.15))
    
    NGRID     = _get_int("NGRID"            , 384)
    
    
    print(f"seed = {SEED}  CE_TRIALS = {CE_TRIALS}  noise = {NOISE}  "
          f"BLOCK_LEN = {BLOCK_LEN}  pad = {PAD}  ngrid = {NGRID}")

    
    # -------------------- unify config (NO hard-coded overrides) --------------------
    
    # Canonical run parameters (upper-case)
    # Keep these names as the source of truth everywhere.
    # If later code expects lower-case names, provide aliases (below).
    # -------------------------------------------------------------------------------
    # SEED, TRIALS, N_SPIN, BLOCK_LEN, NOISE, PAD, NGRID already set above
    
    # Aliases for legacy code blocks (no need to refactor everything at once)
    seed_sur = SEED
    trials   = TRIALS
    noise    = NOISE
    pad      = PAD
    ngrid    = NGRID
    
    # If you have a distinct "surrogate count" used in some sections, define it explicitly.
    # Otherwise, keep it tied to  TRIALS.
    NS_SUR   = _get_int("NS_SUR", TRIALS)
    ns_sur   = NS_SUR
    
    can_imfs_raw = None
    can_imfs_det = None
    
    
    # ---------- small pretty-print helpers (align '=') ------------------- #
    def _aline(pairs, kpad = None, sep = "  "):
        """
        pairs: list of (key, value_str)
        returns a single line with keys left-padded so all '=' align.
        """
        w        = kpad or max(len(str(k)) for k, _ in pairs)
        
        return sep.join(f"{str(k):<{w}} = {v}" for k, v in pairs)


    def _print2(label_left, val_left, label_right, val_right):
        """
        Two aligned key=value items on one line with same key width.
        """
        w         = max(len(label_left), len(label_right))
        
        print(_aline([(label_left, val_left), (label_right, val_right)], kpad = w, sep = "   "))

    DEFAULT_P_RAW = Path("/Users/boyde/.spyder-py3/plamb_runs/hht_bridge/FR_log1pz_hht_export_V5/sn_residuals_enriched.csv")
    DEFAULT_P_DET = Path("/Users/boyde/.spyder-py3/plamb_runs/hht_bridge/FR_log1pz_hht_export_V5/sn_residuals_detrended.csv")
    
    # env overrides
    p_raw = Path(os.environ.get("HHT_P_RAW", str(DEFAULT_P_RAW))).expanduser()
    p_det = Path(os.environ.get("HHT_P_DET", str(DEFAULT_P_DET))).expanduser()
    
    # target-script autodetect overrides env/defaults (if present)
    if found:
        p_raw = Path(found.get("sn_residuals_enriched.csv", p_raw)).expanduser()
        p_det = Path(found.get("sn_residuals_detrended.csv", p_det)).expanduser()

    
    # --- output folder: prefer explicit, else use det folder, else script folder ---
    out_dir = Path(os.environ.get("HHT_OUTDIR", str(p_det.parent if p_det.exists() else Path(__file__).resolve().parent))).expanduser()
    out_dir.mkdir(parents = True, exist_ok = True)

    _ping(f"loaded files: raw={p_raw.name}, det={p_det.name}")

    print(_aline([
                  ("seed",      str(seed_sur)),
                  ("trials",    str(trials)),
                  ("noise",     f"{noise}"),
                  ("BLOCK_LEN", str(BLOCK_LEN)),
                  ("pad",       f"{pad}"),
                  ("ngrid",     str(ngrid)),
                 ]))


    # ----------- Load ----------------------------------------------------- #
    df_raw              = pd.read_csv(p_raw)
    df_det              = pd.read_csv(p_det)
    
    # === calendar setup =======================================================
    # Ensure solar λ available as a column for convenience (radians)
    try:
        lam_raw_col, _  = ecliptic_angles_from_radec(df_raw["RA"], df_raw["DEC"])
        df_raw["lam"]   = np.mod(lam_raw_col, _TWOPI)
        
    except Exception:
        pass
    
    try:
        lam_det_col, _  = ecliptic_angles_from_radec(df_det["RA"], df_det["DEC"])
        df_det["lam"]   = np.mod(lam_det_col, _TWOPI)
    
    except Exception:
        pass
    
    # Add lunar columns (no-op if MJD missing)
    add_lunar_columns(df_raw, mjd_col = MJD_COL, lam_sun_col = ("lam" if "lam" in df_raw.columns else None))
    add_lunar_columns(df_det, mjd_col = MJD_COL, lam_sun_col = ("lam" if "lam" in df_det.columns else None))

    
    if WRITE_LUNAR_COLS:
       df_raw.to_csv(p_raw.with_name(f"{p_raw.stem}_enriched_lunar.csv"),  index = False)
       df_det.to_csv(p_det.with_name(f"{p_det.stem}_detrended_lunar.csv"), index = False)

    # Calendar keys for keyed block nulls
    try:
        CAL_KEY_RAW     = calendar_key_from_df(df_raw, mode = CALENDAR_MODE) if CALENDAR_MODE != "none" else None
        CAL_KEY_DET     = calendar_key_from_df(df_det, mode = CALENDAR_MODE) if CALENDAR_MODE != "none" else None
        
    except Exception as _e:
        print(f"[lunar] warning: could not build calendar key ({_e}). Falling back to unkeyed.", flush = True)
        CAL_KEY_RAW     = CAL_KEY_DET = None
    
    print(f"[calendar] mode = {CALENDAR_MODE}   BLOCK_LEN = {BLOCK_LEN}")


    # Sanity: ensure required columns
    y_det_col           = pick_det_ycol(df_det)
    need_raw            = {"z", "resid", "RA",  "DEC"}
    need_det            = {"z", "RA",    "DEC", y_det_col}

    if not need_raw.issubset(df_raw.columns):
        raise ValueError(f"{p_raw} must contain {need_raw}; has {list(df_raw.columns)}")
    
    if not need_det.issubset(df_det.columns):
        raise ValueError(f"{p_det} must contain {need_det}; has {list(df_det.columns)}")

    # ----------- Build uniform L-grids ----------------------------------- #
    Lr, yr              = to_uniform_L(df_raw["z"], 
                                       df_raw["resid"],   
                                       pad         = pad, 
                                       ngrid       = ngrid)
    
    Ld, yd              = to_uniform_L(df_det["z"], 
                                       df_det[y_det_col], 
                                       pad         = pad, 
                                       ngrid       = ngrid)

    # ----------- E2 surrogate tests -------------------------------------- #
    s_raw               = e2_surrogate_test(Lr, yr, 
                                            trials = CE_TRIALS, 
                                            noise  = noise, 
                                            ns     = N_SUR_E2, 
                                            seed   = seed_sur)
    
    s_det               = e2_surrogate_test(Ld, yd, 
                                            trials = CE_TRIALS, 
                                            noise  = noise, 
                                            ns     = N_SUR_E2, 
                                            seed   = seed_sur)

    print("\n[E2 surrogate test]")
    print(_aline([
                  ("raw.E_obs ", f"{s_raw['E_obs']:.6f}"),
                  ("raw.p",     f"{s_raw['p']:.4g}"),
                  ("raw.S[q05,q50,q95]", f"[{s_raw['s_q05']:.6f}, {s_raw['s_q50']:.6f}, {s_raw['s_q95']:.6f}]"),
                  ("raw.n_valid", f"{s_raw['n_valid']}"),
                 ], kpad = 18))
    
    print(_aline([
                  ("det.E_obs ", f"{s_det['E_obs']:.6f}"),
                  ("det.p",     f"{s_det['p']:.4g}"),
                  ("det.S[q05,q50,q95]", f"[{s_det['s_q05']:.6f}, {s_det['s_q50']:.6f}, {s_det['s_q95']:.6f}]"),
                  ("det.n_valid", f"{s_det['n_valid']}"),
                 ], kpad = 18))

    # ----------- IMF2 phases --------------------------------------------- #
    imfs_r                  = ceemdan_imfs(Lr, yr, 
                                           trials = CE_TRIALS, 
                                           noise  = noise, 
                                           seed   = seed_sur)
    
    imfs_d                  = ceemdan_imfs(Ld, yd, 
                                           trials = CE_TRIALS, 
                                           noise  = noise, 
                                           seed   = seed_sur)

    phi_r                   = inst_phase_from_imf2(imfs_r)
    phi_d                   = inst_phase_from_imf2(imfs_d)

    # Ecliptic λ on the two grids (proper interpolation of angles)
    lam_raw, _              = ecliptic_angles_from_radec(df_raw["RA"], df_raw["DEC"])
    L_raw_pts               = np.log1p(df_raw["z"].to_numpy(float))
    o_raw                   = np.argsort(L_raw_pts)
    lam_r_on_Lr             = target_angle_on_L(df_raw, 
                                                Lr, 
                                                mode = CALENDAR_MODE)
    
    lam_det, _              = ecliptic_angles_from_radec(df_det["RA"], df_det["DEC"])
    L_det_pts               = np.log1p(df_det["z"].to_numpy(float))
    o_det                   = np.argsort(L_det_pts)
    lam_d_on_Ld             = target_angle_on_L(df_det, 
                                                Ld, 
                                                mode = CALENDAR_MODE)
    
    # ----------- PLV vs ecliptic + permutation p ------------------------- #
    plv_r_lambda            = plv(phi_r,           lam_r_on_Lr)
    plv_d_lambda            = plv(phi_d,           lam_d_on_Ld)
    
    p_plv_raw               = plv_perm(phi_r,      lam_r_on_Lr, n = N_PERM,  seed = seed_sur)
    p_plv_det               = plv_perm(phi_d,      lam_d_on_Ld, n = N_PERM,  seed = seed_sur)
    p_plv_shift_det         = _circ_shift_p(phi_d, lam_d_on_Ld, n = N_SHIFT, seed = seed_sur)
    
    # Seasonal-block λ shuffle p-value (DET)
    lam_block_shufs         = plv_lambda_block_shuffle(lam_d_on_Ld, 
                                                       block_len = BLOCK_LEN, 
                                                       n         = N_SHIFT, 
                                                       seed      = seed_sur)
      
    # ================= IMF-SCAN PLV (max-over-IMFs, corrected nulls) =================
    HHT_IMF_SCAN = os.environ.get("HHT_IMF_SCAN", "1") == "1"                  # set "0" to disable
    SCAN_EFRAC   = float(os.environ.get("HHT_IMF_SCAN_EFRAC", "0.00"))
    SCAN_START   = int(os.environ.get(  "HHT_IMF_SCAN_START", "2"))              # start at IMF2 by default
    
    if HHT_IMF_SCAN:
        # phases for every IMF (DET)
        phi_d_by_imf             = phases_by_imf(imfs_d)
    
        # choose which IMFs to scan
        scan_imfs_det, efrac_det = pick_scan_imfs(imfs_d, 
                                                  start_imf = SCAN_START, 
                                                  efrac_min = SCAN_EFRAC)
    
        # observed scan statistic
        plvs_det_vec, plv_det_max, best_imf_det0 = plv_scan(phi_d_by_imf, 
                                                            lam_d_on_Ld, 
                                                            scan_imfs_det)
        
        best_imf_det                             = (best_imf_det0 + 1) if best_imf_det0 is not None else None  # 1-based label
    
        print("\n[IMF-scan PLV] DET")
        print(f" scan_imfs = {[i+1 for i in scan_imfs_det]}  (start=IMF{SCAN_START}, efrac_min={SCAN_EFRAC})")
        print(f" obs max PLV = {plv_det_max:.6f}   best IMF = {best_imf_det}")
    
        # permutation max-stat p
        obs_perm, best_perm, p_perm_scan, null_perm = plv_scan_perm_p(
                                                                      phi_d_by_imf, 
                                                                      lam_d_on_Ld, 
                                                                      scan_imfs_det, 
                                                                      n    = N_PERM, 
                                                                      seed = seed_sur
                                                                     )
    
        # circular-shift max-stat p
        obs_shift, best_shift, p_shift_scan, null_shift = plv_scan_shift_p(
                                                                           phi_d_by_imf, 
                                                                           lam_d_on_Ld, 
                                                                           scan_imfs_det, 
                                                                           n    = N_SHIFT, 
                                                                           seed = seed_sur
                                                                          )
    
        # block-λ max-stat p (if available)
        if isinstance(lam_block_shufs, np.ndarray) and lam_block_shufs.ndim == 2:
            obs_blk, best_blk, p_block_scan, null_blk = plv_scan_blocklam_p(
                                                                            phi_d_by_imf, 
                                                                            lam_block_shufs, 
                                                                            scan_imfs_det, 
                                                                            obs_theta = lam_d_on_Ld
                                                                           )
        else:
            p_block_scan = np.nan
    
        # spin max-stat p
        sp_scan = spin_test_plv_scan(
                                     df_det,
                                     phi_d_by_imf,
                                     Ld,
                                     scan_imfs_det,
                                     n               = N_SPIN,
                                     seed            = seed_sur,
                                     preserve_L_rank = False
                                    )
        p_spin_scan         = ((sp_scan >= plv_det_max).sum() + 1.0) / (len(sp_scan) + 1.0)
    
        print("[IMF-scan PLV p-values] DET (look-elsewhere corrected by max-stat)")
        print(f" perm   p = {p_perm_scan:.6g}")
        print(f" shift  p = {p_shift_scan:.6g}")
        print(f" blockλ p = {p_block_scan:.6g}")
        print(f" spin   p = {p_spin_scan:.6g}")
    # ================================================================================

    # === IMF energy diagnostics (env-driven) =================================
    HHT_IMF_PLOTS           =     os.environ.get("HHT_IMF_PLOTS", "0") ==  "1"
    HHT_IMF_PLOT_TOPK       = int(os.environ.get("HHT_IMF_PLOT_TOPK",     "12"))        # max IMF rows to show
    HHT_IMF_PLOT_SMOOTH     = int(os.environ.get("HHT_IMF_PLOT_SMOOTH",    "9"))        # moving-average window (points)
    HHT_IMF_PLOT_HILBERT    =     os.environ.get("HHT_IMF_PLOT_HILBERT",   "0") == "1"  # use envelope^2 instead of imf^2
    HHT_IMF_PLOT_DPI        = int(os.environ.get("HHT_IMF_PLOT_DPI",     "160"))
    
    
    # ---- IMF energy diagnostic plots (optional) ----
    if HHT_IMF_PLOTS:
        plot_imf_energy_diagnostics(
                                    "RAW", Lr, yr, imfs_r, 
                                    out_dir         = p_det.parent,
                                    smooth          = HHT_IMF_PLOT_SMOOTH,
                                    topk            = HHT_IMF_PLOT_TOPK,
                                    dpi             = HHT_IMF_PLOT_DPI,
                                    use_hilbert     = HHT_IMF_PLOT_HILBERT
                                )
        
        plot_imf_energy_diagnostics(
                                    "DET", Ld, yd, imfs_d, 
                                    out_dir         = p_det.parent,
                                    smooth          = HHT_IMF_PLOT_SMOOTH,
                                    topk            = HHT_IMF_PLOT_TOPK,
                                    dpi             = HHT_IMF_PLOT_DPI,
                                    use_hilbert     = HHT_IMF_PLOT_HILBERT
                                )

    
    # --- Calendar keys aligned to the L-grids (match y length) ---
    CAL_KEY_RAW_L, CAL_KEY_DET_L = None, None
    mode                         = CALENDAR_MODE
    
    if os.environ.get("HHT_DECOMP_LOG", "0") == "1":
        topk                = int(os.environ.get("HHT_DECOMP_TOPK", "10"))
        
        decomposition_report("RAW", 
                             Lr, 
                             yr, 
                             imfs_r, 
                             target_ang = lam_r_on_Lr, 
                             topk       = topk)
        
        decomposition_report("DET", 
                             Ld, 
                             yd, 
                             imfs_d, 
                             target_ang = lam_d_on_Ld, 
                             topk       = topk)
   
    if mode == "solar":
        CAL_KEY_RAW_L       = lam_r_on_Lr
        CAL_KEY_DET_L       = lam_d_on_Ld
    
    elif mode == "lunar-phase":
        if "phi_moon" in df_raw.columns and "phi_moon" in df_det.columns:
            phi_raw         = pd.to_numeric(df_raw["phi_moon"], errors = "coerce").to_numpy(float)
            phi_det         = pd.to_numeric(df_det["phi_moon"], errors = "coerce").to_numpy(float)
            CAL_KEY_RAW_L   = interpolate_angle_over_L(L_raw_pts[o_raw], phi_raw[o_raw], Lr)
            CAL_KEY_DET_L   = interpolate_angle_over_L(L_det_pts[o_det], phi_det[o_det], Ld)
    
    elif mode == "lunar-lambda":
        if "lam_moon" in df_raw.columns and "lam_moon" in df_det.columns:
            lm_raw          = pd.to_numeric(df_raw["lam_moon"], errors = "coerce").to_numpy(float)
            lm_det          = pd.to_numeric(df_det["lam_moon"], errors = "coerce").to_numpy(float)
            CAL_KEY_RAW_L   = interpolate_angle_over_L(L_raw_pts[o_raw], lm_raw[o_raw], Lr)
            CAL_KEY_DET_L   = interpolate_angle_over_L(L_det_pts[o_det], lm_det[o_det], Ld)
    
    print(f"[calendar] on-grid keys -> raw={CAL_KEY_RAW_L is not None}, det={CAL_KEY_DET_L is not None}")

    
    if isinstance(lam_block_shufs, np.ndarray) and lam_block_shufs.ndim == 2:
        obs_plv             = plv_d_lambda
        ge                  = 0
        
        for j in range(lam_block_shufs.shape[0]):
            if plv(phi_d, lam_block_shufs[j]) >= obs_plv:
                ge += 1
        
        p_plv_blocklam_det  = (ge + 1.0) / (lam_block_shufs.shape[0] + 1.0)
   
    else:
        p_plv_blocklam_det  = np.nan
    
    # Phase stability (raw vs detrended) on Ld
    phi_r_u                 = np.exp(1j * phi_r)
    re_spl                  = CubicSpline(Lr, np.real(phi_r_u), bc_type = "natural")
    im_spl                  = CubicSpline(Lr, np.imag(phi_r_u), bc_type = "natural")
    phi_r_on_Ld             = np.angle(re_spl(Ld) + 1j * im_spl(Ld))
    plv_phase_stab          = plv(phi_r_on_Ld, phi_d)
    
    print("\n[PLV locking]")
    wlab                    = len("PLV(IMF2 vs λ) det")
    
    print(f"{'PLV(IMF2 vs λ) raw                ':<{wlab}} = {plv_r_lambda:.4f}")
    print(f"{'PLV(IMF2 vs λ) det                ':<{wlab}} = {plv_d_lambda:.4f}")
    print(f"{'PLV(IMF2 phase stability, raw↔det)':<{wlab}} = {plv_phase_stab:.4f}")
    
    _print2("PLV perm p (raw)                   ", f"{p_plv_raw:.4g}         ", "PLV perm p (det)", f"{p_plv_det:.4g}")
    print(_aline([("PLV circ-shift p (det)      ", f"{p_plv_shift_det:.4g}   ")],    kpad = wlab))
    print(_aline([("PLV block-λ    p (det)      ", f"{p_plv_blocklam_det:.4g}")], kpad = wlab))


    # ----------- Spin test: random sky rotations (detrended) ------------- #
    print("\n[Spin test on sky rotations] (detrended)")
    sp               = spin_test_plv(df_det, 
                                     phi_d, 
                                     Ld, 
                                     n               = N_SPIN, 
                                     seed            = seed_sur, 
                                     preserve_L_rank = False)

    
    spin_p           = ((sp >= plv_d_lambda).sum() + 1.0) / (len(sp) + 1.0)
    spin_med         = float(np.median(sp))
    spin_q95         = float(np.quantile(sp, 0.95))
    
    print(_aline([
                  ("PLV_det",     f"{plv_d_lambda:.6f}"),
                  ("spin p",      f"{spin_p:.6g}"),
                  ("spin median", f"{spin_med:.6f}"),
                  ("spin q95",    f"{spin_q95:.6f}"),
                 ], kpad=12))

    # ----------- Rayleigh test for phase–angle coupling ------------------ #
    delta_r          = np.mod(phi_r - lam_r_on_Lr + np.pi, 2*np.pi) - np.pi
    delta_d          = np.mod(phi_d - lam_d_on_Ld + np.pi, 2*np.pi) - np.pi
    
    p_ray_r          = float(rayleigh_p(delta_r))
    p_ray_d          = float(rayleigh_p(delta_d))
    p_ray_d_shift    = rayleigh_shift_p(phi_d, 
                                        lam_d_on_Ld, 
                                        n    = N_SHIFT, 
                                        seed = seed_sur)
    
    print(_aline([("det p (analytic)", f"{p_ray_d:.4g}"),
                  ("det p (shift)",    f"{p_ray_d_shift:.4g}")], kpad = 16))

    
    print("\n[Rayleigh test]")
    
    _print2("raw p", f"{p_ray_r:.4g}", "det p", f"{p_ray_d:.4g}")
    
    # ----------- Extended regression (RAW) ------------------------------- #
    print("\n[Extended regression: resid ~ const + sinλ + cosλ + sinβ + z + sinRA + cosRA]")
    
    reg              = extended_regression(df_raw.rename(columns={"resid": "resid"}))
    se               = reg.get("se_phi_deg")
    se_txt           = f"{se:.1f}" if (se is not None and np.isfinite(se)) else "NA"
    
    print(_aline([
                  ("dipole_dir_raw (λ°)", f"{reg['phi_ecl_deg']:.1f}"),
                  ("±",                   se_txt),
                  ("amp",                 f"{reg['amp_ecl']:.4f}"),
                 ], kpad = 22))
    
    print(_aline([
                  ("amp_ecl", f"{reg['amp_ecl']:.6f}"),
                  ("p(sinλ)", f"{reg.get('p_sin_lam', reg.get('p_sinλ', np.nan)):.3g}"),
                  ("p(cosλ)", f"{reg.get('p_cos_lam', reg.get('p_cosλ', np.nan)):.3g}"),
                  ("Wald F",  f"{reg['wald_F']:.3g}"),
                  ("Wald p",  f"{reg['wald_p']:.3g}"),
                 ], kpad = 10))

    # ----------- E(imf2+imf3) surrogates -------------------------------- #
    print("\n[E(imf2+imf3) with IAAFT surrogates]")
    
    r_iaaft = e_band_surrogate_test(Lr, yr, 
                                    band      = "imf23", 
                                    mode      = "iaaft",
                                    ns        = N_SUR_E23, 
                                    trials    = CE_TRIALS, 
                                    noise     = noise, 
                                    seed      = seed_sur)

    
    d_iaaft = e_band_surrogate_test(Ld, yd, 
                                    band      = "imf23", 
                                    mode      = "iaaft", 
                                    ns        = N_SUR_E23,
                                    trials    = CE_TRIALS, 
                                    noise     = noise, 
                                    seed      = seed_sur)
    
    _print2("raw.E_obs ", f"{r_iaaft['E_obs']:.6f}", "raw.p", f"{r_iaaft['p']:.4g}")
    print(_aline([("raw.S50", f"{r_iaaft['s_q50']:.6f}")], kpad = 10))
   
    _print2("det.E_obs ", f"{d_iaaft['E_obs']:.6f}", "det.p", f"{d_iaaft['p']:.4g}")
    print(_aline([("det.S50", f"{d_iaaft['s_q50']:.6f}")], kpad = 10))

    print(f"\n[E(imf2+imf3) with BLOCK surrogates]  block_len = {BLOCK_LEN}")
    
    # BLOCK
    r_block = e_band_surrogate_test(Lr, yr, 
                                    band      = "imf23", 
                                    mode      = "block",
                                    ns        = N_SUR_E23, 
                                    trials    = CE_TRIALS, 
                                    noise     = noise, 
                                    seed      = seed_sur,
                                    block_len = BLOCK_LEN, 
                                    cal_key   = CAL_KEY_RAW_L, 
                                    tag       = "raw")
    
    d_block = e_band_surrogate_test(Ld, yd, 
                                    band      = "imf23", 
                                    mode      = "block",
                                    ns        = N_SUR_E23, 
                                    trials    = CE_TRIALS, 
                                    noise     = noise, 
                                    seed      = seed_sur,
                                    block_len = BLOCK_LEN, 
                                    cal_key   = CAL_KEY_DET_L, 
                                    tag       = "det")
    
    
    _print2("raw.E_obs ", f"{r_block['E_obs']:.6f}", "raw.p", f"{r_block['p']:.4g}")
    print(_aline([("raw.S50", f"{r_block['s_q50']:.6f}")], kpad = 10))
    
    _print2("det.E_obs ", f"{d_block['E_obs']:.6f}", "det.p", f"{d_block['p']:.4g}")
    print(_aline([("det.S50", f"{d_block['s_q50']:.6f}")], kpad = 10))

    print(f"\n[E(imf2+imf3) with AR1+block surrogates]  block_len = {BLOCK_LEN}")
    
    # AR1+BLOCK
    r_ar1 = e_band_surrogate_test(Lr, yr, 
                                  band        = "imf23", 
                                  mode        = "ar1",
                                  ns          = N_SUR_E23, 
                                  trials      = CE_TRIALS, 
                                  noise       = noise, 
                                  seed        = seed_sur,
                                  block_len   = BLOCK_LEN, 
                                  cal_key     = CAL_KEY_RAW_L, 
                                  tag         = "raw")
    
    d_ar1 = e_band_surrogate_test(Ld, yd, 
                                  band        = "imf23", 
                                  mode        = "ar1",
                                  ns          = N_SUR_E23, 
                                  trials      = CE_TRIALS, 
                                  noise       = noise, 
                                  seed        = seed_sur,
                                  block_len   = BLOCK_LEN, 
                                  cal_key     = CAL_KEY_DET_L, 
                                  tag         = "det")

    
    _print2("raw.E_obs ", f"{r_ar1['E_obs']:.6f}", "raw.p", f"{r_ar1['p']:.4g}")
    print(_aline([("raw.S50", f"{r_ar1['s_q50']:.6f}")], kpad = 10))
    
    _print2("det.E_obs ", f"{d_ar1['E_obs']:.6f}", "det.p", f"{d_ar1['p']:.4g}")
    print(_aline([("det.S50", f"{d_ar1['s_q50']:.6f}")], kpad = 10))

    # -------------- Extended regression (DETRENDED) ------------------------- #
    print("\n[Extended regression on DETRENDED: resid_detrended ~ const + sinλ + cosλ + sinβ + z + sinRA + cosRA]")
    
    cols = ["RA", "DEC", "z", "w_mu", y_det_col]
    for extra in ["sin_phi_moon", "cos_phi_moon", "sin_lam_moon", "cos_lam_moon"]:
        if extra in df_det.columns:
            cols.append(extra)
    
    df_det_reg = pd.DataFrame({
                               "resid": pd.to_numeric(df_det[y_det_col], errors = "coerce"),
                                **{c: df_det[c] for c in cols if c != y_det_col}
                             })
    
    reg_det = extended_regression(df_det_reg)
    
    print(_aline([
                  ("dipole_dir_det (λ°)", f"{reg_det['phi_ecl_deg']:.1f}"),
                  ("±",                   f"{reg_det['se_phi_deg']:.1f}"),
                  ("amp",                 f"{reg_det['amp_ecl']:.4f}"),
                 ], kpad = 22))
    
    _print2("amp_ecl", f"{reg_det['amp_ecl']:.6f}", "Wald p", f"{reg_det['wald_p']:.3g}")

    # ------------  Direction in Galactic coordinates  ----------------------- #
    try:
        l_deg, b_deg    = ecliptic_to_galactic(reg_det["phi_ecl_deg"], 0.0)    # β≈0 for λ-only dipole
        print(_aline([("Direction (Galactic)", f"l={l_deg:.1f}°, b={b_deg:.1f}°")], kpad = 22))
    
    except Exception as e:
        l_deg, b_deg    = (np.nan, np.nan)
        print(_aline([("Direction (Galactic)", f"skipped: {e}")], kpad = 22))

    # -----------  Optional sliding-L diagnosis  ---------------------------- #
    runL_center         = runL_amp = runL_p = np.nan
    
    try:
        sweep           = running_ecl_amp(df_det_reg)
        
        if isinstance(sweep, pd.DataFrame) and len(sweep):
            i           = int(np.nanargmax(sweep["amp_ecl"]))
            runL_center = float(sweep.iloc[i]["L_center"])
            runL_amp    = float(sweep.iloc[i]["amp_ecl"])
            runL_p      = float(sweep.iloc[i]["wald_p"])
            
            
            m           = len(sweep)
            pmin        = float(np.nanmin(sweep["wald_p"]))
            p_bonf      = min(1.0, pmin * m)
            qvals       = fdr_bh(sweep["wald_p"].to_numpy(float))
            qmin        = float(np.nanmin(qvals))
            
            print(_aline([
                          ("Running-L windows", f"{m}"),
                          ("best p (nominal)",  f"{pmin:.3g}"),
                          ("best p (Bonf)",     f"{p_bonf:.3g}"),
                          ("best q (BH-FDR)",   f"{qmin:.3g}"),
                         ], kpad=22))

    except Exception as e:
        print(_aline([("Running-L", f"skipped: {e}")], kpad = 22))
        
    # --- Running-λ (ecliptic sector scan) global spin p (DET) ----------------
    # Uses object-level residuals and ecliptic λ; scans λ centers and
    # corrects for look-elsewhere by spinning λ uniformly on [0,2π).
    
    try:
        resid_det_arr  = pd.to_numeric(df_det_reg["resid"],           errors = "coerce").to_numpy(float)
        lam_det_obj, _ = ecliptic_angles_from_radec(df_det_reg["RA"], df_det_reg["DEC"])
        w_mu_arr       = pd.to_numeric(df_det_reg.get("w_mu", 1.0),   errors = "coerce").fillna(1.0).to_numpy(float)
    
        rng_runl       = np.random.default_rng(SPIN_SEED)

        obs_amp_det, obs_L_det, p_spin_det, q50_det, q95_det = runningL_spin_global_p(
                                                                                      resid_det_arr,
                                                                                      lam_det_obj,
                                                                                      base_weights = w_mu_arr,
                                                                                      N            = N_SPIN_RUNL,
                                                                                      rng          = rng_runl,     # preferred: explicit RNG
                                                                                      centers      = None,
                                                                                      half_width   = np.deg2rad(20),
                                                                                     )


    
        print(f"[Running-λ] obs_max_amp = {obs_amp_det:.6f}  obs_L_deg = {(np.degrees(obs_L_det) % 360):.2f}")
        print(f"[Running-λ spin null] global p = {p_spin_det:.6f}  spin median = {q50_det:.3f}  spin q95 = {q95_det:.3f}")
    
    except Exception as e:
        obs_amp_det = obs_L_det = p_spin_det = q50_det = q95_det = np.nan
        print(f"[Running-λ spin null] skipped: {e}")


    # ----------- Bootstrap dipole direction CI ------------------------------ 
    # --- Dipole direction bootstrap (detrended) ----------------------------- 
    try:
        n_boot                   = _get_int("N_BOOT", 1000)  # configurable (driver can set N_BOOT)
        mu_deg, (phi_lo, phi_hi) = bootstrap_phi(df_det_reg, 
                                                 n    = n_boot, 
                                                 seed = SEED)
    
        print(_aline([
                      ("Dipole dir bootstrap (det) φ", f"{mu_deg:.1f}°"),
                      ("68% CI",                      f"[{phi_lo:.1f}, {phi_hi:.1f}]"),
                      ("n_boot",                      f"{n_boot}"),
                     ], kpad = 28))
    
    except Exception as e:
        mu_deg, phi_lo, phi_hi   = (np.nan, np.nan, np.nan)
        _ping(f"[warn] bootstrap_phi failed: {type(e).__name__}: {e}")

        print(_aline([("Dipole dir bootstrap (det)", f"skipped: {e}")], kpad = 28))

    # ----------- Summary CSV (append) ------------------------------------ #
    out_dir      = p_det.parent

    row = {
            # E2
            "E2_raw"                        : s_raw["E_obs"], 
            "E2_det"                        : s_det["E_obs"],
            "p_raw"                         : s_raw["p"],      
            "p_det"                         : s_det["p"],
    
            # PLV & locks
            "PLV_raw_lambda"                : plv_r_lambda, 
            "PLV_det_lambda"                : plv_d_lambda, 
            "PLV_phase_stab"                : plv_phase_stab,
            "PLV_p_raw"                     : p_plv_raw, 
            "PLV_p_det"                     : p_plv_det,
            "PLV_p_det_spin"                : float(((sp >= plv_d_lambda).sum() + 1.0) / (len(sp) + 1.0)),
            "PLV_p_det_shift"               : p_plv_shift_det,
            "PLV_p_det_block_lambda"        : p_plv_blocklam_det,
            "PLV_spin_med_det"              : spin_med,
            "PLV_spin_q95_det"              : spin_q95,
    
            # Rayleigh
            "Rayleigh_p_raw"                : p_ray_r, 
            "Rayleigh_p_det"                : p_ray_d,
            "Rayleigh_p_det_shift"          : p_ray_d_shift,
            
            # Regression RAW
            "amp_ecl"                       : reg["amp_ecl"], 
            "wald_F"                        : reg["wald_F"], 
            "wald_p"                        : reg["wald_p"],
            "phi_ecl_raw_deg"               : reg["phi_ecl_deg"], 
            "se_phi_raw_deg"                : reg["se_phi_deg"],
    
            # IAAFT/BLOCK/AR1
            "E23_raw_iaaft"                 : r_iaaft["E_obs"], 
            "p23_raw_iaaft"                 : r_iaaft["p"],
            "E23_det_iaaft"                 : d_iaaft["E_obs"], 
            "p23_det_iaaft"                 : d_iaaft["p"],
            "E23_raw_block"                 : r_block["E_obs"], 
            "p23_raw_block"                 : r_block["p"],
            "E23_det_block"                 : d_block["E_obs"], 
            "p23_det_block"                 : d_block["p"],
            "E23_raw_ar1"                   : r_ar1["E_obs"],     
            "p23_raw_ar1"                   : r_ar1["p"],
            "E23_det_ar1"                   : d_ar1["E_obs"],     
            "p23_det_ar1"                   : d_ar1["p"],
    
            # Regression DET
            "amp_ecl_det"                   : reg_det["amp_ecl"],
            "wald_F_det"                    : reg_det.get("wald_F", np.nan),
            "wald_p_det"                    : reg_det["wald_p"],
            "phi_ecl_det_deg"               : reg_det["phi_ecl_deg"],
            "se_phi_det_deg"                : reg_det["se_phi_deg"],
    
            # Coords + diagnostics
            "phi_det_gal_l_deg"             : l_deg, 
            "phi_det_gal_b_deg"             : b_deg,
            "runL_best_center"              : runL_center, 
            "runL_best_amp"                 : runL_amp, 
            "runL_best_p"                   : runL_p,
            "phi_det_boot_deg"              : mu_deg, 
            "phi_det_ci68_lo"               : phi_lo, 
            "phi_det_ci68_hi"               : phi_hi,
    
            # Null sizes / controls / provenance
            "n_perm"                        : N_PERM, 
            "n_spin"                        : N_SPIN, 
            "n_shift"                       : N_SHIFT, 
            "n_iaaft"                       : N_SUR_E23,
            "BLOCK_LEN"                     : BLOCK_LEN, 
            "seed"                          : seed_sur, 
            "trials"                        : CE_TRIALS, 
            "noise"                         : noise,
            "data_raw"                      : str(p_raw), 
            "data_det"                      : str(p_det),
            "provenance_hash"               : _sha256_of_paths([p_raw, p_det]),
            
            "calendar"                      : CALENDAR_MODE,
            
            # Running-λ (DET) spin-global
            "runningL_obs_amp_det"          : float(obs_amp_det),
            "runningL_obs_L_deg_det"        : float(np.degrees(obs_L_det) % 360) if np.isfinite(obs_L_det) else np.nan,
            "runningL_spin_global_p_det"    : float(p_spin_det),
            "runningL_spin_q50_det"         : float(q50_det),
            "runningL_spin_q95_det"         : float(q95_det),
          }

    write_summary_V5(row, out_dir   = out_dir)      # prints a [write] ... (schema=4) line
    

    # ----------- E23 Sanity check ---------------------------------------- #
    e23_frac_raw                    = e23_fraction(yr, imfs_r)
    e23_frac_det                    = e23_fraction(yd, imfs_d)
    
    print(_aline([
                  ("E23 fraction (raw)", f"{e23_frac_raw:.4f}"),
                  ("E23 fraction (det)", f"{e23_frac_det:.4f}"),
                 ], kpad = 20))
    
    print(f"null depths: N_PERM = {N_PERM} ({N_PERM_SRC})  N_SHIFT = {N_SHIFT} ({N_SHIFT_SRC})  "
          f"N_SPIN = {N_SPIN} ({N_SPIN_SRC})  N_SUR_E23 = {N_SUR_E23} ({N_SUR_E23_SRC})  N_SUR_E2 = {N_SUR_E2} ({N_SUR_E2_SRC})")

    print(f"min p: perm = {1/(N_PERM+1):.6g} shift = {1/(N_SHIFT+1):.6g} spin = {1/(N_SPIN+1):.6g} "
          f"e2 = {1/(N_SUR_E2+1):.6g} e23 = {1/(N_SUR_E23+1):.6g}")


    print(f"n_sur_e2   = {N_SUR_E2} ({N_SUR_E2_SRC})  n_sur_e23 = {N_SUR_E23} ({N_SUR_E23_SRC})", flush = True)
    print(f"n_perm     = {N_PERM} ({N_PERM_SRC})  n_shift   = {N_SHIFT} ({N_SHIFT_SRC})  n_spin = {N_SPIN} ({N_SPIN_SRC})", flush = True)
    print(f"CE_TRIALS  = {CE_TRIALS} ({CE_TRIALS_SRC})", flush = True)

    
