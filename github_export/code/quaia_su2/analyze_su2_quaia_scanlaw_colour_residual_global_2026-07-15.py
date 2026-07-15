from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy.io import fits


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
CODE_DIR = ROOT / "github_export" / "code" / "quaia_su2"
OUTPUTS = Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs")
DEFAULT_QUAIA = ROOT / "quaia_G20.0.fits"
DEFAULT_SCAN = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_scan" / "su2_quaia_scan.csv"
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_scanlaw_colour_residual_global_20260715"
EXT_SCRIPT = "analyze_su2_quaia_external_dust_gate_2026-07-15.py"

CMB_L_DEG = 264.021
CMB_B_DEG = 48.253


@dataclass(frozen=True)
class Window:
    window_index: int
    tag: str
    bcut_deg: float
    z_lo: float
    z_hi: float
    z_center: float
    scan_N: int
    raw_amp: float
    raw_snr: float
    raw_priority: float
    raw_l_deg: float
    raw_b_deg: float
    raw_sep_cmb_deg: float
    max_abs_su2_delta_mu: float


def load_module(name: str, candidates: list[Path]) -> Any:
    for path in candidates:
        if path.exists():
            spec = importlib.util.spec_from_file_location(name, path)
            if spec is None or spec.loader is None:
                continue
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            return module
    raise FileNotFoundError(f"Could not find module source for {name}")


EXT = load_module("su2_external_dust_gate", [OUTPUTS / EXT_SCRIPT, CODE_DIR / EXT_SCRIPT])


def parse_float(value: Any, default: float = float("nan")) -> float:
    if value in ("", None):
        return default
    try:
        return float(value)
    except Exception:
        return default


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.6g}" if math.isfinite(value) else ""
    return str(value)


def markdown_table(rows: list[dict[str, Any]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines)


def p_ge(values: np.ndarray, observed: float) -> float:
    values = values[np.isfinite(values)]
    if len(values) == 0 or not math.isfinite(observed):
        return float("nan")
    return float((1 + np.sum(values >= observed)) / (1 + len(values)))


def qtile(values: np.ndarray, prob: float) -> float:
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return float("nan")
    return float(np.quantile(values, prob))


def finite_mask(*arrays: np.ndarray) -> np.ndarray:
    mask = np.ones(len(arrays[0]), dtype=bool)
    for arr in arrays:
        mask &= np.isfinite(arr)
    return mask


def standardize_global(values: np.ndarray, mask: np.ndarray) -> tuple[np.ndarray, dict[str, float]]:
    center = float(np.nanmedian(values[mask]))
    scale = float(np.nanstd(values[mask]))
    if not math.isfinite(scale) or scale <= 0.0:
        scale = 1.0
    return (values - center) / scale, {"center": center, "scale": scale}


def standardize_window(values: np.ndarray) -> tuple[np.ndarray, float, float]:
    center = float(np.nanmedian(values))
    scale = float(np.nanstd(values))
    if not math.isfinite(scale) or scale <= 0.0:
        scale = 1.0
    return (values - center) / scale, center, scale


def gal_unit_from_lb(l_deg: float, b_deg: float) -> np.ndarray:
    l = math.radians(l_deg)
    b = math.radians(b_deg)
    cb = math.cos(b)
    return np.array([cb * math.cos(l), cb * math.sin(l), math.sin(b)], dtype=float)


def unit_to_lb(vec: np.ndarray) -> tuple[float, float]:
    norm = float(np.linalg.norm(vec))
    if norm <= 0.0 or not math.isfinite(norm):
        return float("nan"), float("nan")
    unit = vec / norm
    l_deg = math.degrees(math.atan2(float(unit[1]), float(unit[0]))) % 360.0
    b_deg = math.degrees(math.asin(float(np.clip(unit[2], -1.0, 1.0))))
    return l_deg, b_deg


def angular_sep_deg(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 0.0:
        return float("nan")
    dot = float(np.clip(np.dot(a, b) / denom, -1.0, 1.0))
    return math.degrees(math.acos(dot))


def read_windows(scan_csv: Path, min_count: int, max_windows: int | None = None) -> list[Window]:
    df = pd.read_csv(scan_csv)
    df = df[df["tag"].astype(str).str.startswith("z")].copy()
    df["N"] = pd.to_numeric(df["N"], errors="coerce")
    df = df[df["N"] >= min_count].sort_values(["bcut_deg", "z_lo", "z_hi", "tag"]).reset_index(drop=True)
    if max_windows is not None:
        df = df.head(max_windows).copy()
    windows: list[Window] = []
    for i, row in df.iterrows():
        max_delta = parse_float(row.get("max_abs_su2_delta_mu"), 0.0)
        raw_snr = parse_float(row.get("formal_snr"), 0.0)
        raw_priority = parse_float(row.get("su2_quaia_priority_score"), raw_snr * max_delta)
        windows.append(
            Window(
                window_index=i + 1,
                tag=str(row["tag"]),
                bcut_deg=parse_float(row.get("bcut_deg")),
                z_lo=parse_float(row.get("z_lo")),
                z_hi=parse_float(row.get("z_hi")),
                z_center=parse_float(row.get("z_center")),
                scan_N=int(parse_float(row.get("N"), 0.0)),
                raw_amp=parse_float(row.get("amp")),
                raw_snr=raw_snr,
                raw_priority=raw_priority,
                raw_l_deg=parse_float(row.get("l_deg")),
                raw_b_deg=parse_float(row.get("b_deg")),
                raw_sep_cmb_deg=parse_float(row.get("sep_cmb_deg")),
                max_abs_su2_delta_mu=max_delta,
            )
        )
    return windows


def load_quaia(args: argparse.Namespace, windows: list[Window]) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    quaia_path = args.quaia
    if not quaia_path.exists() and hasattr(EXT, "DEFAULT_QUAIA"):
        quaia_path = Path(EXT.DEFAULT_QUAIA)
    gaiascanlaw_data_dir = args.gaiascanlaw_data_dir or EXT.resolve_gaiascanlaw_data_dir()
    with fits.open(quaia_path, memmap=True) as hdul:
        data = hdul[1].data
        q = {
            "z": np.asarray(data["redshift_quaia"], dtype=float),
            "ra": np.asarray(data["ra"], dtype=float),
            "dec": np.asarray(data["dec"], dtype=float),
            "l": np.asarray(data["l"], dtype=float),
            "b": np.asarray(data["b"], dtype=float),
            "bp_rp": np.asarray(data["phot_bp_mean_mag"], dtype=float) - np.asarray(data["phot_rp_mean_mag"], dtype=float),
            "w1_w2": np.asarray(data["mag_w1_vg"], dtype=float) - np.asarray(data["mag_w2_vg"], dtype=float),
        }
    q["vec"] = EXT.gal_unit_array(q["l"], q["b"])
    q.update(EXT.load_gaia_scanlaw_templates(q["ra"], q["dec"], gaiascanlaw_data_dir))

    min_bcut = min(window.bcut_deg for window in windows)
    min_z = min(window.z_lo for window in windows)
    max_z = max(window.z_hi for window in windows)
    base_mask = (
        (q["z"] >= min_z)
        & (q["z"] < max_z)
        & (np.abs(q["b"]) >= min_bcut)
        & finite_mask(q["bp_rp"], q["w1_w2"])
    )
    q["bp_rp_z"], bp_meta = standardize_global(q["bp_rp"], base_mask)
    q["w1_w2_z"], w_meta = standardize_global(q["w1_w2"], base_mask)
    cross_raw = q["bp_rp_z"] * q["w1_w2_z"]
    q["bp_rp_x_w1_w2_z"], cross_meta = standardize_global(cross_raw, base_mask)

    control_names = control_columns()
    controls_finite = np.column_stack([np.isfinite(q[name]) for name in control_names])
    meta = {
        "quaia": str(quaia_path),
        "gaiascanlaw_data_dir": str(gaiascanlaw_data_dir) if gaiascanlaw_data_dir else None,
        "n_rows": int(len(q["z"])),
        "base_standardization_mask_N": int(np.sum(base_mask)),
        "scanlaw_templates_present": bool(np.all(controls_finite[base_mask, :4])),
        "global_colour_standardization": {
            "bp_rp_z": bp_meta,
            "w1_w2_z": w_meta,
            "bp_rp_x_w1_w2_z": cross_meta,
        },
    }
    return q, meta


def control_columns() -> list[str]:
    return [
        "gaia_scan_count_log1p_dr3",
        "gaia_scan_angle_cos2_mean_dr3",
        "gaia_scan_angle_sin2_mean_dr3",
        "gaia_scan_angle_resultant_dr3",
        "bp_rp_z",
        "w1_w2_z",
        "bp_rp_x_w1_w2_z",
    ]


def fit_controls(z: np.ndarray, controls: np.ndarray) -> tuple[np.ndarray, dict[str, Any]]:
    columns = [np.ones(len(z), dtype=float)]
    control_centers: list[float] = []
    control_scales: list[float] = []
    for j in range(controls.shape[1]):
        col, center, scale = standardize_window(controls[:, j])
        columns.append(col)
        control_centers.append(center)
        control_scales.append(scale)
    x = np.column_stack(columns)
    beta, *_ = np.linalg.lstsq(x, z, rcond=None)
    fitted = x @ beta
    resid = z - fitted
    rss = float(np.dot(resid, resid))
    tss = float(np.dot(z - np.mean(z), z - np.mean(z)))
    meta = {
        "control_r2": 1.0 - rss / tss if tss > 0.0 else float("nan"),
        "control_rms_resid": math.sqrt(rss / max(len(z) - x.shape[1], 1)),
        "control_k": int(x.shape[1]),
        "control_beta_intercept": float(beta[0]),
        "control_betas": ";".join(f"{value:.10g}" for value in beta[1:]),
        "control_centers": ";".join(f"{value:.10g}" for value in control_centers),
        "control_scales": ";".join(f"{value:.10g}" for value in control_scales),
    }
    return resid, meta


def cached_dipole_inputs(vec: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    x = np.column_stack([np.ones(len(vec), dtype=float), vec[:, 0], vec[:, 1], vec[:, 2]])
    xtx = x.T @ x
    pinv_xtx = np.linalg.pinv(xtx)
    info_d = np.linalg.pinv(pinv_xtx[1:4, 1:4])
    dof = max(len(vec) - 4, 1)
    return x, pinv_xtx, info_d, dof


def fit_dipole_cached(y: np.ndarray, x: np.ndarray, pinv_xtx: np.ndarray, info_d: np.ndarray, dof: int, yty: float | None = None) -> dict[str, float]:
    if yty is None:
        yty = float(np.dot(y, y))
    xty = x.T @ y
    beta = pinv_xtx @ xty
    sse = max(yty - float(beta @ xty), 0.0)
    sigma2 = sse / dof
    dvec = beta[1:4]
    amp = float(np.linalg.norm(dvec))
    l_deg, b_deg = unit_to_lb(dvec)
    cmb = gal_unit_from_lb(CMB_L_DEG, CMB_B_DEG)
    sep = angular_sep_deg(dvec, cmb)
    amp_par = float(np.dot(dvec, cmb))
    frac_par = amp_par / amp if amp > 0.0 else float("nan")
    try:
        snr = float(math.sqrt(max(float(dvec @ info_d @ dvec) / max(sigma2, 1e-300), 0.0)))
    except Exception:
        snr = float("nan")
    return {
        "amp": amp,
        "l_deg": l_deg,
        "b_deg": b_deg,
        "sep_cmb_deg": sep,
        "amp_par_cmb": amp_par,
        "frac_par_cmb": frac_par,
        "formal_snr": snr,
        "rms_resid": float(math.sqrt(max(sigma2, 0.0))),
    }


def strata_from_score(score: np.ndarray, n_strata: int) -> np.ndarray:
    finite = score[np.isfinite(score)]
    if len(finite) == 0:
        return np.zeros(len(score), dtype=np.int16)
    edges = np.unique(np.quantile(finite, np.linspace(0.0, 1.0, n_strata + 1)))
    if len(edges) <= 2:
        return np.zeros(len(score), dtype=np.int16)
    return np.searchsorted(edges[1:-1], score, side="right").astype(np.int16)


def prepare_windows(args: argparse.Namespace, q: dict[str, np.ndarray], windows: list[Window]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    names = control_columns()
    prepared: list[dict[str, Any]] = []
    observed_rows: list[dict[str, Any]] = []
    for window in windows:
        idx0 = np.flatnonzero(
            (q["z"] >= window.z_lo)
            & (q["z"] < window.z_hi)
            & (np.abs(q["b"]) >= window.bcut_deg)
            & finite_mask(q["z"], q["vec"][:, 0], q["vec"][:, 1], q["vec"][:, 2])
        )
        if len(idx0) < args.min_count:
            continue
        controls0 = np.column_stack([q[name][idx0] for name in names])
        finite = np.all(np.isfinite(controls0), axis=1)
        idx = idx0[finite]
        if len(idx) < args.min_count:
            continue
        z = q["z"][idx]
        vec = q["vec"][idx]
        controls = np.column_stack([q[name][idx] for name in names])
        residual, control_meta = fit_controls(z, controls)
        x, pinv_xtx, info_d, dof = cached_dipole_inputs(vec)
        yty = float(np.dot(residual, residual))
        fit = fit_dipole_cached(residual, x, pinv_xtx, info_d, dof, yty=yty)
        priority = fit["formal_snr"] * window.max_abs_su2_delta_mu
        strata_score = z - residual
        strata = strata_from_score(strata_score, args.n_strata)
        stratum_sizes = np.bincount(strata, minlength=int(np.max(strata)) + 1)
        observed_rows.append(
            {
                "window_index": window.window_index,
                "tag": window.tag,
                "bcut_deg": window.bcut_deg,
                "z_lo": window.z_lo,
                "z_hi": window.z_hi,
                "z_center": window.z_center,
                "scan_N": window.scan_N,
                "N": int(len(idx)),
                "n_controls": len(names),
                "control_names": ";".join(names),
                "n_strata": int(len(stratum_sizes)),
                "min_stratum_N": int(np.min(stratum_sizes)) if len(stratum_sizes) else 0,
                "raw_amp": window.raw_amp,
                "raw_formal_snr": window.raw_snr,
                "raw_priority": window.raw_priority,
                "raw_l_deg": window.raw_l_deg,
                "raw_b_deg": window.raw_b_deg,
                "raw_sep_cmb_deg": window.raw_sep_cmb_deg,
                "residual_amp": fit["amp"],
                "residual_formal_snr": fit["formal_snr"],
                "residual_priority": priority,
                "residual_l_deg": fit["l_deg"],
                "residual_b_deg": fit["b_deg"],
                "residual_sep_cmb_deg": fit["sep_cmb_deg"],
                "residual_amp_par_cmb": fit["amp_par_cmb"],
                "residual_frac_par_cmb": fit["frac_par_cmb"],
                "residual_rms": fit["rms_resid"],
                "amp_ratio_residual_vs_raw": fit["amp"] / window.raw_amp if window.raw_amp > 0.0 else float("nan"),
                "snr_delta_residual_minus_raw": fit["formal_snr"] - window.raw_snr,
                "max_abs_su2_delta_mu": window.max_abs_su2_delta_mu,
                **control_meta,
            }
        )
        prepared.append(
            {
                "window": window,
                "residual": residual.astype(np.float64, copy=False),
                "x": x.astype(np.float64, copy=False),
                "pinv_xtx": pinv_xtx,
                "info_d": info_d,
                "dof": dof,
                "yty": yty,
                "strata": strata,
            }
        )
        if args.progress_every_windows and (len(prepared) == 1 or len(prepared) % args.progress_every_windows == 0):
            print(
                f"prepared {len(prepared)}/{len(windows)} windows: {window.tag} |b|>{window.bcut_deg:g}, "
                f"N={len(idx)}, residual_snr={fit['formal_snr']:.4g}",
                flush=True,
            )
    return prepared, observed_rows


def initial_best() -> dict[str, Any]:
    return {
        "window_index": -1,
        "tag": "",
        "bcut_deg": float("nan"),
        "value": -np.inf,
        "snr": float("nan"),
        "priority": float("nan"),
        "amp": float("nan"),
        "l_deg": float("nan"),
        "b_deg": float("nan"),
        "sep_cmb_deg": float("nan"),
    }


def update_best(best: dict[str, Any], value: float, window: Window, fit: dict[str, float], priority: float) -> dict[str, Any]:
    if math.isfinite(value) and value > float(best["value"]):
        return {
            "window_index": window.window_index,
            "tag": window.tag,
            "bcut_deg": window.bcut_deg,
            "value": value,
            "snr": fit["formal_snr"],
            "priority": priority,
            "amp": fit["amp"],
            "l_deg": fit["l_deg"],
            "b_deg": fit["b_deg"],
            "sep_cmb_deg": fit["sep_cmb_deg"],
        }
    return best


def shuffle_within_strata(values: np.ndarray, strata: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    out = values.copy()
    for sid in np.unique(strata):
        sub = np.flatnonzero(strata == sid)
        if len(sub) > 1:
            out[sub] = values[sub][rng.permutation(len(sub))]
    return out


def run_mocks(args: argparse.Namespace, prepared: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rng = np.random.default_rng(args.seed)
    rows: list[dict[str, Any]] = []
    partial_path = args.out_dir / "su2_quaia_scanlaw_colour_residual_global_mock_maxima.partial.csv"
    for mock_id in range(args.n_mocks):
        best_snr = initial_best()
        best_priority = initial_best()
        best_abs_par = initial_best()
        best_z1p00_snr = initial_best()
        best_z1p00_priority = initial_best()
        for item in prepared:
            window: Window = item["window"]
            y_mock = shuffle_within_strata(item["residual"], item["strata"], rng)
            fit = fit_dipole_cached(y_mock, item["x"], item["pinv_xtx"], item["info_d"], item["dof"], yty=item["yty"])
            priority = fit["formal_snr"] * window.max_abs_su2_delta_mu
            abs_par = abs(fit["amp_par_cmb"])
            best_snr = update_best(best_snr, fit["formal_snr"], window, fit, priority)
            best_priority = update_best(best_priority, priority, window, fit, priority)
            best_abs_par = update_best(best_abs_par, abs_par, window, fit, priority)
            if window.tag == "z1p00_1p50":
                best_z1p00_snr = update_best(best_z1p00_snr, fit["formal_snr"], window, fit, priority)
                best_z1p00_priority = update_best(best_z1p00_priority, priority, window, fit, priority)

        row = {
            "mock_id": mock_id,
            "max_snr": best_snr["value"],
            "max_snr_window_index": best_snr["window_index"],
            "max_snr_tag": best_snr["tag"],
            "max_snr_bcut_deg": best_snr["bcut_deg"],
            "max_snr_amp": best_snr["amp"],
            "max_snr_l_deg": best_snr["l_deg"],
            "max_snr_b_deg": best_snr["b_deg"],
            "max_snr_sep_cmb_deg": best_snr["sep_cmb_deg"],
            "max_priority": best_priority["value"],
            "max_priority_window_index": best_priority["window_index"],
            "max_priority_tag": best_priority["tag"],
            "max_priority_bcut_deg": best_priority["bcut_deg"],
            "max_priority_snr": best_priority["snr"],
            "max_priority_amp": best_priority["amp"],
            "max_priority_l_deg": best_priority["l_deg"],
            "max_priority_b_deg": best_priority["b_deg"],
            "max_priority_sep_cmb_deg": best_priority["sep_cmb_deg"],
            "max_abs_par_cmb": best_abs_par["value"],
            "max_abs_par_cmb_window_index": best_abs_par["window_index"],
            "max_abs_par_cmb_tag": best_abs_par["tag"],
            "max_abs_par_cmb_bcut_deg": best_abs_par["bcut_deg"],
            "max_abs_par_cmb_snr": best_abs_par["snr"],
            "max_z1p00_snr": best_z1p00_snr["value"],
            "max_z1p00_snr_bcut_deg": best_z1p00_snr["bcut_deg"],
            "max_z1p00_priority": best_z1p00_priority["value"],
            "max_z1p00_priority_bcut_deg": best_z1p00_priority["bcut_deg"],
        }
        rows.append(row)
        if args.progress_every and (mock_id == 0 or (mock_id + 1) % args.progress_every == 0):
            print(
                f"mock {mock_id + 1}/{args.n_mocks}: max_snr={row['max_snr']:.4g} "
                f"({row['max_snr_tag']} |b|>{row['max_snr_bcut_deg']:g}), "
                f"max_priority={row['max_priority']:.5g} ({row['max_priority_tag']} |b|>{row['max_priority_bcut_deg']:g})",
                flush=True,
            )
        if args.checkpoint_every and ((mock_id + 1) % args.checkpoint_every == 0 or mock_id == 0):
            write_csv(partial_path, rows)
            if args.copy_dir is not None:
                args.copy_dir.mkdir(parents=True, exist_ok=True)
                shutil.copy2(partial_path, args.copy_dir / partial_path.name)
    return rows


def threshold_rows(observed_rows: list[dict[str, Any]], mock_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    max_snr = np.array([row["max_snr"] for row in mock_rows], dtype=float)
    max_priority = np.array([row["max_priority"] for row in mock_rows], dtype=float)
    max_abs_par = np.array([row["max_abs_par_cmb"] for row in mock_rows], dtype=float)
    max_z1p00_snr = np.array([row["max_z1p00_snr"] for row in mock_rows], dtype=float)
    max_z1p00_priority = np.array([row["max_z1p00_priority"] for row in mock_rows], dtype=float)

    targets: list[tuple[str, dict[str, Any]]] = []
    targets.append(("observed_max_residual_snr", max(observed_rows, key=lambda row: row["residual_formal_snr"])))
    targets.append(("observed_max_residual_priority", max(observed_rows, key=lambda row: row["residual_priority"])))
    z_rows = [row for row in observed_rows if row["tag"] == "z1p00_1p50"]
    if z_rows:
        targets.append(("z1p00_1p50_family_max_residual_snr", max(z_rows, key=lambda row: row["residual_formal_snr"])))
        targets.append(("z1p00_1p50_family_max_residual_priority", max(z_rows, key=lambda row: row["residual_priority"])))
        for row in sorted(z_rows, key=lambda value: value["bcut_deg"]):
            targets.append((f"z1p00_1p50_bcut_{row['bcut_deg']:g}", row))

    out: list[dict[str, Any]] = []
    seen: set[tuple[str, float, str]] = set()
    for role, target in targets:
        key = (target["tag"], float(target["bcut_deg"]), role)
        if key in seen:
            continue
        seen.add(key)
        out.append(
            {
                "role": role,
                "window_index": target["window_index"],
                "tag": target["tag"],
                "bcut_deg": target["bcut_deg"],
                "z_center": target["z_center"],
                "N": target["N"],
                "raw_formal_snr": target["raw_formal_snr"],
                "raw_priority": target["raw_priority"],
                "control_r2": target["control_r2"],
                "residual_amp": target["residual_amp"],
                "residual_formal_snr": target["residual_formal_snr"],
                "residual_priority": target["residual_priority"],
                "residual_l_deg": target["residual_l_deg"],
                "residual_b_deg": target["residual_b_deg"],
                "residual_sep_cmb_deg": target["residual_sep_cmb_deg"],
                "residual_amp_par_cmb": target["residual_amp_par_cmb"],
                "max_abs_su2_delta_mu": target["max_abs_su2_delta_mu"],
                "global_p_any_window_snr_ge_observed": p_ge(max_snr, target["residual_formal_snr"]),
                "global_p_any_window_priority_ge_observed": p_ge(max_priority, target["residual_priority"]),
                "global_p_any_window_abs_par_cmb_ge_observed": p_ge(max_abs_par, abs(target["residual_amp_par_cmb"])),
                "z1p00_family_p_any_bcut_snr_ge_observed": p_ge(max_z1p00_snr, target["residual_formal_snr"]) if target["tag"] == "z1p00_1p50" else float("nan"),
                "z1p00_family_p_any_bcut_priority_ge_observed": p_ge(max_z1p00_priority, target["residual_priority"]) if target["tag"] == "z1p00_1p50" else float("nan"),
            }
        )
    return out


def summarize_gate(thresholds: list[dict[str, Any]]) -> str:
    finite_global = [
        row["global_p_any_window_snr_ge_observed"]
        for row in thresholds
        if math.isfinite(float(row.get("global_p_any_window_snr_ge_observed", float("nan"))))
    ]
    finite_z = [
        row["z1p00_family_p_any_bcut_snr_ge_observed"]
        for row in thresholds
        if row.get("tag") == "z1p00_1p50" and math.isfinite(float(row.get("z1p00_family_p_any_bcut_snr_ge_observed", float("nan"))))
    ]
    min_global = min(finite_global) if finite_global else float("nan")
    min_z = min(finite_z) if finite_z else float("nan")
    if math.isfinite(min_global) and min_global <= 0.01:
        return "A residual angular mode is globally unusual at the 1% scan gate and should be treated as a follow-up candidate only."
    if math.isfinite(min_z) and min_z <= 0.01:
        return "The z1p00_1p50 residual family is unusual within its family gate but does not by itself constitute a discovery claim."
    return "No scan-law/colour residual mode survives the global promotion gate; the scan-law/colour explanation remains favoured."


def write_reports(args: argparse.Namespace, observed_rows: list[dict[str, Any]], mock_rows: list[dict[str, Any]], thresholds: list[dict[str, Any]], config: dict[str, Any]) -> None:
    max_snr = np.array([row["max_snr"] for row in mock_rows], dtype=float)
    max_priority = np.array([row["max_priority"] for row in mock_rows], dtype=float)
    top = sorted(observed_rows, key=lambda row: row["residual_priority"], reverse=True)[:15]
    bottom = summarize_gate(thresholds)
    report = [
        "# SU2 / Quaia Scan-Law + Colour Residual Global Look-Elsewhere Search",
        "",
        "Date: July 15, 2026",
        "",
        "## Purpose",
        "",
        "This promotion gate asks whether any scanned redshift-window / Galactic-latitude-cut angular mode remains globally unusual after regressing Quaia redshift against Gaia scan-law geometry and catalogue-colour controls.",
        "",
        "## Residual Model",
        "",
        "```text",
        "z_i      = alpha + S_i beta_s + C_i beta_c + epsilon_i",
        "r_i      = z_i - fitted(z_i | S_i, C_i)",
        "r_i      = a0 + d . n_i + eta_i",
        "priority = residual_SNR_d * max_abs_delta_mu_SU2(window)",
        "```",
        "",
        "`S_i` contains Gaia scan count log1p, scan-angle cos2 mean, scan-angle sin2 mean, and scan-angle resultant. `C_i` contains Gaia BP-RP z-score, WISE W1-W2 z-score, and their cross term.",
        "",
        "## Null",
        "",
        "For each window/cut, residuals are shuffled within quantiles of the scan-law/colour fitted redshift. Each mock records the maximum statistic over all prepared windows.",
        "",
        "## Bottom Line",
        "",
        bottom,
        "",
        "## Configuration",
        "",
        f"- prepared windows: `{len(observed_rows)}`",
        f"- residual-shuffle mocks: `{len(mock_rows)}`",
        f"- strata per window: `{args.n_strata}`",
        f"- seed: `{args.seed}`",
        "",
        "## Mock Maxima",
        "",
        f"- max residual SNR mean / p95 / p99: `{np.nanmean(max_snr):.6g}` / `{qtile(max_snr, 0.95):.6g}` / `{qtile(max_snr, 0.99):.6g}`",
        f"- max residual priority mean / p95 / p99: `{np.nanmean(max_priority):.6g}` / `{qtile(max_priority, 0.95):.6g}` / `{qtile(max_priority, 0.99):.6g}`",
        "",
        "## Observed Thresholds",
        "",
        markdown_table(
            thresholds,
            [
                "role",
                "tag",
                "bcut_deg",
                "N",
                "residual_formal_snr",
                "residual_priority",
                "global_p_any_window_snr_ge_observed",
                "global_p_any_window_priority_ge_observed",
                "z1p00_family_p_any_bcut_snr_ge_observed",
                "z1p00_family_p_any_bcut_priority_ge_observed",
            ],
        ),
        "",
        "## Top Observed Residual Windows",
        "",
        markdown_table(
            top,
            [
                "tag",
                "bcut_deg",
                "N",
                "control_r2",
                "raw_formal_snr",
                "residual_formal_snr",
                "residual_priority",
                "residual_l_deg",
                "residual_b_deg",
                "residual_sep_cmb_deg",
            ],
        ),
        "",
        "## Claim Boundary",
        "",
        "This is not a discovery run. Ordinary residual global p-values favour the scan-law/colour explanation. A surviving residual mode is a follow-up candidate only.",
        "",
        "## Files",
        "",
        f"- observed windows: `{args.out_dir / 'su2_quaia_scanlaw_colour_residual_global_observed_windows.csv'}`",
        f"- mock maxima: `{args.out_dir / 'su2_quaia_scanlaw_colour_residual_global_mock_maxima.csv'}`",
        f"- thresholds: `{args.out_dir / 'su2_quaia_scanlaw_colour_residual_global_thresholds.csv'}`",
    ]
    (args.out_dir / "su2_quaia_scanlaw_colour_residual_global_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    readout = [
        "# SU2 / Quaia Scan-Law + Colour Residual Global Readout",
        "",
        "Date: July 15, 2026",
        "",
        "## Bottom Line",
        "",
        bottom,
        "",
        "## Thresholds",
        "",
        markdown_table(
            thresholds,
            [
                "role",
                "tag",
                "bcut_deg",
                "residual_formal_snr",
                "residual_priority",
                "global_p_any_window_snr_ge_observed",
                "z1p00_family_p_any_bcut_snr_ge_observed",
            ],
        ),
        "",
    ]
    (args.out_dir / "su2_quaia_scanlaw_colour_residual_global_readout.md").write_text("\n".join(readout) + "\n", encoding="utf-8")
    (args.out_dir / "su2_quaia_scanlaw_colour_residual_global_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")


def copy_user_facing(args: argparse.Namespace) -> None:
    if args.copy_dir is None:
        return
    args.copy_dir.mkdir(parents=True, exist_ok=True)
    for name in [
        "su2_quaia_scanlaw_colour_residual_global_observed_windows.csv",
        "su2_quaia_scanlaw_colour_residual_global_mock_maxima.csv",
        "su2_quaia_scanlaw_colour_residual_global_thresholds.csv",
        "su2_quaia_scanlaw_colour_residual_global_report.md",
        "su2_quaia_scanlaw_colour_residual_global_readout.md",
        "su2_quaia_scanlaw_colour_residual_global_config.json",
        "su2_quaia_scanlaw_colour_residual_global_manifest.json",
    ]:
        src = args.out_dir / name
        if src.exists():
            shutil.copy2(src, args.copy_dir / name)


def write_manifest(args: argparse.Namespace, config: dict[str, Any], observed_rows: list[dict[str, Any]], mock_rows: list[dict[str, Any]], thresholds: list[dict[str, Any]]) -> None:
    manifest = {
        "date": "2026-07-15",
        "analysis": "su2_quaia_scanlaw_colour_residual_global",
        "claim_boundary": "Promotion gate only; surviving residual modes are follow-up candidates, not discoveries.",
        "config": config,
        "counts": {
            "observed_windows": len(observed_rows),
            "mocks": len(mock_rows),
            "threshold_rows": len(thresholds),
        },
        "primary_files": [
            str(args.out_dir / "su2_quaia_scanlaw_colour_residual_global_observed_windows.csv"),
            str(args.out_dir / "su2_quaia_scanlaw_colour_residual_global_mock_maxima.csv"),
            str(args.out_dir / "su2_quaia_scanlaw_colour_residual_global_thresholds.csv"),
            str(args.out_dir / "su2_quaia_scanlaw_colour_residual_global_report.md"),
            str(args.out_dir / "su2_quaia_scanlaw_colour_residual_global_readout.md"),
        ],
    }
    (args.out_dir / "su2_quaia_scanlaw_colour_residual_global_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scan-law/colour residual global look-elsewhere gate for SU2/Quaia.")
    parser.add_argument("--quaia", type=Path, default=DEFAULT_QUAIA)
    parser.add_argument("--scan-csv", type=Path, default=DEFAULT_SCAN)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--copy-dir", type=Path, default=None)
    parser.add_argument("--gaiascanlaw-data-dir", type=Path, default=None)
    parser.add_argument("--n-mocks", type=int, default=1000)
    parser.add_argument("--n-strata", type=int, default=10)
    parser.add_argument("--min-count", type=int, default=50)
    parser.add_argument("--seed", type=int, default=190715)
    parser.add_argument("--progress-every", type=int, default=10)
    parser.add_argument("--progress-every-windows", type=int, default=25)
    parser.add_argument("--checkpoint-every", type=int, default=50)
    parser.add_argument("--max-windows", type=int, default=None, help="Debug/smoke limit on the number of scan windows.")
    parser.add_argument("--prepare-only", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    windows = read_windows(args.scan_csv, args.min_count, args.max_windows)
    if not windows:
        raise RuntimeError("No scan windows passed the input N threshold.")
    q, data_meta = load_quaia(args, windows)
    print(f"Loaded Quaia rows: {len(q['z'])}", flush=True)
    print(f"Candidate scan windows: {len(windows)}", flush=True)
    print(f"Control columns: {', '.join(control_columns())}", flush=True)
    prepared, observed_rows = prepare_windows(args, q, windows)
    print(f"Prepared residual windows: {len(prepared)}", flush=True)
    if not prepared:
        raise RuntimeError("No windows remained after finite scan-law/colour control filtering.")

    config = {
        "date": "2026-07-15",
        "analysis": "scanlaw_colour_residual_global_lookelsewhere",
        "quaia": data_meta["quaia"],
        "scan_csv": str(args.scan_csv),
        "out_dir": str(args.out_dir),
        "n_mocks": args.n_mocks,
        "n_strata": args.n_strata,
        "min_count": args.min_count,
        "seed": args.seed,
        "max_windows": args.max_windows,
        "checkpoint_every": args.checkpoint_every,
        "controls": control_columns(),
        "null": "Residuals from z~scanlaw+colour are shuffled within quantiles of the fitted scanlaw+colour redshift for each window/cut.",
        "data_meta": data_meta,
    }
    write_csv(args.out_dir / "su2_quaia_scanlaw_colour_residual_global_observed_windows.csv", observed_rows)
    (args.out_dir / "su2_quaia_scanlaw_colour_residual_global_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    if args.prepare_only:
        mock_rows: list[dict[str, Any]] = []
        thresholds: list[dict[str, Any]] = []
    else:
        mock_rows = run_mocks(args, prepared)
        thresholds = threshold_rows(observed_rows, mock_rows)

    write_csv(args.out_dir / "su2_quaia_scanlaw_colour_residual_global_observed_windows.csv", observed_rows)
    if mock_rows:
        write_csv(args.out_dir / "su2_quaia_scanlaw_colour_residual_global_mock_maxima.csv", mock_rows)
        write_csv(args.out_dir / "su2_quaia_scanlaw_colour_residual_global_thresholds.csv", thresholds)
        write_reports(args, observed_rows, mock_rows, thresholds, config)
    else:
        (args.out_dir / "su2_quaia_scanlaw_colour_residual_global_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    write_manifest(args, config, observed_rows, mock_rows, thresholds)
    copy_user_facing(args)
    print(f"Saved observed windows: {args.out_dir / 'su2_quaia_scanlaw_colour_residual_global_observed_windows.csv'}", flush=True)
    if mock_rows:
        print(f"Saved mock maxima: {args.out_dir / 'su2_quaia_scanlaw_colour_residual_global_mock_maxima.csv'}", flush=True)
        print(f"Saved thresholds: {args.out_dir / 'su2_quaia_scanlaw_colour_residual_global_thresholds.csv'}", flush=True)
        print(f"Saved report: {args.out_dir / 'su2_quaia_scanlaw_colour_residual_global_report.md'}", flush=True)
        print(summarize_gate(thresholds), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
