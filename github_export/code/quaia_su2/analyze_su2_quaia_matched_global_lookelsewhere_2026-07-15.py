from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy.io import fits
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
DEFAULT_QUAIA = ROOT / "quaia_G20.0.fits"
DEFAULT_SCAN = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_scan" / "su2_quaia_scan.csv"
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_matched_global_lookelsewhere_20260715"


@dataclass(frozen=True)
class Window:
    window_index: int
    tag: str
    bcut_deg: float
    z_lo: float
    z_hi: float
    z_center: float
    max_abs_su2_delta_mu: float


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


def gal_unit_array(l_deg: np.ndarray, b_deg: np.ndarray) -> np.ndarray:
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    cb = np.cos(b)
    return np.column_stack([cb * np.cos(l), cb * np.sin(l), np.sin(b)])


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


def load_quaia(path: Path) -> dict[str, np.ndarray]:
    with fits.open(path, memmap=True) as hdul:
        data = hdul[1].data
        q = {
            "z": np.asarray(data["redshift_quaia"], dtype=float),
            "zerr": np.asarray(data["redshift_quaia_err"], dtype=float),
            "l": np.asarray(data["l"], dtype=float),
            "b": np.asarray(data["b"], dtype=float),
            "g": np.asarray(data["phot_g_mean_mag"], dtype=float),
            "bp_rp": np.asarray(data["phot_bp_mean_mag"], dtype=float) - np.asarray(data["phot_rp_mean_mag"], dtype=float),
            "w1": np.asarray(data["mag_w1_vg"], dtype=float),
            "w2": np.asarray(data["mag_w2_vg"], dtype=float),
            "w1_w2": np.asarray(data["mag_w1_vg"], dtype=float) - np.asarray(data["mag_w2_vg"], dtype=float),
            "pmra_error": np.asarray(data["pmra_error"], dtype=float),
            "pmdec_error": np.asarray(data["pmdec_error"], dtype=float),
        }
    q["vec"] = gal_unit_array(q["l"], q["b"])
    return q


def read_windows(scan_csv: Path, min_count: int) -> list[Window]:
    df = pd.read_csv(scan_csv)
    df = df[df["tag"].astype(str).str.startswith("z")].copy()
    df["N"] = pd.to_numeric(df["N"], errors="coerce")
    df = df[df["N"] >= min_count].sort_values(["bcut_deg", "z_lo", "z_hi", "tag"]).reset_index(drop=True)
    windows: list[Window] = []
    for i, row in df.iterrows():
        windows.append(
            Window(
                window_index=i + 1,
                tag=str(row["tag"]),
                bcut_deg=float(row["bcut_deg"]),
                z_lo=float(row["z_lo"]),
                z_hi=float(row["z_hi"]),
                z_center=float(row["z_center"]),
                max_abs_su2_delta_mu=float(row.get("max_abs_su2_delta_mu", 0.0) or 0.0),
            )
        )
    return windows


def fit_weighted_dipole(z: np.ndarray, vec: np.ndarray, weights: np.ndarray | None = None) -> dict[str, float]:
    n = len(z)
    x = np.column_stack([np.ones(n), vec[:, 0], vec[:, 1], vec[:, 2]])
    if weights is None:
        weights = np.ones(n, dtype=float)
    weights = np.asarray(weights, dtype=float)
    weights = np.where(np.isfinite(weights) & (weights > 0.0), weights, 0.0)
    if int(np.sum(weights > 0.0)) < 5:
        return {"amp": float("nan"), "l_deg": float("nan"), "b_deg": float("nan"), "formal_snr": float("nan"), "rms_resid": float("nan")}
    w = weights / np.mean(weights[weights > 0.0])
    xtwx = x.T @ (x * w[:, None])
    xtwz = x.T @ (w * z)
    pinv = np.linalg.pinv(xtwx)
    beta = pinv @ xtwz
    resid = z - x @ beta
    ess = float((np.sum(w) ** 2) / np.sum(w * w))
    dof = max(ess - x.shape[1], 1.0)
    sigma2 = float(np.sum(w * resid * resid) / dof)
    dvec = beta[1:4]
    amp = float(np.linalg.norm(dvec))
    l_deg, b_deg = unit_to_lb(dvec)
    try:
        cov_d = sigma2 * pinv[1:4, 1:4]
        snr = float(math.sqrt(max(float(dvec @ np.linalg.pinv(cov_d) @ dvec), 0.0)))
    except Exception:
        snr = float("nan")
    return {"amp": amp, "l_deg": l_deg, "b_deg": b_deg, "formal_snr": snr, "rms_resid": float(math.sqrt(max(sigma2, 0.0)))}


def weighted_fit_with_cache(z: np.ndarray, vec: np.ndarray, weights: np.ndarray, pinv_xtwx: np.ndarray) -> dict[str, float]:
    n = len(z)
    x = np.column_stack([np.ones(n), vec[:, 0], vec[:, 1], vec[:, 2]])
    w = weights / np.mean(weights[weights > 0.0])
    xtwz = x.T @ (w * z)
    beta = pinv_xtwx @ xtwz
    resid = z - x @ beta
    ess = float((np.sum(w) ** 2) / np.sum(w * w))
    dof = max(ess - x.shape[1], 1.0)
    sigma2 = float(np.sum(w * resid * resid) / dof)
    dvec = beta[1:4]
    amp = float(np.linalg.norm(dvec))
    l_deg, b_deg = unit_to_lb(dvec)
    try:
        cov_d = sigma2 * pinv_xtwx[1:4, 1:4]
        snr = float(math.sqrt(max(float(dvec @ np.linalg.pinv(cov_d) @ dvec), 0.0)))
    except Exception:
        snr = float("nan")
    return {"amp": amp, "l_deg": l_deg, "b_deg": b_deg, "formal_snr": snr, "rms_resid": float(math.sqrt(max(sigma2, 0.0)))}


def build_features(q: dict[str, np.ndarray], idx: np.ndarray, variables: list[str]) -> np.ndarray:
    return np.column_stack([q[name][idx] for name in variables]).astype(float)


def overlap_weights(features: np.ndarray, hemisphere: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray, float]:
    model = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=600, C=1.0, solver="lbfgs", random_state=seed),
    )
    model.fit(features, hemisphere.astype(int))
    p = np.clip(model.predict_proba(features)[:, 1], 0.02, 0.98)
    weights = np.where(hemisphere, 1.0 - p, p)
    weights /= np.mean(weights)
    ess = float((np.sum(weights) ** 2) / np.sum(weights * weights))
    return weights, p, ess


def propensity_strata(propensity: np.ndarray, n_strata: int) -> np.ndarray:
    finite = propensity[np.isfinite(propensity)]
    edges = np.unique(np.quantile(finite, np.linspace(0.0, 1.0, n_strata + 1)))
    if len(edges) <= 2:
        return np.zeros(len(propensity), dtype=np.int64)
    return np.searchsorted(edges[1:-1], propensity, side="right").astype(np.int64)


def prepare_windows(args: argparse.Namespace, q: dict[str, np.ndarray], windows: list[Window]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    variables = ["zerr", "g", "bp_rp", "w1", "w2", "w1_w2", "pmra_error", "pmdec_error"]
    prepared: list[dict[str, Any]] = []
    observed_rows: list[dict[str, Any]] = []
    for window in windows:
        idx = np.flatnonzero(
            (q["z"] >= window.z_lo)
            & (q["z"] < window.z_hi)
            & (np.abs(q["b"]) >= window.bcut_deg)
            & np.isfinite(q["z"])
            & np.all(np.column_stack([np.isfinite(q[name]) for name in variables]), axis=1)
        )
        if len(idx) < args.min_count:
            continue
        z = q["z"][idx]
        vec = q["vec"][idx]
        baseline = fit_weighted_dipole(z, vec)
        axis = np.array(
            [
                math.cos(math.radians(baseline["b_deg"])) * math.cos(math.radians(baseline["l_deg"])),
                math.cos(math.radians(baseline["b_deg"])) * math.sin(math.radians(baseline["l_deg"])),
                math.sin(math.radians(baseline["b_deg"])),
            ],
            dtype=float,
        )
        hemisphere = (vec @ axis) >= 0.0
        if np.sum(hemisphere) < args.min_hemisphere or np.sum(~hemisphere) < args.min_hemisphere:
            continue
        weights, propensity, ess = overlap_weights(build_features(q, idx, variables), hemisphere, args.seed + window.window_index)
        strata = propensity_strata(propensity, args.n_strata)
        weighted = fit_weighted_dipole(z, vec, weights)
        priority = weighted["formal_snr"] * window.max_abs_su2_delta_mu
        x = np.column_stack([np.ones(len(z)), vec[:, 0], vec[:, 1], vec[:, 2]])
        w = weights / np.mean(weights[weights > 0.0])
        pinv_xtwx = np.linalg.pinv(x.T @ (x * w[:, None]))
        prepared.append(
            {
                "window": window,
                "idx": idx,
                "z": z,
                "vec": vec,
                "weights": weights,
                "strata": strata,
                "pinv_xtwx": pinv_xtwx,
            }
        )
        observed_rows.append(
            {
                "window_index": window.window_index,
                "tag": window.tag,
                "bcut_deg": window.bcut_deg,
                "z_center": window.z_center,
                "N": int(len(idx)),
                "ess": ess,
                "baseline_snr": baseline["formal_snr"],
                "baseline_amp": baseline["amp"],
                "weighted_snr": weighted["formal_snr"],
                "weighted_amp": weighted["amp"],
                "weighted_l_deg": weighted["l_deg"],
                "weighted_b_deg": weighted["b_deg"],
                "weighted_priority": priority,
                "amp_ratio_weighted_vs_baseline": weighted["amp"] / baseline["amp"] if baseline["amp"] > 0.0 else float("nan"),
                "direction_sep_weighted_vs_baseline_deg": angular_sep_deg(
                    np.array(
                        [
                            weighted["amp"] * math.cos(math.radians(weighted["b_deg"])) * math.cos(math.radians(weighted["l_deg"])),
                            weighted["amp"] * math.cos(math.radians(weighted["b_deg"])) * math.sin(math.radians(weighted["l_deg"])),
                            weighted["amp"] * math.sin(math.radians(weighted["b_deg"])),
                        ],
                        dtype=float,
                    ),
                    axis,
                ),
            }
        )
        if args.progress_every_windows and len(observed_rows) % args.progress_every_windows == 0:
            print(f"prepared {len(observed_rows)} matched windows", flush=True)
    return prepared, observed_rows


def best_template() -> dict[str, Any]:
    return {"window_index": -1, "tag": "", "bcut_deg": float("nan"), "value": -np.inf, "snr": float("nan"), "amp": float("nan")}


def update_best(best: dict[str, Any], value: float, window: Window, fit: dict[str, float]) -> dict[str, Any]:
    if math.isfinite(value) and value > float(best["value"]):
        return {
            "window_index": window.window_index,
            "tag": window.tag,
            "bcut_deg": window.bcut_deg,
            "value": value,
            "snr": fit["formal_snr"],
            "amp": fit["amp"],
            "l_deg": fit["l_deg"],
            "b_deg": fit["b_deg"],
        }
    return best


def run_mocks(args: argparse.Namespace, prepared: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rng = np.random.default_rng(args.seed)
    rows: list[dict[str, Any]] = []
    for mock_id in range(args.n_mocks):
        best_snr = best_template()
        best_priority = best_template()
        best_z1p00_snr = best_template()
        best_z1p00_priority = best_template()
        for item in prepared:
            window: Window = item["window"]
            z_mock = item["z"].copy()
            strata = item["strata"]
            for sid in np.unique(strata):
                sub = np.flatnonzero(strata == sid)
                if len(sub) > 1:
                    z_mock[sub] = rng.permutation(z_mock[sub])
            fit = weighted_fit_with_cache(z_mock, item["vec"], item["weights"], item["pinv_xtwx"])
            priority = fit["formal_snr"] * window.max_abs_su2_delta_mu
            best_snr = update_best(best_snr, fit["formal_snr"], window, fit)
            best_priority = update_best(best_priority, priority, window, fit)
            if window.tag == "z1p00_1p50":
                best_z1p00_snr = update_best(best_z1p00_snr, fit["formal_snr"], window, fit)
                best_z1p00_priority = update_best(best_z1p00_priority, priority, window, fit)
        row = {
            "mock_id": mock_id,
            "max_snr": best_snr["value"],
            "max_snr_tag": best_snr["tag"],
            "max_snr_bcut_deg": best_snr["bcut_deg"],
            "max_snr_amp": best_snr["amp"],
            "max_priority": best_priority["value"],
            "max_priority_tag": best_priority["tag"],
            "max_priority_bcut_deg": best_priority["bcut_deg"],
            "max_priority_snr": best_priority["snr"],
            "max_z1p00_snr": best_z1p00_snr["value"],
            "max_z1p00_snr_bcut_deg": best_z1p00_snr["bcut_deg"],
            "max_z1p00_priority": best_z1p00_priority["value"],
            "max_z1p00_priority_bcut_deg": best_z1p00_priority["bcut_deg"],
        }
        rows.append(row)
        if args.progress_every and ((mock_id + 1) % args.progress_every == 0 or mock_id == 0):
            print(
                f"mock {mock_id + 1}/{args.n_mocks}: max_snr={row['max_snr']:.4g} "
                f"({row['max_snr_tag']} |b|>{row['max_snr_bcut_deg']:g}), "
                f"max_priority={row['max_priority']:.5g}",
                flush=True,
            )
    return rows


def threshold_rows(observed_rows: list[dict[str, Any]], mock_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    max_snr = np.array([row["max_snr"] for row in mock_rows], dtype=float)
    max_priority = np.array([row["max_priority"] for row in mock_rows], dtype=float)
    max_z1p00_snr = np.array([row["max_z1p00_snr"] for row in mock_rows], dtype=float)
    max_z1p00_priority = np.array([row["max_z1p00_priority"] for row in mock_rows], dtype=float)

    targets: list[dict[str, Any]] = []
    targets.append(max(observed_rows, key=lambda row: row["weighted_snr"]))
    targets.append(max(observed_rows, key=lambda row: row["weighted_priority"]))
    targets.extend(
        sorted(
            [row for row in observed_rows if row["tag"] == "z1p00_1p50" and row["bcut_deg"] in (10.0, 15.0, 25.0, 35.0)],
            key=lambda row: row["bcut_deg"],
        )
    )
    out: list[dict[str, Any]] = []
    seen: set[tuple[str, float, str]] = set()
    for target in targets:
        role = "target"
        if target is targets[0]:
            role = "matched_observed_max_snr"
        if target is targets[1]:
            role = "matched_observed_max_priority" if role == "target" else role + "+matched_observed_max_priority"
        key = (target["tag"], target["bcut_deg"], role)
        if key in seen:
            continue
        seen.add(key)
        out.append(
            {
                "role": role,
                "tag": target["tag"],
                "bcut_deg": target["bcut_deg"],
                "N": target["N"],
                "ess": target["ess"],
                "weighted_snr": target["weighted_snr"],
                "weighted_priority": target["weighted_priority"],
                "weighted_amp": target["weighted_amp"],
                "amp_ratio_weighted_vs_baseline": target["amp_ratio_weighted_vs_baseline"],
                "direction_sep_weighted_vs_baseline_deg": target["direction_sep_weighted_vs_baseline_deg"],
                "global_p_any_window_snr_ge_observed": p_ge(max_snr, target["weighted_snr"]),
                "global_p_any_window_priority_ge_observed": p_ge(max_priority, target["weighted_priority"]),
                "z1p00_family_p_any_bcut_snr_ge_observed": p_ge(max_z1p00_snr, target["weighted_snr"]) if target["tag"] == "z1p00_1p50" else float("nan"),
                "z1p00_family_p_any_bcut_priority_ge_observed": p_ge(max_z1p00_priority, target["weighted_priority"]) if target["tag"] == "z1p00_1p50" else float("nan"),
            }
        )
    return out


def write_reports(args: argparse.Namespace, observed_rows: list[dict[str, Any]], mock_rows: list[dict[str, Any]], thresholds: list[dict[str, Any]]) -> None:
    max_snr = np.array([row["max_snr"] for row in mock_rows], dtype=float)
    max_priority = np.array([row["max_priority"] for row in mock_rows], dtype=float)
    top = sorted(observed_rows, key=lambda row: row["weighted_priority"], reverse=True)[:12]
    config = {
        "date": "2026-07-15",
        "analysis": "matched_quality_global_lookelsewhere",
        "quaia_fits": str(args.quaia_fits),
        "scan_csv": str(args.scan_csv),
        "n_mocks": args.n_mocks,
        "n_strata": args.n_strata,
        "min_count": args.min_count,
        "seed": args.seed,
        "quality_variables": ["zerr", "g", "bp_rp", "w1", "w2", "w1_w2", "pmra_error", "pmdec_error"],
        "method": "For each scanned redshift-window/latitude-cut row, fit the baseline dipole, construct all-catalogue-quality overlap weights across the baseline hemispheres, refit the weighted dipole, and compare scan maxima to propensity-stratified redshift shuffles.",
    }
    (args.out_dir / "su2_quaia_matched_global_lookelsewhere_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    report = [
        "# SU2 / Quaia Matched-Quality Global Look-Elsewhere Audit",
        "",
        "Date: July 15, 2026",
        "",
        "## Purpose",
        "",
        "This test repeats the all-window SU2/Quaia look-elsewhere scan after all-catalogue-quality overlap weighting. It asks whether the targeted z ~ 1-1.5 mode remains globally unusual once the catalogue-quality mechanism gate is imposed across the whole scan.",
        "",
        "## Model",
        "",
        "```text",
        "z_i      = c0 + d . n_i + epsilon_i",
        "w_i      = overlap propensity weight for all catalogue-quality variables",
        "SNR_d    = sqrt(d^T Cov(d)^(-1) d)",
        "priority = SNR_d * max_abs_delta_mu_SU2(window)",
        "```",
        "",
        "## Configuration",
        "",
        f"- scan windows prepared: `{len(observed_rows)}`",
        f"- redshift-shuffle mocks: `{len(mock_rows)}`",
        f"- propensity strata: `{args.n_strata}`",
        "",
        "## Mock Maxima",
        "",
        f"- max matched SNR mean / p95 / p99: `{np.nanmean(max_snr):.6g}` / `{qtile(max_snr, 0.95):.6g}` / `{qtile(max_snr, 0.99):.6g}`",
        f"- max matched priority mean / p95 / p99: `{np.nanmean(max_priority):.6g}` / `{qtile(max_priority, 0.95):.6g}` / `{qtile(max_priority, 0.99):.6g}`",
        "",
        "## Observed Matched Thresholds",
        "",
        markdown_table(
            thresholds,
            [
                "role",
                "tag",
                "bcut_deg",
                "N",
                "ess",
                "weighted_snr",
                "weighted_priority",
                "global_p_any_window_snr_ge_observed",
                "global_p_any_window_priority_ge_observed",
                "z1p00_family_p_any_bcut_snr_ge_observed",
            ],
        ),
        "",
        "## Top Matched Observed Windows",
        "",
        markdown_table(
            top,
            [
                "tag",
                "bcut_deg",
                "N",
                "ess",
                "weighted_snr",
                "weighted_priority",
                "amp_ratio_weighted_vs_baseline",
                "direction_sep_weighted_vs_baseline_deg",
            ],
        ),
        "",
        "## Interpretation Gate",
        "",
        "- Passing condition: the locked z ~ 1-1.5 family remains small after global look-elsewhere correction under the matched-quality null.",
        "- Failing condition: global p-values are ordinary or the matched weighted mode is no longer the scan maximum.",
        "- This is a promotion gate, not a discovery test by itself.",
    ]
    (args.out_dir / "su2_quaia_matched_global_lookelsewhere_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    readout = [
        "# SU2 / Quaia Matched-Quality Global Look-Elsewhere Readout",
        "",
        "Date: July 15, 2026",
        "",
        "## Bottom Line",
        "",
    ]
    zrows = [row for row in thresholds if row["tag"] == "z1p00_1p50"]
    if zrows:
        best_z = min(zrows, key=lambda row: row["global_p_any_window_snr_ge_observed"])
        readout.append(
            f"The best matched z1p00_1p50 threshold has global SNR p = `{fmt(best_z['global_p_any_window_snr_ge_observed'])}` "
            f"and z1p00-family p = `{fmt(best_z['z1p00_family_p_any_bcut_snr_ge_observed'])}`."
        )
    else:
        readout.append("No z1p00_1p50 matched threshold row was available.")
    readout.extend(["", "## Thresholds", "", markdown_table(thresholds, ["role", "tag", "bcut_deg", "weighted_snr", "weighted_priority", "global_p_any_window_snr_ge_observed", "z1p00_family_p_any_bcut_snr_ge_observed"])])
    (args.out_dir / "su2_quaia_matched_global_lookelsewhere_readout.md").write_text("\n".join(readout) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Matched-quality global look-elsewhere audit for SU2/Quaia.")
    parser.add_argument("--quaia-fits", type=Path, default=DEFAULT_QUAIA)
    parser.add_argument("--scan-csv", type=Path, default=DEFAULT_SCAN)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--n-mocks", type=int, default=200)
    parser.add_argument("--n-strata", type=int, default=10)
    parser.add_argument("--min-count", type=int, default=50)
    parser.add_argument("--min-hemisphere", type=int, default=100)
    parser.add_argument("--seed", type=int, default=180715)
    parser.add_argument("--progress-every", type=int, default=10)
    parser.add_argument("--progress-every-windows", type=int, default=25)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    q = load_quaia(args.quaia_fits)
    windows = read_windows(args.scan_csv, args.min_count)
    print(f"Loaded Quaia rows: {len(q['z'])}", flush=True)
    print(f"Candidate scan windows: {len(windows)}", flush=True)
    prepared, observed_rows = prepare_windows(args, q, windows)
    print(f"Prepared matched-quality windows: {len(prepared)}", flush=True)
    if not prepared:
        raise RuntimeError("No matched-quality windows prepared.")
    mock_rows = run_mocks(args, prepared)
    thresholds = threshold_rows(observed_rows, mock_rows)
    write_csv(args.out_dir / "su2_quaia_matched_global_lookelsewhere_observed_windows.csv", observed_rows)
    write_csv(args.out_dir / "su2_quaia_matched_global_lookelsewhere_mock_maxima.csv", mock_rows)
    write_csv(args.out_dir / "su2_quaia_matched_global_lookelsewhere_thresholds.csv", thresholds)
    write_reports(args, observed_rows, mock_rows, thresholds)
    print(f"Saved report: {args.out_dir / 'su2_quaia_matched_global_lookelsewhere_report.md'}", flush=True)
    print(f"Saved readout: {args.out_dir / 'su2_quaia_matched_global_lookelsewhere_readout.md'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
