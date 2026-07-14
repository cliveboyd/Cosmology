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

from quaia_redshift_window_full_mocks import (
    DEFAULT_QUAIA,
    DEFAULT_RANDOMS,
    DEFAULT_SCAN,
    CMB_L_DEG,
    CMB_B_DEG,
    fit_redshift_dipole,
    gal_unit_from_lb,
    load_quaia,
    load_randoms,
)


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_global_lookelsewhere_20260714"


@dataclass(frozen=True)
class ScanWindow:
    window_index: int
    bcut_deg: float
    tag: str
    z_lo: float
    z_hi: float
    z_center: float
    observed_N: int
    observed_amp: float
    observed_abs_par: float
    observed_formal_snr: float
    observed_priority: float
    max_abs_su2_delta_mu: float
    observed_l_deg: float
    observed_b_deg: float
    observed_sep_cmb_deg: float


def parse_float(value: Any, default: float = float("nan")) -> float:
    if value in ("", None):
        return default
    try:
        return float(value)
    except Exception:
        return default


def p_ge(values: np.ndarray, observed: float) -> float:
    values = values[np.isfinite(values)]
    if len(values) == 0 or not math.isfinite(observed):
        return float("nan")
    return float((1 + np.sum(values >= observed)) / (1 + len(values)))


def q(values: np.ndarray, prob: float) -> float:
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return float("nan")
    return float(np.quantile(values, prob))


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


def read_scan_windows(scan_csv: Path, min_count: int) -> list[ScanWindow]:
    df = pd.read_csv(scan_csv)
    df = df[df["tag"].astype(str).str.startswith("z")].copy()
    df["N"] = pd.to_numeric(df["N"], errors="coerce")
    df = df[df["N"] >= min_count].copy()
    df = df.sort_values(["bcut_deg", "z_lo", "z_hi", "tag"]).reset_index(drop=True)
    windows: list[ScanWindow] = []
    for i, row in df.iterrows():
        delta = parse_float(row.get("max_abs_su2_delta_mu"), 0.0)
        snr = parse_float(row.get("formal_snr"), 0.0)
        priority = parse_float(row.get("su2_quaia_priority_score"), snr * delta)
        windows.append(
            ScanWindow(
                window_index=i + 1,
                bcut_deg=parse_float(row.get("bcut_deg")),
                tag=str(row.get("tag")),
                z_lo=parse_float(row.get("z_lo")),
                z_hi=parse_float(row.get("z_hi")),
                z_center=parse_float(row.get("z_center")),
                observed_N=int(parse_float(row.get("N"), 0.0)),
                observed_amp=parse_float(row.get("amp")),
                observed_abs_par=abs(parse_float(row.get("amp_par_cmb"))),
                observed_formal_snr=snr,
                observed_priority=priority,
                max_abs_su2_delta_mu=delta,
                observed_l_deg=parse_float(row.get("l_deg")),
                observed_b_deg=parse_float(row.get("b_deg")),
                observed_sep_cmb_deg=parse_float(row.get("sep_cmb_deg")),
            )
        )
    return windows


def prepare_window_samples(
    windows: list[ScanWindow],
    quaia: dict[str, np.ndarray],
    randoms: dict[str, np.ndarray],
    min_count: int,
) -> list[dict[str, Any]]:
    pools_by_bcut: dict[float, np.ndarray] = {}
    prepared: list[dict[str, Any]] = []
    for window in windows:
        obs_mask = (
            (quaia["z"] >= window.z_lo)
            & (quaia["z"] < window.z_hi)
            & (np.abs(quaia["b_gal"]) >= window.bcut_deg)
        )
        obs_idx = np.flatnonzero(obs_mask)
        if len(obs_idx) < min_count:
            continue
        if window.bcut_deg not in pools_by_bcut:
            pools_by_bcut[window.bcut_deg] = np.flatnonzero(np.abs(randoms["b_gal"]) >= window.bcut_deg)
        rand_pool = pools_by_bcut[window.bcut_deg]
        if len(rand_pool) == 0:
            continue
        prepared.append(
            {
                "window": window,
                "z": np.asarray(quaia["z"][obs_idx], dtype=float),
                "rand_pool": rand_pool,
                "replace_randoms": len(rand_pool) < len(obs_idx),
            }
        )
    return prepared


def initial_max_row() -> dict[str, Any]:
    return {
        "window_index": -1,
        "tag": "",
        "bcut_deg": float("nan"),
        "z_center": float("nan"),
        "N": 0,
        "value": -np.inf,
        "amp": float("nan"),
        "formal_snr": float("nan"),
        "priority": float("nan"),
        "abs_par": float("nan"),
        "l_deg": float("nan"),
        "b_deg": float("nan"),
        "sep_cmb_deg": float("nan"),
    }


def update_best(best: dict[str, Any], value: float, window: ScanWindow, stats: dict[str, float], priority: float) -> dict[str, Any]:
    if math.isfinite(value) and value > float(best["value"]):
        return {
            "window_index": window.window_index,
            "tag": window.tag,
            "bcut_deg": window.bcut_deg,
            "z_center": window.z_center,
            "N": window.observed_N,
            "value": value,
            "amp": stats["amp"],
            "formal_snr": stats["formal_snr"],
            "priority": priority,
            "abs_par": abs(stats["amp_par_cmb"]),
            "l_deg": stats["l_deg"],
            "b_deg": stats["b_deg"],
            "sep_cmb_deg": stats["sep_cmb_deg"],
        }
    return best


def run_global_mocks(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[ScanWindow]]:
    args.out.mkdir(parents=True, exist_ok=True)
    quaia = load_quaia(args.quaia_csv)
    cache_path = args.out / "quaia_randoms_galactic_cache.npz" if args.cache_randoms else None
    randoms = load_randoms(args.randoms_fits, cache_path)
    windows = read_scan_windows(args.scan_csv, args.min_count)
    prepared = prepare_window_samples(windows, quaia, randoms, args.min_count)
    if not prepared:
        raise RuntimeError("No scan windows passed the data/mask filters.")

    print(f"Loaded Quaia rows: {len(quaia['z'])}", flush=True)
    print(f"Prepared scan windows: {len(prepared)}", flush=True)
    print(f"Mock count: {args.n_mocks}", flush=True)

    rng = np.random.default_rng(args.seed)
    rows: list[dict[str, Any]] = []
    cmb = gal_unit_from_lb(CMB_L_DEG, CMB_B_DEG)

    for mock_id in range(args.n_mocks):
        best_snr = initial_max_row()
        best_priority = initial_max_row()
        best_abs_par = initial_max_row()
        best_z1p00_snr = initial_max_row()
        best_z1p00_priority = initial_max_row()

        for item in prepared:
            window: ScanWindow = item["window"]
            z = item["z"]
            pool = item["rand_pool"]
            rand_idx = rng.choice(pool, size=len(z), replace=item["replace_randoms"])
            stats = fit_redshift_dipole(z, randoms["vec"][rand_idx])
            priority = float(stats["formal_snr"] * window.max_abs_su2_delta_mu)
            abs_par = abs(float(np.dot(
                np.array(
                    [
                        stats["amp"] * math.cos(math.radians(stats["b_deg"])) * math.cos(math.radians(stats["l_deg"])),
                        stats["amp"] * math.cos(math.radians(stats["b_deg"])) * math.sin(math.radians(stats["l_deg"])),
                        stats["amp"] * math.sin(math.radians(stats["b_deg"])),
                    ],
                    dtype=float,
                ),
                cmb,
            )))

            best_snr = update_best(best_snr, stats["formal_snr"], window, stats, priority)
            best_priority = update_best(best_priority, priority, window, stats, priority)
            best_abs_par = update_best(best_abs_par, abs_par, window, stats, priority)
            if window.tag == "z1p00_1p50":
                best_z1p00_snr = update_best(best_z1p00_snr, stats["formal_snr"], window, stats, priority)
                best_z1p00_priority = update_best(best_z1p00_priority, priority, window, stats, priority)

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
            "max_priority_snr": best_priority["formal_snr"],
            "max_priority_amp": best_priority["amp"],
            "max_priority_l_deg": best_priority["l_deg"],
            "max_priority_b_deg": best_priority["b_deg"],
            "max_priority_sep_cmb_deg": best_priority["sep_cmb_deg"],
            "max_abs_par": best_abs_par["value"],
            "max_abs_par_window_index": best_abs_par["window_index"],
            "max_abs_par_tag": best_abs_par["tag"],
            "max_abs_par_bcut_deg": best_abs_par["bcut_deg"],
            "max_abs_par_snr": best_abs_par["formal_snr"],
            "max_z1p00_snr": best_z1p00_snr["value"],
            "max_z1p00_snr_bcut_deg": best_z1p00_snr["bcut_deg"],
            "max_z1p00_priority": best_z1p00_priority["value"],
            "max_z1p00_priority_bcut_deg": best_z1p00_priority["bcut_deg"],
        }
        rows.append(row)
        if (mock_id + 1) % args.progress_every == 0 or mock_id == 0:
            print(
                f"mock {mock_id + 1}/{args.n_mocks}: "
                f"max_snr={row['max_snr']:.4g} ({row['max_snr_tag']} |b|>{row['max_snr_bcut_deg']:g}), "
                f"max_priority={row['max_priority']:.5g} ({row['max_priority_tag']} |b|>{row['max_priority_bcut_deg']:g})",
                flush=True,
            )
    return rows, [item["window"] for item in prepared]


def observed_threshold_rows(mock_rows: list[dict[str, Any]], windows: list[ScanWindow]) -> list[dict[str, Any]]:
    mock_max_snr = np.array([row["max_snr"] for row in mock_rows], dtype=float)
    mock_max_priority = np.array([row["max_priority"] for row in mock_rows], dtype=float)
    mock_max_abs_par = np.array([row["max_abs_par"] for row in mock_rows], dtype=float)
    mock_max_z1p00_snr = np.array([row["max_z1p00_snr"] for row in mock_rows], dtype=float)
    mock_max_z1p00_priority = np.array([row["max_z1p00_priority"] for row in mock_rows], dtype=float)

    targets: list[ScanWindow] = []
    targets.append(max(windows, key=lambda w: w.observed_formal_snr))
    targets.append(max(windows, key=lambda w: w.observed_priority))
    targets.extend(
        sorted(
            [w for w in windows if w.tag == "z1p00_1p50" and w.bcut_deg in (10.0, 15.0, 25.0, 35.0)],
            key=lambda w: (w.bcut_deg, w.tag),
        )
    )

    out: list[dict[str, Any]] = []
    seen: set[tuple[str, float, str]] = set()
    for target in targets:
        role = "target"
        if target == max(windows, key=lambda w: w.observed_formal_snr):
            role = "observed_max_snr"
        if target == max(windows, key=lambda w: w.observed_priority):
            role = "observed_max_priority" if role == "target" else role + "+observed_max_priority"
        key = (target.tag, target.bcut_deg, role)
        if key in seen:
            continue
        seen.add(key)
        out.append(
            {
                "role": role,
                "window_index": target.window_index,
                "tag": target.tag,
                "bcut_deg": target.bcut_deg,
                "z_center": target.z_center,
                "N": target.observed_N,
                "observed_amp": target.observed_amp,
                "observed_abs_par_cmb": target.observed_abs_par,
                "observed_formal_snr": target.observed_formal_snr,
                "observed_priority": target.observed_priority,
                "max_abs_su2_delta_mu": target.max_abs_su2_delta_mu,
                "observed_l_deg": target.observed_l_deg,
                "observed_b_deg": target.observed_b_deg,
                "observed_sep_cmb_deg": target.observed_sep_cmb_deg,
                "global_p_any_window_snr_ge_observed": p_ge(mock_max_snr, target.observed_formal_snr),
                "global_p_any_window_priority_ge_observed": p_ge(mock_max_priority, target.observed_priority),
                "global_p_any_window_abs_par_ge_observed": p_ge(mock_max_abs_par, target.observed_abs_par),
                "z1p00_family_p_any_bcut_snr_ge_observed": p_ge(mock_max_z1p00_snr, target.observed_formal_snr) if target.tag == "z1p00_1p50" else float("nan"),
                "z1p00_family_p_any_bcut_priority_ge_observed": p_ge(mock_max_z1p00_priority, target.observed_priority) if target.tag == "z1p00_1p50" else float("nan"),
            }
        )
    return out


def markdown_table(rows: list[dict[str, Any]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        vals: list[str] = []
        for col in columns:
            value = row.get(col, "")
            if isinstance(value, float):
                vals.append(f"{value:.6g}" if math.isfinite(value) else "")
            else:
                vals.append(str(value))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def write_report(path: Path, mock_rows: list[dict[str, Any]], threshold_rows: list[dict[str, Any]], args: argparse.Namespace) -> None:
    max_snr = np.array([row["max_snr"] for row in mock_rows], dtype=float)
    max_priority = np.array([row["max_priority"] for row in mock_rows], dtype=float)
    max_abs_par = np.array([row["max_abs_par"] for row in mock_rows], dtype=float)
    lines = [
        "# SU2 / Quaia Global Look-Elsewhere Mock Audit",
        "",
        "Date: July 14, 2026",
        "",
        "## Purpose",
        "",
        "This run tests whether the strongest SU2/Quaia redshift-window targets remain unusual after accounting for the fact that the scan searched many redshift windows and Galactic latitude cuts.",
        "",
        "For each mock, every non-full redshift-window row from the original SU2/Quaia scan is regenerated using the observed window redshift distribution and random Quaia sky positions under the same latitude cut. The recorded null statistic is the maximum over the scan.",
        "",
        "## Configuration",
        "",
        f"- scan CSV: `{args.scan_csv}`",
        f"- Quaia table: `{args.quaia_csv}`",
        f"- random sky catalog: `{args.randoms_fits}`",
        f"- n_mocks: `{args.n_mocks}`",
        f"- seed: `{args.seed}`",
        f"- min_count: `{args.min_count}`",
        "",
        "## Mock Maxima",
        "",
        f"- max formal SNR mean / p95 / p99: `{np.nanmean(max_snr):.6g}` / `{q(max_snr, 0.95):.6g}` / `{q(max_snr, 0.99):.6g}`",
        f"- max SU2-priority mean / p95 / p99: `{np.nanmean(max_priority):.6g}` / `{q(max_priority, 0.95):.6g}` / `{q(max_priority, 0.99):.6g}`",
        f"- max CMB-projected component mean / p95 / p99: `{np.nanmean(max_abs_par):.6g}` / `{q(max_abs_par, 0.95):.6g}` / `{q(max_abs_par, 0.99):.6g}`",
        "",
        "## Observed Thresholds",
        "",
        markdown_table(
            threshold_rows,
            [
                "role",
                "tag",
                "bcut_deg",
                "N",
                "observed_formal_snr",
                "observed_priority",
                "global_p_any_window_snr_ge_observed",
                "global_p_any_window_priority_ge_observed",
                "z1p00_family_p_any_bcut_snr_ge_observed",
            ],
        ),
        "",
        "## Interpretation Gate",
        "",
        "- If the global p-values are ordinary, the SU2/Quaia result remains a targeting diagnostic, not evidence for anisotropy.",
        "- If a z1p00_1p50 target stays small under the global scan correction, the next gate is a selection-function/template regression.",
        "- CMB-projected components remain controls unless their global p-values become small under the same protocol.",
        "",
        "## Outputs",
        "",
        f"- per-mock maxima: `{args.out / 'su2_quaia_global_lookelsewhere_mock_maxima.csv'}`",
        f"- observed thresholds: `{args.out / 'su2_quaia_global_lookelsewhere_observed_thresholds.csv'}`",
        f"- config: `{args.out / 'su2_quaia_global_lookelsewhere_config.json'}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Global look-elsewhere mock audit for the SU2/Quaia scan.")
    parser.add_argument("--quaia-csv", type=Path, default=DEFAULT_QUAIA)
    parser.add_argument("--randoms-fits", type=Path, default=DEFAULT_RANDOMS)
    parser.add_argument("--scan-csv", type=Path, default=DEFAULT_SCAN)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--n-mocks", type=int, default=500)
    parser.add_argument("--min-count", type=int, default=50)
    parser.add_argument("--seed", type=int, default=280714)
    parser.add_argument("--progress-every", type=int, default=10)
    parser.add_argument("--cache-randoms", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    mock_rows, windows = run_global_mocks(args)
    threshold_rows = observed_threshold_rows(mock_rows, windows)
    write_csv(args.out / "su2_quaia_global_lookelsewhere_mock_maxima.csv", mock_rows)
    write_csv(args.out / "su2_quaia_global_lookelsewhere_observed_thresholds.csv", threshold_rows)
    (args.out / "su2_quaia_global_lookelsewhere_config.json").write_text(json.dumps(vars(args), indent=2, default=str), encoding="utf-8")
    write_report(args.out / "su2_quaia_global_lookelsewhere_report.md", mock_rows, threshold_rows, args)
    print(f"Saved per-mock maxima: {args.out / 'su2_quaia_global_lookelsewhere_mock_maxima.csv'}", flush=True)
    print(f"Saved observed thresholds: {args.out / 'su2_quaia_global_lookelsewhere_observed_thresholds.csv'}", flush=True)
    print(f"Saved report: {args.out / 'su2_quaia_global_lookelsewhere_report.md'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
