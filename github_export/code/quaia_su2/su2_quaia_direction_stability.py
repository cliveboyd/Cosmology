from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
DEFAULT_QUAIA = ROOT / "quaia_G20p0_basic_gal.csv"
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_direction_stability_20260714"


def gal_unit_array(l_deg: np.ndarray, b_deg: np.ndarray) -> np.ndarray:
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    cb = np.cos(b)
    return np.column_stack([cb * np.cos(l), cb * np.sin(l), np.sin(b)])


def unit_from_lb(l_deg: float, b_deg: float) -> np.ndarray:
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


def fit_redshift_dipole(z: np.ndarray, vec: np.ndarray) -> dict[str, float]:
    n = len(z)
    sx = np.sum(vec, axis=0)
    xtx = np.empty((4, 4), dtype=float)
    xtx[0, 0] = float(n)
    xtx[0, 1:] = sx
    xtx[1:, 0] = sx
    xtx[1:, 1:] = vec.T @ vec
    xty = np.empty(4, dtype=float)
    xty[0] = float(np.sum(z))
    xty[1:] = vec.T @ z
    beta = np.linalg.pinv(xtx) @ xty
    yty = float(np.dot(z, z))
    sse = max(yty - float(beta @ xty), 0.0)
    dof = max(n - 4, 1)
    sigma2 = sse / dof
    dvec = beta[1:4]
    amp = float(np.linalg.norm(dvec))
    l_deg, b_deg = unit_to_lb(dvec)
    try:
        cov_d = sigma2 * np.linalg.pinv(xtx)[1:, 1:]
        snr = float(math.sqrt(max(float(dvec @ np.linalg.pinv(cov_d) @ dvec), 0.0)))
    except Exception:
        snr = float("nan")
    return {
        "a0": float(beta[0]),
        "amp": amp,
        "l_deg": l_deg,
        "b_deg": b_deg,
        "formal_snr": snr,
        "rms_resid": float(math.sqrt(sigma2)),
    }


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


def markdown_table(rows: list[dict[str, Any]], columns: list[str], max_rows: int | None = None) -> str:
    rows = rows[:max_rows] if max_rows is not None else rows
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


def window_grid() -> list[dict[str, Any]]:
    windows: list[dict[str, Any]] = []
    fixed = [
        (0.90, 1.40, "shift_m0p10_width0p50"),
        (0.95, 1.45, "shift_m0p05_width0p50"),
        (1.00, 1.50, "baseline_width0p50"),
        (1.05, 1.55, "shift_p0p05_width0p50"),
        (1.10, 1.60, "shift_p0p10_width0p50"),
    ]
    for z_lo, z_hi, tag in fixed:
        windows.append({"tag": tag, "z_lo": z_lo, "z_hi": z_hi, "z_center": 0.5 * (z_lo + z_hi), "width": z_hi - z_lo, "family": "width0p50"})
    for width in (0.20, 0.30):
        for center in np.arange(1.05, 1.56, 0.05):
            z_lo = round(float(center - 0.5 * width), 4)
            z_hi = round(float(center + 0.5 * width), 4)
            tag = f"slide_c{center:.2f}_w{width:.2f}".replace(".", "p")
            windows.append({"tag": tag, "z_lo": z_lo, "z_hi": z_hi, "z_center": float(center), "width": width, "family": f"width{width:.2f}".replace(".", "p")})
    return windows


def run(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    df = pd.read_csv(args.quaia_csv)
    cols = {c.upper(): c for c in df.columns}
    z = df[cols["Z"]].to_numpy(float)
    l = df[cols["L_GAL"]].to_numpy(float)
    b = df[cols["B_GAL"]].to_numpy(float)
    vec = gal_unit_array(l, b)
    ok = np.isfinite(z) & np.isfinite(l) & np.isfinite(b) & (z > 0.0)
    z = z[ok]
    b = b[ok]
    vec = vec[ok]
    windows = window_grid()
    bcuts = [float(x) for x in args.bcuts]

    ref_mask = (z >= 1.0) & (z < 1.5) & (np.abs(b) >= 35.0)
    ref_stats = fit_redshift_dipole(z[ref_mask], vec[ref_mask])
    ref_vec = unit_from_lb(ref_stats["l_deg"], ref_stats["b_deg"])
    cmb_vec = unit_from_lb(264.021, 48.253)

    fit_rows: list[dict[str, Any]] = []
    for window in windows:
        for bcut in bcuts:
            mask = (z >= window["z_lo"]) & (z < window["z_hi"]) & (np.abs(b) >= bcut)
            n = int(np.sum(mask))
            if n < args.min_count:
                continue
            stats = fit_redshift_dipole(z[mask], vec[mask])
            dvec = unit_from_lb(stats["l_deg"], stats["b_deg"])
            fit_rows.append(
                {
                    **window,
                    "bcut_deg": bcut,
                    "N": n,
                    **stats,
                    "sep_from_baseline_b35_deg": angular_sep_deg(dvec, ref_vec),
                    "sep_from_cmb_deg": angular_sep_deg(dvec, cmb_vec),
                    "amp_ratio_vs_baseline_b35": stats["amp"] / ref_stats["amp"] if ref_stats["amp"] > 0 else float("nan"),
                    "reference_l_deg": ref_stats["l_deg"],
                    "reference_b_deg": ref_stats["b_deg"],
                    "reference_amp": ref_stats["amp"],
                    "reference_snr": ref_stats["formal_snr"],
                }
            )
            print(
                f"{window['tag']} |b|>{bcut:g}: N={n} amp={stats['amp']:.4g} "
                f"snr={stats['formal_snr']:.3g} sep_ref={fit_rows[-1]['sep_from_baseline_b35_deg']:.1f}",
                flush=True,
            )

    summary_rows: list[dict[str, Any]] = []
    fit_df = pd.DataFrame(fit_rows)
    for (tag, family), group in fit_df.groupby(["tag", "family"], sort=False):
        amps = group["amp"].to_numpy(float)
        snrs = group["formal_snr"].to_numpy(float)
        seps = group["sep_from_baseline_b35_deg"].to_numpy(float)
        dirs = np.array([unit_from_lb(float(row.l_deg), float(row.b_deg)) for row in group.itertuples()], dtype=float)
        max_pair_sep = 0.0
        for i in range(len(dirs)):
            for j in range(i + 1, len(dirs)):
                max_pair_sep = max(max_pair_sep, angular_sep_deg(dirs[i], dirs[j]))
        mean_dir = np.sum(dirs, axis=0)
        mean_resultant = float(np.linalg.norm(mean_dir) / max(len(dirs), 1))
        summary_rows.append(
            {
                "tag": tag,
                "family": family,
                "z_lo": float(group["z_lo"].iloc[0]),
                "z_hi": float(group["z_hi"].iloc[0]),
                "z_center": float(group["z_center"].iloc[0]),
                "width": float(group["width"].iloc[0]),
                "n_bcuts": int(len(group)),
                "N_min": int(group["N"].min()),
                "N_max": int(group["N"].max()),
                "amp_min": float(np.nanmin(amps)),
                "amp_max": float(np.nanmax(amps)),
                "amp_mean": float(np.nanmean(amps)),
                "amp_cv": float(np.nanstd(amps) / np.nanmean(amps)) if np.nanmean(amps) > 0 else float("nan"),
                "snr_max": float(np.nanmax(snrs)),
                "sep_ref_mean_deg": float(np.nanmean(seps)),
                "sep_ref_max_deg": float(np.nanmax(seps)),
                "max_pair_direction_sep_deg": float(max_pair_sep),
                "mean_resultant_length": mean_resultant,
                "passes_30deg_ref_gate": bool(np.nanmax(seps) <= 30.0),
                "passes_45deg_pair_gate": bool(max_pair_sep <= 45.0),
            }
        )

    bcut_rows: list[dict[str, Any]] = []
    for bcut, group in fit_df.groupby("bcut_deg", sort=True):
        base = group[group["tag"] == "baseline_width0p50"]
        if base.empty:
            continue
        base_amp = float(base["amp"].iloc[0])
        dirs = np.array([unit_from_lb(float(row.l_deg), float(row.b_deg)) for row in group.itertuples()], dtype=float)
        base_dir = unit_from_lb(float(base["l_deg"].iloc[0]), float(base["b_deg"].iloc[0]))
        seps_base = np.array([angular_sep_deg(d, base_dir) for d in dirs], dtype=float)
        bcut_rows.append(
            {
                "bcut_deg": float(bcut),
                "n_windows": int(len(group)),
                "baseline_amp": base_amp,
                "baseline_snr": float(base["formal_snr"].iloc[0]),
                "amp_min": float(group["amp"].min()),
                "amp_max": float(group["amp"].max()),
                "amp_ratio_min_vs_bcut_baseline": float(group["amp"].min() / base_amp) if base_amp > 0 else float("nan"),
                "amp_ratio_max_vs_bcut_baseline": float(group["amp"].max() / base_amp) if base_amp > 0 else float("nan"),
                "sep_from_bcut_baseline_mean_deg": float(np.nanmean(seps_base)),
                "sep_from_bcut_baseline_max_deg": float(np.nanmax(seps_base)),
                "passes_30deg_bcut_baseline_gate": bool(np.nanmax(seps_base) <= 30.0),
            }
        )
    return fit_rows, summary_rows, bcut_rows


def write_report(path: Path, fit_rows: list[dict[str, Any]], summary_rows: list[dict[str, Any]], bcut_rows: list[dict[str, Any]], args: argparse.Namespace) -> None:
    fixed = [row for row in summary_rows if row["family"] == "width0p50"]
    most_stable = sorted(summary_rows, key=lambda row: (row["sep_ref_max_deg"], row["amp_cv"]))[:10]
    bcut_sorted = sorted(bcut_rows, key=lambda row: row["bcut_deg"])
    lines = [
        "# SU2 / Quaia Direction-Amplitude Stability Grid",
        "",
        "Date: July 14, 2026",
        "",
        "## Purpose",
        "",
        "This gate tests whether the remaining `z ~ 1.0-1.5` SU2/Quaia angular hook behaves like a stable angular mode under nearby redshift-window and Galactic-latitude perturbations.",
        "",
        "The reference direction is the strongest original family row: `1.0 <= z < 1.5`, `|b| > 35`.",
        "",
        "## Fixed Width 0.50 Window Shifts",
        "",
        markdown_table(
            fixed,
            [
                "tag",
                "z_lo",
                "z_hi",
                "amp_mean",
                "amp_cv",
                "snr_max",
                "sep_ref_max_deg",
                "max_pair_direction_sep_deg",
                "passes_30deg_ref_gate",
                "passes_45deg_pair_gate",
            ],
        ),
        "",
        "## Most Stable Rows By Reference Separation",
        "",
        markdown_table(
            most_stable,
            [
                "tag",
                "family",
                "z_lo",
                "z_hi",
                "amp_mean",
                "amp_cv",
                "snr_max",
                "sep_ref_max_deg",
                "max_pair_direction_sep_deg",
            ],
            max_rows=10,
        ),
        "",
        "## Stability By Latitude Cut",
        "",
        markdown_table(
            bcut_sorted,
            [
                "bcut_deg",
                "n_windows",
                "baseline_amp",
                "baseline_snr",
                "amp_ratio_min_vs_bcut_baseline",
                "amp_ratio_max_vs_bcut_baseline",
                "sep_from_bcut_baseline_max_deg",
                "passes_30deg_bcut_baseline_gate",
            ],
        ),
        "",
        "## Interpretation Gate",
        "",
        "- A stable mode should keep direction clustered under small redshift-window shifts and latitude-cut changes.",
        "- Smooth amplitude changes are acceptable; sharp direction flips or large window sensitivity are not.",
        "- This observed-grid test is diagnostic only. Any candidate-stable region still needs mock stability percentiles and external selection templates.",
        "",
        "## Outputs",
        "",
        f"- fit grid: `{args.out / 'su2_quaia_direction_stability_fit_grid.csv'}`",
        f"- window summary: `{args.out / 'su2_quaia_direction_stability_window_summary.csv'}`",
        f"- bcut summary: `{args.out / 'su2_quaia_direction_stability_bcut_summary.csv'}`",
        f"- report: `{args.out / 'su2_quaia_direction_stability_report.md'}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Observed direction/amplitude stability grid for SU2/Quaia z~1.0-1.5.")
    parser.add_argument("--quaia-csv", type=Path, default=DEFAULT_QUAIA)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--bcuts", type=float, nargs="+", default=[10, 15, 20, 25, 30, 35])
    parser.add_argument("--min-count", type=int, default=1000)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    fit_rows, summary_rows, bcut_rows = run(args)
    write_csv(args.out / "su2_quaia_direction_stability_fit_grid.csv", fit_rows)
    write_csv(args.out / "su2_quaia_direction_stability_window_summary.csv", summary_rows)
    write_csv(args.out / "su2_quaia_direction_stability_bcut_summary.csv", bcut_rows)
    write_report(args.out / "su2_quaia_direction_stability_report.md", fit_rows, summary_rows, bcut_rows, args)
    print(f"Saved fit grid: {args.out / 'su2_quaia_direction_stability_fit_grid.csv'}", flush=True)
    print(f"Saved window summary: {args.out / 'su2_quaia_direction_stability_window_summary.csv'}", flush=True)
    print(f"Saved bcut summary: {args.out / 'su2_quaia_direction_stability_bcut_summary.csv'}", flush=True)
    print(f"Saved report: {args.out / 'su2_quaia_direction_stability_report.md'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
