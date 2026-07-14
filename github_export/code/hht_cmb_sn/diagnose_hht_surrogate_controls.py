r"""
Surrogate controls for HHT residual diagnostics.

This script tests whether an HHT IMF result is stronger than shuffled-order
controls. It preserves the x grid and the one-point distribution of the chosen
y column, randomly permutes y, runs the same EMD/HHT summary, and compares the
original IMF energy metrics to the surrogate distribution.

Outputs:
    plamb_runs/diagnostics/hht_resonance/hht_surrogate_<label>_summary.csv
    plamb_runs/diagnostics/hht_resonance/hht_surrogate_<label>_summary.md
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Any

import numpy as np

import diagnose_hht_resonance as hht


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "hht_resonance"


def safe_slug(value: str) -> str:
    return hht.safe_slug(value)


def read_series(path: Path, x_col: str, y_col: str) -> tuple[np.ndarray, np.ndarray]:
    cols = hht.read_csv_columns(path)
    if x_col not in cols:
        raise ValueError(f"x column {x_col!r} not found")
    if y_col not in cols:
        raise ValueError(f"y column {y_col!r} not found")
    x = hht.as_float_array(cols[x_col], x_col)
    y = hht.as_float_array(cols[y_col], y_col)
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    unique = np.concatenate(([True], np.diff(x) > 0))
    return x[unique], y[unique]


def imf_metrics(x: np.ndarray, y: np.ndarray, max_imfs: int, sift_max_iter: int, sift_sd: float) -> dict[str, Any]:
    imfs, residue = hht.emd(y, max_imfs=max_imfs, sift_max_iter=sift_max_iter, sift_sd=sift_sd)
    rows = hht.summarize_imfs("series", x, y, imfs, residue)
    imf_rows = [r for r in rows if str(r["component"]).startswith("IMF")]
    if not imf_rows:
        return {
            "n_imfs": 0,
            "max_imf_energy_fraction": 0.0,
            "sum_imf_energy_fraction": 0.0,
            "dominant_imf_frequency": math.nan,
            "dominant_imf_index": math.nan,
        }
    dominant = max(imf_rows, key=lambda r: float(r["energy_fraction"]))
    return {
        "n_imfs": len(imf_rows),
        "max_imf_energy_fraction": float(dominant["energy_fraction"]),
        "sum_imf_energy_fraction": float(sum(float(r["energy_fraction"]) for r in imf_rows)),
        "dominant_imf_frequency": float(dominant["amplitude_weighted_frequency"]),
        "dominant_imf_index": int(dominant["component_index"]),
    }


def percentile_rank(value: float, samples: np.ndarray) -> float:
    if samples.size == 0:
        return math.nan
    return float(np.mean(samples <= value))


def upper_tail_p(value: float, samples: np.ndarray) -> float:
    if samples.size == 0:
        return math.nan
    return float((np.count_nonzero(samples >= value) + 1) / (samples.size + 1))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, args, original: dict[str, Any], summary: dict[str, Any]) -> None:
    lines = [
        "# HHT surrogate control",
        "",
        f"Input: `{args.input_csv}`",
        f"x column: `{args.x_col}`",
        f"y column: `{args.y_col}`",
        f"surrogates: {args.n_surrogates}",
        "",
        "## Original metrics",
        "",
        "| metric | value |",
        "|---|---:|",
        f"| n IMFs | {original['n_imfs']} |",
        f"| max IMF energy fraction | {original['max_imf_energy_fraction']:.6g} |",
        f"| sum IMF energy fraction | {original['sum_imf_energy_fraction']:.6g} |",
        f"| dominant IMF frequency | {original['dominant_imf_frequency']:.6g} |",
        "",
        "## Shuffled-control comparison",
        "",
        "| metric | original | surrogate median | surrogate p90 | upper-tail p | percentile |",
        "|---|---:|---:|---:|---:|---:|",
        (
            f"| max IMF energy fraction | {original['max_imf_energy_fraction']:.6g} | "
            f"{summary['max_energy_median']:.6g} | {summary['max_energy_p90']:.6g} | "
            f"{summary['max_energy_upper_tail_p']:.6g} | {summary['max_energy_percentile']:.6g} |"
        ),
        (
            f"| sum IMF energy fraction | {original['sum_imf_energy_fraction']:.6g} | "
            f"{summary['sum_energy_median']:.6g} | {summary['sum_energy_p90']:.6g} | "
            f"{summary['sum_energy_upper_tail_p']:.6g} | {summary['sum_energy_percentile']:.6g} |"
        ),
        "",
        "## Interpretation",
        "",
    ]
    if summary["max_energy_upper_tail_p"] <= 0.05 or summary["sum_energy_upper_tail_p"] <= 0.05:
        lines.append("The original IMF energy is unusually high relative to shuffled-order controls. Treat this as a candidate feature requiring probe-by-probe checks.")
    else:
        lines.append("The original IMF energy is not unusual relative to shuffled-order controls. Treat the HHT IMFs as likely sampling/probe-mixing structure rather than resonance evidence.")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run shuffled-order surrogate controls for an HHT CSV series.")
    parser.add_argument("--input-csv", type=Path, required=True)
    parser.add_argument("--x-col", required=True)
    parser.add_argument("--y-col", required=True)
    parser.add_argument("--label", default="")
    parser.add_argument("--n-surrogates", type=int, default=100)
    parser.add_argument("--seed", type=int, default=260520333)
    parser.add_argument("--max-imfs", type=int, default=6)
    parser.add_argument("--sift-max-iter", type=int, default=80)
    parser.add_argument("--sift-sd", type=float, default=0.08)
    args = parser.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    x, y = read_series(args.input_csv, args.x_col, args.y_col)
    original = imf_metrics(x, y, args.max_imfs, args.sift_max_iter, args.sift_sd)

    rng = np.random.default_rng(args.seed)
    rows: list[dict[str, Any]] = []
    for i in range(args.n_surrogates):
        y_perm = rng.permutation(y)
        metrics = imf_metrics(x, y_perm, args.max_imfs, args.sift_max_iter, args.sift_sd)
        row = {"kind": "surrogate", "index": i}
        row.update(metrics)
        rows.append(row)

    max_energy = np.asarray([float(r["max_imf_energy_fraction"]) for r in rows], dtype=float)
    sum_energy = np.asarray([float(r["sum_imf_energy_fraction"]) for r in rows], dtype=float)
    summary = {
        "max_energy_median": float(np.median(max_energy)),
        "max_energy_p90": float(np.percentile(max_energy, 90)),
        "max_energy_upper_tail_p": upper_tail_p(float(original["max_imf_energy_fraction"]), max_energy),
        "max_energy_percentile": percentile_rank(float(original["max_imf_energy_fraction"]), max_energy),
        "sum_energy_median": float(np.median(sum_energy)),
        "sum_energy_p90": float(np.percentile(sum_energy, 90)),
        "sum_energy_upper_tail_p": upper_tail_p(float(original["sum_imf_energy_fraction"]), sum_energy),
        "sum_energy_percentile": percentile_rank(float(original["sum_imf_energy_fraction"]), sum_energy),
    }

    label = safe_slug(args.label or args.y_col)
    summary_row = {"kind": "original", "index": -1}
    summary_row.update(original)
    all_rows = [summary_row] + rows
    csv_path = OUTDIR / f"hht_surrogate_{label}_summary.csv"
    md_path = OUTDIR / f"hht_surrogate_{label}_summary.md"
    write_csv(csv_path, all_rows)
    write_markdown(md_path, args, original, summary)

    print(f"Saved surrogate CSV: {csv_path}")
    print(f"Saved surrogate summary: {md_path}")
    print(
        "Original max IMF energy="
        f"{float(original['max_imf_energy_fraction']):.4g}, "
        f"shuffle p={summary['max_energy_upper_tail_p']:.4g}; "
        "sum IMF energy="
        f"{float(original['sum_imf_energy_fraction']):.4g}, "
        f"shuffle p={summary['sum_energy_upper_tail_p']:.4g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
