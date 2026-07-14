#!/usr/bin/env python3
r"""
Shuffle-control HHT checks for the SU2-Quaia diagnostic series.

This script reuses the HHT/EMD helpers from diagnose_hht_resonance.py and
compares the observed dominant IMF against redshift-shuffled surrogates. The
control preserves the set of Quaia window amplitudes but destroys their ordered
relationship with redshift.

Outputs:
    plamb_runs/diagnostics/su2_quaia_hht_controls/su2_quaia_hht_control_summary.csv
    plamb_runs/diagnostics/su2_quaia_hht_controls/su2_quaia_hht_control_draws.csv
    plamb_runs/diagnostics/su2_quaia_hht_controls/su2_quaia_hht_control_report.md

Run:
    python diagnose_su2_quaia_hht_controls.py
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

from diagnose_hht_resonance import decompose_series


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_scan" / "su2_quaia_hht_series.csv"
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_hht_controls"


def read_csv_columns(path: Path) -> dict[str, list[str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        columns = {name: [] for name in (reader.fieldnames or [])}
        for row in reader:
            for name in columns:
                columns[name].append(row.get(name, ""))
    return columns


def as_float(values: list[str], label: str) -> np.ndarray:
    out: list[float] = []
    for value in values:
        try:
            out.append(float(value))
        except ValueError as exc:
            raise ValueError(f"Column {label!r} contains nonnumeric value {value!r}") from exc
    return np.asarray(out, dtype=float)


def imf_rows_only(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [row for row in rows if str(row.get("component", "")).startswith("IMF")]


def top_imf_metrics(t: np.ndarray, y: np.ndarray, label: str, args: argparse.Namespace) -> dict[str, float | str]:
    _, _, rows = decompose_series(
        t,
        {label: y},
        max_imfs=args.max_imfs,
        sift_max_iter=args.sift_max_iter,
        sift_sd=args.sift_sd,
    )
    imfs = imf_rows_only(rows)
    if not imfs:
        return {
            "component": "",
            "energy_fraction": float("nan"),
            "weighted_frequency": float("nan"),
            "median_positive_frequency": float("nan"),
            "zero_crossings": float("nan"),
        }
    top = max(imfs, key=lambda row: float(row.get("energy_fraction", float("nan"))))
    return {
        "component": str(top.get("component", "")),
        "energy_fraction": float(top.get("energy_fraction", float("nan"))),
        "weighted_frequency": float(top.get("amplitude_weighted_frequency", float("nan"))),
        "median_positive_frequency": float(top.get("median_positive_frequency", float("nan"))),
        "zero_crossings": float(top.get("zero_crossings", float("nan"))),
    }


def run_column(t: np.ndarray, y: np.ndarray, column: str, args: argparse.Namespace) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rng = np.random.default_rng(args.seed + abs(hash(column)) % 1_000_000)
    observed = top_imf_metrics(t, y, column, args)
    observed_energy = float(observed["energy_fraction"])
    observed_freq = float(observed["weighted_frequency"])

    draws: list[dict[str, Any]] = []
    shuffled_energy: list[float] = []
    nearest_band_energy: list[float] = []
    for draw in range(args.shuffles):
        y_shuffled = np.array(y, copy=True)
        rng.shuffle(y_shuffled)
        metrics = top_imf_metrics(t, y_shuffled, f"{column}_shuffle", args)
        energy = float(metrics["energy_fraction"])
        freq = float(metrics["weighted_frequency"])
        shuffled_energy.append(energy)
        in_band = (
            math.isfinite(freq)
            and math.isfinite(observed_freq)
            and freq > 0.0
            and observed_freq > 0.0
            and abs(math.log(freq / observed_freq)) <= math.log(args.freq_ratio)
        )
        nearest_band_energy.append(energy if in_band else 0.0)
        draws.append(
            {
                "column": column,
                "draw": draw,
                "shuffle_top_energy_fraction": energy,
                "shuffle_top_weighted_frequency": freq,
                "shuffle_top_component": metrics["component"],
                "within_observed_frequency_band": int(in_band),
            }
        )

    shuffled_arr = np.asarray(shuffled_energy, dtype=float)
    band_arr = np.asarray(nearest_band_energy, dtype=float)
    p_max = float((1 + np.sum(shuffled_arr >= observed_energy)) / (1 + len(shuffled_arr)))
    p_band = float((1 + np.sum(band_arr >= observed_energy)) / (1 + len(band_arr)))
    summary: dict[str, Any] = {
        "column": column,
        "N": len(y),
        "shuffles": args.shuffles,
        "observed_component": observed["component"],
        "observed_top_energy_fraction": observed_energy,
        "observed_weighted_frequency": observed_freq,
        "observed_median_positive_frequency": observed["median_positive_frequency"],
        "observed_zero_crossings": observed["zero_crossings"],
        "shuffle_mean_top_energy_fraction": float(np.mean(shuffled_arr)),
        "shuffle_median_top_energy_fraction": float(np.median(shuffled_arr)),
        "shuffle_p95_top_energy_fraction": float(np.quantile(shuffled_arr, 0.95)),
        "shuffle_p99_top_energy_fraction": float(np.quantile(shuffled_arr, 0.99)),
        "p_value_top_energy": p_max,
        "p_value_same_frequency_band": p_band,
        "frequency_band_ratio": args.freq_ratio,
    }
    return summary, draws


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, Any]], columns: list[str]) -> list[str]:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        values: list[str] = []
        for col in columns:
            value = row.get(col, "")
            if isinstance(value, float):
                values.append(f"{value:.6g}" if math.isfinite(value) else "")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return lines


def write_report(path: Path, summaries: list[dict[str, Any]], args: argparse.Namespace) -> None:
    lines = [
        "# SU2-Quaia HHT Shuffle Controls",
        "",
        "## Headline",
        "",
        f"- Input series: `{args.input_csv}`",
        f"- Shuffles per column: `{args.shuffles}`",
        "- Control: shuffle y-values against redshift while preserving the same values and same x grid.",
        "- Low p-values would indicate ordered redshift structure beyond this simple shuffle control.",
        "",
        "## Results",
        "",
    ]
    lines.extend(
        markdown_table(
            summaries,
            [
                "column",
                "N",
                "observed_component",
                "observed_top_energy_fraction",
                "observed_weighted_frequency",
                "shuffle_p95_top_energy_fraction",
                "p_value_top_energy",
                "p_value_same_frequency_band",
            ],
        )
    )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- `p_value_top_energy` compares the observed dominant IMF energy fraction against the strongest IMF in each shuffled series.",
            "- `p_value_same_frequency_band` additionally asks whether shuffled series produce a comparably strong IMF within a factor of the observed frequency.",
            "- This is a first-pass null test. It does not replace Quaia selection-function mocks, sky masks, or survey-footprint controls.",
            "",
            "## Next Step",
            "",
            "Use the existing Quaia random/mock products to repeat the full dipole scan, not just shuffle the derived series. That will test mask and sky-selection effects directly.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out", type=Path, default=OUTDIR)
    parser.add_argument("--x-col", default="x")
    parser.add_argument("--columns", nargs="+", default=["amp_par_cmb", "frac_par_cmb", "dipole_amp"])
    parser.add_argument("--shuffles", type=int, default=500)
    parser.add_argument("--seed", type=int, default=271828)
    parser.add_argument("--freq-ratio", type=float, default=1.35)
    parser.add_argument("--max-imfs", type=int, default=6)
    parser.add_argument("--sift-max-iter", type=int, default=80)
    parser.add_argument("--sift-sd", type=float, default=0.08)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if not args.input_csv.exists():
        print(f"Missing input CSV: {args.input_csv}", file=sys.stderr)
        return 2
    columns = read_csv_columns(args.input_csv)
    if args.x_col not in columns:
        print(f"Missing x column: {args.x_col}", file=sys.stderr)
        return 2
    t = as_float(columns[args.x_col], args.x_col)
    order = np.argsort(t)
    t = t[order]

    summaries: list[dict[str, Any]] = []
    all_draws: list[dict[str, Any]] = []
    for column in args.columns:
        if column not in columns:
            print(f"Missing y column: {column}", file=sys.stderr)
            return 2
        y = as_float(columns[column], column)[order]
        print(f"[controls] {column}: N={len(y)} shuffles={args.shuffles}", flush=True)
        summary, draws = run_column(t, y, column, args)
        summaries.append(summary)
        all_draws.extend(draws)

    summary_csv = args.out / "su2_quaia_hht_control_summary.csv"
    draws_csv = args.out / "su2_quaia_hht_control_draws.csv"
    report_md = args.out / "su2_quaia_hht_control_report.md"
    write_csv(summary_csv, summaries)
    write_csv(draws_csv, all_draws)
    write_report(report_md, summaries, args)
    print(f"Saved summary: {summary_csv}")
    print(f"Saved draws: {draws_csv}")
    print(f"Saved report: {report_md}")
    best = min(summaries, key=lambda row: float(row["p_value_top_energy"]))
    print(
        "Most unusual column by shuffle control: "
        f"{best['column']} p={best['p_value_top_energy']:.4g} "
        f"obs_energy={best['observed_top_energy_fraction']:.4g}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
