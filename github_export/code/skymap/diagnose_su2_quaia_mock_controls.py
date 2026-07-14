#!/usr/bin/env python3
r"""
Mock-control comparison for the SU2-Quaia angular diagnostic scan.

This script compares the object-level full-sample Quaia dipoles from
diagnose_su2_quaia_scan.py against the existing Quaia mock summary products:

    * quaia_G20p0_random_mocks_stats.csv
    * quaia_G20p0_random_mocks_b30_stats.csv
    * quaia_G20p0_random_mocks_zdipole_quad.csv

The available mock products are summary-level dipole tables, not full
redshift-window mock catalogs. Therefore this script can test full-sample
dipole amplitudes and CMB-projected components, but it cannot yet validate the
top redshift-window HHT targets.

Outputs:
    plamb_runs/diagnostics/su2_quaia_mock_controls/su2_quaia_mock_control_summary.csv
    plamb_runs/diagnostics/su2_quaia_mock_controls/su2_quaia_mock_control_report.md
    plamb_runs/diagnostics/su2_quaia_mock_controls/su2_quaia_mock_control_amp.png

Run:
    python diagnose_su2_quaia_mock_controls.py
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover - optional plotting
    plt = None
    MATPLOTLIB_IMPORT_ERROR = exc
else:
    MATPLOTLIB_IMPORT_ERROR = None


ROOT = Path(__file__).resolve().parent
DEFAULT_SCAN = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_scan" / "su2_quaia_scan.csv"
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_mock_controls"
CMB_L_DEG = 264.021
CMB_B_DEG = 48.253


def unit_from_lb(l_deg: float, b_deg: float) -> np.ndarray:
    l = math.radians(l_deg)
    b = math.radians(b_deg)
    cb = math.cos(b)
    return np.array([cb * math.cos(l), cb * math.sin(l), math.sin(b)], dtype=float)


def angular_sep_deg(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 0.0:
        return float("nan")
    dot = float(np.clip(np.dot(a, b) / denom, -1.0, 1.0))
    return math.degrees(math.acos(dot))


def parse_float(value: Any, default: float = float("nan")) -> float:
    if value in ("", None):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def load_observed_full_rows(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in read_csv(path):
        if row.get("tag") != "full":
            continue
        rows.append(
            {
                "bcut_deg": parse_float(row.get("bcut_deg")),
                "N": int(parse_float(row.get("N"), 0.0)),
                "amp": parse_float(row.get("amp")),
                "l_deg": parse_float(row.get("l_deg")),
                "b_deg": parse_float(row.get("b_deg")),
                "sep_cmb_deg": parse_float(row.get("sep_cmb_deg")),
                "amp_par_cmb": parse_float(row.get("amp_par_cmb")),
                "frac_par_cmb": parse_float(row.get("frac_par_cmb")),
                "formal_snr": parse_float(row.get("formal_snr")),
                "max_abs_su2_delta_mu": parse_float(row.get("max_abs_su2_delta_mu")),
            }
        )
    rows.sort(key=lambda row: row["bcut_deg"])
    return rows


def load_mock_rows(path: Path) -> list[dict[str, Any]]:
    cmb = unit_from_lb(CMB_L_DEG, CMB_B_DEG)
    out: list[dict[str, Any]] = []
    for row in read_csv(path):
        amp = parse_float(row.get("dipole_amp"))
        l_deg = parse_float(row.get("l_dip_deg"))
        b_deg = parse_float(row.get("b_dip_deg"))
        if not all(math.isfinite(x) for x in [amp, l_deg, b_deg]):
            continue
        vec = unit_from_lb(l_deg, b_deg)
        sep = angular_sep_deg(vec, cmb)
        cos_sep = float(np.dot(vec, cmb))
        amp_par = amp * cos_sep
        out.append(
            {
                "mock_id": row.get("mock_id", ""),
                "z_global_mean": parse_float(row.get("z_global_mean")),
                "amp": amp,
                "frac": parse_float(row.get("dipole_frac")),
                "l_deg": l_deg,
                "b_deg": b_deg,
                "sep_cmb_deg": sep,
                "amp_par_cmb": amp_par,
                "abs_amp_par_cmb": abs(amp_par),
                "C2_residual": parse_float(row.get("C2_residual")),
                "D2_residual": parse_float(row.get("D2_residual")),
            }
        )
    return out


def p_ge(values: np.ndarray, observed: float) -> float:
    values = values[np.isfinite(values)]
    if len(values) == 0 or not math.isfinite(observed):
        return float("nan")
    return float((1 + np.sum(values >= observed)) / (1 + len(values)))


def p_le(values: np.ndarray, observed: float) -> float:
    values = values[np.isfinite(values)]
    if len(values) == 0 or not math.isfinite(observed):
        return float("nan")
    return float((1 + np.sum(values <= observed)) / (1 + len(values)))


def quantiles(values: np.ndarray) -> dict[str, float]:
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return {"mean": float("nan"), "p50": float("nan"), "p95": float("nan"), "p99": float("nan"), "max": float("nan")}
    return {
        "mean": float(np.mean(values)),
        "p50": float(np.quantile(values, 0.50)),
        "p95": float(np.quantile(values, 0.95)),
        "p99": float(np.quantile(values, 0.99)),
        "max": float(np.max(values)),
    }


def compare(observed_rows: list[dict[str, Any]], mock_name: str, mock_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    amp = np.array([row["amp"] for row in mock_rows], dtype=float)
    abs_par = np.array([row["abs_amp_par_cmb"] for row in mock_rows], dtype=float)
    sep = np.array([row["sep_cmb_deg"] for row in mock_rows], dtype=float)
    amp_q = quantiles(amp)
    par_q = quantiles(abs_par)
    sep_q = quantiles(sep)
    rows: list[dict[str, Any]] = []
    for obs in observed_rows:
        obs_abs_par = abs(float(obs["amp_par_cmb"]))
        p_joint = float(
            (1 + np.sum((amp >= obs["amp"]) & (abs_par >= obs_abs_par)))
            / (1 + len(mock_rows))
        )
        rows.append(
            {
                "mock_set": mock_name,
                "mock_N": len(mock_rows),
                "bcut_deg": obs["bcut_deg"],
                "observed_N": obs["N"],
                "observed_amp": obs["amp"],
                "observed_abs_amp_par_cmb": obs_abs_par,
                "observed_sep_cmb_deg": obs["sep_cmb_deg"],
                "observed_l_deg": obs["l_deg"],
                "observed_b_deg": obs["b_deg"],
                "observed_formal_snr": obs["formal_snr"],
                "mock_amp_mean": amp_q["mean"],
                "mock_amp_p95": amp_q["p95"],
                "mock_amp_p99": amp_q["p99"],
                "mock_amp_max": amp_q["max"],
                "p_amp_ge_observed": p_ge(amp, obs["amp"]),
                "mock_abs_par_mean": par_q["mean"],
                "mock_abs_par_p95": par_q["p95"],
                "mock_abs_par_p99": par_q["p99"],
                "mock_abs_par_max": par_q["max"],
                "p_abs_par_ge_observed": p_ge(abs_par, obs_abs_par),
                "mock_sep_p50": sep_q["p50"],
                "mock_sep_p05": float(np.quantile(sep[np.isfinite(sep)], 0.05)) if np.any(np.isfinite(sep)) else float("nan"),
                "p_sep_le_observed": p_le(sep, obs["sep_cmb_deg"]),
                "p_joint_amp_and_abs_par": p_joint,
                "control_note": "summary-level mock control; no redshift-window mock scan",
            }
        )
    return rows


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


def write_plot(path: Path, observed_rows: list[dict[str, Any]], mock_sets: dict[str, list[dict[str, Any]]]) -> str:
    if plt is None:
        return f"Plot skipped because matplotlib is unavailable: {MATPLOTLIB_IMPORT_ERROR}"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), dpi=150)
    for name, rows in mock_sets.items():
        amp = np.array([row["amp"] for row in rows], dtype=float)
        axes[0].hist(amp, bins=28, alpha=0.35, density=True, label=name)
        abs_par = np.array([row["abs_amp_par_cmb"] for row in rows], dtype=float)
        axes[1].hist(abs_par, bins=28, alpha=0.35, density=True, label=name)
    for obs in observed_rows:
        if float(obs["bcut_deg"]) in (10.0, 25.0, 30.0, 35.0):
            axes[0].axvline(obs["amp"], lw=1.2, ls="--", label=f"obs |b|>{obs['bcut_deg']:g}")
            axes[1].axvline(abs(obs["amp_par_cmb"]), lw=1.2, ls="--")
    axes[0].set_xlabel("Dipole amplitude")
    axes[0].set_ylabel("Mock density")
    axes[0].set_title("Full-sample dipole amplitude")
    axes[1].set_xlabel("|CMB-projected dipole component|")
    axes[1].set_title("CMB-projected component")
    for ax in axes:
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return f"Saved plot: {path}"


def write_report(path: Path, observed_rows: list[dict[str, Any]], comparison_rows: list[dict[str, Any]], args: argparse.Namespace, plot_status: str) -> None:
    focus = [
        row
        for row in comparison_rows
        if row["mock_set"] in {"random_mocks_stats", "random_mocks_zdipole_quad"}
        and float(row["bcut_deg"]) in {25.0, 30.0, 35.0}
    ]
    strongest = sorted(comparison_rows, key=lambda row: float(row["p_amp_ge_observed"]))[:8]
    lines = [
        "# SU2-Quaia Mock-Control Comparison",
        "",
        "## Headline",
        "",
        f"- Observed scan: `{args.scan_csv}`",
        "- Existing mock products are summary-level dipole tables, not full redshift-window mock catalogs.",
        "- This control tests full-sample dipole amplitude and CMB-projected components only.",
        "",
        "## Focus Rows",
        "",
    ]
    lines.extend(
        markdown_table(
            focus,
            [
                "mock_set",
                "bcut_deg",
                "observed_amp",
                "mock_amp_p95",
                "mock_amp_p99",
                "p_amp_ge_observed",
                "observed_abs_amp_par_cmb",
                "mock_abs_par_p95",
                "p_abs_par_ge_observed",
                "p_sep_le_observed",
            ],
        )
    )
    lines.extend(
        [
            "",
            "## Lowest Amplitude p-values",
            "",
        ]
    )
    lines.extend(
        markdown_table(
            strongest,
            [
                "mock_set",
                "bcut_deg",
                "observed_amp",
                "mock_amp_p95",
                "mock_amp_p99",
                "p_amp_ge_observed",
                "p_joint_amp_and_abs_par",
            ],
        )
    )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- Against the 50-row random-mock summary products, the low-latitude full-sample observed dipoles can look large, but that mock set is small and appears to be a legacy summary product.",
            "- Against the 500-row `zdipole_quad` mock product, the observed full-sample amplitudes are not unusually large.",
            "- The current top SU2-Quaia redshift windows cannot be validated by these summary mocks; they require full mock catalogs or rerunning the old mock-generation pipeline per redshift bin.",
            "- Treat this as a cautionary control: the full-sample angular signal is not yet robust evidence for SU2 physics.",
            "",
            "## Local Outputs",
            "",
            f"- Summary CSV: `{OUTDIR / 'su2_quaia_mock_control_summary.csv'}`",
            f"- Report: `{OUTDIR / 'su2_quaia_mock_control_report.md'}`",
            f"- Plot: `{OUTDIR / 'su2_quaia_mock_control_amp.png'}`",
            f"- Plot status: {plot_status}",
            "",
            "## Next Step",
            "",
            "Regenerate Quaia mocks at the object/redshift-window level, or recover the full mock catalog generation script, then rerun `diagnose_su2_quaia_scan.py` on each mock to produce p-values for the top redshift windows.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-csv", type=Path, default=DEFAULT_SCAN)
    parser.add_argument("--out", type=Path, default=OUTDIR)
    parser.add_argument(
        "--mock-files",
        type=Path,
        nargs="+",
        default=[
            ROOT / "quaia_G20p0_random_mocks_stats.csv",
            ROOT / "quaia_G20p0_random_mocks_b30_stats.csv",
            ROOT / "quaia_G20p0_random_mocks_zdipole_quad.csv",
        ],
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if not args.scan_csv.exists():
        print(f"Missing scan CSV: {args.scan_csv}", file=sys.stderr)
        return 2
    observed_rows = load_observed_full_rows(args.scan_csv)
    if not observed_rows:
        print("No full-sample observed rows found.", file=sys.stderr)
        return 2

    mock_sets: dict[str, list[dict[str, Any]]] = {}
    for path in args.mock_files:
        if not path.exists():
            print(f"Warning: missing mock file {path}", file=sys.stderr)
            continue
        name = path.stem.replace("quaia_G20p0_", "")
        mock_sets[name] = load_mock_rows(path)
        print(f"Loaded {len(mock_sets[name])} mock rows from {path}")
    if not mock_sets:
        print("No mock rows loaded.", file=sys.stderr)
        return 2

    comparison_rows: list[dict[str, Any]] = []
    for name, rows in mock_sets.items():
        comparison_rows.extend(compare(observed_rows, name, rows))

    summary_csv = args.out / "su2_quaia_mock_control_summary.csv"
    report_md = args.out / "su2_quaia_mock_control_report.md"
    plot_png = args.out / "su2_quaia_mock_control_amp.png"
    write_csv(summary_csv, comparison_rows)
    plot_status = write_plot(plot_png, observed_rows, mock_sets)
    write_report(report_md, observed_rows, comparison_rows, args, plot_status)
    print(f"Saved summary: {summary_csv}")
    print(f"Saved report: {report_md}")
    print(plot_status)
    best = min(comparison_rows, key=lambda row: float(row["p_amp_ge_observed"]))
    print(
        "Smallest full-sample amplitude p-value: "
        f"mock={best['mock_set']} bcut={best['bcut_deg']} "
        f"p={best['p_amp_ge_observed']:.4g} obs_amp={best['observed_amp']:.4g}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
