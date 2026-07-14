#!/usr/bin/env python3
r"""
High-leverage jackknife for SPARC PLAMB/RAR comparison.

Uses the high-leverage galaxies identified in the SPARC distance-method audit,
removes the top N by absolute PLAMB-vs-RAR chi2 difference, and recomputes
fixed-exponent PLAMB vs RAR on comparable splits.

Run:
    python jackknife_sparc_high_leverage.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np

try:
    from scipy.optimize import minimize
except Exception as exc:  # pragma: no cover
    minimize = None
    SCIPY_IMPORT_ERROR = exc
else:
    SCIPY_IMPORT_ERROR = None

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover
    plt = None
    MATPLOTLIB_IMPORT_ERROR = exc
else:
    MATPLOTLIB_IMPORT_ERROR = None

from diagnose_plamb_sparc_fixed_exponent_grid import (
    chi2,
    predicted_velocity_fixed_plamb,
    predicted_velocity_rar,
    subset_dataset,
)
from diagnose_plamb_sparc_h0_nuisance import rescale_for_h0
from diagnose_plamb_sparc_rotation import DEFAULT_POINTS, DEFAULT_SAMPLE, Dataset, load_dataset


ROOT = Path(__file__).resolve().parent
LEVERAGE = ROOT / "plamb_runs" / "diagnostics" / "sparc_distance_method_audit" / "sparc_distance_method_leverage_galaxies.csv"
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "sparc_high_leverage_jackknife"


SPLIT_CONFIGS = [
    {
        "split": "all_Q2",
        "quality_max": 2,
        "distance_method": "all",
        "plamb_h0": 65.0,
        "plamb_n": 0.425,
        "plamb_s": 1.25,
        "rar_h0": 65.0,
    },
    {
        "split": "all_Q1",
        "quality_max": 1,
        "distance_method": "all",
        "plamb_h0": 70.0,
        "plamb_n": 0.425,
        "plamb_s": 1.5,
        "rar_h0": 62.5,
    },
    {
        "split": "hubble_flow_Q2",
        "quality_max": 2,
        "distance_method": "hubble_flow",
        "plamb_h0": 70.0,
        "plamb_n": 0.4,
        "plamb_s": 1.25,
        "rar_h0": 62.5,
    },
    {
        "split": "non_hubble_Q2",
        "quality_max": 2,
        "distance_method": "non_hubble",
        "plamb_h0": 70.0,
        "plamb_n": 0.5,
        "plamb_s": 1.25,
        "rar_h0": 62.5,
    },
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_int_list(value: str) -> list[int]:
    vals = [int(float(token.strip())) for token in value.split(",") if token.strip()]
    if not vals:
        raise ValueError("No removal stages supplied.")
    return vals


def top_leverage_galaxies(path: Path) -> list[str]:
    rows = read_csv(path)
    rows.sort(key=lambda row: abs(float(row.get("delta_chi2_PLAMB_minus_RAR", "0"))), reverse=True)
    return [row["Galaxy"] for row in rows]


def remove_galaxies(data: Dataset, removed: set[str]) -> Dataset:
    mask = np.array([str(galaxy) not in removed for galaxy in data.galaxy], dtype=bool)
    return subset_dataset(data, mask)


def fit_plamb_ml(data: Dataset, h0: float, inertia_index: float, screen_index: float) -> tuple[np.ndarray, float]:
    if minimize is None:
        raise RuntimeError(f"SciPy optimizer import failed: {SCIPY_IMPORT_ERROR}")

    def objective(theta: np.ndarray) -> float:
        ydisk, ybul = theta
        pred = predicted_velocity_fixed_plamb(data, ydisk, ybul, h0, inertia_index, screen_index)
        return chi2(pred, data)

    best_x = None
    best_fun = 1.0e100
    for start in ([0.5, 0.7], [0.35, 0.6], [0.8, 0.9], [1.0, 1.0]):
        result = minimize(
            objective,
            np.asarray(start, dtype=float),
            method="L-BFGS-B",
            bounds=[(0.05, 1.25), (0.05, 1.60)],
            options={"maxiter": 800, "ftol": 1e-10},
        )
        if result.success and float(result.fun) < best_fun:
            best_x = np.asarray(result.x, dtype=float)
            best_fun = float(result.fun)
    if best_x is None:
        raise RuntimeError("PLAMB jackknife fit failed.")
    return best_x, best_fun


def fit_rar(data: Dataset) -> tuple[np.ndarray, float]:
    if minimize is None:
        raise RuntimeError(f"SciPy optimizer import failed: {SCIPY_IMPORT_ERROR}")

    def objective(theta: np.ndarray) -> float:
        ydisk, ybul, log10_gdagger = theta
        pred = predicted_velocity_rar(data, ydisk, ybul, log10_gdagger)
        return chi2(pred, data)

    best_x = None
    best_fun = 1.0e100
    for start in ([0.5, 0.7, -10.0], [0.35, 0.6, -10.1], [0.75, 0.9, -9.9]):
        result = minimize(
            objective,
            np.asarray(start, dtype=float),
            method="L-BFGS-B",
            bounds=[(0.05, 1.25), (0.05, 1.60), (-11.6, -9.2)],
            options={"maxiter": 1000, "ftol": 1e-10},
        )
        if result.success and float(result.fun) < best_fun:
            best_x = np.asarray(result.x, dtype=float)
            best_fun = float(result.fun)
    if best_x is None:
        raise RuntimeError("RAR jackknife fit failed.")
    return best_x, best_fun


def score(n: int, k: int, chi2_value: float) -> tuple[float, float, float]:
    aic = chi2_value + 2.0 * k
    bic = chi2_value + k * math.log(max(n, 1))
    return chi2_value / max(n - k, 1), aic, bic


def make_plot(rows: list[dict[str, Any]], path: Path) -> None:
    if plt is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 5.7))
    for split in sorted({row["split"] for row in rows}):
        split_rows = sorted([row for row in rows if row["split"] == split], key=lambda row: int(row["removed_count"]))
        ax.plot(
            [int(row["removed_count"]) for row in split_rows],
            [float(row["delta_AIC_PLAMB_minus_RAR"]) for row in split_rows],
            marker="o",
            linewidth=1.8,
            label=split,
        )
    ax.axhline(0.0, color="#555555", linewidth=1.0)
    ax.set_xlabel("Top high-leverage galaxies removed")
    ax.set_ylabel("Delta AIC: PLAMB fixed exponent minus RAR")
    ax.set_title("SPARC High-Leverage Jackknife")
    ax.grid(True, alpha=0.22)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def write_report(rows: list[dict[str, Any]], stages: list[int], leverage: list[str]) -> None:
    lines = [
        "# SPARC High-Leverage Jackknife",
        "",
        "## Executive Summary",
        "",
        "This removes the top high-leverage galaxies from the distance-method audit and refits fixed-exponent PLAMB and RAR within each comparable split.",
        "",
        f"- Removal stages: {', '.join(str(stage) for stage in stages)}",
        f"- Top removed order starts: {', '.join(leverage[:10])}",
        "",
        "| Split | Removed | N galaxies | N points | Delta AIC PLAMB-RAR | PLAMB chi2/dof | RAR chi2/dof |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['split']} | {row['removed_count']} | {row['N_galaxies']} | {row['N_points']} | "
            f"{row['delta_AIC_PLAMB_minus_RAR']:.3f} | {row['plamb_chi2_dof']:.3f} | {row['rar_chi2_dof']:.3f} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Guide",
            "",
            "- Negative Delta AIC means fixed-exponent PLAMB beats RAR.",
            "- Positive Delta AIC means RAR beats fixed-exponent PLAMB.",
            "- If the sign flips after removing a small number of galaxies, the preference is high-leverage sensitive.",
            "- If the sign remains stable across removal stages, the preference is more robust.",
        ]
    )
    (OUTDIR / "sparc_high_leverage_jackknife.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (OUTDIR / "sparc_high_leverage_jackknife.json").write_text(
        json.dumps({"stages": stages, "leverage_order": leverage}, indent=2),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="SPARC high-leverage jackknife.")
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--points", type=Path, default=DEFAULT_POINTS)
    parser.add_argument("--leverage", type=Path, default=LEVERAGE)
    parser.add_argument("--remove-counts", default="0,5,10,15,20")
    parser.add_argument("--err-floor-kms", type=float, default=3.0)
    args = parser.parse_args()

    stages = parse_int_list(args.remove_counts)
    leverage_order = top_leverage_galaxies(args.leverage)
    rows: list[dict[str, Any]] = []
    for stage in stages:
        removed = set(leverage_order[:stage])
        for config in SPLIT_CONFIGS:
            base = load_dataset(
                args.sample,
                args.points,
                quality_max=int(config["quality_max"]),
                distance_method=str(config["distance_method"]),
                err_floor=args.err_floor_kms,
            )
            base = remove_galaxies(base, removed)
            plamb_data = rescale_for_h0(base, float(config["plamb_h0"]), 73.0, "hubble_flow_only")
            rar_data = rescale_for_h0(base, float(config["rar_h0"]), 73.0, "hubble_flow_only")

            plamb_theta, plamb_chi2 = fit_plamb_ml(
                plamb_data,
                float(config["plamb_h0"]),
                float(config["plamb_n"]),
                float(config["plamb_s"]),
            )
            rar_theta, rar_chi2 = fit_rar(rar_data)
            n = len(base.vobs)
            plamb_chi2_dof, plamb_aic, plamb_bic = score(n, 2, plamb_chi2)
            rar_chi2_dof, rar_aic, rar_bic = score(n, 3, rar_chi2)
            row = {
                "split": config["split"],
                "removed_count": stage,
                "removed_galaxies": ";".join(leverage_order[:stage]),
                "N_points": n,
                "N_galaxies": len(set(base.galaxy.tolist())),
                "plamb_H0": config["plamb_h0"],
                "plamb_n": config["plamb_n"],
                "plamb_s": config["plamb_s"],
                "plamb_ydisk": plamb_theta[0],
                "plamb_ybul": plamb_theta[1],
                "plamb_chi2": plamb_chi2,
                "plamb_chi2_dof": plamb_chi2_dof,
                "plamb_AIC": plamb_aic,
                "plamb_BIC": plamb_bic,
                "rar_H0": config["rar_h0"],
                "rar_ydisk": rar_theta[0],
                "rar_ybul": rar_theta[1],
                "rar_log10_gdagger": rar_theta[2],
                "rar_chi2": rar_chi2,
                "rar_chi2_dof": rar_chi2_dof,
                "rar_AIC": rar_aic,
                "rar_BIC": rar_bic,
                "delta_AIC_PLAMB_minus_RAR": plamb_aic - rar_aic,
                "delta_BIC_PLAMB_minus_RAR": plamb_bic - rar_bic,
            }
            rows.append(row)
            print(
                f"remove={stage:2d} {config['split']:16s} "
                f"dAIC={row['delta_AIC_PLAMB_minus_RAR']:9.3f} Ngal={row['N_galaxies']}"
            )

    fields = [
        "split",
        "removed_count",
        "removed_galaxies",
        "N_points",
        "N_galaxies",
        "plamb_H0",
        "plamb_n",
        "plamb_s",
        "plamb_ydisk",
        "plamb_ybul",
        "plamb_chi2",
        "plamb_chi2_dof",
        "plamb_AIC",
        "plamb_BIC",
        "rar_H0",
        "rar_ydisk",
        "rar_ybul",
        "rar_log10_gdagger",
        "rar_chi2",
        "rar_chi2_dof",
        "rar_AIC",
        "rar_BIC",
        "delta_AIC_PLAMB_minus_RAR",
        "delta_BIC_PLAMB_minus_RAR",
    ]
    OUTDIR.mkdir(parents=True, exist_ok=True)
    write_csv(OUTDIR / "sparc_high_leverage_jackknife.csv", rows, fields)
    make_plot(rows, OUTDIR / "sparc_high_leverage_jackknife_delta_aic.png")
    write_report(rows, stages, leverage_order)
    print(f"Saved report: {OUTDIR / 'sparc_high_leverage_jackknife.md'}")


if __name__ == "__main__":
    main()
