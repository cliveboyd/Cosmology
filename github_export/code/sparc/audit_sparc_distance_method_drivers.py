#!/usr/bin/env python3
r"""
Audit why SPARC distance-method splits change PLAMB/RAR preference.

This script compares SPARC Hubble-flow-distance galaxies (f_D=1) with
non-Hubble-distance galaxies and summarizes whether model preference changes
are plausibly driven by sample composition, distance scale, quality, or a few
high-leverage galaxies.

Inputs:
    external_datasets/current_cosmology_sources/SPARC/sparc_galaxy_sample.csv
    external_datasets/current_cosmology_sources/SPARC/sparc_rotation_points.csv
    plamb_runs/diagnostics/sparc_rotation_fixed_exponent_grid/sparc_fixed_exponent_grid.csv

Run:
    python audit_sparc_distance_method_drivers.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

try:
    from scipy.stats import mannwhitneyu, spearmanr
except Exception as exc:  # pragma: no cover - reported at runtime
    mannwhitneyu = None
    spearmanr = None
    SCIPY_IMPORT_ERROR = exc
else:
    SCIPY_IMPORT_ERROR = None

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover - plotting is optional
    plt = None
    MATPLOTLIB_IMPORT_ERROR = exc
else:
    MATPLOTLIB_IMPORT_ERROR = None

from diagnose_plamb_sparc_h0_nuisance import rescale_for_h0
from diagnose_plamb_sparc_fixed_exponent_grid import (
    predicted_velocity_fixed_plamb,
    predicted_velocity_rar,
    subset_dataset,
)
from diagnose_plamb_sparc_rotation import DEFAULT_POINTS, DEFAULT_SAMPLE, Dataset, load_dataset


ROOT = Path(__file__).resolve().parent
GRID_PATH = ROOT / "plamb_runs" / "diagnostics" / "sparc_rotation_fixed_exponent_grid" / "sparc_fixed_exponent_grid.csv"
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "sparc_distance_method_audit"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_float(value: Any, default: float = float("nan")) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def as_int(value: Any, default: int = -1) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def median_or_nan(values: list[float] | np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return float("nan")
    return float(np.median(arr))


def mean_or_nan(values: list[float] | np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return float("nan")
    return float(np.mean(arr))


def mann_p(a: np.ndarray, b: np.ndarray) -> float:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if mannwhitneyu is None or len(a) < 3 or len(b) < 3:
        return float("nan")
    result = mannwhitneyu(a, b, alternative="two-sided")
    return float(result.pvalue)


def load_sample_rows(sample_path: Path) -> dict[str, dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    for row in read_csv(sample_path):
        galaxy = row["Galaxy"]
        rows[galaxy] = {
            "Galaxy": galaxy,
            "T": as_int(row.get("T")),
            "D_Mpc": as_float(row.get("D_Mpc")),
            "e_D_Mpc": as_float(row.get("e_D_Mpc")),
            "f_D": as_int(row.get("f_D")),
            "Inc_deg": as_float(row.get("Inc_deg")),
            "e_Inc_deg": as_float(row.get("e_Inc_deg")),
            "L36_1e9_Lsun": as_float(row.get("L36_1e9_Lsun")),
            "Reff_kpc": as_float(row.get("Reff_kpc")),
            "SBeff_Lsun_pc2": as_float(row.get("SBeff_Lsun_pc2")),
            "Rdisk_kpc": as_float(row.get("Rdisk_kpc")),
            "SBdisk_Lsun_pc2": as_float(row.get("SBdisk_Lsun_pc2")),
            "MHI_1e9_Msun": as_float(row.get("MHI_1e9_Msun")),
            "RHI_kpc": as_float(row.get("RHI_kpc")),
            "Vflat_kms": as_float(row.get("Vflat_kms")),
            "Q": as_int(row.get("Q")),
        }
    return rows


def group_label(f_d: int) -> str:
    if f_d == 1:
        return "hubble_flow"
    if f_d in {2, 3, 4, 5}:
        return "non_hubble"
    return "other"


def sample_property_rows(sample: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    properties = [
        "D_Mpc",
        "e_D_Mpc",
        "Inc_deg",
        "e_Inc_deg",
        "T",
        "L36_1e9_Lsun",
        "Reff_kpc",
        "SBeff_Lsun_pc2",
        "Rdisk_kpc",
        "SBdisk_Lsun_pc2",
        "MHI_1e9_Msun",
        "RHI_kpc",
        "Vflat_kms",
        "Q",
    ]
    rows = []
    hubble = [row for row in sample.values() if row["f_D"] == 1]
    non_hubble = [row for row in sample.values() if row["f_D"] != 1]
    for prop in properties:
        h_vals = np.asarray([row[prop] for row in hubble], dtype=float)
        nh_vals = np.asarray([row[prop] for row in non_hubble], dtype=float)
        rows.append(
            {
                "property": prop,
                "hubble_N": int(np.sum(np.isfinite(h_vals))),
                "non_hubble_N": int(np.sum(np.isfinite(nh_vals))),
                "hubble_median": median_or_nan(h_vals),
                "non_hubble_median": median_or_nan(nh_vals),
                "hubble_mean": mean_or_nan(h_vals),
                "non_hubble_mean": mean_or_nan(nh_vals),
                "mann_whitney_p": mann_p(h_vals, nh_vals),
            }
        )
    return rows


def best_grid_rows(grid_path: Path) -> list[dict[str, Any]]:
    rows = read_csv(grid_path)
    keep: list[dict[str, Any]] = []
    for split in sorted({row["split"] for row in rows if not row["model"].endswith("_HOLDOUT")}):
        split_rows = [row for row in rows if row["split"] == split and not row["model"].endswith("_HOLDOUT")]
        plamb = [row for row in split_rows if row["model"] == "PLAMB_FIXED_EXPONENT"]
        rar = [row for row in split_rows if row["model"] == "RAR_FREE_GDAGGER"]
        if not plamb or not rar:
            continue
        best_plamb = min(plamb, key=lambda row: as_float(row.get("AIC")))
        best_rar = min(rar, key=lambda row: as_float(row.get("AIC")))
        keep.append(
            {
                "split": split,
                "best_PLAMB_H0": as_float(best_plamb.get("H0_trial")),
                "best_PLAMB_n": as_float(best_plamb.get("inertia_index")),
                "best_PLAMB_s": as_float(best_plamb.get("screen_index")),
                "best_PLAMB_AIC": as_float(best_plamb.get("AIC")),
                "best_PLAMB_chi2_dof": as_float(best_plamb.get("chi2_dof")),
                "best_RAR_H0": as_float(best_rar.get("H0_trial")),
                "best_RAR_AIC": as_float(best_rar.get("AIC")),
                "best_RAR_chi2_dof": as_float(best_rar.get("chi2_dof")),
                "delta_AIC_PLAMB_minus_RAR": as_float(best_plamb.get("AIC")) - as_float(best_rar.get("AIC")),
                "N_points": as_int(best_plamb.get("N_points")),
            }
        )
    return keep


def chi2_by_galaxy(data: Dataset, pred: np.ndarray, label: str, h0: float, split: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for galaxy in sorted(set(data.galaxy.tolist())):
        mask = data.galaxy == galaxy
        sub = subset_dataset(data, mask)
        residual = (sub.vobs - pred[mask]) / sub.sigma_v
        chi2 = float(np.sum(residual**2))
        rows.append(
            {
                "Galaxy": galaxy,
                "split": split,
                "model": label,
                "H0_trial": h0,
                "N_points": int(np.sum(mask)),
                "chi2": chi2,
                "chi2_per_point": chi2 / max(int(np.sum(mask)), 1),
                "Q": int(sub.quality[0]),
                "f_D": int(sub.distance_method[0]),
                "distance_group": group_label(int(sub.distance_method[0])),
                "D_Mpc": float(sub.d_mpc[0]),
                "Inc_deg": float(sub.inc_deg[0]),
                "Vobs_outer_kms": float(sub.vobs[np.argmax(sub.rad_kpc)]),
                "Rad_outer_kpc": float(np.max(sub.rad_kpc)),
            }
        )
    return rows


def residual_leverage_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    data = load_dataset(args.sample, args.points, quality_max=2, distance_method="all", err_floor=args.err_floor_kms)

    # Use the all_Q2 fixed-exponent best and the all_Q2 RAR benchmark so the
    # same galaxies and same N are compared directly.
    h0_plamb = 65.0
    h0_rar = 65.0
    data_plamb = rescale_for_h0(data, h0_plamb, 73.0, "hubble_flow_only")
    data_rar = rescale_for_h0(data, h0_rar, 73.0, "hubble_flow_only")

    # Values from the fixed-exponent all_Q2 summary and RAR grid row.
    pred_plamb = predicted_velocity_fixed_plamb(data_plamb, ydisk=0.54903, ybul=0.68539, h0=h0_plamb, inertia_index=0.425, screen_index=1.25)
    pred_rar = predicted_velocity_rar(data_rar, ydisk=0.53271, ybul=0.65518, log10_gdagger=-10.029)

    plamb_rows = chi2_by_galaxy(data_plamb, pred_plamb, "PLAMB_FIXED_allQ2_best", h0_plamb, "all_Q2")
    rar_rows = chi2_by_galaxy(data_rar, pred_rar, "RAR_allQ2_best", h0_rar, "all_Q2")
    by_key = {(row["Galaxy"], row["model"]): row for row in plamb_rows + rar_rows}
    for galaxy in sorted({row["Galaxy"] for row in plamb_rows}):
        p = by_key[(galaxy, "PLAMB_FIXED_allQ2_best")]
        r = by_key[(galaxy, "RAR_allQ2_best")]
        merged = dict(p)
        merged["chi2_PLAMB"] = p["chi2"]
        merged["chi2_RAR"] = r["chi2"]
        merged["delta_chi2_PLAMB_minus_RAR"] = p["chi2"] - r["chi2"]
        rows.append(merged)
    return sorted(rows, key=lambda row: abs(float(row["delta_chi2_PLAMB_minus_RAR"])), reverse=True)


def make_property_plot(rows: list[dict[str, Any]], path: Path) -> None:
    if plt is None:
        return
    selected = [row for row in rows if row["property"] in {"D_Mpc", "Inc_deg", "Vflat_kms", "SBdisk_Lsun_pc2", "MHI_1e9_Msun", "Q"}]
    labels = [row["property"] for row in selected]
    hubble = [row["hubble_median"] for row in selected]
    non = [row["non_hubble_median"] for row in selected]
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(9, 5.6))
    ax.bar(x - 0.18, hubble, width=0.36, label="Hubble-flow f_D=1")
    ax.bar(x + 0.18, non, width=0.36, label="Non-Hubble f_D!=1")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_ylabel("Median value")
    ax.set_title("SPARC Distance-Method Sample Differences")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def write_report(property_rows: list[dict[str, Any]], grid_rows: list[dict[str, Any]], leverage_rows: list[dict[str, Any]]) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    notable_props = sorted(property_rows, key=lambda row: as_float(row["mann_whitney_p"]))[:8]
    top_leverage = leverage_rows[:15]
    lines = [
        "# SPARC Distance-Method Driver Audit",
        "",
        "## Executive Summary",
        "",
        "- The fixed-exponent PLAMB branch beats RAR on all-sample and Hubble-flow-only splits.",
        "- RAR beats fixed-exponent PLAMB on non-Hubble-only splits.",
        "- Therefore, the current SPARC result is distance-method dependent and should not be presented as distance-calibration independent yet.",
        "",
        "## Comparable Split Summary",
        "",
        "| Split | N | Best PLAMB `(H0,n,s)` | Delta AIC PLAMB-RAR |",
        "|---|---:|---|---:|",
    ]
    for row in grid_rows:
        lines.append(
            f"| {row['split']} | {row['N_points']} | "
            f"({row['best_PLAMB_H0']:.1f}, {row['best_PLAMB_n']:.4g}, {row['best_PLAMB_s']:.4g}) | "
            f"{row['delta_AIC_PLAMB_minus_RAR']:.3f} |"
        )
    lines.extend(
        [
            "",
            "## Largest Hubble vs Non-Hubble Sample Differences",
            "",
            "| Property | Hubble median | Non-Hubble median | Mann-Whitney p |",
            "|---|---:|---:|---:|",
        ]
    )
    for row in notable_props:
        lines.append(
            f"| {row['property']} | {row['hubble_median']:.5g} | {row['non_hubble_median']:.5g} | {row['mann_whitney_p']:.3g} |"
        )
    lines.extend(
        [
            "",
            "## Highest-Leverage Galaxies In all_Q2 PLAMB-vs-RAR",
            "",
            "| Galaxy | group | Q | chi2_PLAMB | chi2_RAR | delta chi2 |",
            "|---|---|---:|---:|---:|---:|",
        ]
    )
    for row in top_leverage:
        lines.append(
            f"| {row['Galaxy']} | {row['distance_group']} | {row['Q']} | "
            f"{row['chi2_PLAMB']:.3f} | {row['chi2_RAR']:.3f} | {row['delta_chi2_PLAMB_minus_RAR']:.3f} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The preferred fixed exponent is not universal across distance groups. Full/all and Hubble-flow subsets prefer approximately `n=0.4-0.425`, while non-Hubble subsets prefer `n=0.5` and stronger screening. This could mean either:",
            "",
            "1. the PLAMB rescaling is picking up a distance-calibration effect in Hubble-flow galaxies,",
            "2. the non-Hubble galaxies have different morphology/surface-brightness/velocity distributions,",
            "3. the fixed exponent grid is still too coarse, or",
            "4. a few high-leverage galaxies dominate the model preference.",
            "",
            "Before making a strong claim, repeat this with galaxy-level distance uncertainty propagation and either hierarchical M/L priors or jackknife removal of high-leverage systems.",
        ]
    )
    (OUTDIR / "sparc_distance_method_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    metadata = {
        "n_property_rows": len(property_rows),
        "n_grid_rows": len(grid_rows),
        "n_leverage_rows": len(leverage_rows),
    }
    (OUTDIR / "sparc_distance_method_audit.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit SPARC distance-method drivers.")
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--points", type=Path, default=DEFAULT_POINTS)
    parser.add_argument("--grid", type=Path, default=GRID_PATH)
    parser.add_argument("--err-floor-kms", type=float, default=3.0)
    args = parser.parse_args()

    sample = load_sample_rows(args.sample)
    property_rows = sample_property_rows(sample)
    grid_rows = best_grid_rows(args.grid)
    leverage_rows = residual_leverage_rows(args)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    write_csv(
        OUTDIR / "sparc_distance_method_sample_properties.csv",
        property_rows,
        [
            "property",
            "hubble_N",
            "non_hubble_N",
            "hubble_median",
            "non_hubble_median",
            "hubble_mean",
            "non_hubble_mean",
            "mann_whitney_p",
        ],
    )
    write_csv(
        OUTDIR / "sparc_distance_method_model_preference.csv",
        grid_rows,
        [
            "split",
            "N_points",
            "best_PLAMB_H0",
            "best_PLAMB_n",
            "best_PLAMB_s",
            "best_PLAMB_AIC",
            "best_PLAMB_chi2_dof",
            "best_RAR_H0",
            "best_RAR_AIC",
            "best_RAR_chi2_dof",
            "delta_AIC_PLAMB_minus_RAR",
        ],
    )
    write_csv(
        OUTDIR / "sparc_distance_method_leverage_galaxies.csv",
        leverage_rows,
        [
            "Galaxy",
            "distance_group",
            "Q",
            "f_D",
            "D_Mpc",
            "Inc_deg",
            "Vobs_outer_kms",
            "Rad_outer_kpc",
            "N_points",
            "chi2_PLAMB",
            "chi2_RAR",
            "delta_chi2_PLAMB_minus_RAR",
        ],
    )
    make_property_plot(property_rows, OUTDIR / "sparc_distance_method_property_medians.png")
    write_report(property_rows, grid_rows, leverage_rows)

    print(f"Saved report: {OUTDIR / 'sparc_distance_method_audit.md'}")
    for row in grid_rows:
        print(
            f"{row['split']:16s} dAIC={row['delta_AIC_PLAMB_minus_RAR']:9.3f} "
            f"PLAMB=({row['best_PLAMB_H0']:.1f}, {row['best_PLAMB_n']:.4g}, {row['best_PLAMB_s']:.4g})"
        )


if __name__ == "__main__":
    main()
