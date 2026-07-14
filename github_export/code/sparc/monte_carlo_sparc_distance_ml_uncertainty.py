#!/usr/bin/env python3
r"""
Monte Carlo SPARC uncertainty propagation for fixed PLAMB vs RAR.

This test propagates galaxy-level distance and stellar mass-to-light uncertainty
instead of relying on a single pointwise AIC. For each draw it perturbs each
galaxy distance, applies a galaxy-level stellar M/L multiplier, and compares
fixed-exponent PLAMB with RAR.

The implementation is intentionally a diagnostic approximation:

  - observed velocities are left fixed,
  - radii scale as D_new / D_old,
  - gas/disk/bulge model velocities scale as sqrt(D_new / D_old),
  - disk/bulge velocities additionally scale as sqrt(eta_ML),
  - RAR and PLAMB retain their selected shape/scale settings.

Run:
    python monte_carlo_sparc_distance_ml_uncertainty.py
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
)
from diagnose_plamb_sparc_rotation import DEFAULT_POINTS, DEFAULT_SAMPLE, Dataset, load_dataset


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "sparc_distance_ml_monte_carlo"


SPLIT_CONFIGS = [
    {
        "split": "all_Q2",
        "quality_max": 2,
        "distance_method": "all",
        "plamb_h0": 65.0,
        "plamb_n": 0.425,
        "plamb_s": 1.25,
        "rar_h0": 65.0,
        "rar_log10_gdagger": -10.029,
    },
    {
        "split": "hubble_flow_Q2",
        "quality_max": 2,
        "distance_method": "hubble_flow",
        "plamb_h0": 70.0,
        "plamb_n": 0.4,
        "plamb_s": 1.25,
        "rar_h0": 62.5,
        "rar_log10_gdagger": -10.0,
    },
    {
        "split": "non_hubble_Q2",
        "quality_max": 2,
        "distance_method": "non_hubble",
        "plamb_h0": 70.0,
        "plamb_n": 0.5,
        "plamb_s": 1.25,
        "rar_h0": 62.5,
        "rar_log10_gdagger": -10.0,
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


def parse_float(value: Any, default: float = float("nan")) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def sample_distance_errors(sample_path: Path) -> dict[str, tuple[float, float]]:
    rows = {}
    for row in read_csv(sample_path):
        rows[row["Galaxy"]] = (parse_float(row.get("D_Mpc")), parse_float(row.get("e_D_Mpc")))
    return rows


def apply_draw(
    data: Dataset,
    distance_info: dict[str, tuple[float, float]],
    rng: np.random.Generator,
    distance_sigma_scale: float,
    sigma_ln_ml: float,
    distance_floor_frac: float,
) -> Dataset:
    galaxies = sorted(set(data.galaxy.tolist()))
    distance_scale_by_galaxy: dict[str, float] = {}
    ml_eta_by_galaxy: dict[str, float] = {}
    for galaxy in galaxies:
        d_old, e_d = distance_info.get(str(galaxy), (float("nan"), float("nan")))
        if not math.isfinite(d_old) or d_old <= 0.0:
            d_old = float(np.nanmedian(data.d_mpc[data.galaxy == galaxy]))
        if not math.isfinite(e_d) or e_d <= 0.0:
            e_d = distance_floor_frac * d_old
        sigma_d = max(e_d * distance_sigma_scale, distance_floor_frac * d_old)
        d_new = rng.normal(d_old, sigma_d)
        d_new = max(d_new, 0.05 * d_old, 0.01)
        distance_scale_by_galaxy[str(galaxy)] = d_new / d_old
        ml_eta_by_galaxy[str(galaxy)] = float(np.exp(rng.normal(0.0, sigma_ln_ml)))

    scale = np.array([distance_scale_by_galaxy[str(galaxy)] for galaxy in data.galaxy], dtype=float)
    eta = np.array([ml_eta_by_galaxy[str(galaxy)] for galaxy in data.galaxy], dtype=float)
    velocity_distance_scale = np.sqrt(scale)
    velocity_stellar_scale = np.sqrt(eta)
    return replace(
        data,
        rad_kpc=data.rad_kpc * scale,
        d_mpc=data.d_mpc * scale,
        vgas=data.vgas * velocity_distance_scale,
        vdisk=data.vdisk * velocity_distance_scale * velocity_stellar_scale,
        vbul=data.vbul * velocity_distance_scale * velocity_stellar_scale,
    )


def score_models(data: Dataset, config: dict[str, Any]) -> dict[str, float]:
    if minimize is None:
        raise RuntimeError(f"SciPy optimizer import failed: {SCIPY_IMPORT_ERROR}")

    def fit_plamb() -> tuple[np.ndarray, float]:
        def objective(theta: np.ndarray) -> float:
            pred = predicted_velocity_fixed_plamb(
                data,
                float(theta[0]),
                float(theta[1]),
                float(config["plamb_h0"]),
                float(config["plamb_n"]),
                float(config["plamb_s"]),
            )
            return chi2(pred, data)

        return fit_with_starts(objective, [0.5, 0.7], [(0.05, 1.25), (0.05, 1.60)])

    def fit_rar() -> tuple[np.ndarray, float]:
        def objective(theta: np.ndarray) -> float:
            pred = predicted_velocity_rar(
                data,
                float(theta[0]),
                float(theta[1]),
                float(theta[2]),
            )
            return chi2(pred, data)

        return fit_with_starts(objective, [0.5, 0.7, float(config["rar_log10_gdagger"])], [(0.05, 1.25), (0.05, 1.60), (-11.6, -9.2)])

    plamb_theta, plamb_chi2 = fit_plamb()
    rar_theta, rar_chi2 = fit_rar()
    n = len(data.vobs)
    plamb_aic = plamb_chi2 + 4.0
    rar_aic = rar_chi2 + 6.0
    return {
        "N_points": n,
        "N_galaxies": len(set(data.galaxy.tolist())),
        "plamb_ydisk": float(plamb_theta[0]),
        "plamb_ybul": float(plamb_theta[1]),
        "rar_ydisk": float(rar_theta[0]),
        "rar_ybul": float(rar_theta[1]),
        "rar_log10_gdagger_fit": float(rar_theta[2]),
        "plamb_chi2": plamb_chi2,
        "rar_chi2": rar_chi2,
        "plamb_chi2_dof": plamb_chi2 / max(n - 2, 1),
        "rar_chi2_dof": rar_chi2 / max(n - 3, 1),
        "plamb_AIC": plamb_aic,
        "rar_AIC": rar_aic,
        "delta_AIC_PLAMB_minus_RAR": plamb_aic - rar_aic,
    }


def fit_with_starts(objective, default_start: list[float], bounds: list[tuple[float, float]]) -> tuple[np.ndarray, float]:
    starts = [default_start]
    if len(default_start) == 2:
        starts.extend([[0.35, 0.55], [0.8, 0.9], [1.0, 1.0]])
    elif len(default_start) == 3:
        starts.extend([[0.35, 0.55, -10.1], [0.8, 0.9, -9.9], [1.0, 1.0, -10.0]])
    best_x: np.ndarray | None = None
    best_fun = 1.0e100
    for start in starts:
        result = minimize(
            objective,
            np.asarray(start, dtype=float),
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 600, "ftol": 1e-9},
        )
        if result.success and float(result.fun) < best_fun:
            best_x = np.asarray(result.x, dtype=float)
            best_fun = float(result.fun)
    if best_x is None:
        raise RuntimeError("Monte Carlo fit failed.")
    return best_x, best_fun


def summarize(draw_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for split in sorted({row["split"] for row in draw_rows}):
        vals = np.asarray([float(row["delta_AIC_PLAMB_minus_RAR"]) for row in draw_rows if row["split"] == split], dtype=float)
        rows.append(
            {
                "split": split,
                "draws": len(vals),
                "frac_PLAMB_wins": float(np.mean(vals < 0.0)),
                "delta_AIC_mean": float(np.mean(vals)),
                "delta_AIC_median": float(np.median(vals)),
                "delta_AIC_p16": float(np.percentile(vals, 16)),
                "delta_AIC_p84": float(np.percentile(vals, 84)),
                "delta_AIC_min": float(np.min(vals)),
                "delta_AIC_max": float(np.max(vals)),
            }
        )
    return rows


def make_plot(draw_rows: list[dict[str, Any]], path: Path) -> None:
    if plt is None:
        return
    splits = sorted({row["split"] for row in draw_rows})
    data = [[float(row["delta_AIC_PLAMB_minus_RAR"]) for row in draw_rows if row["split"] == split] for split in splits]
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 5.8))
    ax.boxplot(data, labels=splits, showfliers=False)
    ax.axhline(0.0, color="#555555", linewidth=1.0)
    ax.set_ylabel("Delta AIC: PLAMB fixed exponent minus RAR")
    ax.set_title("SPARC Distance/M-L Monte Carlo")
    ax.tick_params(axis="x", rotation=20)
    ax.grid(axis="y", alpha=0.22)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def write_report(summary_rows: list[dict[str, Any]], args: argparse.Namespace) -> None:
    lines = [
        "# SPARC Distance and M/L Monte Carlo",
        "",
        "## Executive Summary",
        "",
        "Negative Delta AIC means fixed-exponent PLAMB wins; positive means RAR wins.",
        "",
        "| Split | Draws | PLAMB win fraction | Median Delta AIC | 16-84% interval |",
        "|---|---:|---:|---:|---|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['split']} | {row['draws']} | {row['frac_PLAMB_wins']:.3f} | "
            f"{row['delta_AIC_median']:.3f} | [{row['delta_AIC_p16']:.3f}, {row['delta_AIC_p84']:.3f}] |"
        )
    lines.extend(
        [
            "",
            "## Assumptions",
            "",
            f"- Distance draws use SPARC `D_Mpc` and `e_D_Mpc`, scaled by `{args.distance_sigma_scale}`.",
            f"- A minimum fractional distance floor of `{args.distance_floor_frac}` is applied.",
            f"- Stellar M/L eta is lognormal with sigma_ln = `{args.sigma_ln_ml}`.",
            "- Radius and model baryonic velocities are rescaled consistently with distance.",
            "- Each draw refits global PLAMB disk/bulge M/L and RAR disk/bulge/g_dagger parameters after perturbation.",
            "- This is a diagnostic approximation, not a full hierarchical posterior.",
        ]
    )
    (OUTDIR / "sparc_distance_ml_monte_carlo_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="SPARC distance/M-L uncertainty Monte Carlo.")
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--points", type=Path, default=DEFAULT_POINTS)
    parser.add_argument("--draws", type=int, default=500)
    parser.add_argument("--seed", type=int, default=260713)
    parser.add_argument("--err-floor-kms", type=float, default=3.0)
    parser.add_argument("--distance-sigma-scale", type=float, default=1.0)
    parser.add_argument("--distance-floor-frac", type=float, default=0.03)
    parser.add_argument("--sigma-ln-ml", type=float, default=0.25)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    distance_info = sample_distance_errors(args.sample)
    draw_rows: list[dict[str, Any]] = []
    for config in SPLIT_CONFIGS:
        base = load_dataset(
            args.sample,
            args.points,
            quality_max=int(config["quality_max"]),
            distance_method=str(config["distance_method"]),
            err_floor=args.err_floor_kms,
        )
        print(f"=== {config['split']} N={len(base.vobs)} galaxies={len(set(base.galaxy.tolist()))} ===")
        for draw in range(args.draws):
            drawn = apply_draw(
                base,
                distance_info,
                rng,
                distance_sigma_scale=args.distance_sigma_scale,
                sigma_ln_ml=args.sigma_ln_ml,
                distance_floor_frac=args.distance_floor_frac,
            )
            scores = score_models(drawn, config)
            row: dict[str, Any] = {
                "split": config["split"],
                "draw": draw,
                "plamb_H0": config["plamb_h0"],
                "plamb_n": config["plamb_n"],
                "plamb_s": config["plamb_s"],
                "rar_H0": config["rar_h0"],
                "rar_log10_gdagger": config["rar_log10_gdagger"],
            }
            row.update(scores)
            draw_rows.append(row)
        vals = [row["delta_AIC_PLAMB_minus_RAR"] for row in draw_rows if row["split"] == config["split"]]
        print(f"    PLAMB win fraction={np.mean(np.asarray(vals) < 0.0):.3f}, median dAIC={np.median(vals):.3f}")

    summary_rows = summarize(draw_rows)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    write_csv(
        OUTDIR / "sparc_distance_ml_monte_carlo_draws.csv",
        draw_rows,
        [
            "split",
            "draw",
            "N_points",
            "N_galaxies",
            "plamb_H0",
            "plamb_n",
            "plamb_s",
            "rar_H0",
            "rar_log10_gdagger",
            "plamb_ydisk",
            "plamb_ybul",
            "rar_ydisk",
            "rar_ybul",
            "rar_log10_gdagger_fit",
            "plamb_chi2",
            "rar_chi2",
            "plamb_chi2_dof",
            "rar_chi2_dof",
            "plamb_AIC",
            "rar_AIC",
            "delta_AIC_PLAMB_minus_RAR",
        ],
    )
    write_csv(
        OUTDIR / "sparc_distance_ml_monte_carlo_summary.csv",
        summary_rows,
        [
            "split",
            "draws",
            "frac_PLAMB_wins",
            "delta_AIC_mean",
            "delta_AIC_median",
            "delta_AIC_p16",
            "delta_AIC_p84",
            "delta_AIC_min",
            "delta_AIC_max",
        ],
    )
    make_plot(draw_rows, OUTDIR / "sparc_distance_ml_monte_carlo_delta_aic.png")
    write_report(summary_rows, args)
    (OUTDIR / "sparc_distance_ml_monte_carlo_metadata.json").write_text(
        json.dumps({"args": vars(args)}, indent=2, default=str),
        encoding="utf-8",
    )
    print(f"Saved report: {OUTDIR / 'sparc_distance_ml_monte_carlo_report.md'}")


if __name__ == "__main__":
    main()
