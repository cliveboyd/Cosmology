#!/usr/bin/env python3
r"""
Fixed-exponent PLAMB SPARC grid.

This script stress-tests the promising screened PLAMB rotation form by removing
most shape freedom. For each H0/distance/quality split and each fixed exponent
pair, it fits only global disk/bulge stellar mass-to-light ratios:

    g_pred = g_bar * [1 + (g0/g_bar)^n / (1 + (g_bar/g0)^s)]

with

    g0 = c H0 / (2 pi)

The goal is to discover whether simple, physically discussable exponent pairs
survive compared with the empirical RAR benchmark.

Run:
    python diagnose_plamb_sparc_fixed_exponent_grid.py
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
    from scipy.stats import spearmanr
except Exception as exc:  # pragma: no cover - reported at runtime
    minimize = None
    spearmanr = None
    SCIPY_IMPORT_ERROR = exc
else:
    SCIPY_IMPORT_ERROR = None

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover - optional plotting
    plt = None
    MATPLOTLIB_IMPORT_ERROR = exc
else:
    MATPLOTLIB_IMPORT_ERROR = None

from diagnose_plamb_sparc_h0_nuisance import (
    acceleration_cH0_over_2pi,
    rescale_for_h0,
)
from diagnose_plamb_sparc_rotation import (
    DEFAULT_POINTS,
    DEFAULT_SAMPLE,
    Dataset,
    acceleration_from_v2,
    baryon_v2,
    load_dataset,
    velocity_from_acceleration,
)


ROOT = Path(__file__).resolve().parents[3]
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "sparc_rotation_fixed_exponent_grid"


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_float_list(value: str) -> list[float]:
    vals = [float(token.strip()) for token in value.split(",") if token.strip()]
    if not vals:
        raise ValueError("Expected at least one numeric value.")
    return vals


def parse_int_list(value: str) -> list[int]:
    vals = [int(float(token.strip())) for token in value.split(",") if token.strip()]
    if not vals:
        raise ValueError("Expected at least one integer value.")
    return vals


def subset_dataset(data: Dataset, mask: np.ndarray) -> Dataset:
    return replace(
        data,
        galaxy=data.galaxy[mask],
        rad_kpc=data.rad_kpc[mask],
        vobs=data.vobs[mask],
        sigma_v=data.sigma_v[mask],
        vgas=data.vgas[mask],
        vdisk=data.vdisk[mask],
        vbul=data.vbul[mask],
        quality=data.quality[mask],
        distance_method=data.distance_method[mask],
        inc_deg=data.inc_deg[mask],
        d_mpc=data.d_mpc[mask],
    )


def dataset_for_split(
    sample: Path,
    points: Path,
    quality_max: int,
    distance_method: str,
    err_floor_kms: float,
) -> Dataset:
    return load_dataset(sample, points, quality_max, distance_method, err_floor_kms)


def predicted_velocity_fixed_plamb(data: Dataset, ydisk: float, ybul: float, h0: float, inertia_index: float, screen_index: float) -> np.ndarray:
    g0 = acceleration_cH0_over_2pi(h0)
    gbar = acceleration_from_v2(baryon_v2(data, ydisk, ybul), data.rad_kpc)
    ratio = np.clip(g0 / np.maximum(gbar, 1e-30), 1e-8, 1e8)
    boost = 1.0 + ratio**inertia_index / (1.0 + np.clip(1.0 / ratio, 1e-8, 1e8) ** screen_index)
    return velocity_from_acceleration(gbar * boost, data.rad_kpc)


def predicted_velocity_rar(data: Dataset, ydisk: float, ybul: float, log10_gdagger: float) -> np.ndarray:
    gbar = acceleration_from_v2(baryon_v2(data, ydisk, ybul), data.rad_kpc)
    gdagger = 10.0**log10_gdagger
    x = np.sqrt(np.maximum(gbar / gdagger, 1e-30))
    gpred = gbar / np.maximum(1.0 - np.exp(-x), 1e-12)
    return velocity_from_acceleration(gpred, data.rad_kpc)


def predicted_velocity_baryon(data: Dataset, ydisk: float, ybul: float) -> np.ndarray:
    return np.sqrt(np.maximum(baryon_v2(data, ydisk, ybul), 0.0))


def chi2(pred: np.ndarray, data: Dataset) -> float:
    if pred.shape != data.vobs.shape or not np.all(np.isfinite(pred)):
        return 1.0e100
    return float(np.sum(((data.vobs - pred) / data.sigma_v) ** 2))


def fit_two_ml(data: Dataset, predictor) -> tuple[np.ndarray, float]:
    if minimize is None:
        raise RuntimeError(f"SciPy optimizer import failed: {SCIPY_IMPORT_ERROR}")

    def objective(theta: np.ndarray) -> float:
        ydisk, ybul = theta
        pred = predictor(data, ydisk, ybul)
        return chi2(pred, data)

    starts = ([0.5, 0.7], [0.3, 0.5], [0.8, 0.9], [1.0, 1.0])
    best_x: np.ndarray | None = None
    best_fun = 1.0e100
    for start in starts:
        result = minimize(
            objective,
            np.asarray(start, dtype=float),
            method="L-BFGS-B",
            bounds=[(0.05, 1.25), (0.05, 1.60)],
            options={"maxiter": 800, "ftol": 1e-10},
        )
        if result.success and float(result.fun) < best_fun:
            best_fun = float(result.fun)
            best_x = np.asarray(result.x, dtype=float)
    if best_x is None:
        raise RuntimeError("M/L fit did not converge.")
    return best_x, best_fun


def fit_rar(data: Dataset) -> tuple[np.ndarray, float]:
    if minimize is None:
        raise RuntimeError(f"SciPy optimizer import failed: {SCIPY_IMPORT_ERROR}")

    def objective(theta: np.ndarray) -> float:
        ydisk, ybul, log10_gdagger = theta
        pred = predicted_velocity_rar(data, ydisk, ybul, log10_gdagger)
        return chi2(pred, data)

    starts = ([0.5, 0.7, -10.0], [0.35, 0.6, -10.1], [0.75, 0.9, -9.9])
    best_x: np.ndarray | None = None
    best_fun = 1.0e100
    for start in starts:
        result = minimize(
            objective,
            np.asarray(start, dtype=float),
            method="L-BFGS-B",
            bounds=[(0.05, 1.25), (0.05, 1.60), (-11.6, -9.2)],
            options={"maxiter": 1000, "ftol": 1e-10},
        )
        if result.success and float(result.fun) < best_fun:
            best_fun = float(result.fun)
            best_x = np.asarray(result.x, dtype=float)
    if best_x is None:
        raise RuntimeError("RAR fit did not converge.")
    return best_x, best_fun


def score_row(
    model: str,
    split_label: str,
    h0: float,
    data: Dataset,
    k: int,
    chi2_value: float,
    params: dict[str, Any],
) -> dict[str, Any]:
    n = len(data.vobs)
    aic = chi2_value + 2.0 * k
    bic = chi2_value + k * math.log(max(n, 1))
    row = {
        "model": model,
        "split": split_label,
        "H0_trial": h0,
        "N_points": n,
        "N_galaxies": len(set(data.galaxy.tolist())),
        "k": k,
        "chi2": chi2_value,
        "dof": n - k,
        "chi2_dof": chi2_value / max(n - k, 1),
        "AIC": aic,
        "BIC": bic,
    }
    row.update(params)
    return row


def residual_correlations(data: Dataset, pred: np.ndarray, label: str, split: str, h0: float) -> list[dict[str, Any]]:
    residual = data.vobs - pred
    abs_resid = np.abs(residual)
    gbar = acceleration_from_v2(baryon_v2(data, 0.5, 0.7), data.rad_kpc)
    variables = {
        "radius_kpc": data.rad_kpc,
        "vobs_kms": data.vobs,
        "sigma_v_kms": data.sigma_v,
        "inclination_deg": data.inc_deg,
        "distance_mpc": data.d_mpc,
        "quality_flag": data.quality.astype(float),
        "distance_method": data.distance_method.astype(float),
        "log10_gbar_fixed_ml": np.log10(np.maximum(gbar, 1e-30)),
    }
    rows: list[dict[str, Any]] = []
    for name, values in variables.items():
        mask = np.isfinite(values) & np.isfinite(abs_resid)
        if spearmanr is None or int(np.sum(mask)) < 8:
            rho = float("nan")
            pval = float("nan")
        else:
            result = spearmanr(values[mask], abs_resid[mask])
            rho = float(result.statistic)
            pval = float(result.pvalue)
        rows.append(
            {
                "model": label,
                "split": split,
                "H0_trial": h0,
                "variable": name,
                "spearman_rho_abs_resid": rho,
                "p_value": pval,
                "N": int(np.sum(mask)),
            }
        )
    return rows


def stable_hash(value: str) -> int:
    total = 0
    for char in value:
        total = (total * 131 + ord(char)) % 1_000_003
    return total


def holdout_masks(data: Dataset, fold: int) -> tuple[np.ndarray, np.ndarray]:
    galaxies = np.array([str(galaxy) for galaxy in data.galaxy], dtype=object)
    test_galaxies = {galaxy for galaxy in sorted(set(galaxies.tolist())) if stable_hash(galaxy) % 5 == fold}
    test_mask = np.array([galaxy in test_galaxies for galaxy in galaxies], dtype=bool)
    train_mask = ~test_mask
    return train_mask, test_mask


def write_report(rows: list[dict[str, Any]], corr_rows: list[dict[str, Any]], args: argparse.Namespace) -> None:
    global_rows = [row for row in rows if not str(row["model"]).endswith("_HOLDOUT")]
    holdout_rows = [row for row in rows if str(row["model"]).endswith("_HOLDOUT")]
    best = min(global_rows, key=lambda row: float(row["AIC"]))
    best_rar = min([row for row in global_rows if row["model"] == "RAR_FREE_GDAGGER"], key=lambda row: float(row["AIC"]))
    best_plamb = min([row for row in global_rows if row["model"] == "PLAMB_FIXED_EXPONENT"], key=lambda row: float(row["AIC"]))
    lines = [
        "# Fixed-Exponent PLAMB SPARC Grid",
        "",
        "## Executive Summary",
        "",
        f"- H0 values: {args.h0_values}",
        f"- inertia exponents: {args.inertia_indices}",
        f"- screen exponents: {args.screen_indices}",
        f"- distance splits: {args.distance_splits}",
        f"- quality cuts: {args.quality_max_values}",
        f"- Best global row: {best['model']} split={best['split']} H0={float(best['H0_trial']):.3f} AIC={float(best['AIC']):.3f}",
        f"- Best fixed PLAMB row: split={best_plamb['split']} H0={float(best_plamb['H0_trial']):.3f} n={float(best_plamb['inertia_index']):.4g} s={float(best_plamb['screen_index']):.4g} AIC={float(best_plamb['AIC']):.3f}",
        f"- Best RAR row: split={best_rar['split']} H0={float(best_rar['H0_trial']):.3f} AIC={float(best_rar['AIC']):.3f}",
        f"- Delta AIC fixed PLAMB minus RAR: {float(best_plamb['AIC']) - float(best_rar['AIC']):.3f}",
        "",
        "## Top Rows By AIC",
        "",
        "| Model | Split | H0 | n | s | chi2/dof | AIC | BIC |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in sorted(global_rows, key=lambda item: float(item["AIC"]))[:25]:
        lines.append(
            f"| {row['model']} | {row['split']} | {float(row['H0_trial']):.3f} | "
            f"{float(row.get('inertia_index', float('nan'))):.4g} | "
            f"{float(row.get('screen_index', float('nan'))):.4g} | "
            f"{float(row['chi2_dof']):.3f} | {float(row['AIC']):.3f} | {float(row['BIC']):.3f} |"
        )
    if holdout_rows:
        lines.extend(
            [
                "",
                "## Holdout Rows",
                "",
                "| Model | Split | H0 | n | s | train chi2/dof | test chi2/dof |",
                "|---|---|---:|---:|---:|---:|---:|",
            ]
        )
        for row in sorted(holdout_rows, key=lambda item: float(item["test_chi2_dof"]))[:20]:
            lines.append(
                f"| {row['model']} | {row['split']} | {float(row['H0_trial']):.3f} | "
                f"{float(row.get('inertia_index', float('nan'))):.4g} | "
                f"{float(row.get('screen_index', float('nan'))):.4g} | "
                f"{float(row['train_chi2_dof']):.3f} | {float(row['test_chi2_dof']):.3f} |"
            )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- Grid CSV: `{OUTDIR / 'sparc_fixed_exponent_grid.csv'}`",
            f"- Top CSV: `{OUTDIR / 'sparc_fixed_exponent_top.csv'}`",
            f"- Residual correlations: `{OUTDIR / 'sparc_fixed_exponent_residual_correlations.csv'}`",
            f"- Plot: `{OUTDIR / 'sparc_fixed_exponent_aic_heatmap.png'}`",
        ]
    )
    (OUTDIR / "sparc_fixed_exponent_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    meta = {
        "args": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "n_rows": len(rows),
        "n_corr_rows": len(corr_rows),
        "best_row": best,
        "best_plamb": best_plamb,
        "best_rar": best_rar,
    }
    (OUTDIR / "sparc_fixed_exponent_metadata.json").write_text(json.dumps(meta, indent=2, sort_keys=True), encoding="utf-8")


def make_heatmap(rows: list[dict[str, Any]], path: Path) -> None:
    if plt is None:
        return
    plamb_rows = [row for row in rows if row["model"] == "PLAMB_FIXED_EXPONENT" and row["split"] == "all_Q2"]
    if not plamb_rows:
        return
    best_by_pair: dict[tuple[float, float], float] = {}
    for row in plamb_rows:
        pair = (float(row["inertia_index"]), float(row["screen_index"]))
        best_by_pair[pair] = min(best_by_pair.get(pair, float("inf")), float(row["AIC"]))
    ns = sorted({pair[0] for pair in best_by_pair})
    ss = sorted({pair[1] for pair in best_by_pair})
    matrix = np.full((len(ss), len(ns)), np.nan)
    for i_s, s_value in enumerate(ss):
        for i_n, n_value in enumerate(ns):
            matrix[i_s, i_n] = best_by_pair.get((n_value, s_value), np.nan)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    im = ax.imshow(matrix, origin="lower", aspect="auto", cmap="viridis_r")
    ax.set_xticks(range(len(ns)))
    ax.set_xticklabels([f"{value:.3g}" for value in ns])
    ax.set_yticks(range(len(ss)))
    ax.set_yticklabels([f"{value:.3g}" for value in ss])
    ax.set_xlabel("inertia exponent n")
    ax.set_ylabel("screen exponent s")
    ax.set_title("Best AIC Across H0 for Fixed PLAMB Exponents, all_Q2")
    fig.colorbar(im, ax=ax, label="AIC")
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Fixed-exponent PLAMB SPARC grid.")
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--points", type=Path, default=DEFAULT_POINTS)
    parser.add_argument("--h0-values", default="62.5,65,67.5,70")
    parser.add_argument("--inertia-indices", default="0.333333,0.4,0.425,0.5")
    parser.add_argument("--screen-indices", default="1.0,1.25,1.5,2.0")
    parser.add_argument("--distance-splits", default="all,non_hubble,hubble_flow")
    parser.add_argument("--quality-max-values", default="1,2")
    parser.add_argument("--distance-scaling-mode", choices=["hubble_flow_only", "all", "none"], default="hubble_flow_only")
    parser.add_argument("--err-floor-kms", type=float, default=3.0)
    parser.add_argument("--holdout", action="store_true")
    parser.add_argument("--holdout-fold", type=int, default=0)
    parser.add_argument("--output-label", default="")
    args = parser.parse_args()

    global OUTDIR
    if args.output_label:
        safe = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in args.output_label)
        OUTDIR = OUTDIR / safe
    OUTDIR.mkdir(parents=True, exist_ok=True)

    h0_values = parse_float_list(args.h0_values)
    inertia_indices = parse_float_list(args.inertia_indices)
    screen_indices = parse_float_list(args.screen_indices)
    distance_splits = [token.strip() for token in args.distance_splits.split(",") if token.strip()]
    quality_values = parse_int_list(args.quality_max_values)

    rows: list[dict[str, Any]] = []
    corr_rows: list[dict[str, Any]] = []
    for quality_max in quality_values:
        for distance_split in distance_splits:
            data_base = dataset_for_split(args.sample, args.points, quality_max, distance_split, args.err_floor_kms)
            split_label = f"{distance_split}_Q{quality_max}"
            for h0 in h0_values:
                data = rescale_for_h0(data_base, h0, 73.0, args.distance_scaling_mode)
                print(f"=== split={split_label} H0={h0:.3f} N={len(data.vobs)} ===")

                bary_theta, bary_chi2 = fit_two_ml(data, predicted_velocity_baryon)
                rows.append(
                    score_row(
                        "BARYON_FREE_ML",
                        split_label,
                        h0,
                        data,
                        2,
                        bary_chi2,
                        {"ydisk": bary_theta[0], "ybul": bary_theta[1]},
                    )
                )

                rar_theta, rar_chi2 = fit_rar(data)
                rar_pred = predicted_velocity_rar(data, rar_theta[0], rar_theta[1], rar_theta[2])
                rows.append(
                    score_row(
                        "RAR_FREE_GDAGGER",
                        split_label,
                        h0,
                        data,
                        3,
                        rar_chi2,
                        {"ydisk": rar_theta[0], "ybul": rar_theta[1], "log10_gdagger": rar_theta[2]},
                    )
                )
                corr_rows.extend(residual_correlations(data, rar_pred, "RAR_FREE_GDAGGER", split_label, h0))

                for inertia_index in inertia_indices:
                    for screen_index in screen_indices:
                        predictor = (
                            lambda dataset, ydisk, ybul, n=inertia_index, s=screen_index, h=h0:
                            predicted_velocity_fixed_plamb(dataset, ydisk, ybul, h, n, s)
                        )
                        theta, fit_chi2 = fit_two_ml(data, predictor)
                        pred = predictor(data, theta[0], theta[1])
                        row = score_row(
                            "PLAMB_FIXED_EXPONENT",
                            split_label,
                            h0,
                            data,
                            2,
                            fit_chi2,
                            {
                                "ydisk": theta[0],
                                "ybul": theta[1],
                                "log10_g0_cosmo": math.log10(acceleration_cH0_over_2pi(h0)),
                                "inertia_index": inertia_index,
                                "screen_index": screen_index,
                            },
                        )
                        rows.append(row)
                        corr_rows.extend(residual_correlations(data, pred, "PLAMB_FIXED_EXPONENT", split_label, h0))
                        print(
                            f"    n={inertia_index:.4g} s={screen_index:.4g} "
                            f"chi2/dof={row['chi2_dof']:.4g} AIC={row['AIC']:.4g}"
                        )

                if args.holdout and split_label == "all_Q2":
                    train_mask, test_mask = holdout_masks(data, args.holdout_fold)
                    train = subset_dataset(data, train_mask)
                    test = subset_dataset(data, test_mask)
                    for inertia_index in inertia_indices:
                        for screen_index in screen_indices:
                            predictor = (
                                lambda dataset, ydisk, ybul, n=inertia_index, s=screen_index, h=h0:
                                predicted_velocity_fixed_plamb(dataset, ydisk, ybul, h, n, s)
                            )
                            theta, train_chi2 = fit_two_ml(train, predictor)
                            test_chi2 = chi2(predictor(test, theta[0], theta[1]), test)
                            rows.append(
                                {
                                    "model": "PLAMB_FIXED_EXPONENT_HOLDOUT",
                                    "split": split_label,
                                    "H0_trial": h0,
                                    "N_points": len(data.vobs),
                                    "N_train": len(train.vobs),
                                    "N_test": len(test.vobs),
                                    "N_galaxies": len(set(data.galaxy.tolist())),
                                    "k": 2,
                                    "train_chi2": train_chi2,
                                    "test_chi2": test_chi2,
                                    "train_chi2_dof": train_chi2 / max(len(train.vobs) - 2, 1),
                                    "test_chi2_dof": test_chi2 / max(len(test.vobs) - 2, 1),
                                    "chi2": train_chi2,
                                    "dof": len(train.vobs) - 2,
                                    "chi2_dof": train_chi2 / max(len(train.vobs) - 2, 1),
                                    "AIC": train_chi2 + 4.0,
                                    "BIC": train_chi2 + 2.0 * math.log(max(len(train.vobs), 1)),
                                    "ydisk": theta[0],
                                    "ybul": theta[1],
                                    "log10_g0_cosmo": math.log10(acceleration_cH0_over_2pi(h0)),
                                    "inertia_index": inertia_index,
                                    "screen_index": screen_index,
                                }
                            )

    fields = [
        "model",
        "split",
        "H0_trial",
        "N_points",
        "N_train",
        "N_test",
        "N_galaxies",
        "k",
        "chi2",
        "dof",
        "chi2_dof",
        "AIC",
        "BIC",
        "train_chi2",
        "test_chi2",
        "train_chi2_dof",
        "test_chi2_dof",
        "ydisk",
        "ybul",
        "log10_gdagger",
        "log10_g0_cosmo",
        "inertia_index",
        "screen_index",
    ]
    rows_sorted = sorted(rows, key=lambda row: float(row.get("AIC", float("inf"))))
    write_csv(OUTDIR / "sparc_fixed_exponent_grid.csv", rows_sorted, fields)
    write_csv(OUTDIR / "sparc_fixed_exponent_top.csv", rows_sorted[:80], fields)
    write_csv(
        OUTDIR / "sparc_fixed_exponent_residual_correlations.csv",
        sorted(corr_rows, key=lambda row: (str(row["model"]), str(row["split"]), float(row["H0_trial"]), str(row["variable"]))),
        ["model", "split", "H0_trial", "variable", "spearman_rho_abs_resid", "p_value", "N"],
    )
    make_heatmap(rows, OUTDIR / "sparc_fixed_exponent_aic_heatmap.png")
    write_report(rows, corr_rows, args)

    best_global = min([row for row in rows if not str(row["model"]).endswith("_HOLDOUT")], key=lambda row: float(row["AIC"]))
    print(f"Saved grid: {OUTDIR / 'sparc_fixed_exponent_grid.csv'}")
    print(f"Saved report: {OUTDIR / 'sparc_fixed_exponent_report.md'}")
    print(
        f"Best: {best_global['model']} split={best_global['split']} H0={best_global['H0_trial']} "
        f"AIC={best_global['AIC']:.4g} chi2/dof={best_global['chi2_dof']:.4g}"
    )


if __name__ == "__main__":
    main()
