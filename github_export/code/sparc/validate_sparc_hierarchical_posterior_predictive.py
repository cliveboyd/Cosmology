#!/usr/bin/env python3
r"""
Approximate posterior-predictive validation for the hierarchical SPARC models.

This is the follow-up to the hierarchical MAP diagnostic. It asks a stricter
question than "which model optimizes the full data best?": if a galaxy is held
out, how well does each model predict its rotation curve after marginalizing
over galaxy-level distance and stellar M/L nuisance priors?

The global-parameter posterior used here is a local Gaussian approximation
around the MAP values from fit_sparc_hierarchical_map.py. It is intentionally
lighter than a full high-dimensional MCMC run, but it avoids refitting the
held-out galaxy's nuisance parameters to its own data.

Run:
    python validate_sparc_hierarchical_posterior_predictive.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover
    plt = None
    MATPLOTLIB_IMPORT_ERROR = exc
else:
    MATPLOTLIB_IMPORT_ERROR = None

from diagnose_plamb_sparc_h0_nuisance import acceleration_cH0_over_2pi
from diagnose_plamb_sparc_rotation import DEFAULT_POINTS, DEFAULT_SAMPLE, KPC_M, Dataset
from fit_sparc_hierarchical_map import OUTDIR as MAP_BASE_OUTDIR
from fit_sparc_hierarchical_map import prepare, split_configs


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "sparc_posterior_predictive_loo"
DEFAULT_MAP_DIR = MAP_BASE_OUTDIR / "optical_depth_hierarchical_20260714"


MODEL_ORDER = [
    "RAR",
    "PLAMB_OPTICAL_DEPTH",
    "PLAMB_OPTICAL_DEPTH_KAPPA",
    "PLAMB_OPTICAL_DEPTH_KAPPA_P",
]


@dataclass
class GlobalDraws:
    ydisk: np.ndarray
    ybul: np.ndarray
    log10_gdagger: np.ndarray
    log10_kappa: np.ndarray
    bridge_exponent: np.ndarray


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
        if value in ("", None):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def logmeanexp(values: np.ndarray) -> float:
    vmax = float(np.max(values))
    return vmax + math.log(float(np.mean(np.exp(values - vmax))))


def load_map_rows(map_dir: Path) -> dict[tuple[str, str], dict[str, str]]:
    path = map_dir / "sparc_hierarchical_map_summary.csv"
    rows = read_csv(path)
    return {(row["split"], row["model"]): row for row in rows}


def draw_global_params(row: dict[str, str], model: str, n_draws: int, rng: np.random.Generator, args: argparse.Namespace) -> GlobalDraws:
    ydisk0 = parse_float(row["ydisk"])
    ybul0 = parse_float(row["ybul"])
    log_ydisk = rng.normal(math.log(ydisk0), args.global_log_ml_sd, size=n_draws)
    log_ybul = rng.normal(math.log(ybul0), args.global_log_ml_sd, size=n_draws)
    ydisk = np.clip(np.exp(log_ydisk), 0.05, 1.50)
    ybul = np.clip(np.exp(log_ybul), 0.05, 1.80)

    log10_gdagger = np.full(n_draws, np.nan, dtype=float)
    log10_kappa = np.full(n_draws, np.nan, dtype=float)
    bridge_exponent = np.full(n_draws, np.nan, dtype=float)
    if model == "RAR":
        center = parse_float(row["log10_gdagger"])
        log10_gdagger = np.clip(rng.normal(center, args.log10_scale_sd, size=n_draws), -11.6, -9.2)
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA":
        center = parse_float(row["log10_kappa"])
        log10_kappa = np.clip(rng.normal(center, args.log10_scale_sd, size=n_draws), -0.50, 0.50)
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA_P":
        center_k = parse_float(row["log10_kappa"])
        center_p = parse_float(row["bridge_exponent"])
        log10_kappa = np.clip(rng.normal(center_k, args.log10_scale_sd, size=n_draws), -0.50, 0.50)
        bridge_exponent = np.clip(rng.normal(center_p, args.p_sd, size=n_draws), 0.15, 1.20)
    return GlobalDraws(ydisk=ydisk, ybul=ybul, log10_gdagger=log10_gdagger, log10_kappa=log10_kappa, bridge_exponent=bridge_exponent)


def subset_dataset(data: Dataset, mask: np.ndarray) -> Dataset:
    return Dataset(
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


def predict_samples(
    data: Dataset,
    model: str,
    h0: float,
    draws: GlobalDraws,
    logd: np.ndarray,
    logeta: np.ndarray,
) -> np.ndarray:
    n_global = len(draws.ydisk)
    n_nuis = len(logd)
    ydisk = np.repeat(draws.ydisk, n_nuis)
    ybul = np.repeat(draws.ybul, n_nuis)
    log10_gdagger = np.repeat(draws.log10_gdagger, n_nuis)
    log10_kappa = np.repeat(draws.log10_kappa, n_nuis)
    bridge_exponent = np.repeat(draws.bridge_exponent, n_nuis)
    logd_all = np.tile(logd, n_global)
    logeta_all = np.tile(logeta, n_global)

    scale = np.exp(logd_all)[:, None]
    eta = np.exp(logeta_all)[:, None]
    root_scale = np.sqrt(scale)
    root_eta = np.sqrt(eta)
    rad = data.rad_kpc[None, :] * scale
    vgas = data.vgas[None, :] * root_scale
    vdisk = data.vdisk[None, :] * root_scale * root_eta
    vbul = data.vbul[None, :] * root_scale * root_eta

    v2_bar = vgas * np.abs(vgas) + ydisk[:, None] * vdisk**2 + ybul[:, None] * vbul**2
    gbar = np.maximum(v2_bar, 1e-18) * 1.0e6 / (rad * KPC_M)
    if model == "RAR":
        gdagger = (10.0**log10_gdagger)[:, None]
        tau = np.sqrt(np.maximum(gbar / gdagger, 1e-30))
    elif model == "PLAMB_OPTICAL_DEPTH":
        g0 = acceleration_cH0_over_2pi(h0)
        tau = np.sqrt(np.maximum(gbar / g0, 1e-30))
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA":
        g0 = (10.0**log10_kappa)[:, None] * acceleration_cH0_over_2pi(h0)
        tau = np.sqrt(np.maximum(gbar / g0, 1e-30))
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA_P":
        g0 = (10.0**log10_kappa)[:, None] * acceleration_cH0_over_2pi(h0)
        tau = np.maximum(gbar / g0, 1e-30) ** bridge_exponent[:, None]
    else:
        raise ValueError(model)
    gpred = gbar / np.maximum(1.0 - np.exp(-tau), 1e-12)
    return np.sqrt(np.maximum(gpred, 0.0) * rad * KPC_M) / 1000.0


def galaxy_log_predictive(
    data: Dataset,
    model: str,
    h0: float,
    draws: GlobalDraws,
    mu_logd: float,
    sigma_logd: float,
    sigma_ln_ml: float,
    n_nuisance_draws: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    logd = rng.normal(mu_logd, sigma_logd, size=n_nuisance_draws)
    logeta = rng.normal(0.0, sigma_ln_ml, size=n_nuisance_draws)
    pred = predict_samples(data, model, h0, draws, logd, logeta)
    resid = (data.vobs[None, :] - pred) / data.sigma_v[None, :]
    log_norm = np.log(2.0 * math.pi * data.sigma_v[None, :] ** 2)
    sample_loglike = -0.5 * np.sum(resid**2 + log_norm, axis=1)
    return logmeanexp(sample_loglike), float(np.mean(np.sum(resid**2, axis=1)))


def score_split_model(
    split: str,
    model: str,
    prep: Any,
    map_row: dict[str, str],
    args: argparse.Namespace,
    rng: np.random.Generator,
) -> list[dict[str, Any]]:
    h0 = parse_float(map_row["H0_trial"])
    draws = draw_global_params(map_row, model, args.global_draws, rng, args)
    scores: list[dict[str, Any]] = []
    for galaxy_index, galaxy in enumerate(prep.galaxies):
        mask = prep.data.galaxy == galaxy
        gdata = subset_dataset(prep.data, mask)
        mu_logd = float(prep.mu_logd_by_h0[h0][galaxy_index])
        sigma_logd = float(prep.sigma_logd[galaxy_index])
        lp, mean_chi2 = galaxy_log_predictive(
            gdata,
            model,
            h0,
            draws,
            mu_logd,
            sigma_logd,
            prep.sigma_ln_ml,
            args.nuisance_draws,
            rng,
        )
        scores.append(
            {
                "split": split,
                "model": model,
                "galaxy": galaxy,
                "N_points": int(len(gdata.vobs)),
                "log_predictive": lp,
                "mean_predictive_chi2": mean_chi2,
                "log_predictive_per_point": lp / max(len(gdata.vobs), 1),
                "H0_trial": h0,
            }
        )
    return scores


def summarize_scores(galaxy_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    splits = sorted({str(row["split"]) for row in galaxy_rows})
    for split in splits:
        split_rows = [row for row in galaxy_rows if row["split"] == split]
        rar_by_gal = {
            str(row["galaxy"]): row
            for row in split_rows
            if row["model"] == "RAR"
        }
        for model in MODEL_ORDER:
            model_rows = [row for row in split_rows if row["model"] == model]
            if not model_rows:
                continue
            total_lp = float(sum(float(row["log_predictive"]) for row in model_rows))
            total_points = int(sum(int(row["N_points"]) for row in model_rows))
            total_rar = float(sum(float(rar_by_gal[str(row["galaxy"])]["log_predictive"]) for row in model_rows))
            galaxy_deltas = [
                float(row["log_predictive"]) - float(rar_by_gal[str(row["galaxy"])]["log_predictive"])
                for row in model_rows
                if str(row["galaxy"]) in rar_by_gal
            ]
            rows.append(
                {
                    "split": split,
                    "model": model,
                    "N_galaxies": len(model_rows),
                    "N_points": total_points,
                    "total_log_predictive": total_lp,
                    "delta_log_predictive_vs_RAR": total_lp - total_rar,
                    "mean_log_predictive_per_point": total_lp / max(total_points, 1),
                    "median_galaxy_delta_vs_RAR": float(np.median(galaxy_deltas)) if galaxy_deltas else 0.0,
                    "galaxies_better_than_RAR": int(sum(delta > 0.0 for delta in galaxy_deltas)),
                    "galaxies_worse_than_RAR": int(sum(delta < 0.0 for delta in galaxy_deltas)),
                }
            )
    return rows


def add_galaxy_deltas(galaxy_rows: list[dict[str, Any]]) -> None:
    for split in sorted({str(row["split"]) for row in galaxy_rows}):
        rar_by_gal = {
            str(row["galaxy"]): float(row["log_predictive"])
            for row in galaxy_rows
            if row["split"] == split and row["model"] == "RAR"
        }
        for row in galaxy_rows:
            if row["split"] == split:
                row["delta_log_predictive_vs_RAR"] = float(row["log_predictive"]) - rar_by_gal.get(str(row["galaxy"]), float("nan"))


def make_plot(summary_rows: list[dict[str, Any]], path: Path) -> None:
    if plt is None:
        return
    splits = sorted({str(row["split"]) for row in summary_rows})
    models = [model for model in MODEL_ORDER if model != "RAR"]
    x = np.arange(len(splits), dtype=float)
    offsets = np.linspace(-0.24, 0.24, len(models))
    width = 0.18
    colors = {
        "PLAMB_OPTICAL_DEPTH": "#2471a3",
        "PLAMB_OPTICAL_DEPTH_KAPPA": "#229954",
        "PLAMB_OPTICAL_DEPTH_KAPPA_P": "#d68910",
    }
    fig, ax = plt.subplots(figsize=(10.4, 5.8))
    for offset, model in zip(offsets, models):
        vals = []
        for split in splits:
            row = next((r for r in summary_rows if r["split"] == split and r["model"] == model), None)
            vals.append(float(row["delta_log_predictive_vs_RAR"]) if row else float("nan"))
        ax.bar(x + offset, vals, width=width, label=model, color=colors[model])
    ax.axhline(0.0, color="#333333", linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(splits, rotation=20, ha="right")
    ax.set_ylabel("Delta held-out log predictive density vs RAR")
    ax.set_title("SPARC Approximate Posterior-Predictive Galaxy Holdout")
    ax.legend(fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def markdown_table(rows: list[dict[str, Any]], fields: list[str]) -> list[str]:
    lines = ["| " + " | ".join(fields) + " |", "|" + "|".join(["---"] * len(fields)) + "|"]
    for row in rows:
        cells = []
        for field in fields:
            value = row.get(field, "")
            if isinstance(value, float):
                cells.append(f"{value:.6g}")
            else:
                cells.append(str(value))
        lines.append("| " + " | ".join(cells) + " |")
    return lines


def write_report(summary_rows: list[dict[str, Any]], args: argparse.Namespace) -> None:
    lines = [
        "# SPARC Posterior-Predictive Galaxy Holdout",
        "",
        "## Executive Summary",
        "",
        "This diagnostic scores each galaxy as held out, using the full-data hierarchical MAP solution as a local posterior center and marginalizing over new-galaxy distance and stellar M/L nuisance priors. Positive delta log predictive density means the PLAMB branch predicts held-out galaxies better than RAR under this approximation.",
        "",
        "| Split | Model | Delta log predictive vs RAR | Galaxies better | Galaxies worse | Mean log predictive / point |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        if row["model"] == "RAR":
            continue
        lines.append(
            f"| {row['split']} | {row['model']} | {float(row['delta_log_predictive_vs_RAR']):.3f} | "
            f"{row['galaxies_better_than_RAR']} | {row['galaxies_worse_than_RAR']} | "
            f"{float(row['mean_log_predictive_per_point']):.5f} |"
        )
    lines.extend(
        [
            "",
            "## Method Notes",
            "",
            f"- Global draws per model: `{args.global_draws}`",
            f"- New-galaxy nuisance draws per global draw: `{args.nuisance_draws}`",
            f"- Global log-M/L draw scale: `{args.global_log_ml_sd}`",
            f"- log10 scale draw SD: `{args.log10_scale_sd}`",
            f"- bridge exponent draw SD: `{args.p_sd}`",
            "- This is not a full refit-leave-one-galaxy-out MCMC. It is a fast local-posterior predictive test that avoids optimizing the held-out galaxy nuisance parameters on its own data.",
            "",
            "## Outputs",
            "",
            f"- Summary CSV: `{OUTDIR / 'sparc_posterior_predictive_summary.csv'}`",
            f"- Galaxy scores CSV: `{OUTDIR / 'sparc_posterior_predictive_galaxy_scores.csv'}`",
            f"- Plot: `{OUTDIR / 'sparc_posterior_predictive_delta_logp.png'}`",
        ]
    )
    (OUTDIR / "sparc_posterior_predictive_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Approximate posterior-predictive SPARC galaxy holdout validation.")
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--points", type=Path, default=DEFAULT_POINTS)
    parser.add_argument("--map-dir", type=Path, default=DEFAULT_MAP_DIR)
    parser.add_argument("--output-label", default="posterior_predictive_20260714")
    parser.add_argument("--global-draws", type=int, default=128)
    parser.add_argument("--nuisance-draws", type=int, default=32)
    parser.add_argument("--global-log-ml-sd", type=float, default=0.055)
    parser.add_argument("--log10-scale-sd", type=float, default=0.030)
    parser.add_argument("--p-sd", type=float, default=0.030)
    parser.add_argument("--err-floor-kms", type=float, default=3.0)
    parser.add_argument("--sigma-ln-ml", type=float, default=0.25)
    parser.add_argument("--distance-floor-frac", type=float, default=0.03)
    parser.add_argument("--hubble-prior-mode", choices=["model_h0_rescale", "published"], default="model_h0_rescale")
    parser.add_argument("--seed", type=int, default=260714)
    parser.add_argument("--only", default="", help="Optional split name.")
    args = parser.parse_args()

    if args.output_label:
        safe = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in args.output_label)
        globals()["OUTDIR"] = OUTDIR / safe

    rng = np.random.default_rng(args.seed)
    map_rows = load_map_rows(args.map_dir)
    galaxy_rows: list[dict[str, Any]] = []
    for config in split_configs():
        split = str(config["split"])
        if args.only and split != args.only:
            continue
        h0_values = sorted(
            {
                parse_float(map_rows[(split, model)]["H0_trial"])
                for model in MODEL_ORDER
                if (split, model) in map_rows
            }
        )
        prep = prepare(
            args.sample,
            args.points,
            int(config["quality_max"]),
            str(config["distance_method"]),
            args.err_floor_kms,
            h0_values,
            args.sigma_ln_ml,
            args.distance_floor_frac,
            args.hubble_prior_mode,
        )
        print(f"=== {split} galaxies={len(prep.galaxies)} points={len(prep.data.vobs)} ===")
        for model in MODEL_ORDER:
            row = map_rows[(split, model)]
            scores = score_split_model(split, model, prep, row, args, rng)
            galaxy_rows.extend(scores)
            print(f"{split} {model} total_logp={sum(s['log_predictive'] for s in scores):.4g}")

    add_galaxy_deltas(galaxy_rows)
    summary_rows = summarize_scores(galaxy_rows)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    write_csv(
        OUTDIR / "sparc_posterior_predictive_galaxy_scores.csv",
        galaxy_rows,
        [
            "split",
            "model",
            "galaxy",
            "N_points",
            "log_predictive",
            "delta_log_predictive_vs_RAR",
            "mean_predictive_chi2",
            "log_predictive_per_point",
            "H0_trial",
        ],
    )
    write_csv(
        OUTDIR / "sparc_posterior_predictive_summary.csv",
        summary_rows,
        [
            "split",
            "model",
            "N_galaxies",
            "N_points",
            "total_log_predictive",
            "delta_log_predictive_vs_RAR",
            "mean_log_predictive_per_point",
            "median_galaxy_delta_vs_RAR",
            "galaxies_better_than_RAR",
            "galaxies_worse_than_RAR",
        ],
    )
    make_plot(summary_rows, OUTDIR / "sparc_posterior_predictive_delta_logp.png")
    write_report(summary_rows, args)
    (OUTDIR / "sparc_posterior_predictive_metadata.json").write_text(json.dumps({"args": vars(args)}, indent=2, default=str), encoding="utf-8")
    print(f"Saved report: {OUTDIR / 'sparc_posterior_predictive_report.md'}")


if __name__ == "__main__":
    main()
