#!/usr/bin/env python3
r"""
Lightweight null controls and split stress tests for SPARC PLAMB diagnostics.

This script is intentionally cheap enough to run while the overnight posterior
sampler is active. It reuses the MAP-global posterior-predictive scoring idea:
global parameters are fixed at the hierarchical MAP rows, while held-out galaxy
distance and stellar M/L nuisance priors are marginalized by Monte Carlo draws.

Controls:
  - scramble galaxy labels for distance-prior centers/sigmas,
  - blind Hubble-flow prior metadata by using published-distance priors,
  - shuffle distance-method labels and rescore pseudo Hubble/non-Hubble splits.

Stress tests:
  - gas-dominated galaxies,
  - bulge-dominated galaxies removed,
  - low-acceleration outer points,
  - high- and low-inclination galaxies,
  - high-quality Q1 only.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

from diagnose_plamb_sparc_h0_nuisance import acceleration_cH0_over_2pi
from diagnose_plamb_sparc_rotation import DEFAULT_POINTS, DEFAULT_SAMPLE, KPC_M, Dataset
from fit_sparc_hierarchical_map import prepare, split_configs
from validate_sparc_hierarchical_posterior_predictive import (
    MODEL_ORDER,
    draw_global_params,
    galaxy_log_predictive,
    load_map_rows,
    parse_float,
    subset_dataset,
    write_csv,
)


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "sparc_null_stress_controls_20260714"
DEFAULT_MAP_DIR = ROOT / "plamb_runs" / "diagnostics" / "sparc_hierarchical_map" / "optical_depth_hierarchical_20260714"
MAIN_MODELS = ["RAR", "PLAMB_OPTICAL_DEPTH_KAPPA_P"]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_sample_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def split_config(name: str) -> dict[str, Any]:
    return next(item for item in split_configs() if item["split"] == name)


def h0_values_for_models(map_rows: dict[tuple[str, str], dict[str, str]], split: str, models: list[str]) -> list[float]:
    return sorted({parse_float(map_rows[(split, model)]["H0_trial"]) for model in models})


def baryon_components(data: Dataset, ydisk: float = 0.5, ybul: float = 0.7) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    gas = data.vgas * np.abs(data.vgas)
    disk = ydisk * data.vdisk**2
    bulge = ybul * data.vbul**2
    total = np.maximum(gas + disk + bulge, 1e-18)
    return gas, disk, bulge, total


def gbar_fiducial(data: Dataset) -> np.ndarray:
    gas, disk, bulge, _total = baryon_components(data)
    v2 = gas + disk + bulge
    return np.maximum(v2, 1e-18) * 1.0e6 / (data.rad_kpc * KPC_M)


def galaxy_metric(data: Dataset, galaxies: list[str], values: np.ndarray, reducer: str = "median") -> dict[str, float]:
    result: dict[str, float] = {}
    for galaxy in galaxies:
        vals = values[data.galaxy == galaxy]
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            result[galaxy] = float("nan")
        elif reducer == "max":
            result[galaxy] = float(np.max(vals))
        else:
            result[galaxy] = float(np.median(vals))
    return result


def filter_prepared(prep: Any, galaxy_keep: set[str] | None = None, point_mask: np.ndarray | None = None) -> Any:
    mask = np.ones(len(prep.data.vobs), dtype=bool)
    if galaxy_keep is not None:
        mask &= np.asarray([str(galaxy) in galaxy_keep for galaxy in prep.data.galaxy], dtype=bool)
    if point_mask is not None:
        mask &= point_mask
    data = subset_dataset(prep.data, mask)
    kept = [galaxy for galaxy in prep.galaxies if np.any(data.galaxy == galaxy)]
    idx = {galaxy: i for i, galaxy in enumerate(prep.galaxies)}
    kept_idx = np.asarray([idx[galaxy] for galaxy in kept], dtype=int)
    return SimpleNamespace(
        data=data,
        galaxies=kept,
        sigma_logd=prep.sigma_logd[kept_idx],
        mu_logd_by_h0={h0: values[kept_idx] for h0, values in prep.mu_logd_by_h0.items()},
        sigma_ln_ml=prep.sigma_ln_ml,
    )


def make_shuffled_distance_sample(sample: Path, outdir: Path, rng: np.random.Generator) -> Path:
    rows = read_csv(sample)
    fields = list(rows[0].keys())
    labels = [row["f_D"] for row in rows]
    shuffled = labels[:]
    rng.shuffle(shuffled)
    new_rows = [dict(row) for row in rows]
    for row, label in zip(new_rows, shuffled):
        row["f_D"] = label
    path = outdir / "scratch" / "sparc_galaxy_sample_distance_labels_shuffled.csv"
    write_sample_csv(path, new_rows, fields)
    return path


def score_case(
    case: str,
    split_for_map: str,
    prep: Any,
    map_rows: dict[tuple[str, str], dict[str, str]],
    models: list[str],
    args: argparse.Namespace,
    rng: np.random.Generator,
    prior_perm: np.ndarray | None = None,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if len(prep.galaxies) == 0 or len(prep.data.vobs) == 0:
        return rows
    for model in models:
        map_row = map_rows[(split_for_map, model)]
        h0 = parse_float(map_row["H0_trial"])
        draws = draw_global_params(map_row, model, args.global_draws, rng, args)
        for galaxy_index, galaxy in enumerate(prep.galaxies):
            mask = prep.data.galaxy == galaxy
            gdata = subset_dataset(prep.data, mask)
            prior_index = int(prior_perm[galaxy_index]) if prior_perm is not None else galaxy_index
            mu_logd = float(prep.mu_logd_by_h0[h0][prior_index])
            sigma_logd = float(prep.sigma_logd[prior_index])
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
            rows.append(
                {
                    "case": case,
                    "split_for_map": split_for_map,
                    "model": model,
                    "galaxy": galaxy,
                    "N_points": int(len(gdata.vobs)),
                    "log_predictive": lp,
                    "mean_predictive_chi2": mean_chi2,
                    "log_predictive_per_point": lp / max(len(gdata.vobs), 1),
                }
            )
    return rows


def add_case_deltas(rows: list[dict[str, Any]]) -> None:
    for case in sorted({row["case"] for row in rows}):
        case_rows = [row for row in rows if row["case"] == case]
        rar = {row["galaxy"]: float(row["log_predictive"]) for row in case_rows if row["model"] == "RAR"}
        for row in case_rows:
            row["delta_log_predictive_vs_RAR"] = float(row["log_predictive"]) - rar.get(row["galaxy"], float("nan"))


def summarize(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    summary: list[dict[str, Any]] = []
    for case in sorted({row["case"] for row in rows}):
        case_rows = [row for row in rows if row["case"] == case]
        rar_total = sum(float(row["log_predictive"]) for row in case_rows if row["model"] == "RAR")
        for model in sorted({row["model"] for row in case_rows}, key=lambda name: MODEL_ORDER.index(name) if name in MODEL_ORDER else 999):
            model_rows = [row for row in case_rows if row["model"] == model]
            total_lp = float(sum(float(row["log_predictive"]) for row in model_rows))
            deltas = [float(row.get("delta_log_predictive_vs_RAR", 0.0)) for row in model_rows if model != "RAR"]
            summary.append(
                {
                    "case": case,
                    "model": model,
                    "N_galaxies": len(model_rows),
                    "N_points": int(sum(int(row["N_points"]) for row in model_rows)),
                    "total_log_predictive": total_lp,
                    "delta_log_predictive_vs_RAR": total_lp - rar_total,
                    "median_galaxy_delta_vs_RAR": float(np.median(deltas)) if deltas else 0.0,
                    "galaxies_better_than_RAR": int(sum(delta > 0.0 for delta in deltas)),
                    "galaxies_worse_than_RAR": int(sum(delta < 0.0 for delta in deltas)),
                }
            )
    return summary


def run_controls(map_rows: dict[tuple[str, str], dict[str, str]], args: argparse.Namespace, rng: np.random.Generator) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    all_rows: list[dict[str, Any]] = []
    metadata: list[dict[str, Any]] = []
    models = MAIN_MODELS
    for split in ["all_Q2", "hubble_flow_Q2", "non_hubble_Q2"]:
        cfg = split_config(split)
        prep = prepare(
            args.sample,
            args.points,
            int(cfg["quality_max"]),
            str(cfg["distance_method"]),
            args.err_floor_kms,
            h0_values_for_models(map_rows, split, models),
            args.sigma_ln_ml,
            args.distance_floor_frac,
            "model_h0_rescale",
        )
        all_rows.extend(score_case(f"baseline_{split}", split, prep, map_rows, models, args, rng))
        perm = rng.permutation(len(prep.galaxies))
        all_rows.extend(score_case(f"scrambled_nuisance_priors_{split}", split, prep, map_rows, models, args, rng, prior_perm=perm))
        metadata.append({"case": f"scrambled_nuisance_priors_{split}", "note": "distance-prior centers and sigmas permuted among galaxies"})

        prep_blind = prepare(
            args.sample,
            args.points,
            int(cfg["quality_max"]),
            str(cfg["distance_method"]),
            args.err_floor_kms,
            h0_values_for_models(map_rows, split, models),
            args.sigma_ln_ml,
            args.distance_floor_frac,
            "published",
        )
        all_rows.extend(score_case(f"hubble_prior_blinded_{split}", split, prep_blind, map_rows, models, args, rng))
        metadata.append({"case": f"hubble_prior_blinded_{split}", "note": "Hubble-flow H0 rescale removed from distance priors"})

    shuffled_sample = make_shuffled_distance_sample(args.sample, OUTDIR, rng)
    for split in ["hubble_flow_Q2", "non_hubble_Q2"]:
        cfg = split_config(split)
        prep = prepare(
            shuffled_sample,
            args.points,
            int(cfg["quality_max"]),
            str(cfg["distance_method"]),
            args.err_floor_kms,
            h0_values_for_models(map_rows, split, models),
            args.sigma_ln_ml,
            args.distance_floor_frac,
            "model_h0_rescale",
        )
        all_rows.extend(score_case(f"distance_labels_shuffled_{split}", split, prep, map_rows, models, args, rng))
        metadata.append({"case": f"distance_labels_shuffled_{split}", "note": "f_D labels shuffled before selecting split"})
    return all_rows, metadata


def run_stress_tests(map_rows: dict[tuple[str, str], dict[str, str]], args: argparse.Namespace, rng: np.random.Generator) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    metadata: list[dict[str, Any]] = []
    models = MAIN_MODELS
    split = "all_Q2"
    cfg = split_config(split)
    prep = prepare(
        args.sample,
        args.points,
        int(cfg["quality_max"]),
        str(cfg["distance_method"]),
        args.err_floor_kms,
        h0_values_for_models(map_rows, split, models),
        args.sigma_ln_ml,
        args.distance_floor_frac,
        "model_h0_rescale",
    )
    gas, disk, bulge, total = baryon_components(prep.data)
    gas_frac = gas / total
    bulge_frac = bulge / total
    gas_by_gal = galaxy_metric(prep.data, prep.galaxies, gas_frac)
    bulge_by_gal = galaxy_metric(prep.data, prep.galaxies, bulge_frac)
    gas_keep = {galaxy for galaxy, value in gas_by_gal.items() if math.isfinite(value) and value >= args.gas_fraction_min}
    nobulge_keep = {galaxy for galaxy, value in bulge_by_gal.items() if not math.isfinite(value) or value <= args.bulge_fraction_max}
    rows.extend(score_case("stress_gas_dominated_all_Q2", split, filter_prepared(prep, gas_keep), map_rows, models, args, rng))
    rows.extend(score_case("stress_bulge_removed_all_Q2", split, filter_prepared(prep, nobulge_keep), map_rows, models, args, rng))
    metadata.extend(
        [
            {"case": "stress_gas_dominated_all_Q2", "note": f"median gas baryon fraction >= {args.gas_fraction_min}"},
            {"case": "stress_bulge_removed_all_Q2", "note": f"median bulge baryon fraction <= {args.bulge_fraction_max}"},
        ]
    )

    gbar = gbar_fiducial(prep.data)
    g0 = acceleration_cH0_over_2pi(65.0)
    median_rad = galaxy_metric(prep.data, prep.galaxies, prep.data.rad_kpc)
    outer = np.asarray([prep.data.rad_kpc[i] >= median_rad[str(galaxy)] for i, galaxy in enumerate(prep.data.galaxy)], dtype=bool)
    low_accel = gbar < g0
    low_outer_mask = outer & low_accel
    low_outer_prep = filter_prepared(prep, point_mask=low_outer_mask)
    rows.extend(score_case("stress_low_accel_outer_points_all_Q2", split, low_outer_prep, map_rows, models, args, rng))
    metadata.append({"case": "stress_low_accel_outer_points_all_Q2", "note": "point-level subset: radius above galaxy median and gbar < cH0/2pi at H0=65"})

    inc_by_gal = galaxy_metric(prep.data, prep.galaxies, prep.data.inc_deg)
    high_inc = {galaxy for galaxy, value in inc_by_gal.items() if math.isfinite(value) and value >= args.high_inc_min}
    low_inc = {galaxy for galaxy, value in inc_by_gal.items() if math.isfinite(value) and value < args.high_inc_min}
    rows.extend(score_case("stress_high_inclination_all_Q2", split, filter_prepared(prep, high_inc), map_rows, models, args, rng))
    rows.extend(score_case("stress_low_inclination_all_Q2", split, filter_prepared(prep, low_inc), map_rows, models, args, rng))
    metadata.extend(
        [
            {"case": "stress_high_inclination_all_Q2", "note": f"median inclination >= {args.high_inc_min} deg"},
            {"case": "stress_low_inclination_all_Q2", "note": f"median inclination < {args.high_inc_min} deg"},
        ]
    )

    prep_q1 = prepare(
        args.sample,
        args.points,
        1,
        "all",
        args.err_floor_kms,
        h0_values_for_models(map_rows, split, models),
        args.sigma_ln_ml,
        args.distance_floor_frac,
        "model_h0_rescale",
    )
    rows.extend(score_case("stress_high_quality_Q1_all", split, prep_q1, map_rows, models, args, rng))
    metadata.append({"case": "stress_high_quality_Q1_all", "note": "Q=1 galaxies only; all distance methods"})
    return rows, metadata


def write_report(summary_rows: list[dict[str, Any]], metadata: list[dict[str, Any]], args: argparse.Namespace) -> None:
    by_case_note = {row["case"]: row["note"] for row in metadata}
    lines = [
        "# SPARC Null Controls and Stress Tests",
        "",
        "Positive delta log predictive density means PLAMB_OPTICAL_DEPTH_KAPPA_P beats RAR under MAP-global predictive scoring.",
        "",
        "| Case | N galaxies | N points | Delta log predictive | Median galaxy delta | Better / worse | Note |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in summary_rows:
        if row["model"] != "PLAMB_OPTICAL_DEPTH_KAPPA_P":
            continue
        lines.append(
            f"| {row['case']} | {row['N_galaxies']} | {row['N_points']} | "
            f"{float(row['delta_log_predictive_vs_RAR']):.3f} | {float(row['median_galaxy_delta_vs_RAR']):.3f} | "
            f"{row['galaxies_better_than_RAR']} / {row['galaxies_worse_than_RAR']} | {by_case_note.get(row['case'], '')} |"
        )
    lines.extend(
        [
            "",
            "## Method",
            "",
            f"- nuisance draws per galaxy: `{args.nuisance_draws}`",
            "- global parameters fixed at the hierarchical MAP values",
            "- models compared: `RAR` and `PLAMB_OPTICAL_DEPTH_KAPPA_P`",
            "",
            "## Outputs",
            "",
            f"- summary CSV: `{OUTDIR / 'sparc_null_stress_summary.csv'}`",
            f"- galaxy score CSV: `{OUTDIR / 'sparc_null_stress_galaxy_scores.csv'}`",
        ]
    )
    (OUTDIR / "sparc_null_stress_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run lightweight SPARC null controls and split stress tests.")
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--points", type=Path, default=DEFAULT_POINTS)
    parser.add_argument("--map-dir", type=Path, default=DEFAULT_MAP_DIR)
    parser.add_argument("--nuisance-draws", type=int, default=384)
    parser.add_argument("--global-draws", type=int, default=1)
    parser.add_argument("--global-log-ml-sd", type=float, default=0.0)
    parser.add_argument("--log10-scale-sd", type=float, default=0.0)
    parser.add_argument("--p-sd", type=float, default=0.0)
    parser.add_argument("--err-floor-kms", type=float, default=3.0)
    parser.add_argument("--sigma-ln-ml", type=float, default=0.25)
    parser.add_argument("--distance-floor-frac", type=float, default=0.03)
    parser.add_argument("--gas-fraction-min", type=float, default=0.5)
    parser.add_argument("--bulge-fraction-max", type=float, default=0.2)
    parser.add_argument("--high-inc-min", type=float, default=60.0)
    parser.add_argument("--seed", type=int, default=260714)
    args = parser.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    map_rows = load_map_rows(args.map_dir)
    control_rows, control_meta = run_controls(map_rows, args, rng)
    stress_rows, stress_meta = run_stress_tests(map_rows, args, rng)
    rows = control_rows + stress_rows
    add_case_deltas(rows)
    summary_rows = summarize(rows)
    write_csv(
        OUTDIR / "sparc_null_stress_galaxy_scores.csv",
        rows,
        [
            "case",
            "split_for_map",
            "model",
            "galaxy",
            "N_points",
            "log_predictive",
            "delta_log_predictive_vs_RAR",
            "mean_predictive_chi2",
            "log_predictive_per_point",
        ],
    )
    write_csv(
        OUTDIR / "sparc_null_stress_summary.csv",
        summary_rows,
        [
            "case",
            "model",
            "N_galaxies",
            "N_points",
            "total_log_predictive",
            "delta_log_predictive_vs_RAR",
            "median_galaxy_delta_vs_RAR",
            "galaxies_better_than_RAR",
            "galaxies_worse_than_RAR",
        ],
    )
    metadata = control_meta + stress_meta
    (OUTDIR / "sparc_null_stress_metadata.json").write_text(json.dumps({"args": vars(args), "case_notes": metadata}, indent=2, default=str), encoding="utf-8")
    write_report(summary_rows, metadata, args)
    print(f"Saved report: {OUTDIR / 'sparc_null_stress_report.md'}")


if __name__ == "__main__":
    main()
