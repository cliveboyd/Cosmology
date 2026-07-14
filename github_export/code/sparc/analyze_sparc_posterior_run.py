#!/usr/bin/env python3
r"""
Analyze completed SPARC hierarchical posterior sampler runs.

The overnight sampler writes one case directory per split/model:

    <run_dir>/<split>/<model>/posterior_samples.npz

This analyzer turns those chains into immediately usable science products:

  - chain diagnostics and Rhat for global parameters,
  - posterior checks for PLAMB kappa and bridge exponent p,
  - nuisance-pull summaries for galaxy distance and stellar M/L shifts,
  - posterior-predictive RAR-vs-PLAMB scores using actual chain samples,
  - stress-subset predictive scores for Q1, gas-dominated, bulge-removed,
    low-acceleration outer points, and inclination splits.

It is designed to run after the overnight job completes, but the smoke-test
run can be used to validate the analyzer itself.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
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
from fit_sparc_hierarchical_map import prepare, split_configs
from sample_sparc_hierarchical_posterior import DEFAULT_MAP_DIR, OUTDIR as POSTERIOR_BASE_OUTDIR
from validate_sparc_hierarchical_posterior_predictive import load_map_rows, parse_float


DEFAULT_RUN_DIR = POSTERIOR_BASE_OUTDIR / "posterior_overnight_20260714"
MODEL_ORDER = ["RAR", "PLAMB_OPTICAL_DEPTH", "PLAMB_OPTICAL_DEPTH_KAPPA", "PLAMB_OPTICAL_DEPTH_KAPPA_P"]


@dataclass
class Case:
    split: str
    model: str
    path: Path
    global_samples: np.ndarray
    nuisance_samples: np.ndarray
    logpost_samples: np.ndarray
    param_names: list[str]
    galaxies: list[str]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def logmeanexp(values: np.ndarray) -> float:
    vmax = float(np.max(values))
    return vmax + math.log(float(np.mean(np.exp(values - vmax))))


def split_config(split: str) -> dict[str, Any]:
    return next(item for item in split_configs() if item["split"] == split)


def discover_cases(run_dir: Path) -> list[Case]:
    cases: list[Case] = []
    for path in sorted(run_dir.glob("*/*/posterior_samples.npz")):
        split = path.parent.parent.name
        model = path.parent.name
        z = np.load(path, allow_pickle=True)
        cases.append(
            Case(
                split=split,
                model=model,
                path=path.parent,
                global_samples=np.asarray(z["global_samples"], dtype=float),
                nuisance_samples=np.asarray(z["nuisance_samples"], dtype=float),
                logpost_samples=np.asarray(z["logpost_samples"], dtype=float),
                param_names=[str(item) for item in z["param_names"].tolist()],
                galaxies=[str(item) for item in z["galaxies"].tolist()],
            )
        )
    return cases


def summarize_vector(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return {"mean": np.nan, "sd": np.nan, "q05": np.nan, "q16": np.nan, "median": np.nan, "q84": np.nan, "q95": np.nan}
    return {
        "mean": float(np.mean(values)),
        "sd": float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
        "q05": float(np.quantile(values, 0.05)),
        "q16": float(np.quantile(values, 0.16)),
        "median": float(np.quantile(values, 0.50)),
        "q84": float(np.quantile(values, 0.84)),
        "q95": float(np.quantile(values, 0.95)),
    }


def rhat(chains: np.ndarray) -> np.ndarray:
    if chains.ndim != 3 or chains.shape[0] < 2 or chains.shape[1] < 2:
        return np.full(chains.shape[-1] if chains.ndim == 3 else 0, np.nan)
    m, n, _p = chains.shape
    chain_means = np.mean(chains, axis=1)
    chain_vars = np.var(chains, axis=1, ddof=1)
    w = np.mean(chain_vars, axis=0)
    b = n * np.var(chain_means, axis=0, ddof=1)
    var_hat = ((n - 1) / n) * w + b / n
    return np.sqrt(np.divide(var_hat, w, out=np.full_like(w, np.nan), where=w > 0))


def load_chain_diagnostics(case: Case) -> list[dict[str, Any]]:
    path = case.path / "chain_diagnostics.csv"
    if not path.exists():
        return []
    return read_csv(path)


def global_parameter_rows(case: Case) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    flat = case.global_samples.reshape(-1, case.global_samples.shape[-1])
    rh = rhat(case.global_samples)
    for i, name in enumerate(case.param_names):
        stats = summarize_vector(flat[:, i])
        rows.append({"split": case.split, "model": case.model, "parameter": name, "rhat": float(rh[i]) if len(rh) else "", **stats})
        if name == "log10_kappa":
            kappa_stats = summarize_vector(10.0 ** flat[:, i])
            rows.append({"split": case.split, "model": case.model, "parameter": "kappa", "rhat": "", **kappa_stats})
        if name == "log_ydisk":
            rows.append({"split": case.split, "model": case.model, "parameter": "ydisk", "rhat": "", **summarize_vector(np.exp(flat[:, i]))})
        if name == "log_ybul":
            rows.append({"split": case.split, "model": case.model, "parameter": "ybul", "rhat": "", **summarize_vector(np.exp(flat[:, i]))})
    logpost_stats = summarize_vector(case.logpost_samples.reshape(-1))
    rows.append({"split": case.split, "model": case.model, "parameter": "logpost", "rhat": "", **logpost_stats})
    return rows


def case_diagnostic_row(case: Case, param_rows: list[dict[str, Any]]) -> dict[str, Any]:
    diag = load_chain_diagnostics(case)
    global_accept = [parse_float(row.get("accept_global")) for row in diag]
    galaxy_accept = [parse_float(row.get("accept_galaxy")) for row in diag]
    rhat_values = [float(row["rhat"]) for row in param_rows if row["split"] == case.split and row["model"] == case.model and row.get("rhat", "") not in ("", None) and math.isfinite(float(row["rhat"]))]
    return {
        "split": case.split,
        "model": case.model,
        "n_chains": int(case.global_samples.shape[0]),
        "n_saved_per_chain": int(case.global_samples.shape[1]),
        "n_global_params": int(case.global_samples.shape[2]),
        "n_galaxies": len(case.galaxies),
        "mean_accept_global": float(np.nanmean(global_accept)) if global_accept else "",
        "mean_accept_galaxy": float(np.nanmean(galaxy_accept)) if galaxy_accept else "",
        "max_global_rhat": float(np.nanmax(rhat_values)) if rhat_values else "",
        "median_global_rhat": float(np.nanmedian(rhat_values)) if rhat_values else "",
        "complete": bool(case.global_samples.shape[1] > 0),
    }


def posterior_checks(case: Case, param_rows: list[dict[str, Any]]) -> dict[str, Any]:
    by_param = {
        row["parameter"]: row
        for row in param_rows
        if row["split"] == case.split and row["model"] == case.model
    }
    row: dict[str, Any] = {
        "split": case.split,
        "model": case.model,
        "p_median": "",
        "p_q16": "",
        "p_q84": "",
        "p_q05": "",
        "p_q95": "",
        "p_half_inside_68": "",
        "p_half_inside_90": "",
        "kappa_median": "",
        "kappa_q16": "",
        "kappa_q84": "",
        "kappa_q05": "",
        "kappa_q95": "",
        "kappa_one_inside_68": "",
        "kappa_one_inside_90": "",
    }
    if "bridge_exponent" in by_param:
        p = by_param["bridge_exponent"]
        row.update(
            {
                "p_median": p["median"],
                "p_q16": p["q16"],
                "p_q84": p["q84"],
                "p_q05": p["q05"],
                "p_q95": p["q95"],
                "p_half_inside_68": float(p["q16"]) <= 0.5 <= float(p["q84"]),
                "p_half_inside_90": float(p["q05"]) <= 0.5 <= float(p["q95"]),
            }
        )
    if "kappa" in by_param:
        k = by_param["kappa"]
        row.update(
            {
                "kappa_median": k["median"],
                "kappa_q16": k["q16"],
                "kappa_q84": k["q84"],
                "kappa_q05": k["q05"],
                "kappa_q95": k["q95"],
                "kappa_one_inside_68": float(k["q16"]) <= 1.0 <= float(k["q84"]),
                "kappa_one_inside_90": float(k["q05"]) <= 1.0 <= float(k["q95"]),
            }
        )
    return row


def prepare_for_case(case: Case, map_rows: dict[tuple[str, str], dict[str, str]], sample: Path, points: Path, args: argparse.Namespace) -> Any:
    cfg = split_config(case.split)
    h0 = parse_float(map_rows[(case.split, case.model)]["H0_trial"])
    return prepare(
        sample,
        points,
        int(cfg["quality_max"]),
        str(cfg["distance_method"]),
        args.err_floor_kms,
        [h0],
        args.sigma_ln_ml,
        args.distance_floor_frac,
        args.hubble_prior_mode,
    )


def nuisance_pull_rows(case: Case, prep: Any, map_rows: dict[tuple[str, str], dict[str, str]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    h0 = parse_float(map_rows[(case.split, case.model)]["H0_trial"])
    nuisance = case.nuisance_samples.reshape(-1, case.nuisance_samples.shape[-1])
    n_gal = len(case.galaxies)
    mu_logd = prep.mu_logd_by_h0[h0]
    rows: list[dict[str, Any]] = []
    d_pulls: list[float] = []
    ml_pulls: list[float] = []
    for i, galaxy in enumerate(case.galaxies):
        logd = nuisance[:, i]
        logeta = nuisance[:, n_gal + i]
        d_med = float(np.median(logd))
        e_med = float(np.median(logeta))
        d_pull = (d_med - float(mu_logd[i])) / float(prep.sigma_logd[i])
        ml_pull = e_med / float(prep.sigma_ln_ml)
        d_pulls.append(abs(d_pull))
        ml_pulls.append(abs(ml_pull))
        rows.append(
            {
                "split": case.split,
                "model": case.model,
                "galaxy": galaxy,
                "distance_pull_median": d_pull,
                "stellar_ml_pull_median": ml_pull,
                "distance_multiplier_median": math.exp(d_med),
                "stellar_ml_multiplier_median": math.exp(e_med),
            }
        )
    summary = {
        "split": case.split,
        "model": case.model,
        "n_galaxies": n_gal,
        "max_abs_distance_pull": float(np.max(d_pulls)) if d_pulls else "",
        "median_abs_distance_pull": float(np.median(d_pulls)) if d_pulls else "",
        "distance_pull_gt2": int(sum(value > 2 for value in d_pulls)),
        "distance_pull_gt3": int(sum(value > 3 for value in d_pulls)),
        "max_abs_stellar_ml_pull": float(np.max(ml_pulls)) if ml_pulls else "",
        "median_abs_stellar_ml_pull": float(np.median(ml_pulls)) if ml_pulls else "",
        "stellar_ml_pull_gt2": int(sum(value > 2 for value in ml_pulls)),
        "stellar_ml_pull_gt3": int(sum(value > 3 for value in ml_pulls)),
    }
    return rows, summary


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


def galaxy_metric(data: Dataset, galaxies: list[str], values: np.ndarray) -> dict[str, float]:
    out: dict[str, float] = {}
    for galaxy in galaxies:
        vals = values[data.galaxy == galaxy]
        vals = vals[np.isfinite(vals)]
        out[galaxy] = float(np.median(vals)) if len(vals) else float("nan")
    return out


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


def model_globals_from_samples(case: Case, max_draws: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
    flat = case.global_samples.reshape(-1, case.global_samples.shape[-1])
    if len(flat) > max_draws:
        idx = rng.choice(len(flat), size=max_draws, replace=False)
        flat = flat[idx]
    values: dict[str, np.ndarray] = {name: flat[:, i] for i, name in enumerate(case.param_names)}
    values["ydisk"] = np.exp(values["log_ydisk"])
    values["ybul"] = np.exp(values["log_ybul"])
    return values


def predict_from_global_arrays(
    data: Dataset,
    model: str,
    h0: float,
    glob: dict[str, np.ndarray],
    logd: np.ndarray,
    logeta: np.ndarray,
) -> np.ndarray:
    n_global = len(glob["ydisk"])
    n_nuis = len(logd)
    ydisk = np.repeat(glob["ydisk"], n_nuis)
    ybul = np.repeat(glob["ybul"], n_nuis)
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
        gdagger = np.repeat(10.0 ** glob["log10_gdagger"], n_nuis)[:, None]
        tau = np.sqrt(np.maximum(gbar / gdagger, 1e-30))
    elif model == "PLAMB_OPTICAL_DEPTH":
        tau = np.sqrt(np.maximum(gbar / acceleration_cH0_over_2pi(h0), 1e-30))
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA":
        g0 = np.repeat(10.0 ** glob["log10_kappa"], n_nuis)[:, None] * acceleration_cH0_over_2pi(h0)
        tau = np.sqrt(np.maximum(gbar / g0, 1e-30))
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA_P":
        g0 = np.repeat(10.0 ** glob["log10_kappa"], n_nuis)[:, None] * acceleration_cH0_over_2pi(h0)
        p = np.repeat(glob["bridge_exponent"], n_nuis)[:, None]
        tau = np.maximum(gbar / g0, 1e-30) ** p
    else:
        raise ValueError(model)
    gpred = gbar / np.maximum(1.0 - np.exp(-tau), 1e-12)
    return np.sqrt(np.maximum(gpred, 0.0) * rad * KPC_M) / 1000.0


def score_galaxy_predictive(
    data: Dataset,
    model: str,
    h0: float,
    glob: dict[str, np.ndarray],
    mu_logd: float,
    sigma_logd: float,
    sigma_ln_ml: float,
    nuisance_draws: int,
    rng: np.random.Generator,
) -> float:
    logd = rng.normal(mu_logd, sigma_logd, size=nuisance_draws)
    logeta = rng.normal(0.0, sigma_ln_ml, size=nuisance_draws)
    pred = predict_from_global_arrays(data, model, h0, glob, logd, logeta)
    resid = (data.vobs[None, :] - pred) / data.sigma_v[None, :]
    log_norm = np.log(2.0 * math.pi * data.sigma_v[None, :] ** 2)
    sample_loglike = -0.5 * np.sum(resid**2 + log_norm, axis=1)
    return logmeanexp(sample_loglike)


def predictive_rows_for_case(
    score_case: str,
    case: Case,
    prep: Any,
    map_rows: dict[tuple[str, str], dict[str, str]],
    args: argparse.Namespace,
    rng: np.random.Generator,
) -> list[dict[str, Any]]:
    h0 = parse_float(map_rows[(case.split, case.model)]["H0_trial"])
    glob = model_globals_from_samples(case, args.max_global_draws, rng)
    rows: list[dict[str, Any]] = []
    for galaxy_index, galaxy in enumerate(prep.galaxies):
        mask = prep.data.galaxy == galaxy
        data = subset_dataset(prep.data, mask)
        lp = score_galaxy_predictive(
            data,
            case.model,
            h0,
            glob,
            float(prep.mu_logd_by_h0[h0][galaxy_index]),
            float(prep.sigma_logd[galaxy_index]),
            float(prep.sigma_ln_ml),
            args.nuisance_draws,
            rng,
        )
        rows.append(
            {
                "score_case": score_case,
                "split": case.split,
                "model": case.model,
                "galaxy": galaxy,
                "N_points": int(len(data.vobs)),
                "log_predictive": lp,
            }
        )
    return rows


def add_predictive_deltas(rows: list[dict[str, Any]]) -> None:
    for score_case in sorted({row["score_case"] for row in rows}):
        case_rows = [row for row in rows if row["score_case"] == score_case]
        rar = {row["galaxy"]: float(row["log_predictive"]) for row in case_rows if row["model"] == "RAR"}
        for row in case_rows:
            row["delta_log_predictive_vs_RAR"] = float(row["log_predictive"]) - rar.get(row["galaxy"], float("nan"))


def summarize_predictive(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for score_case in sorted({row["score_case"] for row in rows}):
        case_rows = [row for row in rows if row["score_case"] == score_case]
        rar_total = sum(float(row["log_predictive"]) for row in case_rows if row["model"] == "RAR")
        for model in sorted({row["model"] for row in case_rows}, key=lambda x: MODEL_ORDER.index(x) if x in MODEL_ORDER else 999):
            mrows = [row for row in case_rows if row["model"] == model]
            total = sum(float(row["log_predictive"]) for row in mrows)
            deltas = [float(row.get("delta_log_predictive_vs_RAR", np.nan)) for row in mrows if model != "RAR"]
            deltas = [delta for delta in deltas if math.isfinite(delta)]
            out.append(
                {
                    "score_case": score_case,
                    "model": model,
                    "N_galaxies": len(mrows),
                    "N_points": int(sum(int(row["N_points"]) for row in mrows)),
                    "total_log_predictive": total,
                    "delta_log_predictive_vs_RAR": total - rar_total,
                    "median_galaxy_delta_vs_RAR": float(np.median(deltas)) if deltas else 0.0,
                    "galaxies_better_than_RAR": int(sum(delta > 0 for delta in deltas)),
                    "galaxies_worse_than_RAR": int(sum(delta < 0 for delta in deltas)),
                }
            )
    return out


def build_predictive_cases(
    cases: list[Case],
    map_rows: dict[tuple[str, str], dict[str, str]],
    args: argparse.Namespace,
    rng: np.random.Generator,
) -> list[dict[str, Any]]:
    by_key = {(case.split, case.model): case for case in cases}
    rows: list[dict[str, Any]] = []
    for split in sorted({case.split for case in cases}):
        if (split, "RAR") not in by_key:
            continue
        cfg = split_config(split)
        h0_values = [parse_float(map_rows[(split, case.model)]["H0_trial"]) for case in cases if case.split == split]
        prep = prepare(args.sample, args.points, int(cfg["quality_max"]), str(cfg["distance_method"]), args.err_floor_kms, sorted(set(h0_values)), args.sigma_ln_ml, args.distance_floor_frac, args.hubble_prior_mode)
        for model in MODEL_ORDER:
            key = (split, model)
            if key in by_key:
                rows.extend(predictive_rows_for_case(f"posterior_baseline_{split}", by_key[key], prep, map_rows, args, rng))

    if args.include_stress and ("all_Q2", "RAR") in by_key:
        split = "all_Q2"
        cfg = split_config(split)
        h0_values = [parse_float(map_rows[(split, case.model)]["H0_trial"]) for case in cases if case.split == split]
        prep = prepare(args.sample, args.points, int(cfg["quality_max"]), str(cfg["distance_method"]), args.err_floor_kms, sorted(set(h0_values)), args.sigma_ln_ml, args.distance_floor_frac, args.hubble_prior_mode)
        gas, _disk, bulge, total = baryon_components(prep.data)
        gas_frac = gas / total
        bulge_frac = bulge / total
        gas_by_gal = galaxy_metric(prep.data, prep.galaxies, gas_frac)
        bulge_by_gal = galaxy_metric(prep.data, prep.galaxies, bulge_frac)
        stress_preps = {
            "posterior_stress_gas_dominated_all_Q2": filter_prepared(prep, {g for g, v in gas_by_gal.items() if math.isfinite(v) and v >= args.gas_fraction_min}),
            "posterior_stress_bulge_removed_all_Q2": filter_prepared(prep, {g for g, v in bulge_by_gal.items() if not math.isfinite(v) or v <= args.bulge_fraction_max}),
        }
        gbar = gbar_fiducial(prep.data)
        g0 = acceleration_cH0_over_2pi(65.0)
        med_rad = galaxy_metric(prep.data, prep.galaxies, prep.data.rad_kpc)
        outer = np.asarray([prep.data.rad_kpc[i] >= med_rad[str(g)] for i, g in enumerate(prep.data.galaxy)], dtype=bool)
        stress_preps["posterior_stress_low_accel_outer_points_all_Q2"] = filter_prepared(prep, point_mask=outer & (gbar < g0))
        inc_by_gal = galaxy_metric(prep.data, prep.galaxies, prep.data.inc_deg)
        stress_preps["posterior_stress_high_inclination_all_Q2"] = filter_prepared(prep, {g for g, v in inc_by_gal.items() if math.isfinite(v) and v >= args.high_inc_min})
        stress_preps["posterior_stress_low_inclination_all_Q2"] = filter_prepared(prep, {g for g, v in inc_by_gal.items() if math.isfinite(v) and v < args.high_inc_min})
        prep_q1 = prepare(args.sample, args.points, 1, "all", args.err_floor_kms, sorted(set(h0_values)), args.sigma_ln_ml, args.distance_floor_frac, args.hubble_prior_mode)
        stress_preps["posterior_stress_high_quality_Q1_all"] = prep_q1
        for score_case, sprep in stress_preps.items():
            for model in MODEL_ORDER:
                key = (split, model)
                if key in by_key:
                    rows.extend(predictive_rows_for_case(score_case, by_key[key], sprep, map_rows, args, rng))
    add_predictive_deltas(rows)
    return rows


def make_predictive_plot(summary: list[dict[str, Any]], path: Path) -> None:
    if plt is None:
        return
    plamb_rows = [row for row in summary if row["model"] != "RAR"]
    if not plamb_rows:
        return
    cases = [row["score_case"] for row in plamb_rows]
    values = [float(row["delta_log_predictive_vs_RAR"]) for row in plamb_rows]
    fig_h = max(4.5, 0.32 * len(cases))
    fig, ax = plt.subplots(figsize=(10.5, fig_h))
    y = np.arange(len(cases))
    colors = ["#2e7d32" if v >= 0 else "#b23b3b" for v in values]
    ax.barh(y, values, color=colors)
    ax.axvline(0.0, color="#333333", linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(cases, fontsize=7)
    ax.set_xlabel("Delta posterior-predictive log score vs RAR")
    ax.set_title("SPARC Chain-Based Predictive Scores")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def write_report(
    outdir: Path,
    run_dir: Path,
    cases: list[Case],
    diag_rows: list[dict[str, Any]],
    check_rows: list[dict[str, Any]],
    nuisance_summary: list[dict[str, Any]],
    predictive_summary: list[dict[str, Any]],
    args: argparse.Namespace,
) -> None:
    lines = [
        "# SPARC Posterior Run Analysis",
        "",
        f"Run directory: `{run_dir}`",
        "",
        "## Chain Diagnostics",
        "",
        "| Split | Model | Chains | Saved / chain | Max Rhat | Global accept | Galaxy accept |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in diag_rows:
        lines.append(
            f"| {row['split']} | {row['model']} | {row['n_chains']} | {row['n_saved_per_chain']} | "
            f"{float(row['max_global_rhat']):.4g} | {float(row['mean_accept_global']):.3f} | {float(row['mean_accept_galaxy']):.3f} |"
        )
    lines.extend(["", "## PLAMB Checks", "", "| Split | Model | p median [16,84] | p=0.5 in 68%? | kappa median [16,84] | kappa=1 in 68%? |", "|---|---|---:|---:|---:|---:|"])
    for row in check_rows:
        if row["model"] == "RAR":
            continue
        p_text = ""
        if row["p_median"] != "":
            p_text = f"{float(row['p_median']):.5g} [{float(row['p_q16']):.5g}, {float(row['p_q84']):.5g}]"
        k_text = ""
        if row["kappa_median"] != "":
            k_text = f"{float(row['kappa_median']):.5g} [{float(row['kappa_q16']):.5g}, {float(row['kappa_q84']):.5g}]"
        lines.append(f"| {row['split']} | {row['model']} | {p_text} | {row['p_half_inside_68']} | {k_text} | {row['kappa_one_inside_68']} |")
    lines.extend(["", "## Nuisance Pulls", "", "| Split | Model | Distance pull >3 | M/L pull >3 | Max |d pull| | Max |M/L pull| |", "|---|---|---:|---:|---:|---:|"])
    for row in nuisance_summary:
        lines.append(
            f"| {row['split']} | {row['model']} | {row['distance_pull_gt3']} | {row['stellar_ml_pull_gt3']} | "
            f"{float(row['max_abs_distance_pull']):.3f} | {float(row['max_abs_stellar_ml_pull']):.3f} |"
        )
    if predictive_summary:
        lines.extend(["", "## Chain-Based Predictive Scores", "", "| Case | Model | Delta log score vs RAR | Better / worse galaxies |", "|---|---|---:|---:|"])
        for row in predictive_summary:
            if row["model"] == "RAR":
                continue
            lines.append(
                f"| {row['score_case']} | {row['model']} | {float(row['delta_log_predictive_vs_RAR']):.3f} | "
                f"{row['galaxies_better_than_RAR']} / {row['galaxies_worse_than_RAR']} |"
            )
    lines.extend(
        [
            "",
            "## Output Files",
            "",
            f"- case diagnostics: `{outdir / 'posterior_case_diagnostics.csv'}`",
            f"- global parameter summary: `{outdir / 'posterior_global_parameter_summary.csv'}`",
            f"- PLAMB checks: `{outdir / 'posterior_plamb_checks.csv'}`",
            f"- nuisance summary: `{outdir / 'posterior_nuisance_pull_summary.csv'}`",
            f"- predictive summary: `{outdir / 'posterior_predictive_summary.csv'}`",
        ]
    )
    (outdir / "posterior_run_analysis_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze a completed SPARC hierarchical posterior sampler run.")
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--map-dir", type=Path, default=DEFAULT_MAP_DIR)
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--points", type=Path, default=DEFAULT_POINTS)
    parser.add_argument("--output-label", default="")
    parser.add_argument("--max-global-draws", type=int, default=600)
    parser.add_argument("--nuisance-draws", type=int, default=8)
    parser.add_argument("--skip-predictive", action="store_true")
    parser.add_argument("--include-stress", action="store_true", default=True)
    parser.add_argument("--err-floor-kms", type=float, default=3.0)
    parser.add_argument("--sigma-ln-ml", type=float, default=0.25)
    parser.add_argument("--distance-floor-frac", type=float, default=0.03)
    parser.add_argument("--hubble-prior-mode", choices=["model_h0_rescale", "published"], default="model_h0_rescale")
    parser.add_argument("--gas-fraction-min", type=float, default=0.5)
    parser.add_argument("--bulge-fraction-max", type=float, default=0.2)
    parser.add_argument("--high-inc-min", type=float, default=60.0)
    parser.add_argument("--seed", type=int, default=260714)
    args = parser.parse_args()

    run_dir = args.run_dir
    cases = discover_cases(run_dir)
    if not cases:
        raise FileNotFoundError(f"No posterior_samples.npz files found under {run_dir}")
    label = args.output_label or f"{run_dir.name}_analysis"
    outdir = run_dir / label
    outdir.mkdir(parents=True, exist_ok=True)
    map_rows = load_map_rows(args.map_dir)
    rng = np.random.default_rng(args.seed)

    global_rows: list[dict[str, Any]] = []
    diag_rows: list[dict[str, Any]] = []
    check_rows: list[dict[str, Any]] = []
    nuisance_rows: list[dict[str, Any]] = []
    nuisance_summary: list[dict[str, Any]] = []
    for case in cases:
        prows = global_parameter_rows(case)
        global_rows.extend(prows)
        diag_rows.append(case_diagnostic_row(case, prows))
        check_rows.append(posterior_checks(case, prows))
        prep = prepare_for_case(case, map_rows, args.sample, args.points, args)
        nrows, nsummary = nuisance_pull_rows(case, prep, map_rows)
        nuisance_rows.extend(nrows)
        nuisance_summary.append(nsummary)

    predictive_rows: list[dict[str, Any]] = []
    predictive_summary: list[dict[str, Any]] = []
    if not args.skip_predictive:
        predictive_rows = build_predictive_cases(cases, map_rows, args, rng)
        predictive_summary = summarize_predictive(predictive_rows)

    write_csv(outdir / "posterior_case_diagnostics.csv", diag_rows, ["split", "model", "n_chains", "n_saved_per_chain", "n_global_params", "n_galaxies", "mean_accept_global", "mean_accept_galaxy", "max_global_rhat", "median_global_rhat", "complete"])
    write_csv(outdir / "posterior_global_parameter_summary.csv", global_rows, ["split", "model", "parameter", "mean", "sd", "q05", "q16", "median", "q84", "q95", "rhat"])
    write_csv(outdir / "posterior_plamb_checks.csv", check_rows, ["split", "model", "p_median", "p_q16", "p_q84", "p_q05", "p_q95", "p_half_inside_68", "p_half_inside_90", "kappa_median", "kappa_q16", "kappa_q84", "kappa_q05", "kappa_q95", "kappa_one_inside_68", "kappa_one_inside_90"])
    write_csv(outdir / "posterior_galaxy_nuisance_pulls.csv", nuisance_rows, ["split", "model", "galaxy", "distance_pull_median", "stellar_ml_pull_median", "distance_multiplier_median", "stellar_ml_multiplier_median"])
    write_csv(outdir / "posterior_nuisance_pull_summary.csv", nuisance_summary, ["split", "model", "n_galaxies", "max_abs_distance_pull", "median_abs_distance_pull", "distance_pull_gt2", "distance_pull_gt3", "max_abs_stellar_ml_pull", "median_abs_stellar_ml_pull", "stellar_ml_pull_gt2", "stellar_ml_pull_gt3"])
    if predictive_rows:
        write_csv(outdir / "posterior_predictive_galaxy_scores.csv", predictive_rows, ["score_case", "split", "model", "galaxy", "N_points", "log_predictive", "delta_log_predictive_vs_RAR"])
        write_csv(outdir / "posterior_predictive_summary.csv", predictive_summary, ["score_case", "model", "N_galaxies", "N_points", "total_log_predictive", "delta_log_predictive_vs_RAR", "median_galaxy_delta_vs_RAR", "galaxies_better_than_RAR", "galaxies_worse_than_RAR"])
        make_predictive_plot(predictive_summary, outdir / "posterior_predictive_delta_logscore.png")
    (outdir / "posterior_run_analysis_metadata.json").write_text(json.dumps({"args": vars(args), "n_cases": len(cases)}, indent=2, default=str), encoding="utf-8")
    write_report(outdir, run_dir, cases, diag_rows, check_rows, nuisance_summary, predictive_summary, args)
    print(f"Saved analysis report: {outdir / 'posterior_run_analysis_report.md'}")


if __name__ == "__main__":
    main()
