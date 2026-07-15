#!/usr/bin/env python3
r"""
Blocked posterior sampler for the hierarchical SPARC optical-depth models.

This is the overnight-grade follow-up to the MAP and fast posterior-predictive
diagnostics. It samples the full hierarchical parameter vector:

  - global disk and bulge stellar M/L scales,
  - optional RAR or PLAMB scale/exponent parameters,
  - one distance shift per galaxy,
  - one stellar M/L shift per galaxy.

It uses a blocked random-walk Metropolis sampler rather than emcee because the
full parameter vector has hundreds of galaxy-level nuisance dimensions. Global
parameters are proposed as one block; galaxy nuisance pairs are updated as local
blocks using only that galaxy's likelihood and priors.

Run the default overnight comparison with:
    run_sparc_posterior_sampler_overnight.cmd

Or directly:
    python -u sample_sparc_hierarchical_posterior.py --splits all_Q2,hubble_flow_Q2,non_hubble_Q2 --models RAR,PLAMB_OPTICAL_DEPTH_KAPPA_P
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


ROOT = Path(__file__).resolve().parents[3]
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "sparc_hierarchical_posterior"
DEFAULT_MAP_DIR = MAP_BASE_OUTDIR / "optical_depth_hierarchical_20260714"

SUPPORTED_MODELS = {
    "RAR",
    "PLAMB_OPTICAL_DEPTH",
    "PLAMB_OPTICAL_DEPTH_KAPPA",
    "PLAMB_OPTICAL_DEPTH_KAPPA_P",
}


@dataclass
class ChainResult:
    global_samples: np.ndarray
    nuisance_samples: np.ndarray
    logpost_samples: np.ndarray
    loglike_samples: np.ndarray
    logprior_samples: np.ndarray
    accept_global: float
    accept_galaxy: float
    final_state: np.ndarray
    final_logpost: float
    final_local_loglike: np.ndarray
    final_local_logprior: np.ndarray


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


def parse_csv_list(value: str) -> list[str]:
    return [token.strip() for token in value.split(",") if token.strip()]


def load_map_rows(map_dir: Path) -> dict[tuple[str, str], dict[str, str]]:
    rows = read_csv(map_dir / "sparc_hierarchical_map_summary.csv")
    return {(row["split"], row["model"]): row for row in rows}


def load_nuisance_rows(map_dir: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    path = map_dir / "sparc_hierarchical_map_galaxy_nuisance.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Re-run fit_sparc_hierarchical_map.py so the sampler has MAP nuisance starts."
        )
    rows = read_csv(path)
    return {(row["split"], row["model"], row["galaxy"]): row for row in rows}


def model_param_names(model: str) -> list[str]:
    if model == "RAR":
        return ["log_ydisk", "log_ybul", "log10_gdagger"]
    if model == "PLAMB_OPTICAL_DEPTH":
        return ["log_ydisk", "log_ybul"]
    if model == "PLAMB_OPTICAL_DEPTH_KAPPA":
        return ["log_ydisk", "log_ybul", "log10_kappa"]
    if model == "PLAMB_OPTICAL_DEPTH_KAPPA_P":
        return ["log_ydisk", "log_ybul", "log10_kappa", "bridge_exponent"]
    raise ValueError(model)


def global_bounds(model: str) -> list[tuple[float, float]]:
    bounds = [(math.log(0.05), math.log(1.50)), (math.log(0.05), math.log(1.80))]
    if model == "RAR":
        bounds.append((-11.6, -9.2))
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA":
        bounds.append((-0.50, 0.50))
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA_P":
        bounds.extend([(-0.50, 0.50), (0.15, 1.20)])
    elif model != "PLAMB_OPTICAL_DEPTH":
        raise ValueError(model)
    return bounds


def split_data_by_galaxy(data: Dataset, galaxies: list[str]) -> list[Dataset]:
    subsets: list[Dataset] = []
    for galaxy in galaxies:
        mask = data.galaxy == galaxy
        subsets.append(
            Dataset(
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
        )
    return subsets


def initial_state(
    split: str,
    model: str,
    prep: Any,
    map_row: dict[str, str],
    nuisance_by_key: dict[tuple[str, str, str], dict[str, str]],
    rng: np.random.Generator,
    jitter: float,
) -> np.ndarray:
    params: list[float] = [math.log(parse_float(map_row["ydisk"])), math.log(parse_float(map_row["ybul"]))]
    if model == "RAR":
        params.append(parse_float(map_row["log10_gdagger"]))
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA":
        params.append(parse_float(map_row["log10_kappa"]))
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA_P":
        params.extend([parse_float(map_row["log10_kappa"]), parse_float(map_row["bridge_exponent"])])

    logd: list[float] = []
    logeta: list[float] = []
    for galaxy in prep.galaxies:
        row = nuisance_by_key.get((split, model, galaxy))
        if row is None:
            raise KeyError(f"Missing nuisance MAP row for {split} {model} {galaxy}")
        logd.append(parse_float(row["log_distance_shift"]))
        logeta.append(parse_float(row["log_stellar_ml_shift"]))

    state = np.asarray(params + logd + logeta, dtype=float)
    if jitter > 0.0:
        bounds = state_bounds(model, prep)
        noise = rng.normal(0.0, jitter, size=len(state))
        # Keep physical scale/exponent parameters tighter than nuisance starts.
        noise[: len(params)] *= 0.35
        state = np.asarray(
            [min(max(value + eps, lo + 1e-8), hi - 1e-8) for value, eps, (lo, hi) in zip(state, noise, bounds)],
            dtype=float,
        )
    return state


def state_bounds(model: str, prep: Any) -> list[tuple[float, float]]:
    bounds = global_bounds(model)
    for mu, sigma in zip(prep.mu_logd_by_h0[next(iter(prep.mu_logd_by_h0.keys()))], prep.sigma_logd):
        half_width = max(4.0 * float(sigma), 0.08)
        bounds.append((float(mu - half_width), float(mu + half_width)))
    eta_width = max(4.0 * float(prep.sigma_ln_ml), 0.2)
    bounds.extend([(-eta_width, eta_width)] * len(prep.galaxies))
    return bounds


def unpack_global(state: np.ndarray, model: str) -> dict[str, float]:
    values = {name: float(state[i]) for i, name in enumerate(model_param_names(model))}
    values["ydisk"] = math.exp(values["log_ydisk"])
    values["ybul"] = math.exp(values["log_ybul"])
    return values


def predict_galaxy(data: Dataset, model: str, h0: float, glob: dict[str, float], logd: float, logeta: float) -> np.ndarray:
    scale = math.exp(logd)
    eta = math.exp(logeta)
    root_scale = math.sqrt(scale)
    root_eta = math.sqrt(eta)
    rad = data.rad_kpc * scale
    vgas = data.vgas * root_scale
    vdisk = data.vdisk * root_scale * root_eta
    vbul = data.vbul * root_scale * root_eta
    v2_bar = vgas * np.abs(vgas) + glob["ydisk"] * vdisk**2 + glob["ybul"] * vbul**2
    gbar = np.maximum(v2_bar, 1e-18) * 1.0e6 / (rad * KPC_M)
    if model == "RAR":
        tau = np.sqrt(np.maximum(gbar / (10.0 ** glob["log10_gdagger"]), 1e-30))
    elif model == "PLAMB_OPTICAL_DEPTH":
        tau = np.sqrt(np.maximum(gbar / acceleration_cH0_over_2pi(h0), 1e-30))
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA":
        g0 = (10.0 ** glob["log10_kappa"]) * acceleration_cH0_over_2pi(h0)
        tau = np.sqrt(np.maximum(gbar / g0, 1e-30))
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA_P":
        g0 = (10.0 ** glob["log10_kappa"]) * acceleration_cH0_over_2pi(h0)
        tau = np.maximum(gbar / g0, 1e-30) ** glob["bridge_exponent"]
    else:
        raise ValueError(model)
    gpred = gbar / np.maximum(1.0 - np.exp(-tau), 1e-12)
    return np.sqrt(np.maximum(gpred, 0.0) * rad * KPC_M) / 1000.0


def local_loglike(data: Dataset, model: str, h0: float, glob: dict[str, float], logd: float, logeta: float) -> float:
    pred = predict_galaxy(data, model, h0, glob, logd, logeta)
    residual = (data.vobs - pred) / data.sigma_v
    return float(-0.5 * np.sum(residual**2 + np.log(2.0 * math.pi * data.sigma_v**2)))


def local_logprior(logd: float, logeta: float, mu_logd: float, sigma_logd: float, sigma_ln_ml: float) -> float:
    return float(-0.5 * ((logd - mu_logd) / sigma_logd) ** 2 - math.log(max(sigma_logd, 1e-12)) - 0.5 * (logeta / sigma_ln_ml) ** 2)


def global_logprior(glob: dict[str, float], model: str, sigma_global_ml: float) -> float:
    lp = -0.5 * (glob["log_ydisk"] - math.log(0.5)) ** 2 / sigma_global_ml**2
    lp += -0.5 * (glob["log_ybul"] - math.log(0.7)) ** 2 / sigma_global_ml**2
    for value, (lo, hi) in zip([glob[name] for name in model_param_names(model)], global_bounds(model)):
        if not (lo <= value <= hi):
            return -np.inf
    return float(lp)


def evaluate_state(
    state: np.ndarray,
    model: str,
    h0: float,
    prep: Any,
    galaxy_data: list[Dataset],
    mu_logd: np.ndarray,
    sigma_global_ml: float,
) -> tuple[float, np.ndarray, np.ndarray]:
    glob = unpack_global(state, model)
    n_gal = len(prep.galaxies)
    offset = len(model_param_names(model))
    logd = state[offset: offset + n_gal]
    logeta = state[offset + n_gal: offset + 2 * n_gal]
    local_ll = np.empty(n_gal, dtype=float)
    local_lp = np.empty(n_gal, dtype=float)
    for i, data in enumerate(galaxy_data):
        local_ll[i] = local_loglike(data, model, h0, glob, float(logd[i]), float(logeta[i]))
        local_lp[i] = local_logprior(float(logd[i]), float(logeta[i]), float(mu_logd[i]), float(prep.sigma_logd[i]), float(prep.sigma_ln_ml))
    logpost = global_logprior(glob, model, sigma_global_ml) + float(np.sum(local_ll + local_lp))
    return logpost, local_ll, local_lp


def global_step_scales(model: str, args: argparse.Namespace) -> np.ndarray:
    scales: list[float] = [args.global_log_ml_step, args.global_log_ml_step]
    if model == "RAR":
        scales.append(args.log10_scale_step)
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA":
        scales.append(args.log10_scale_step)
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA_P":
        scales.extend([args.log10_scale_step, args.p_step])
    return np.asarray(scales, dtype=float)


def run_chain(
    split: str,
    model: str,
    chain_id: int,
    prep: Any,
    map_row: dict[str, str],
    nuisance_by_key: dict[tuple[str, str, str], dict[str, str]],
    args: argparse.Namespace,
) -> ChainResult:
    rng = np.random.default_rng(args.seed + 1009 * chain_id + 9176 * (abs(hash(split + model)) % 1000))
    h0 = parse_float(map_row["H0_trial"])
    galaxy_data = split_data_by_galaxy(prep.data, prep.galaxies)
    mu_logd = prep.mu_logd_by_h0[h0]
    state = initial_state(split, model, prep, map_row, nuisance_by_key, rng, args.start_jitter)
    bounds = state_bounds(model, prep)
    logpost, local_ll, local_lp = evaluate_state(state, model, h0, prep, galaxy_data, mu_logd, args.sigma_global_ml)
    n_gal = len(prep.galaxies)
    n_global = len(model_param_names(model))
    global_scales = global_step_scales(model, args)
    logd_step = args.logd_step
    logeta_step = args.logeta_step
    save_global: list[np.ndarray] = []
    save_nuisance: list[np.ndarray] = []
    save_logpost: list[float] = []
    save_loglike: list[float] = []
    save_logprior: list[float] = []
    accept_global = 0
    attempt_global = 0
    accept_galaxy = 0
    attempt_galaxy = 0

    for sweep in range(1, args.n_sweeps + 1):
        proposal = state.copy()
        proposal[:n_global] += rng.normal(0.0, global_scales)
        if all(bounds[i][0] <= proposal[i] <= bounds[i][1] for i in range(n_global)):
            p_logpost, p_ll, p_lp = evaluate_state(proposal, model, h0, prep, galaxy_data, mu_logd, args.sigma_global_ml)
            attempt_global += 1
            if math.log(rng.random()) < p_logpost - logpost:
                state = proposal
                logpost = p_logpost
                local_ll = p_ll
                local_lp = p_lp
                accept_global += 1

        n_update = max(1, int(round(args.galaxy_update_fraction * n_gal)))
        for gal_i in rng.choice(n_gal, size=n_update, replace=False):
            offset = n_global
            idx_d = offset + gal_i
            idx_eta = offset + n_gal + gal_i
            proposal_d = float(state[idx_d] + rng.normal(0.0, logd_step))
            proposal_eta = float(state[idx_eta] + rng.normal(0.0, logeta_step))
            if not (bounds[idx_d][0] <= proposal_d <= bounds[idx_d][1] and bounds[idx_eta][0] <= proposal_eta <= bounds[idx_eta][1]):
                attempt_galaxy += 1
                continue
            glob = unpack_global(state, model)
            p_ll_i = local_loglike(galaxy_data[gal_i], model, h0, glob, proposal_d, proposal_eta)
            p_lp_i = local_logprior(proposal_d, proposal_eta, float(mu_logd[gal_i]), float(prep.sigma_logd[gal_i]), float(prep.sigma_ln_ml))
            delta = (p_ll_i + p_lp_i) - (local_ll[gal_i] + local_lp[gal_i])
            attempt_galaxy += 1
            if math.log(rng.random()) < delta:
                state[idx_d] = proposal_d
                state[idx_eta] = proposal_eta
                logpost += delta
                local_ll[gal_i] = p_ll_i
                local_lp[gal_i] = p_lp_i
                accept_galaxy += 1

        if sweep > args.burn_in and (sweep - args.burn_in) % args.thin == 0:
            save_global.append(state[:n_global].copy())
            save_nuisance.append(state[n_global:].copy())
            save_logpost.append(float(logpost))
            save_loglike.append(float(np.sum(local_ll)))
            save_logprior.append(float(global_logprior(unpack_global(state, model), model, args.sigma_global_ml) + np.sum(local_lp)))

        if args.progress_every > 0 and sweep % args.progress_every == 0:
            print(
                f"{split} {model} chain={chain_id} sweep={sweep}/{args.n_sweeps} "
                f"logpost={logpost:.3g} acc_global={accept_global / max(attempt_global, 1):.3f} "
                f"acc_galaxy={accept_galaxy / max(attempt_galaxy, 1):.3f} saved={len(save_global)}",
                flush=True,
            )

    return ChainResult(
        global_samples=np.asarray(save_global, dtype=float),
        nuisance_samples=np.asarray(save_nuisance, dtype=float),
        logpost_samples=np.asarray(save_logpost, dtype=float),
        loglike_samples=np.asarray(save_loglike, dtype=float),
        logprior_samples=np.asarray(save_logprior, dtype=float),
        accept_global=accept_global / max(attempt_global, 1),
        accept_galaxy=accept_galaxy / max(attempt_galaxy, 1),
        final_state=state,
        final_logpost=float(logpost),
        final_local_loglike=local_ll.copy(),
        final_local_logprior=local_lp.copy(),
    )


def summarize_vector(values: np.ndarray) -> dict[str, float]:
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
    m, n, p = chains.shape
    chain_means = np.mean(chains, axis=1)
    chain_vars = np.var(chains, axis=1, ddof=1)
    w = np.mean(chain_vars, axis=0)
    b = n * np.var(chain_means, axis=0, ddof=1)
    var_hat = ((n - 1) / n) * w + b / n
    return np.sqrt(np.divide(var_hat, w, out=np.full(p, np.nan), where=w > 0))


def write_case_outputs(
    case_dir: Path,
    split: str,
    model: str,
    prep: Any,
    map_row: dict[str, str],
    results: list[ChainResult],
    args: argparse.Namespace,
) -> None:
    case_dir.mkdir(parents=True, exist_ok=True)
    param_names = model_param_names(model)
    sample_rows: list[dict[str, Any]] = []
    for chain_id, res in enumerate(results, start=1):
        for draw_id, values in enumerate(res.global_samples):
            row: dict[str, Any] = {
                "split": split,
                "model": model,
                "chain": chain_id,
                "draw": draw_id,
                "logpost": res.logpost_samples[draw_id],
                "loglike": res.loglike_samples[draw_id],
                "logprior": res.logprior_samples[draw_id],
            }
            for name, value in zip(param_names, values):
                row[name] = value
            row["ydisk"] = math.exp(row["log_ydisk"])
            row["ybul"] = math.exp(row["log_ybul"])
            sample_rows.append(row)
    sample_fields = ["split", "model", "chain", "draw", "logpost", "loglike", "logprior"] + param_names + ["ydisk", "ybul"]
    write_csv(case_dir / "global_chain_samples.csv", sample_rows, sample_fields)

    all_global = np.vstack([res.global_samples for res in results if len(res.global_samples)])
    chain_global = np.stack([res.global_samples for res in results if len(res.global_samples)], axis=0)
    rh = rhat(chain_global)
    summary_rows: list[dict[str, Any]] = []
    for i, name in enumerate(param_names):
        stats = summarize_vector(all_global[:, i])
        stats.update({"parameter": name, "rhat": float(rh[i]) if len(rh) else ""})
        summary_rows.append(stats)
    if "log_ydisk" in param_names:
        stats = summarize_vector(np.exp(all_global[:, param_names.index("log_ydisk")]))
        stats.update({"parameter": "ydisk", "rhat": ""})
        summary_rows.append(stats)
    if "log_ybul" in param_names:
        stats = summarize_vector(np.exp(all_global[:, param_names.index("log_ybul")]))
        stats.update({"parameter": "ybul", "rhat": ""})
        summary_rows.append(stats)
    write_csv(case_dir / "posterior_global_summary.csv", summary_rows, ["parameter", "mean", "sd", "q05", "q16", "median", "q84", "q95", "rhat"])

    nuisance = np.vstack([res.nuisance_samples for res in results if len(res.nuisance_samples)])
    n_gal = len(prep.galaxies)
    nuisance_rows: list[dict[str, Any]] = []
    for i, galaxy in enumerate(prep.galaxies):
        logd_stats = summarize_vector(nuisance[:, i])
        logeta_stats = summarize_vector(nuisance[:, n_gal + i])
        nuisance_rows.append(
            {
                "split": split,
                "model": model,
                "galaxy": galaxy,
                "logd_median": logd_stats["median"],
                "logd_q16": logd_stats["q16"],
                "logd_q84": logd_stats["q84"],
                "distance_multiplier_median": math.exp(logd_stats["median"]),
                "logeta_median": logeta_stats["median"],
                "logeta_q16": logeta_stats["q16"],
                "logeta_q84": logeta_stats["q84"],
                "stellar_ml_multiplier_median": math.exp(logeta_stats["median"]),
            }
        )
    write_csv(
        case_dir / "posterior_galaxy_nuisance_summary.csv",
        nuisance_rows,
        [
            "split",
            "model",
            "galaxy",
            "logd_median",
            "logd_q16",
            "logd_q84",
            "distance_multiplier_median",
            "logeta_median",
            "logeta_q16",
            "logeta_q84",
            "stellar_ml_multiplier_median",
        ],
    )

    diagnostics = [
        {
            "split": split,
            "model": model,
            "chain": i + 1,
            "accept_global": res.accept_global,
            "accept_galaxy": res.accept_galaxy,
            "n_saved": len(res.global_samples),
            "final_logpost": res.final_logpost,
        }
        for i, res in enumerate(results)
    ]
    write_csv(case_dir / "chain_diagnostics.csv", diagnostics, ["split", "model", "chain", "accept_global", "accept_galaxy", "n_saved", "final_logpost"])
    np.savez_compressed(
        case_dir / "posterior_samples.npz",
        global_samples=chain_global,
        nuisance_samples=np.stack([res.nuisance_samples for res in results if len(res.nuisance_samples)], axis=0),
        logpost_samples=np.stack([res.logpost_samples for res in results if len(res.logpost_samples)], axis=0),
        param_names=np.asarray(param_names, dtype=object),
        galaxies=np.asarray(prep.galaxies, dtype=object),
    )
    write_case_report(case_dir, split, model, map_row, summary_rows, diagnostics, args)


def write_case_report(
    case_dir: Path,
    split: str,
    model: str,
    map_row: dict[str, str],
    summary_rows: list[dict[str, Any]],
    diagnostics: list[dict[str, Any]],
    args: argparse.Namespace,
) -> None:
    lines = [
        f"# SPARC Hierarchical Posterior Sampler: {split} / {model}",
        "",
        "## Setup",
        "",
        f"- MAP source objective: `{map_row.get('objective', '')}`",
        f"- chains: `{args.chains}`",
        f"- sweeps per chain: `{args.n_sweeps}`",
        f"- burn-in sweeps: `{args.burn_in}`",
        f"- thinning: `{args.thin}`",
        f"- galaxy update fraction per sweep: `{args.galaxy_update_fraction}`",
        "",
        "## Global Posterior Summary",
        "",
        "| Parameter | median | q16 | q84 | Rhat |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        rhat_text = "" if row.get("rhat", "") == "" or not math.isfinite(float(row.get("rhat", float("nan")))) else f"{float(row['rhat']):.4f}"
        lines.append(
            f"| {row['parameter']} | {float(row['median']):.6g} | {float(row['q16']):.6g} | {float(row['q84']):.6g} | {rhat_text} |"
        )
    lines.extend(
        [
            "",
            "## Chain Diagnostics",
            "",
            "| Chain | Global accept | Galaxy accept | Saved draws | Final logpost |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in diagnostics:
        lines.append(
            f"| {row['chain']} | {float(row['accept_global']):.3f} | {float(row['accept_galaxy']):.3f} | "
            f"{row['n_saved']} | {float(row['final_logpost']):.6g} |"
        )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- global samples: `{case_dir / 'global_chain_samples.csv'}`",
            f"- posterior summary: `{case_dir / 'posterior_global_summary.csv'}`",
            f"- galaxy nuisance summary: `{case_dir / 'posterior_galaxy_nuisance_summary.csv'}`",
            f"- compressed samples: `{case_dir / 'posterior_samples.npz'}`",
        ]
    )
    (case_dir / "posterior_sampler_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_overall_report(outdir: Path, case_dirs: list[Path], args: argparse.Namespace) -> None:
    lines = [
        "# SPARC Hierarchical Posterior Overnight Run",
        "",
        "## Run Configuration",
        "",
        f"- splits: `{args.splits}`",
        f"- models: `{args.models}`",
        f"- chains: `{args.chains}`",
        f"- sweeps per chain: `{args.n_sweeps}`",
        f"- burn-in: `{args.burn_in}`",
        f"- thin: `{args.thin}`",
        "",
        "## Case Reports",
        "",
    ]
    for case_dir in case_dirs:
        lines.append(f"- `{case_dir / 'posterior_sampler_report.md'}`")
    (outdir / "posterior_sampler_index.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_case(
    split: str,
    model: str,
    map_rows: dict[tuple[str, str], dict[str, str]],
    nuisance_by_key: dict[tuple[str, str, str], dict[str, str]],
    args: argparse.Namespace,
) -> Path:
    map_row = map_rows[(split, model)]
    h0 = parse_float(map_row["H0_trial"])
    config = next(item for item in split_configs() if item["split"] == split)
    prep = prepare(
        args.sample,
        args.points,
        int(config["quality_max"]),
        str(config["distance_method"]),
        args.err_floor_kms,
        [h0],
        args.sigma_ln_ml,
        args.distance_floor_frac,
        args.hubble_prior_mode,
    )
    case_dir = OUTDIR / split / model
    print(f"=== Sampling {split} {model}: galaxies={len(prep.galaxies)} points={len(prep.data.vobs)} ===", flush=True)
    results = [
        run_chain(split, model, chain_id + 1, prep, map_row, nuisance_by_key, args)
        for chain_id in range(args.chains)
    ]
    write_case_outputs(case_dir, split, model, prep, map_row, results, args)
    print(f"Saved case report: {case_dir / 'posterior_sampler_report.md'}", flush=True)
    return case_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Blocked MCMC posterior sampler for hierarchical SPARC PLAMB/RAR models.")
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--points", type=Path, default=DEFAULT_POINTS)
    parser.add_argument("--map-dir", type=Path, default=DEFAULT_MAP_DIR)
    parser.add_argument("--output-label", default="posterior_overnight_20260714")
    parser.add_argument("--splits", default="all_Q2,hubble_flow_Q2,non_hubble_Q2")
    parser.add_argument("--models", default="RAR,PLAMB_OPTICAL_DEPTH_KAPPA_P")
    parser.add_argument("--chains", type=int, default=3)
    parser.add_argument("--n-sweeps", type=int, default=24000)
    parser.add_argument("--burn-in", type=int, default=6000)
    parser.add_argument("--thin", type=int, default=12)
    parser.add_argument("--galaxy-update-fraction", type=float, default=0.5)
    parser.add_argument("--sigma-ln-ml", type=float, default=0.25)
    parser.add_argument("--sigma-global-ml", type=float, default=0.55)
    parser.add_argument("--distance-floor-frac", type=float, default=0.03)
    parser.add_argument("--hubble-prior-mode", choices=["model_h0_rescale", "published"], default="model_h0_rescale")
    parser.add_argument("--err-floor-kms", type=float, default=3.0)
    parser.add_argument("--global-log-ml-step", type=float, default=0.008)
    parser.add_argument("--log10-scale-step", type=float, default=0.006)
    parser.add_argument("--p-step", type=float, default=0.004)
    parser.add_argument("--logd-step", type=float, default=0.018)
    parser.add_argument("--logeta-step", type=float, default=0.035)
    parser.add_argument("--start-jitter", type=float, default=0.01)
    parser.add_argument("--seed", type=int, default=260714)
    parser.add_argument("--progress-every", type=int, default=1000)
    args = parser.parse_args()

    if args.output_label:
        safe = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in args.output_label)
        globals()["OUTDIR"] = OUTDIR / safe
    OUTDIR.mkdir(parents=True, exist_ok=True)

    splits = parse_csv_list(args.splits)
    models = parse_csv_list(args.models)
    unknown = [model for model in models if model not in SUPPORTED_MODELS]
    if unknown:
        raise ValueError(f"Unsupported model(s): {unknown}")
    map_rows = load_map_rows(args.map_dir)
    nuisance_by_key = load_nuisance_rows(args.map_dir)
    case_dirs: list[Path] = []
    for split in splits:
        for model in models:
            if (split, model) not in map_rows:
                raise KeyError(f"No MAP row for split={split}, model={model}")
            case_dirs.append(run_case(split, model, map_rows, nuisance_by_key, args))
    make_overall_report(OUTDIR, case_dirs, args)
    (OUTDIR / "posterior_sampler_metadata.json").write_text(json.dumps({"args": vars(args)}, indent=2, default=str), encoding="utf-8")
    print(f"Saved index: {OUTDIR / 'posterior_sampler_index.md'}", flush=True)


if __name__ == "__main__":
    main()
