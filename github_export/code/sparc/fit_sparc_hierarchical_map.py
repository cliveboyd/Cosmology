#!/usr/bin/env python3
r"""
Hierarchical MAP fit for SPARC fixed PLAMB variants vs RAR.

This is the first "proper" hierarchical step after the point-estimate,
jackknife, and Monte Carlo diagnostics. It fits galaxy-level nuisance
parameters simultaneously with global model parameters:

  - one distance shift per galaxy with a SPARC distance-error prior,
  - one stellar M/L multiplier per galaxy with a lognormal prior,
  - global disk and bulge M/L scales,
  - RAR has one additional global g_dagger parameter,
  - PLAMB uses fixed exponent choices and g0 = c H0 / (2 pi),
  - PLAMB_OPTICAL_DEPTH uses tau=sqrt(g_bar/g0), with g0 = c H0 / (2 pi),
  - PLAMB_OPTICAL_DEPTH_KAPPA allows g0 = kappa c H0 / (2 pi),
  - PLAMB_OPTICAL_DEPTH_KAPPA_P also allows tau = (g_bar/g0)^p.

It is a MAP optimizer, not a full posterior sampler. Use it to determine
whether the PLAMB/RAR result survives a coherent nuisance fit before spending
overnight time on MCMC.

Run:
    python fit_sparc_hierarchical_map.py
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

from diagnose_plamb_sparc_h0_nuisance import acceleration_cH0_over_2pi
from diagnose_plamb_sparc_rotation import DEFAULT_POINTS, DEFAULT_SAMPLE, KPC_M, Dataset, load_dataset


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "sparc_hierarchical_map"
SPARC_H0_REF = 73.0


@dataclass
class Prepared:
    data: Dataset
    galaxies: list[str]
    gal_idx: np.ndarray
    sigma_logd: np.ndarray
    mu_logd_by_h0: dict[float, np.ndarray]
    sigma_ln_ml: float


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


def sample_distance_info(sample_path: Path) -> dict[str, dict[str, float]]:
    info: dict[str, dict[str, float]] = {}
    for row in read_csv(sample_path):
        galaxy = row["Galaxy"]
        info[galaxy] = {
            "D_Mpc": parse_float(row.get("D_Mpc")),
            "e_D_Mpc": parse_float(row.get("e_D_Mpc")),
            "f_D": parse_float(row.get("f_D")),
        }
    return info


def prepare(
    sample_path: Path,
    points_path: Path,
    quality_max: int,
    distance_method: str,
    err_floor_kms: float,
    h0_values: list[float],
    sigma_ln_ml: float,
    distance_floor_frac: float,
    hubble_prior_mode: str,
) -> Prepared:
    data = load_dataset(sample_path, points_path, quality_max, distance_method, err_floor_kms)
    info = sample_distance_info(sample_path)
    galaxies = sorted(set(str(galaxy) for galaxy in data.galaxy))
    index = {galaxy: i for i, galaxy in enumerate(galaxies)}
    gal_idx = np.asarray([index[str(galaxy)] for galaxy in data.galaxy], dtype=int)
    sigma_logd = np.zeros(len(galaxies), dtype=float)
    f_d = np.zeros(len(galaxies), dtype=float)
    for galaxy, i in index.items():
        item = info.get(galaxy, {})
        d = item.get("D_Mpc", float("nan"))
        e_d = item.get("e_D_Mpc", float("nan"))
        if not math.isfinite(d) or d <= 0.0:
            d = 1.0
        if not math.isfinite(e_d) or e_d <= 0.0:
            e_d = distance_floor_frac * d
        sigma_logd[i] = max(e_d / d, distance_floor_frac)
        f_d[i] = item.get("f_D", 0.0)

    mu_logd_by_h0: dict[float, np.ndarray] = {}
    for h0 in h0_values:
        mu = np.zeros(len(galaxies), dtype=float)
        if hubble_prior_mode == "model_h0_rescale":
            mu = np.where(f_d == 1, math.log(SPARC_H0_REF / h0), 0.0)
        elif hubble_prior_mode == "published":
            mu = np.zeros(len(galaxies), dtype=float)
        else:
            raise ValueError(f"Unknown hubble_prior_mode: {hubble_prior_mode}")
        mu_logd_by_h0[float(h0)] = mu

    return Prepared(
        data=data,
        galaxies=galaxies,
        gal_idx=gal_idx,
        sigma_logd=sigma_logd,
        mu_logd_by_h0=mu_logd_by_h0,
        sigma_ln_ml=sigma_ln_ml,
    )


def unpack(
    theta: np.ndarray,
    model: str,
    n_gal: int,
) -> tuple[float, float, float | None, float | None, float | None, np.ndarray, np.ndarray]:
    log_ydisk = theta[0]
    log_ybul = theta[1]
    offset = 2
    log10_gdagger = None
    log10_kappa = None
    bridge_exponent = None
    if model == "RAR":
        log10_gdagger = float(theta[2])
        offset = 3
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA":
        log10_kappa = float(theta[2])
        offset = 3
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA_P":
        log10_kappa = float(theta[2])
        bridge_exponent = float(theta[3])
        offset = 4
    logd = theta[offset: offset + n_gal]
    logeta = theta[offset + n_gal: offset + 2 * n_gal]
    return float(np.exp(log_ydisk)), float(np.exp(log_ybul)), log10_gdagger, log10_kappa, bridge_exponent, logd, logeta


def scaled_components(prep: Prepared, logd: np.ndarray, logeta: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    data = prep.data
    scale = np.exp(logd[prep.gal_idx])
    eta = np.exp(logeta[prep.gal_idx])
    root_scale = np.sqrt(scale)
    root_eta = np.sqrt(eta)
    rad = data.rad_kpc * scale
    vgas = data.vgas * root_scale
    vdisk = data.vdisk * root_scale * root_eta
    vbul = data.vbul * root_scale * root_eta
    return rad, vgas, vdisk, vbul


def baryon_v2_scaled(vgas: np.ndarray, vdisk: np.ndarray, vbul: np.ndarray, ydisk: float, ybul: float) -> np.ndarray:
    return vgas * np.abs(vgas) + ydisk * vdisk**2 + ybul * vbul**2


def acceleration_from_v2(v2_kms2: np.ndarray, radius_kpc: np.ndarray) -> np.ndarray:
    return np.maximum(v2_kms2, 1e-18) * 1.0e6 / (radius_kpc * KPC_M)


def velocity_from_acceleration(g_si: np.ndarray, radius_kpc: np.ndarray) -> np.ndarray:
    return np.sqrt(np.maximum(g_si, 0.0) * radius_kpc * KPC_M) / 1000.0


def predict(theta: np.ndarray, model: str, prep: Prepared, h0: float, inertia_index: float, screen_index: float) -> tuple[np.ndarray, dict[str, Any]]:
    ydisk, ybul, log10_gdagger, log10_kappa, bridge_exponent, logd, logeta = unpack(theta, model, len(prep.galaxies))
    rad, vgas, vdisk, vbul = scaled_components(prep, logd, logeta)
    v2_bar = baryon_v2_scaled(vgas, vdisk, vbul, ydisk, ybul)
    gbar = acceleration_from_v2(v2_bar, rad)
    if model == "PLAMB":
        g0 = acceleration_cH0_over_2pi(h0)
        ratio = np.clip(g0 / np.maximum(gbar, 1e-30), 1e-8, 1e8)
        boost = 1.0 + ratio**inertia_index / (1.0 + np.clip(1.0 / ratio, 1e-8, 1e8) ** screen_index)
        pred = velocity_from_acceleration(gbar * boost, rad)
    elif model == "PLAMB_OPTICAL_DEPTH":
        g0 = acceleration_cH0_over_2pi(h0)
        tau = np.sqrt(np.maximum(gbar / g0, 1e-30))
        pred = velocity_from_acceleration(gbar / np.maximum(1.0 - np.exp(-tau), 1e-12), rad)
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA":
        assert log10_kappa is not None
        g0 = (10.0**log10_kappa) * acceleration_cH0_over_2pi(h0)
        tau = np.sqrt(np.maximum(gbar / g0, 1e-30))
        pred = velocity_from_acceleration(gbar / np.maximum(1.0 - np.exp(-tau), 1e-12), rad)
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA_P":
        assert log10_kappa is not None and bridge_exponent is not None
        g0 = (10.0**log10_kappa) * acceleration_cH0_over_2pi(h0)
        tau = np.maximum(gbar / g0, 1e-30) ** bridge_exponent
        pred = velocity_from_acceleration(gbar / np.maximum(1.0 - np.exp(-tau), 1e-12), rad)
    elif model == "RAR":
        assert log10_gdagger is not None
        gdagger = 10.0**log10_gdagger
        x = np.sqrt(np.maximum(gbar / gdagger, 1e-30))
        pred = velocity_from_acceleration(gbar / np.maximum(1.0 - np.exp(-x), 1e-12), rad)
    else:
        raise ValueError(model)
    meta = {
        "ydisk": ydisk,
        "ybul": ybul,
        "log10_gdagger": log10_gdagger if log10_gdagger is not None else "",
        "log10_g0_cosmo": math.log10(acceleration_cH0_over_2pi(h0)),
        "log10_kappa": log10_kappa if log10_kappa is not None else "",
        "bridge_exponent": bridge_exponent if bridge_exponent is not None else "",
        "mean_abs_logd": float(np.mean(np.abs(logd))),
        "mean_abs_logeta": float(np.mean(np.abs(logeta))),
        "max_abs_logd": float(np.max(np.abs(logd))),
        "max_abs_logeta": float(np.max(np.abs(logeta))),
    }
    return pred, meta


def objective(theta: np.ndarray, model: str, prep: Prepared, h0: float, inertia_index: float, screen_index: float, sigma_global_ml: float) -> tuple[float, float, float, dict[str, Any]]:
    pred, meta = predict(theta, model, prep, h0, inertia_index, screen_index)
    data = prep.data
    residual = (data.vobs - pred) / data.sigma_v
    chi2_data = float(np.sum(residual**2))
    ydisk, ybul, _g, _kappa, _p, logd, logeta = unpack(theta, model, len(prep.galaxies))
    mu_logd = prep.mu_logd_by_h0[float(h0)]
    chi2_logd = float(np.sum(((logd - mu_logd) / prep.sigma_logd) ** 2))
    chi2_logeta = float(np.sum((logeta / prep.sigma_ln_ml) ** 2))
    chi2_global_ml = ((math.log(ydisk / 0.5) / sigma_global_ml) ** 2) + ((math.log(ybul / 0.7) / sigma_global_ml) ** 2)
    chi2_prior = chi2_logd + chi2_logeta + chi2_global_ml
    total = chi2_data + chi2_prior
    meta.update(
        {
            "chi2_data": chi2_data,
            "chi2_prior": chi2_prior,
            "chi2_logd_prior": chi2_logd,
            "chi2_logeta_prior": chi2_logeta,
            "chi2_global_ml_prior": chi2_global_ml,
            "max_abs_distance_prior_pull": float(np.max(np.abs((logd - mu_logd) / prep.sigma_logd))),
            "max_abs_stellar_ml_prior_pull": float(np.max(np.abs(logeta / prep.sigma_ln_ml))),
        }
    )
    return total, chi2_data, chi2_prior, meta


def initial_and_bounds(model: str, prep: Prepared, h0: float, sigma_global_ml: float) -> tuple[np.ndarray, list[tuple[float, float]]]:
    n_gal = len(prep.galaxies)
    mu_logd = prep.mu_logd_by_h0[float(h0)]
    start = [math.log(0.5), math.log(0.7)]
    bounds = [(math.log(0.05), math.log(1.50)), (math.log(0.05), math.log(1.80))]
    if model == "RAR":
        start.append(-10.0)
        bounds.append((-11.6, -9.2))
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA":
        start.append(0.0)
        bounds.append((-0.50, 0.50))
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA_P":
        start.extend([0.075, 0.464])
        bounds.extend([(-0.50, 0.50), (0.15, 1.20)])
    start.extend(mu_logd.tolist())
    start.extend([0.0] * n_gal)
    for mu, sigma in zip(mu_logd, prep.sigma_logd):
        half_width = max(4.0 * sigma, 0.08)
        bounds.append((float(mu - half_width), float(mu + half_width)))
    eta_width = max(4.0 * prep.sigma_ln_ml, 0.2)
    bounds.extend([(-eta_width, eta_width)] * n_gal)
    return np.asarray(start, dtype=float), bounds


def fit_model(
    model: str,
    prep: Prepared,
    h0: float,
    inertia_index: float,
    screen_index: float,
    sigma_global_ml: float,
    maxiter: int,
    maxfun: int,
) -> dict[str, Any]:
    if minimize is None:
        raise RuntimeError(f"SciPy optimizer import failed: {SCIPY_IMPORT_ERROR}")
    x0, bounds = initial_and_bounds(model, prep, h0, sigma_global_ml)

    def fun(theta: np.ndarray) -> float:
        value, _data, _prior, _meta = objective(theta, model, prep, h0, inertia_index, screen_index, sigma_global_ml)
        return value

    result = minimize(
        fun,
        x0,
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": maxiter, "maxfun": maxfun, "ftol": 1e-8, "maxls": 60},
    )
    value, chi2_data, chi2_prior, meta = objective(result.x, model, prep, h0, inertia_index, screen_index, sigma_global_ml)
    n = len(prep.data.vobs)
    k = len(result.x)
    row: dict[str, Any] = {
        "model": model,
        "H0_trial": h0,
        "inertia_index": inertia_index if model == "PLAMB" else "",
        "screen_index": screen_index if model == "PLAMB" else "",
        "success": bool(result.success),
        "message": str(result.message),
        "N_points": n,
        "N_galaxies": len(prep.galaxies),
        "k_total": k,
        "objective": value,
        "chi2_data": chi2_data,
        "chi2_prior": chi2_prior,
        "chi2_data_dof": chi2_data / max(n - k, 1),
        "AIC_data_naive": chi2_data + 2.0 * k,
        "BIC_data_naive": chi2_data + k * math.log(max(n, 1)),
        "AIC_objective_naive": value + 2.0 * k,
        "BIC_objective_naive": value + k * math.log(max(n, 1)),
    }
    row.update(meta)
    row["_theta"] = result.x
    return row


def galaxy_nuisance_rows(split: str, row: dict[str, Any], prep: Prepared) -> list[dict[str, Any]]:
    theta = row.get("_theta")
    if theta is None:
        return []
    model = str(row["model"])
    h0 = float(row["H0_trial"])
    _ydisk, _ybul, _g, _kappa, _p, logd, logeta = unpack(np.asarray(theta, dtype=float), model, len(prep.galaxies))
    mu_logd = prep.mu_logd_by_h0[h0]
    rows: list[dict[str, Any]] = []
    for i, galaxy in enumerate(prep.galaxies):
        distance_pull = (float(logd[i]) - float(mu_logd[i])) / float(prep.sigma_logd[i])
        ml_pull = float(logeta[i]) / float(prep.sigma_ln_ml)
        rows.append(
            {
                "split": split,
                "model": model,
                "H0_trial": h0,
                "galaxy": galaxy,
                "log_distance_shift": float(logd[i]),
                "distance_multiplier": float(math.exp(logd[i])),
                "distance_prior_center_log": float(mu_logd[i]),
                "sigma_log_distance": float(prep.sigma_logd[i]),
                "distance_prior_pull": distance_pull,
                "log_stellar_ml_shift": float(logeta[i]),
                "stellar_ml_multiplier": float(math.exp(logeta[i])),
                "sigma_ln_stellar_ml": float(prep.sigma_ln_ml),
                "stellar_ml_prior_pull": ml_pull,
            }
        )
    return rows


def split_configs() -> list[dict[str, Any]]:
    return [
        {"split": "all_Q2", "quality_max": 2, "distance_method": "all", "plamb_h0": 65.0, "plamb_n": 0.425, "plamb_s": 1.25, "rar_h0": 65.0},
        {"split": "hubble_flow_Q2", "quality_max": 2, "distance_method": "hubble_flow", "plamb_h0": 70.0, "plamb_n": 0.4, "plamb_s": 1.25, "rar_h0": 62.5},
        {"split": "non_hubble_Q2", "quality_max": 2, "distance_method": "non_hubble", "plamb_h0": 70.0, "plamb_n": 0.5, "plamb_s": 1.25, "rar_h0": 62.5},
    ]


def make_plot(rows: list[dict[str, Any]], path: Path) -> None:
    if plt is None:
        return
    splits = sorted({row["split"] for row in rows})
    model_names = [name for name in ["PLAMB", "PLAMB_OPTICAL_DEPTH", "PLAMB_OPTICAL_DEPTH_KAPPA", "PLAMB_OPTICAL_DEPTH_KAPPA_P"] if any(row["model"] == name for row in rows)]
    offsets = np.linspace(-0.36, 0.36, len(model_names))
    x = np.arange(len(splits), dtype=float)
    fig, ax = plt.subplots(figsize=(10.4, 5.8))
    width = 0.72 / max(len(model_names), 1)
    colors = {
        "PLAMB": "#8e44ad",
        "PLAMB_OPTICAL_DEPTH": "#2471a3",
        "PLAMB_OPTICAL_DEPTH_KAPPA": "#229954",
        "PLAMB_OPTICAL_DEPTH_KAPPA_P": "#d68910",
    }
    for offset, model_name in zip(offsets, model_names):
        values = []
        for split in splits:
            r = next(row for row in rows if row["split"] == split and row["model"] == "RAR")
            m = next((row for row in rows if row["split"] == split and row["model"] == model_name), None)
            values.append(float(m["objective"]) - float(r["objective"]) if m is not None else float("nan"))
        ax.bar(x + offset, values, width=width, color=colors[model_name], label=model_name)
    ax.axhline(0.0, color="#333333", linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(splits, rotation=20, ha="right")
    ax.set_ylabel("Delta MAP objective vs RAR")
    ax.set_title("Hierarchical MAP SPARC Comparison")
    ax.legend()
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def delta_line(rows: list[dict[str, Any]], split: str, model_name: str) -> str:
    m = next(row for row in rows if row["split"] == split and row["model"] == model_name)
    r = next(row for row in rows if row["split"] == split and row["model"] == "RAR")
    return (
        f"| {split} | {model_name} | {float(m['objective']) - float(r['objective']):.3f} | "
        f"{float(m['chi2_data']) - float(r['chi2_data']):.3f} | "
        f"{float(m['chi2_data_dof']):.3f} | {float(r['chi2_data_dof']):.3f} |"
    )


def write_report(rows: list[dict[str, Any]], args: argparse.Namespace) -> None:
    branch_names = [
        name
        for name in ["PLAMB", "PLAMB_OPTICAL_DEPTH", "PLAMB_OPTICAL_DEPTH_KAPPA", "PLAMB_OPTICAL_DEPTH_KAPPA_P"]
        if any(row["model"] == name for row in rows)
    ]
    lines = [
        "# SPARC Hierarchical MAP Fit",
        "",
        "## Executive Summary",
        "",
        "Negative delta means the named PLAMB branch is preferred over RAR for the metric shown.",
        "",
        "| Split | Branch | Delta objective | Delta data chi2 | Branch data chi2/dof | RAR data chi2/dof |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for split in sorted({row["split"] for row in rows}):
        for branch_name in branch_names:
            if any(row["split"] == split and row["model"] == branch_name for row in rows):
                lines.append(delta_line(rows, split, branch_name))
    lines.extend(
        [
            "",
            "## Priors and Nuisance Structure",
            "",
            f"- sigma_ln_ml per galaxy: `{args.sigma_ln_ml}`",
            f"- sigma_global_ml: `{args.sigma_global_ml}`",
            f"- distance floor fraction: `{args.distance_floor_frac}`",
            f"- Hubble-flow distance prior mode: `{args.hubble_prior_mode}`",
            "- Each galaxy has one distance shift and one stellar M/L shift.",
            "- This is MAP optimization, not a sampled posterior.",
            "",
            "## Outputs",
            "",
            f"- Fit CSV: `{OUTDIR / 'sparc_hierarchical_map_summary.csv'}`",
            f"- Plot: `{OUTDIR / 'sparc_hierarchical_map_delta_objective.png'}`",
            f"- Galaxy nuisance MAP CSV: `{OUTDIR / 'sparc_hierarchical_map_galaxy_nuisance.csv'}`",
        ]
    )
    (OUTDIR / "sparc_hierarchical_map_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Hierarchical MAP SPARC PLAMB/RAR fit.")
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--points", type=Path, default=DEFAULT_POINTS)
    parser.add_argument("--err-floor-kms", type=float, default=3.0)
    parser.add_argument("--sigma-ln-ml", type=float, default=0.25)
    parser.add_argument("--sigma-global-ml", type=float, default=0.55)
    parser.add_argument("--distance-floor-frac", type=float, default=0.03)
    parser.add_argument("--hubble-prior-mode", choices=["model_h0_rescale", "published"], default="model_h0_rescale")
    parser.add_argument("--maxiter", type=int, default=900)
    parser.add_argument("--maxfun", type=int, default=180000)
    parser.add_argument("--only", default="", help="Optional split name to run only one split.")
    parser.add_argument("--output-label", default="", help="Optional subdirectory label for outputs.")
    args = parser.parse_args()

    if args.output_label:
        safe = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in args.output_label)
        globals()["OUTDIR"] = OUTDIR / safe

    rows: list[dict[str, Any]] = []
    nuisance_rows: list[dict[str, Any]] = []
    for config in split_configs():
        if args.only and config["split"] != args.only:
            continue
        h0_values = [float(config["plamb_h0"]), float(config["rar_h0"])]
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
        print(f"=== {config['split']} N={len(prep.data.vobs)} galaxies={len(prep.galaxies)} ===")
        plamb = fit_model(
            "PLAMB",
            prep,
            float(config["plamb_h0"]),
            float(config["plamb_n"]),
            float(config["plamb_s"]),
            args.sigma_global_ml,
            args.maxiter,
            args.maxfun,
        )
        plamb["split"] = config["split"]
        rows.append(plamb)
        nuisance_rows.extend(galaxy_nuisance_rows(str(config["split"]), plamb, prep))
        print(f"PLAMB objective={plamb['objective']:.4g} data_chi2={plamb['chi2_data']:.4g} success={plamb['success']}")
        optical_rows: list[dict[str, Any]] = []
        for model_name in ["PLAMB_OPTICAL_DEPTH", "PLAMB_OPTICAL_DEPTH_KAPPA", "PLAMB_OPTICAL_DEPTH_KAPPA_P"]:
            optical = fit_model(
                model_name,
                prep,
                float(config["plamb_h0"]),
                0.0,
                0.0,
                args.sigma_global_ml,
                args.maxiter,
                args.maxfun,
            )
            optical["split"] = config["split"]
            rows.append(optical)
            nuisance_rows.extend(galaxy_nuisance_rows(str(config["split"]), optical, prep))
            optical_rows.append(optical)
            print(
                f"{model_name} "
                f"objective={optical['objective']:.4g} data_chi2={optical['chi2_data']:.4g} "
                f"success={optical['success']}"
            )
        rar = fit_model(
            "RAR",
            prep,
            float(config["rar_h0"]),
            0.0,
            0.0,
            args.sigma_global_ml,
            args.maxiter,
            args.maxfun,
        )
        rar["split"] = config["split"]
        rows.append(rar)
        nuisance_rows.extend(galaxy_nuisance_rows(str(config["split"]), rar, prep))
        print(f"RAR    objective={rar['objective']:.4g} data_chi2={rar['chi2_data']:.4g} success={rar['success']}")
        print(f"delta objective PLAMB-RAR={plamb['objective'] - rar['objective']:.4g}")
        for optical in optical_rows:
            print(f"delta objective {optical['model']}-RAR={optical['objective'] - rar['objective']:.4g}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fields = [
        "split",
        "model",
        "H0_trial",
        "inertia_index",
        "screen_index",
        "success",
        "message",
        "N_points",
        "N_galaxies",
        "k_total",
        "objective",
        "chi2_data",
        "chi2_prior",
        "chi2_logd_prior",
        "chi2_logeta_prior",
        "chi2_global_ml_prior",
        "max_abs_distance_prior_pull",
        "max_abs_stellar_ml_prior_pull",
        "chi2_data_dof",
        "AIC_data_naive",
        "BIC_data_naive",
        "AIC_objective_naive",
        "BIC_objective_naive",
        "ydisk",
        "ybul",
        "log10_gdagger",
        "log10_g0_cosmo",
        "log10_kappa",
        "bridge_exponent",
        "mean_abs_logd",
        "mean_abs_logeta",
        "max_abs_logd",
        "max_abs_logeta",
    ]
    write_csv(OUTDIR / "sparc_hierarchical_map_summary.csv", rows, fields)
    write_csv(
        OUTDIR / "sparc_hierarchical_map_galaxy_nuisance.csv",
        nuisance_rows,
        [
            "split",
            "model",
            "H0_trial",
            "galaxy",
            "log_distance_shift",
            "distance_multiplier",
            "distance_prior_center_log",
            "sigma_log_distance",
            "distance_prior_pull",
            "log_stellar_ml_shift",
            "stellar_ml_multiplier",
            "sigma_ln_stellar_ml",
            "stellar_ml_prior_pull",
        ],
    )
    make_plot(rows, OUTDIR / "sparc_hierarchical_map_delta_objective.png")
    write_report(rows, args)
    (OUTDIR / "sparc_hierarchical_map_metadata.json").write_text(json.dumps({"args": vars(args)}, indent=2, default=str), encoding="utf-8")
    print(f"Saved report: {OUTDIR / 'sparc_hierarchical_map_report.md'}")


if __name__ == "__main__":
    main()
