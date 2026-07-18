#!/usr/bin/env python3
"""Train, sample and exactly validate the hierarchical LINX BBN surrogate."""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import UTC, datetime
from pathlib import Path

import emcee
import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import RidgeCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
OUTDIR = REPO_ROOT / "plamb_runs" / "diagnostics" / "bbn_lithium_hierarchical_surrogate_20260718"
ARCHIVE = REPO_ROOT / "github_export" / "results" / "2026-07-18" / "bbn_lithium"
REACTIONS = [
    "npdg", "dpHe3g", "ddHe3n", "ddtp", "tpag", "tdan",
    "taLi7g", "He3ntp", "He3dap", "He3aBe7g", "Be7nLi7p", "Li7paa",
]
TARGETS = ["D_H_1e5", "Yp_mass", "Li7_H_1e10"]
OBS = {
    "D_H_1e5": (2.533, 0.024),
    "Yp_mass": (0.2458, 0.0013),
    "Li7_H_1e10": (1.45, 0.25),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["train-sample", "finalise"], required=True)
    parser.add_argument("--outdir", type=Path, default=OUTDIR)
    parser.add_argument("--archive-dir", type=Path, default=ARCHIVE)
    parser.add_argument("--walkers", type=int, default=64)
    parser.add_argument("--burn", type=int, default=3000)
    parser.add_argument("--steps", type=int, default=12000)
    parser.add_argument("--thin-save", type=int, default=5)
    parser.add_argument("--validation-per-arm", type=int, default=160)
    parser.add_argument("--seed", type=int, default=20260718)
    return parser.parse_args()


def now() -> str:
    return datetime.now(UTC).astimezone().isoformat(timespec="seconds")


def feature_frame(frame: pd.DataFrame) -> np.ndarray:
    columns = [
        (frame["eta_fac"].to_numpy(float) - 1.0) / (0.04 / 6.12),
        (frame["tau_seconds"].to_numpy(float) - 878.4) / 0.5,
    ]
    columns.extend(frame[f"q_{name}"].to_numpy(float) for name in REACTIONS)
    return np.column_stack(columns)


def train_models(args: argparse.Namespace) -> tuple[dict[str, object], dict[str, object]]:
    path = args.outdir / "exact_linx_design_results.csv"
    frame = pd.read_csv(path)
    frame = frame.loc[frame["status"].eq("ok")].copy()
    if len(frame) < 4700:
        raise RuntimeError(f"Exact design incomplete: only {len(frame)} successful rows")
    train = frame.loc[frame["split"].eq("train")]
    holdout = frame.loc[frame["split"].eq("holdout")]
    x_train, x_holdout = feature_frame(train), feature_frame(holdout)
    models: dict[str, object] = {}
    metrics: dict[str, object] = {"generated": now(), "n_train": len(train), "n_holdout": len(holdout), "targets": {}}
    gates = {target: 0.25 * OBS[target][1] for target in TARGETS}

    for target in TARGETS:
        y_train = train[target].to_numpy(float)
        y_test = holdout[target].to_numpy(float)
        ridge = make_pipeline(
            PolynomialFeatures(degree=2, include_bias=False),
            StandardScaler(),
            RidgeCV(alphas=np.logspace(-8, 1, 16)),
        )
        forest = ExtraTreesRegressor(
            n_estimators=400,
            min_samples_leaf=2,
            max_features=1.0,
            random_state=args.seed,
            n_jobs=-1,
        )
        candidates = {"quadratic_ridge": ridge, "extra_trees": forest}
        scores: dict[str, dict[str, float]] = {}
        for name, model in candidates.items():
            model.fit(x_train, y_train)
            pred = model.predict(x_holdout)
            scores[name] = {
                "rmse": float(np.sqrt(np.mean((pred - y_test) ** 2))),
                "mae": float(np.mean(np.abs(pred - y_test))),
                "max_abs": float(np.max(np.abs(pred - y_test))),
            }
        ridge_passes = scores["quadratic_ridge"]["rmse"] <= gates[target]
        selected = "quadratic_ridge" if ridge_passes else min(scores, key=lambda name: scores[name]["rmse"])
        model = candidates[selected]
        if hasattr(model, "n_jobs"):
            model.n_jobs = 1
        models[target] = model
        metrics["targets"][target] = {
            "selected": selected,
            "observational_sigma": OBS[target][1],
            "gate_rmse": gates[target],
            "passed": scores[selected]["rmse"] <= gates[target],
            "candidates": scores,
        }
        print(f"{target}: selected={selected} rmse={scores[selected]['rmse']:.6g} gate={gates[target]:.6g}", flush=True)

    metrics["passed"] = all(value["passed"] for value in metrics["targets"].values())
    args.outdir.mkdir(parents=True, exist_ok=True)
    joblib.dump(models, args.outdir / "surrogate_models.joblib", compress=3)
    (args.outdir / "surrogate_holdout_metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return models, metrics


def predict(models: dict[str, object], theta14: np.ndarray) -> dict[str, np.ndarray]:
    return {target: np.asarray(models[target].predict(theta14), dtype=float) for target in TARGETS}


def log_probability(models: dict[str, object], theta: np.ndarray, arm: str) -> np.ndarray:
    values = np.atleast_2d(theta)
    base = values[:, :14]
    valid = np.all(np.isfinite(values), axis=1)
    valid &= np.abs(base[:, 0]) <= 5.0
    valid &= np.abs(base[:, 1]) <= 5.0
    valid &= np.all(np.abs(base[:, 2:]) <= 4.0, axis=1)
    if arm == "depletion":
        valid &= (values[:, 14] >= 0.0) & (values[:, 14] <= 1.5)
    result = np.full(len(values), -np.inf)
    if not np.any(valid):
        return result
    selected = values[valid]
    pred = predict(models, selected[:, :14])
    lp = -0.5 * np.sum(selected[:, :14] ** 2, axis=1)
    lp += -0.5 * ((pred["D_H_1e5"] - OBS["D_H_1e5"][0]) / OBS["D_H_1e5"][1]) ** 2
    lp += -0.5 * ((pred["Yp_mass"] - OBS["Yp_mass"][0]) / OBS["Yp_mass"][1]) ** 2
    if arm == "depletion":
        delta = selected[:, 14]
        lp += -0.5 * (delta / 0.30) ** 2
        surface = pred["Li7_H_1e10"] * 10.0 ** (-delta)
        lp += -0.5 * ((surface - OBS["Li7_H_1e10"][0]) / OBS["Li7_H_1e10"][1]) ** 2
    result[valid] = lp
    return result


def split_rhat(chain: np.ndarray) -> float:
    # emcee chain is steps x walkers x dimensions; split each walker in half.
    n_steps, n_walkers, n_dim = chain.shape
    half = n_steps // 2
    pieces = np.concatenate([chain[:half], chain[-half:]], axis=1).transpose(1, 0, 2)
    n = pieces.shape[1]
    means = pieces.mean(axis=1)
    variances = pieces.var(axis=1, ddof=1)
    within = variances.mean(axis=0)
    between = n * means.var(axis=0, ddof=1)
    var_hat = (n - 1.0) / n * within + between / n
    return float(np.nanmax(np.sqrt(var_hat / within)))


def initialise_walkers(models: dict[str, object], arm: str, n_walkers: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n_dim = 15 if arm == "depletion" else 14
    candidates = rng.normal(size=(120000, n_dim))
    candidates[:, 0:2] = np.clip(candidates[:, 0:2], -4.8, 4.8)
    candidates[:, 2:14] = np.clip(candidates[:, 2:14], -3.8, 3.8)
    if arm == "depletion":
        candidates[:, 14] = np.abs(candidates[:, 14]) * 0.30
    score = log_probability(models, candidates, arm)
    pool = np.argpartition(score, -max(1000, n_walkers))[-max(1000, n_walkers):]
    chosen = rng.choice(pool, size=n_walkers, replace=False)
    start = candidates[chosen].copy()
    start += rng.normal(scale=1e-3, size=start.shape)
    if arm == "depletion":
        start[:, 14] = np.maximum(start[:, 14], 1e-4)
    return start


def run_arm(args: argparse.Namespace, models: dict[str, object], arm: str) -> tuple[np.ndarray, dict[str, object]]:
    n_dim = 15 if arm == "depletion" else 14
    if args.walkers < 2 * n_dim:
        raise ValueError(f"Need at least {2 * n_dim} walkers for {arm}")
    start = initialise_walkers(models, arm, args.walkers, args.seed + (1000 if arm == "depletion" else 0))
    sampler = emcee.EnsembleSampler(
        args.walkers,
        n_dim,
        lambda theta: log_probability(models, theta, arm),
        vectorize=True,
    )
    print(f"[{now()}] {arm}: burn-in {args.burn} x {args.walkers}", flush=True)
    state = sampler.run_mcmc(start, args.burn, progress=True)
    sampler.reset()
    print(f"[{now()}] {arm}: production {args.steps} x {args.walkers}", flush=True)
    sampler.run_mcmc(state, args.steps, progress=True)
    chain = sampler.get_chain()
    try:
        tau = sampler.get_autocorr_time(tol=0)
    except Exception:
        tau = np.full(n_dim, np.nan)
    diagnostics = {
        "arm": arm,
        "walkers": args.walkers,
        "steps": args.steps,
        "mean_acceptance": float(np.mean(sampler.acceptance_fraction)),
        "max_split_rhat": split_rhat(chain),
        "autocorr_time": np.asarray(tau).tolist(),
        "min_ess_from_tau": float(np.nanmin(args.walkers * args.steps / tau)) if np.any(np.isfinite(tau)) else math.nan,
    }
    np.savez_compressed(args.outdir / f"{arm}_surrogate_chain.npz", chain=chain[:: args.thin_save])
    (args.outdir / f"{arm}_surrogate_diagnostics.json").write_text(json.dumps(diagnostics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[{now()}] {arm}: {diagnostics}", flush=True)
    return chain, diagnostics


def validation_rows(args: argparse.Namespace, models: dict[str, object], chains: dict[str, np.ndarray]) -> Path:
    rng = np.random.default_rng(args.seed + 9000)
    rows: list[dict[str, object]] = []
    for arm, chain in chains.items():
        flat = chain.reshape(-1, chain.shape[-1])
        indices = rng.choice(len(flat), size=args.validation_per_arm, replace=False)
        values = flat[indices]
        pred = predict(models, values[:, :14])
        for j, theta in enumerate(values):
            row: dict[str, object] = {
                "point_id": f"{arm}_{j:04d}",
                "arm": arm,
                "eta_fac": 1.0 + theta[0] * (0.04 / 6.12),
                "tau_seconds": 878.4 + theta[1] * 0.5,
                "delta_li_dex": theta[14] if arm == "depletion" else 0.0,
                "D_H_1e5_surrogate": pred["D_H_1e5"][j],
                "Yp_mass_surrogate": pred["Yp_mass"][j],
                "Li7_H_1e10_surrogate": pred["Li7_H_1e10"][j],
            }
            row.update({f"q_{name}": theta[k + 2] for k, name in enumerate(REACTIONS)})
            rows.append(row)
    path = args.outdir / "surrogate_posterior_exact_validation_points.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved exact-validation points: {path}", flush=True)
    return path


def train_sample(args: argparse.Namespace) -> None:
    models, metrics = train_models(args)
    if not metrics["passed"]:
        raise RuntimeError("Surrogate holdout gate failed; exact posterior validation was not started")
    chains: dict[str, np.ndarray] = {}
    diagnostics = []
    for arm in ["predictive", "depletion"]:
        chains[arm], diag = run_arm(args, models, arm)
        diagnostics.append(diag)
    validation_rows(args, models, chains)
    gate = {
        "generated": now(),
        "holdout_passed": metrics["passed"],
        "mcmc_diagnostics": diagnostics,
        "mcmc_gate_passed": all(
            d["max_split_rhat"] <= 1.05
            and d["mean_acceptance"] >= 0.15
            and d["mean_acceptance"] <= 0.70
            and d["min_ess_from_tau"] >= 400
            for d in diagnostics
        ),
    }
    (args.outdir / "surrogate_sampling_gate.json").write_text(json.dumps(gate, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def weighted_quantile(values: np.ndarray, weights: np.ndarray, probability: float) -> float:
    order = np.argsort(values)
    cdf = np.cumsum(weights[order])
    return float(np.interp(probability, cdf, values[order]))


def finalise(args: argparse.Namespace) -> None:
    metrics = json.loads((args.outdir / "surrogate_holdout_metrics.json").read_text(encoding="utf-8"))
    sampling = json.loads((args.outdir / "surrogate_sampling_gate.json").read_text(encoding="utf-8"))
    frame = pd.read_csv(args.outdir / "exact_linx_posterior_validation_results.csv")
    frame = frame.loc[frame["status"].eq("ok")].copy()
    rows: list[dict[str, object]] = []
    exact_gate = True
    for arm in ["predictive", "depletion"]:
        subset = frame.loc[frame["arm"].eq(arm)].copy()
        if len(subset) < int(0.95 * args.validation_per_arm):
            exact_gate = False
        logw = np.zeros(len(subset))
        arm_metrics: dict[str, object] = {"arm": arm, "n_exact": len(subset)}
        for target in TARGETS:
            error = subset[target].to_numpy(float) - subset[f"{target}_surrogate"].to_numpy(float)
            q95 = float(np.quantile(np.abs(error), 0.95))
            scaled = q95 / OBS[target][1]
            arm_metrics[f"{target}_q95_abs_error"] = q95
            arm_metrics[f"{target}_q95_abs_error_sigma"] = scaled
            exact_gate &= scaled <= 0.50
        for target in ["D_H_1e5", "Yp_mass"]:
            mean, sigma = OBS[target]
            surrogate = subset[f"{target}_surrogate"].to_numpy(float)
            exact = subset[target].to_numpy(float)
            logw += -0.5 * ((exact - mean) / sigma) ** 2 + 0.5 * ((surrogate - mean) / sigma) ** 2
        if arm == "depletion":
            mean, sigma = OBS["Li7_H_1e10"]
            delta = subset["delta_li_dex"].to_numpy(float)
            surrogate = subset["Li7_H_1e10_surrogate"].to_numpy(float) * 10.0 ** (-delta)
            exact = subset["Li7_H_1e10"].to_numpy(float) * 10.0 ** (-delta)
            logw += -0.5 * ((exact - mean) / sigma) ** 2 + 0.5 * ((surrogate - mean) / sigma) ** 2
        logw -= np.max(logw)
        weights = np.exp(logw)
        weights /= np.sum(weights)
        ess = float(1.0 / np.sum(weights**2))
        arm_metrics["importance_ess"] = ess
        arm_metrics["corrected_li_median"] = weighted_quantile(subset["Li7_H_1e10"].to_numpy(float), weights, 0.5)
        if arm == "depletion":
            arm_metrics["corrected_depletion_dex_median"] = weighted_quantile(subset["delta_li_dex"].to_numpy(float), weights, 0.5)
        exact_gate &= ess >= 100
        rows.append(arm_metrics)

    overall = bool(metrics["passed"] and sampling["mcmc_gate_passed"] and exact_gate)
    payload = {"generated": now(), "promotable": overall, "holdout_gate": metrics["passed"], "mcmc_gate": sampling["mcmc_gate_passed"], "exact_gate": exact_gate, "arms": rows}
    args.archive_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.archive_dir / "bbn_lithium_hierarchical_surrogate_final_gate_2026-07-18.json"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report = [
        "# Hierarchical LINX BBN Surrogate Posterior Readout",
        "",
        f"Generated: {payload['generated']}",
        "",
        f"Overall promotable: **{overall}**",
        f"Holdout gate: `{metrics['passed']}`; MCMC gate: `{sampling['mcmc_gate_passed']}`; exact posterior-draw gate: `{exact_gate}`.",
        "",
        "| arm | exact draws | importance ESS | corrected primordial Li median | corrected depletion median (dex) |",
        "|:--|--:|--:|--:|--:|",
    ]
    for row in rows:
        dep = row.get("corrected_depletion_dex_median", math.nan)
        report.append(f"| {row['arm']} | {row['n_exact']} | {row['importance_ess']:.1f} | {row['corrected_li_median']:.6g} | {dep:.6g} |")
    report.extend([
        "",
        "## Claim Boundary",
        "",
        "This inference is conditional on homogeneous standard BBN and the registered summary likelihoods. It contains no matter-antimatter domain, annihilation or antimatter-galaxy population model.",
    ])
    report_path = args.archive_dir / "bbn_lithium_hierarchical_surrogate_readout_2026-07-18.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Saved final gate: {json_path}", flush=True)
    print(f"Saved final report: {report_path}", flush=True)


def main() -> int:
    args = parse_args()
    args.outdir = args.outdir.resolve()
    args.archive_dir = args.archive_dir.resolve()
    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.mode == "train-sample":
        train_sample(args)
    else:
        finalise(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
