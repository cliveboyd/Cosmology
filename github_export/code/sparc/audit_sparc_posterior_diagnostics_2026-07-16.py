from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import norm, rankdata

ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
import analyze_sparc_posterior_run as base


OUTPUTS = Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs")
DEFAULT_RUN_DIR = (
    ROOT
    / "plamb_runs"
    / "diagnostics"
    / "sparc_hierarchical_posterior"
    / "posterior_minus_persistent_negative_rar_rescue_with_existing_plamb_20260714"
)
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "sparc_posterior_diagnostic_audit_20260716"


def fmt(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.6g}" if math.isfinite(value) else ""
    return str(value)


def markdown_table(rows: list[dict[str, Any]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(column, "")) for column in columns) + " |")
    return "\n".join(lines)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def split_chains(chains: np.ndarray) -> np.ndarray:
    n = chains.shape[1]
    half = n // 2
    if half < 2:
        return chains
    trimmed = chains[:, : 2 * half]
    return np.concatenate([trimmed[:, :half], trimmed[:, half:]], axis=0)


def classic_rhat(chains: np.ndarray) -> np.ndarray:
    m, n, _ = chains.shape
    chain_means = np.mean(chains, axis=1)
    within = np.mean(np.var(chains, axis=1, ddof=1), axis=0)
    between = n * np.var(chain_means, axis=0, ddof=1)
    variance = ((n - 1.0) / n) * within + between / n
    return np.sqrt(np.divide(variance, within, out=np.full_like(within, np.nan), where=within > 0.0))


def rank_normalize(chains: np.ndarray) -> np.ndarray:
    flat = chains.reshape(-1, chains.shape[-1])
    transformed = np.empty_like(flat, dtype=float)
    count = len(flat)
    for column in range(flat.shape[1]):
        ranks = rankdata(flat[:, column], method="average")
        probability = (ranks - 0.375) / (count + 0.25)
        transformed[:, column] = norm.ppf(probability)
    return transformed.reshape(chains.shape)


def rank_split_rhat(chains: np.ndarray) -> np.ndarray:
    split = split_chains(chains)
    ranked = rank_normalize(split)
    folded = np.abs(ranked - np.median(ranked, axis=(0, 1), keepdims=True))
    return np.maximum(classic_rhat(ranked), classic_rhat(folded))


def autocovariance(values: np.ndarray) -> np.ndarray:
    values = values - np.mean(values)
    n = len(values)
    n_fft = 1 << (2 * n - 1).bit_length()
    spectrum = np.fft.rfft(values, n=n_fft)
    return np.fft.irfft(spectrum * np.conjugate(spectrum), n=n_fft)[:n].real / n


def bulk_ess(chains: np.ndarray) -> np.ndarray:
    ranked = rank_normalize(split_chains(chains))
    m, n, p = ranked.shape
    result = np.empty(p, dtype=float)
    for column in range(p):
        autocov = np.vstack([autocovariance(ranked[chain, :, column]) for chain in range(m)])
        within = float(np.mean(autocov[:, 0]))
        between = n * float(np.var(np.mean(ranked[:, :, column], axis=1), ddof=1))
        variance = ((n - 1.0) / n) * within + between / n
        if not math.isfinite(variance) or variance <= 0.0:
            result[column] = float("nan")
            continue
        rho = 1.0 - (within - np.mean(autocov, axis=0)) / variance
        rho[0] = 1.0
        pairs: list[float] = []
        for lag in range(0, n - 1, 2):
            pair = float(rho[lag] + rho[lag + 1])
            if pair < 0.0:
                break
            if pairs and pair > pairs[-1]:
                pair = pairs[-1]
            pairs.append(pair)
        tau = max(-1.0 + 2.0 * sum(pairs), 1.0)
        result[column] = min(float(m * n) / tau, float(m * n))
    return result


def diagnostic_rows(cases: list[base.Case]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for case in cases:
        classic = classic_rhat(case.global_samples)
        robust = rank_split_rhat(case.global_samples)
        ess = bulk_ess(case.global_samples)
        flat = case.global_samples.reshape(-1, case.global_samples.shape[-1])
        for index, name in enumerate(case.param_names):
            sd = float(np.std(flat[:, index], ddof=1))
            rows.append(
                {
                    "split": case.split,
                    "model": case.model,
                    "parameter": name,
                    "chains": case.global_samples.shape[0],
                    "saved_per_chain": case.global_samples.shape[1],
                    "classic_rhat": float(classic[index]),
                    "rank_split_rhat": float(robust[index]),
                    "bulk_ess": float(ess[index]),
                    "mcse_mean_approx": sd / math.sqrt(max(float(ess[index]), 1.0)),
                }
            )
    return rows


def select_paired_draws(case: base.Case, max_draws: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    global_flat = case.global_samples.reshape(-1, case.global_samples.shape[-1])
    nuisance_flat = case.nuisance_samples.reshape(-1, case.nuisance_samples.shape[-1])
    if len(global_flat) != len(nuisance_flat):
        raise ValueError(f"Unpaired posterior shapes for {case.model}")
    if len(global_flat) <= max_draws:
        return global_flat, nuisance_flat
    rng = np.random.default_rng(seed)
    indices = np.sort(rng.choice(len(global_flat), size=max_draws, replace=False))
    return global_flat[indices], nuisance_flat[indices]


def paired_predictions(
    data: base.Dataset,
    model: str,
    h0: float,
    global_draws: np.ndarray,
    param_names: list[str],
    logd: np.ndarray,
    logeta: np.ndarray,
) -> np.ndarray:
    params = {name: global_draws[:, index] for index, name in enumerate(param_names)}
    ydisk = np.exp(params["log_ydisk"])
    ybul = np.exp(params["log_ybul"])
    scale = np.exp(logd)[:, None]
    eta = np.exp(logeta)[:, None]
    root_scale = np.sqrt(scale)
    root_eta = np.sqrt(eta)
    radius = data.rad_kpc[None, :] * scale
    vgas = data.vgas[None, :] * root_scale
    vdisk = data.vdisk[None, :] * root_scale * root_eta
    vbul = data.vbul[None, :] * root_scale * root_eta
    v2_bar = vgas * np.abs(vgas) + ydisk[:, None] * vdisk**2 + ybul[:, None] * vbul**2
    gbar = np.maximum(v2_bar, 1e-18) * 1.0e6 / (radius * base.KPC_M)
    if model == "RAR":
        g0 = 10.0 ** params["log10_gdagger"]
        tau = np.sqrt(np.maximum(gbar / g0[:, None], 1e-30))
    elif model == "PLAMB_OPTICAL_DEPTH_KAPPA_P":
        g0 = (10.0 ** params["log10_kappa"]) * base.acceleration_cH0_over_2pi(h0)
        exponent = params["bridge_exponent"]
        tau = np.maximum(gbar / g0[:, None], 1e-30) ** exponent[:, None]
    else:
        raise ValueError(f"Unsupported audited model: {model}")
    gpred = gbar / np.maximum(1.0 - np.exp(-tau), 1e-12)
    return np.sqrt(np.maximum(gpred, 0.0) * radius * base.KPC_M) / 1000.0


def score_case(
    score_name: str,
    case: base.Case,
    prep: Any,
    map_rows: dict[tuple[str, str], dict[str, str]],
    max_draws: int,
    seed: int,
) -> list[dict[str, Any]]:
    global_draws, nuisance_draws = select_paired_draws(case, max_draws, seed)
    galaxy_lookup = {galaxy: index for index, galaxy in enumerate(case.galaxies)}
    n_galaxies = len(case.galaxies)
    h0 = base.parse_float(map_rows[(case.split, case.model)]["H0_trial"])
    rows: list[dict[str, Any]] = []
    for galaxy in prep.galaxies:
        if galaxy not in galaxy_lookup:
            continue
        galaxy_index = galaxy_lookup[galaxy]
        mask = prep.data.galaxy == galaxy
        data = base.subset_dataset(prep.data, mask)
        prediction = paired_predictions(
            data,
            case.model,
            h0,
            global_draws,
            case.param_names,
            nuisance_draws[:, galaxy_index],
            nuisance_draws[:, n_galaxies + galaxy_index],
        )
        residual = (data.vobs[None, :] - prediction) / data.sigma_v[None, :]
        log_norm = np.log(2.0 * math.pi * data.sigma_v[None, :] ** 2)
        loglike = -0.5 * np.sum(residual**2 + log_norm, axis=1)
        rows.append(
            {
                "score_case": score_name,
                "model": case.model,
                "galaxy": galaxy,
                "N_points": len(data.vobs),
                "posterior_fit_log_density": base.logmeanexp(loglike),
            }
        )
    return rows


def stress_preparations(cases: list[base.Case], map_rows: dict[tuple[str, str], dict[str, str]], args: argparse.Namespace) -> dict[str, Any]:
    cfg = base.split_config("all_Q2")
    h0_values = [base.parse_float(map_rows[("all_Q2", case.model)]["H0_trial"]) for case in cases if case.split == "all_Q2"]
    prep = base.prepare(
        args.sample,
        args.points,
        int(cfg["quality_max"]),
        str(cfg["distance_method"]),
        args.err_floor_kms,
        sorted(set(h0_values)),
        args.sigma_ln_ml,
        args.distance_floor_frac,
        args.hubble_prior_mode,
    )
    gas, _disk, bulge, total = base.baryon_components(prep.data)
    gas_by_gal = base.galaxy_metric(prep.data, prep.galaxies, gas / total)
    bulge_by_gal = base.galaxy_metric(prep.data, prep.galaxies, bulge / total)
    result = {
        "paired_posterior_fit_baseline_all_Q2": prep,
        "paired_posterior_fit_gas_dominated_all_Q2": base.filter_prepared(prep, {g for g, v in gas_by_gal.items() if math.isfinite(v) and v >= args.gas_fraction_min}),
        "paired_posterior_fit_bulge_removed_all_Q2": base.filter_prepared(prep, {g for g, v in bulge_by_gal.items() if not math.isfinite(v) or v <= args.bulge_fraction_max}),
    }
    gbar = base.gbar_fiducial(prep.data)
    g0 = base.acceleration_cH0_over_2pi(65.0)
    median_radius = base.galaxy_metric(prep.data, prep.galaxies, prep.data.rad_kpc)
    outer = np.asarray([prep.data.rad_kpc[index] >= median_radius[str(galaxy)] for index, galaxy in enumerate(prep.data.galaxy)], dtype=bool)
    result["paired_posterior_fit_low_accel_outer_points_all_Q2"] = base.filter_prepared(prep, point_mask=outer & (gbar < g0))
    inclination = base.galaxy_metric(prep.data, prep.galaxies, prep.data.inc_deg)
    result["paired_posterior_fit_high_inclination_all_Q2"] = base.filter_prepared(prep, {g for g, v in inclination.items() if math.isfinite(v) and v >= args.high_inc_min})
    result["paired_posterior_fit_low_inclination_all_Q2"] = base.filter_prepared(prep, {g for g, v in inclination.items() if math.isfinite(v) and v < args.high_inc_min})
    return result


def summarise_scores(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for score_case in sorted({row["score_case"] for row in rows}):
        selected = [row for row in rows if row["score_case"] == score_case]
        rar = {row["galaxy"]: float(row["posterior_fit_log_density"]) for row in selected if row["model"] == "RAR"}
        for model in sorted({row["model"] for row in selected}):
            model_rows = [row for row in selected if row["model"] == model]
            deltas = [float(row["posterior_fit_log_density"]) - rar[row["galaxy"]] for row in model_rows if row["galaxy"] in rar and model != "RAR"]
            total = sum(float(row["posterior_fit_log_density"]) for row in model_rows)
            rar_total = sum(rar.get(row["galaxy"], 0.0) for row in model_rows)
            output.append(
                {
                    "score_case": score_case,
                    "model": model,
                    "N_galaxies": len(model_rows),
                    "N_points": sum(int(row["N_points"]) for row in model_rows),
                    "posterior_fit_log_density": total,
                    "delta_vs_RAR": total - rar_total,
                    "median_galaxy_delta_vs_RAR": float(np.median(deltas)) if deltas else 0.0,
                    "galaxies_better_than_RAR": sum(delta > 0.0 for delta in deltas),
                    "galaxies_worse_than_RAR": sum(delta < 0.0 for delta in deltas),
                }
            )
    return output


def write_report(out_dir: Path, diagnostics: list[dict[str, Any]], scores: list[dict[str, Any]], config: dict[str, Any]) -> tuple[Path, Path]:
    model_summary: list[dict[str, Any]] = []
    for model in sorted({row["model"] for row in diagnostics}):
        selected = [row for row in diagnostics if row["model"] == model]
        model_summary.append(
            {
                "model": model,
                "max_classic_rhat": max(float(row["classic_rhat"]) for row in selected),
                "max_rank_split_rhat": max(float(row["rank_split_rhat"]) for row in selected),
                "min_bulk_ess": min(float(row["bulk_ess"]) for row in selected),
            }
        )
    rar = next(row for row in model_summary if row["model"] == "RAR")
    plamb = next(row for row in model_summary if row["model"] == "PLAMB_OPTICAL_DEPTH_KAPPA_P")
    convergence = "passes" if float(rar["max_rank_split_rhat"]) <= 1.05 else "does not pass"
    report = [
        "# SPARC Posterior Diagnostic Audit",
        "",
        "Date: 2026-07-16",
        "",
        "## Bottom Line",
        "",
        f"The rescued RAR chain {convergence} the rank-normalised split-Rhat <= 1.05 gate. PLAMB has maximum rank-split Rhat {float(plamb['max_rank_split_rhat']):.4g} and is not robustly converged. The older stress table was not a true paired posterior-predictive calculation because galaxy nuisances were redrawn from their priors. The replacement table below uses stored paired global and nuisance posterior draws, but is deliberately labelled an in-sample posterior fit-density check rather than out-of-sample predictive evidence.",
        "",
        "## Robust Convergence Diagnostics",
        "",
        markdown_table(model_summary, ["model", "max_classic_rhat", "max_rank_split_rhat", "min_bulk_ess"]),
        "",
        markdown_table(diagnostics, ["model", "parameter", "classic_rhat", "rank_split_rhat", "bulk_ess", "mcse_mean_approx"]),
        "",
        "## Paired Posterior Fit-Density Checks",
        "",
        markdown_table(scores, ["score_case", "model", "N_galaxies", "N_points", "delta_vs_RAR", "median_galaxy_delta_vs_RAR", "galaxies_better_than_RAR", "galaxies_worse_than_RAR"]),
        "",
        "## Interpretation Boundary",
        "",
        "These scores reuse the fitted observations and therefore measure posterior fit, not generalisation. They can identify severe subset tension, but they cannot promote PLAMB over RAR. Their current positive values are provisional because the PLAMB global chain also fails robust convergence. A publication-grade comparison requires converged chains followed by galaxy-level cross-validation or a genuinely held-out SPARC subset with all selection and nuisance priors frozen in advance.",
        "",
        "The standing SPARC conclusion remains: subset wins, not a full-sample win. Convergence rescue removes one gate failure; any negative high-inclination or gas-dominated paired score remains an explicit gating failure.",
        "",
        "## Configuration",
        "",
        "```json",
        json.dumps(config, indent=2),
        "```",
    ]
    readout = [
        "# SPARC Posterior Diagnostic Audit Readout",
        "",
        "Date: 2026-07-16",
        "",
        f"RAR rank-normalised split-Rhat gate: {convergence}; max `{float(rar['max_rank_split_rhat']):.5g}`, minimum bulk ESS `{float(rar['min_bulk_ess']):.5g}`. PLAMB max rank-split Rhat is `{float(plamb['max_rank_split_rhat']):.5g}` and its minimum bulk ESS is `{float(plamb['min_bulk_ess']):.5g}`.",
        "",
        markdown_table(scores, ["score_case", "model", "delta_vs_RAR", "galaxies_better_than_RAR", "galaxies_worse_than_RAR"]),
        "",
        "Claim boundary: subset wins, not a full-sample win. Both convergence and out-of-sample prediction remain open gates; the paired-posterior table is an in-sample fit check only.",
    ]
    report_path = out_dir / "sparc_posterior_diagnostic_audit_report.md"
    readout_path = out_dir / "sparc_posterior_diagnostic_audit_readout.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    readout_path.write_text("\n".join(readout) + "\n", encoding="utf-8")
    return report_path, readout_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit SPARC convergence and paired posterior fit-density scores.")
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--map-dir", type=Path, default=base.DEFAULT_MAP_DIR)
    parser.add_argument("--sample", type=Path, default=base.DEFAULT_SAMPLE)
    parser.add_argument("--points", type=Path, default=base.DEFAULT_POINTS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--max-paired-draws", type=int, default=4500)
    parser.add_argument("--err-floor-kms", type=float, default=3.0)
    parser.add_argument("--sigma-ln-ml", type=float, default=0.25)
    parser.add_argument("--distance-floor-frac", type=float, default=0.03)
    parser.add_argument("--hubble-prior-mode", choices=["model_h0_rescale", "published"], default="model_h0_rescale")
    parser.add_argument("--gas-fraction-min", type=float, default=0.5)
    parser.add_argument("--bulge-fraction-max", type=float, default=0.2)
    parser.add_argument("--high-inc-min", type=float, default=60.0)
    parser.add_argument("--seed", type=int, default=270716)
    parser.add_argument("--copy-to-outputs", action=argparse.BooleanOptionalAction, default=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    cases = base.discover_cases(args.run_dir)
    if not cases:
        raise FileNotFoundError(f"No posterior cases under {args.run_dir}")
    map_rows = base.load_map_rows(args.map_dir)
    diagnostics = diagnostic_rows(cases)
    preparations = stress_preparations(cases, map_rows, args)
    score_rows: list[dict[str, Any]] = []
    for case_index, case in enumerate(cases):
        for prep_index, (score_name, prep) in enumerate(preparations.items()):
            score_rows.extend(
                score_case(
                    score_name,
                    case,
                    prep,
                    map_rows,
                    args.max_paired_draws,
                    args.seed + case_index * 1009 + prep_index * 9173,
                )
            )
    scores = summarise_scores(score_rows)
    config = {
        "date": "2026-07-16",
        "run_dir": str(args.run_dir),
        "map_dir": str(args.map_dir),
        "max_paired_draws": args.max_paired_draws,
        "seed": args.seed,
        "convergence": "rank-normalised split Rhat plus approximate rank-normalised bulk ESS",
        "fit_score": "log mean likelihood over paired stored global/nuisance posterior draws; in-sample diagnostic only",
    }
    diagnostic_path = args.out_dir / "sparc_posterior_robust_convergence.csv"
    score_path = args.out_dir / "sparc_paired_posterior_fit_summary.csv"
    detail_path = args.out_dir / "sparc_paired_posterior_fit_galaxy_scores.csv"
    config_path = args.out_dir / "sparc_posterior_diagnostic_audit_config.json"
    write_csv(diagnostic_path, diagnostics)
    write_csv(score_path, scores)
    write_csv(detail_path, score_rows)
    config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
    report_path, readout_path = write_report(args.out_dir, diagnostics, scores, config)
    result_files = [diagnostic_path, score_path, detail_path, config_path, report_path, readout_path]
    manifest_path = args.out_dir / "sparc_posterior_diagnostic_audit_manifest.csv"
    write_csv(
        manifest_path,
        [{"file": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in result_files],
    )
    result_files.append(manifest_path)
    if args.copy_to_outputs:
        OUTPUTS.mkdir(parents=True, exist_ok=True)
        for path in result_files:
            shutil.copy2(path, OUTPUTS / f"{path.stem}_2026-07-16{path.suffix}")
    print(f"Saved report: {report_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
