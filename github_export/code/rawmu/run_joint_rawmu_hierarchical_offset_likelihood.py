r"""
Joint DES/Pantheon+/Union3.1 raw-MU hierarchical offset likelihood.

This promotes the July 14 free/profilled offset diagnostic to a random-effect
calibration model with explicit Gaussian priors:

    residual = MU_observed - MU_model(H0, beta)
    residual = X_dataset delta_dataset + X_survey delta_survey + noise
    delta_dataset ~ N(0, sigma_dataset^2)
    delta_survey  ~ N(0, sigma_survey^2)

The covariance treatment is matched to the July 14 joint runner: DES-Dovekie
raw, Pantheon+SH0ES, and Union3.1 are loaded with their available covariance or
precision products and combined as independent block-diagonal precision blocks.

Outputs:

    plamb_runs/diagnostics/joint_rawmu_hierarchical_offset_YYYYMMDD/
        joint_rawmu_hierarchical_summary.csv
        joint_rawmu_hierarchical_offsets.csv
        joint_rawmu_hierarchical_blocks.csv
        joint_rawmu_hierarchical_report.md
        joint_rawmu_hierarchical_config.json
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from fit_plamb_rawmu_nuisance import (
    HAVE_SCIPY,
    PRESETS,
    DatasetBlock,
    ModelChoice,
    combine_blocks,
    load_dataset,
    minimize,
    minimize_scalar,
    mu_model,
    select_block,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTDIR = ROOT / "plamb_runs" / "diagnostics" / f"joint_rawmu_hierarchical_offset_{datetime.now():%Y%m%d}"


FRAME_MAPS = {
    "HD": {
        "PantheonPlusSH0ES": "zHD",
        "DES_Dovekie_raw": "zHD",
        "Union3p1_UNITY1p8": "z",
    },
    "HEL": {
        "PantheonPlusSH0ES": "zHEL",
        "DES_Dovekie_raw": "zHEL",
        "Union3p1_UNITY1p8": "z",
    },
    "CMB_PANTHEON_ONLY": {
        "PantheonPlusSH0ES": "zCMB",
        "DES_Dovekie_raw": "zHD",
        "Union3p1_UNITY1p8": "z",
    },
}


MODELS = [
    ModelChoice("FR_BETA05_H0fixed675", 67.5, 0.5),
    ModelChoice("FR_BETAfree_H0fixed675", 67.5, None),
    ModelChoice("FR_BETA05_H0free", None, 0.5),
    ModelChoice("FR_BETAfree_H0free", None, None),
]


@dataclass(frozen=True)
class PriorConfig:
    name: str
    dataset_sigma: float
    survey_sigma: float


@dataclass(frozen=True)
class RandomDesign:
    mode: str
    x: np.ndarray | None
    labels: list[str]
    component_types: list[str]
    prior_sigmas: np.ndarray
    note: str


@dataclass(frozen=True)
class RandomSolver:
    precision: np.ndarray
    x: np.ndarray | None
    labels: list[str]
    component_types: list[str]
    prior_sigmas: np.ndarray
    xtp: np.ndarray | None
    normal_data: np.ndarray | None
    normal_post: np.ndarray | None
    normal_post_inv: np.ndarray | None
    lambda_inv: np.ndarray | None
    logdet_random_effect: float


def parse_list(raw: str) -> list[str]:
    return [item.strip() for item in str(raw).split(",") if item.strip()]


def parse_prior_configs(raw: str) -> list[PriorConfig]:
    configs: list[PriorConfig] = []
    for item in parse_list(raw):
        parts = item.split(":")
        if len(parts) != 3:
            raise ValueError(f"Prior config must be name:dataset_sigma:survey_sigma, got {item!r}")
        name, dataset_sigma, survey_sigma = parts
        configs.append(PriorConfig(name, float(dataset_sigma), float(survey_sigma)))
    if not configs:
        raise ValueError("At least one prior config is required.")
    return configs


def analysis_args(cli: argparse.Namespace) -> SimpleNamespace:
    return SimpleNamespace(
        mu_col=None,
        calibrator_col=None,
        z_min=cli.z_min,
        max_z=cli.max_z,
        min_n=cli.min_n,
        keep_calibrators=cli.keep_calibrators,
        min_survey_n=cli.min_survey_n,
        allow_precision_submatrix=cli.allow_precision_submatrix,
        h0_min=cli.h0_min,
        h0_max=cli.h0_max,
        beta_min=cli.beta_min,
        beta_max=cli.beta_max,
        grid_steps=cli.grid_steps,
    )


def build_blocks(frame: str, cli: argparse.Namespace) -> tuple[list[DatasetBlock], list[dict[str, object]]]:
    specs = [PRESETS["pantheon"], PRESETS["des_raw"], PRESETS["union31"]]
    z_map = FRAME_MAPS[frame]
    args = analysis_args(cli)
    blocks: list[DatasetBlock] = []
    block_rows: list[dict[str, object]] = []
    for spec in specs:
        loaded = load_dataset(spec, args)
        z_col = z_map[loaded.label]
        block = select_block(loaded, z_col, args)
        blocks.append(block)
        block_rows.append(
            {
                "frame": frame,
                "dataset": block.label,
                "z_col": block.z_col,
                "N_source": block.n_source,
                "N_used": block.n_used,
                "cov_note": block.cov_note,
                "subset_note": block.subset_note,
                "source_path": str(block.source_path),
                "cov_path": str(block.cov_path) if block.cov_path else "",
            }
        )
    return blocks, block_rows


def build_random_design(blocks: list[DatasetBlock], mode: str, prior: PriorConfig) -> RandomDesign:
    n = sum(block.n_used for block in blocks)
    if mode == "none":
        return RandomDesign("none", None, [], [], np.zeros(0), "no calibration offsets")

    row_active: list[list[str]] = []
    labels: list[str] = []
    component_type: dict[str, str] = {}
    sigma_map: dict[str, float] = {}

    def add_label(label: str, typ: str, sigma: float) -> None:
        if label not in component_type:
            labels.append(label)
            component_type[label] = typ
            sigma_map[label] = sigma

    for block in blocks:
        dataset_label = f"{block.label}:DATASET"
        has_survey = any(":IDSURVEY_" in str(value) for value in block.group_values)
        for value in block.group_values:
            active: list[str] = []
            if mode in {"dataset", "dataset+survey"}:
                add_label(dataset_label, "dataset", prior.dataset_sigma)
                active.append(dataset_label)
            if mode in {"survey", "dataset+survey"} and has_survey:
                survey_label = str(value)
                add_label(survey_label, "survey", prior.survey_sigma)
                active.append(survey_label)
            row_active.append(active)

    if not labels:
        return RandomDesign(mode, None, [], [], np.zeros(0), "requested mode produced no calibration columns")

    label_index = {label: idx for idx, label in enumerate(labels)}
    x = np.zeros((n, len(labels)), dtype=float)
    for i, active in enumerate(row_active):
        for label in active:
            x[i, label_index[label]] = 1.0

    component_types = [component_type[label] for label in labels]
    sigmas = np.array([sigma_map[label] for label in labels], dtype=float)
    note = (
        f"{len(labels)} Gaussian calibration offsets; "
        f"dataset sigma={prior.dataset_sigma:g} mag, survey sigma={prior.survey_sigma:g} mag"
    )
    return RandomDesign(mode, x, labels, component_types, sigmas, note)


def prepare_solver(precision: np.ndarray, design: RandomDesign) -> RandomSolver:
    if design.x is None or design.x.shape[1] == 0:
        return RandomSolver(
            precision=precision,
            x=None,
            labels=[],
            component_types=[],
            prior_sigmas=np.zeros(0),
            xtp=None,
            normal_data=None,
            normal_post=None,
            normal_post_inv=None,
            lambda_inv=None,
            logdet_random_effect=0.0,
        )

    x = design.x
    sigmas = design.prior_sigmas
    if np.any(sigmas <= 0.0):
        raise ValueError("All prior sigmas must be positive.")
    xtp = x.T @ precision
    normal_data = xtp @ x
    lambda_inv = np.diag(1.0 / (sigmas * sigmas))
    normal_post = normal_data + lambda_inv
    try:
        normal_post_inv = np.linalg.inv(normal_post)
    except np.linalg.LinAlgError:
        normal_post_inv = np.linalg.pinv(normal_post)
    sign_a, logdet_a = np.linalg.slogdet(normal_post)
    if sign_a <= 0:
        raise np.linalg.LinAlgError("Posterior normal matrix is not positive definite.")
    logdet_lambda = float(np.sum(np.log(sigmas * sigmas)))
    logdet_random_effect = float(logdet_lambda + logdet_a)
    return RandomSolver(
        precision=precision,
        x=x,
        labels=design.labels,
        component_types=design.component_types,
        prior_sigmas=sigmas,
        xtp=xtp,
        normal_data=normal_data,
        normal_post=normal_post,
        normal_post_inv=normal_post_inv,
        lambda_inv=lambda_inv,
        logdet_random_effect=logdet_random_effect,
    )


def score_for_params(
    z: np.ndarray,
    mu: np.ndarray,
    solver: RandomSolver,
    h0: float,
    beta: float,
) -> tuple[dict[str, float], np.ndarray, np.ndarray]:
    model = mu_model(z, h0, beta)
    if not np.all(np.isfinite(model)):
        nan_profiled = np.full_like(mu, np.nan)
        return {
            "raw_chi2": float("inf"),
            "data_chi2_map": float("inf"),
            "prior_chi2": float("inf"),
            "posterior_objective": float("inf"),
            "marginal_score": float("inf"),
            "logdet_random_effect": solver.logdet_random_effect,
        }, np.zeros(0), nan_profiled

    residual = mu - model
    precision = solver.precision
    raw_chi2 = float(residual @ precision @ residual)
    if solver.x is None or solver.xtp is None or solver.normal_post_inv is None or solver.lambda_inv is None:
        return {
            "raw_chi2": raw_chi2,
            "data_chi2_map": raw_chi2,
            "prior_chi2": 0.0,
            "posterior_objective": raw_chi2,
            "marginal_score": raw_chi2,
            "logdet_random_effect": 0.0,
        }, np.zeros(0), residual

    rhs = solver.xtp @ residual
    offsets = solver.normal_post_inv @ rhs
    profiled = residual - solver.x @ offsets
    data_chi2_map = float(profiled @ precision @ profiled)
    prior_chi2 = float(offsets @ solver.lambda_inv @ offsets)
    posterior_objective = float(data_chi2_map + prior_chi2)
    marginal_score = float(posterior_objective + solver.logdet_random_effect)
    return {
        "raw_chi2": raw_chi2,
        "data_chi2_map": data_chi2_map,
        "prior_chi2": prior_chi2,
        "posterior_objective": posterior_objective,
        "marginal_score": marginal_score,
        "logdet_random_effect": solver.logdet_random_effect,
    }, offsets, profiled


def fit_hierarchical_model(
    z: np.ndarray,
    mu: np.ndarray,
    solver: RandomSolver,
    model: ModelChoice,
    args: argparse.Namespace,
) -> tuple[dict[str, object], np.ndarray, np.ndarray]:
    h0_bounds = (args.h0_min, args.h0_max)
    beta_bounds = (args.beta_min, args.beta_max)

    def objective_pair(h0: float, beta: float) -> float:
        scores, _offsets, _profiled = score_for_params(z, mu, solver, h0, beta)
        return scores["posterior_objective"]

    if model.h0_fixed is not None and model.beta_fixed is not None:
        h0 = float(model.h0_fixed)
        beta = float(model.beta_fixed)
        scores, offsets, profiled = score_for_params(z, mu, solver, h0, beta)
        return {"H0": h0, "beta": beta, "method": "fixed", "success": np.isfinite(scores["posterior_objective"]), **scores}, offsets, profiled

    if not HAVE_SCIPY:
        return fit_hierarchical_model_grid(z, mu, solver, model, args)

    if model.h0_fixed is None and model.beta_fixed is not None:
        beta = float(model.beta_fixed)
        opt = minimize_scalar(lambda h0: objective_pair(float(h0), beta), bounds=h0_bounds, method="bounded")
        h0 = float(opt.x)
        scores, offsets, profiled = score_for_params(z, mu, solver, h0, beta)
        return {"H0": h0, "beta": beta, "method": "minimize_H0", "success": bool(opt.success), **scores}, offsets, profiled

    if model.h0_fixed is not None and model.beta_fixed is None:
        h0 = float(model.h0_fixed)
        opt = minimize_scalar(lambda beta: objective_pair(h0, float(beta)), bounds=beta_bounds, method="bounded")
        beta = float(opt.x)
        scores, offsets, profiled = score_for_params(z, mu, solver, h0, beta)
        return {"H0": h0, "beta": beta, "method": "minimize_beta", "success": bool(opt.success), **scores}, offsets, profiled

    starts = [
        np.array([67.5, 0.5]),
        np.array([68.0, 0.54]),
        np.array([71.5, 0.54]),
        np.array([73.0, 0.5]),
    ]
    best = None
    for start in starts:
        opt = minimize(
            lambda theta: objective_pair(float(theta[0]), float(theta[1])),
            start,
            method="L-BFGS-B",
            bounds=[h0_bounds, beta_bounds],
        )
        if best is None or float(opt.fun) < float(best.fun):
            best = opt
    assert best is not None
    h0 = float(best.x[0])
    beta = float(best.x[1])
    scores, offsets, profiled = score_for_params(z, mu, solver, h0, beta)
    return {"H0": h0, "beta": beta, "method": "L-BFGS-B", "success": bool(best.success), **scores}, offsets, profiled


def fit_hierarchical_model_grid(
    z: np.ndarray,
    mu: np.ndarray,
    solver: RandomSolver,
    model: ModelChoice,
    args: argparse.Namespace,
) -> tuple[dict[str, object], np.ndarray, np.ndarray]:
    h0_grid = [model.h0_fixed] if model.h0_fixed is not None else np.linspace(args.h0_min, args.h0_max, args.grid_steps)
    beta_grid = [model.beta_fixed] if model.beta_fixed is not None else np.linspace(args.beta_min, args.beta_max, args.grid_steps)
    best: tuple[float, float, float, dict[str, float], np.ndarray, np.ndarray] | None = None
    for h0 in h0_grid:
        for beta in beta_grid:
            scores, offsets, profiled = score_for_params(z, mu, solver, float(h0), float(beta))
            objective = scores["posterior_objective"]
            if best is None or objective < best[0]:
                best = (objective, float(h0), float(beta), scores, offsets, profiled)
    assert best is not None
    _objective, h0, beta, scores, offsets, profiled = best
    return {"H0": h0, "beta": beta, "method": "grid", "success": np.isfinite(scores["posterior_objective"]), **scores}, offsets, profiled


def row_for_fit(
    run_label: str,
    frame: str,
    blocks: list[DatasetBlock],
    random_mode: str,
    prior: PriorConfig,
    design: RandomDesign,
    model: ModelChoice,
    fit: dict[str, object],
) -> dict[str, object]:
    n = int(sum(block.n_used for block in blocks))
    k_random = len(design.labels)
    k_base = model.k_base
    h0 = float(fit["H0"])
    beta = float(fit["beta"])
    dof_nominal = max(n - k_base - k_random, 0)
    dof_base = max(n - k_base, 0)
    bic_random_marginal_base = float(fit["marginal_score"]) + k_base * math.log(n) if n > 0 else float("nan")
    aic_random_marginal_base = float(fit["marginal_score"]) + 2.0 * k_base
    return {
        "run_label": run_label,
        "frame": frame,
        "frame_z_map": json.dumps(FRAME_MAPS[frame], sort_keys=True),
        "datasets": "+".join(block.label for block in blocks),
        "random_mode": random_mode,
        "prior_config": prior.name,
        "dataset_sigma_mag": prior.dataset_sigma if random_mode != "none" else 0.0,
        "survey_sigma_mag": prior.survey_sigma if "survey" in random_mode else 0.0,
        "model": model.name,
        "N": n,
        "k_base": k_base,
        "k_random": k_random,
        "dof_base": dof_base,
        "dof_nominal": dof_nominal,
        "H0": h0,
        "beta": beta,
        "H0_minus_67p5": h0 - 67.5,
        "beta_minus_0p5": beta - 0.5,
        "raw_chi2_no_offsets": float(fit["raw_chi2"]),
        "data_chi2_map": float(fit["data_chi2_map"]),
        "prior_chi2": float(fit["prior_chi2"]),
        "posterior_objective": float(fit["posterior_objective"]),
        "logdet_random_effect": float(fit["logdet_random_effect"]),
        "marginal_score": float(fit["marginal_score"]),
        "chi2_dof_base": float(fit["data_chi2_map"]) / dof_base if dof_base > 0 else float("nan"),
        "chi2_dof_nominal": float(fit["data_chi2_map"]) / dof_nominal if dof_nominal > 0 else float("nan"),
        "AIC_random_marginal_base": aic_random_marginal_base,
        "BIC_random_marginal_base": bic_random_marginal_base,
        "success": bool(fit["success"]),
        "method": str(fit["method"]),
        "random_effect_note": design.note,
    }


def offset_rows(
    run_label: str,
    frame: str,
    blocks: list[DatasetBlock],
    random_mode: str,
    prior: PriorConfig,
    model_name: str,
    design: RandomDesign,
    offsets: np.ndarray,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if offsets.size == 0:
        return rows
    for i, label in enumerate(design.labels):
        sigma = float(design.prior_sigmas[i])
        offset = float(offsets[i])
        rows.append(
            {
                "run_label": run_label,
                "frame": frame,
                "datasets": "+".join(block.label for block in blocks),
                "random_mode": random_mode,
                "prior_config": prior.name,
                "model": model_name,
                "offset_group": label,
                "component_type": design.component_types[i],
                "prior_sigma_mag": sigma,
                "residual_offset_mag": offset,
                "projection_add_mag": -offset,
                "pull_sigma": offset / sigma if sigma > 0 else float("nan"),
            }
        )
    return rows


def markdown_table(df: pd.DataFrame, cols: list[str], max_rows: int | None = None, floatfmt: str = ".6g") -> str:
    if df.empty:
        return "_No rows._"
    show = df[[c for c in cols if c in df.columns]].copy()
    if max_rows is not None and len(show) > max_rows:
        show = show.head(max_rows)
    return show.to_markdown(index=False, floatfmt=floatfmt)


def write_report(
    path: Path,
    summary: pd.DataFrame,
    offsets: pd.DataFrame,
    blocks: pd.DataFrame,
    cli: argparse.Namespace,
) -> None:
    lines: list[str] = []
    lines.append("# Joint Raw-MU Hierarchical Offset Likelihood")
    lines.append("")
    lines.append(f"Date: {datetime.now():%B %d, %Y}")
    lines.append("")
    lines.append("Datasets: DES-Dovekie raw, Pantheon+SH0ES, and Union3.1 UNITY1.8 compressed MU nodes.")
    lines.append("")
    lines.append("Covariance treatment: each dataset uses its available full covariance or inverse covariance; dataset blocks are combined as independent block-diagonal precision blocks.")
    lines.append("")
    lines.append("Model: `MU = 5 log10((c/H0) z (1 + beta z)) + 25 + dataset_offset + survey_offset`.")
    lines.append("")
    lines.append("Calibration priors are zero-mean Gaussian random effects. Dataset offsets use `sigma_dataset`; nested IDSURVEY offsets use `sigma_survey`. Union3.1 has no IDSURVEY labels, so it receives only a dataset-level offset in `dataset+survey` mode.")
    lines.append("")
    lines.append("Ranking columns:")
    lines.append("")
    lines.append("- `posterior_objective = data_chi2_at_MAP + prior_chi2`.")
    lines.append("- `marginal_score` adds the Gaussian random-effect Occam/log-determinant term, up to a common covariance constant.")
    lines.append("- `BIC_random_marginal_base` additionally counts only the base H0/beta parameters, because calibration offsets have been marginalized with explicit priors.")
    lines.append("")
    lines.append("## Configuration")
    lines.append("")
    lines.append(f"- frames: `{cli.frames}`")
    lines.append(f"- random-effect modes: `{cli.random_modes}`")
    lines.append(f"- prior configs: `{cli.prior_configs}`")
    lines.append(f"- z_min: `{cli.z_min}`")
    lines.append(f"- min IDSURVEY/group size before OTHER pooling: `{cli.min_survey_n}`")
    lines.append(f"- H0/beta bounds: H0 `{cli.h0_min}` to `{cli.h0_max}`, beta `{cli.beta_min}` to `{cli.beta_max}`")
    lines.append("")
    lines.append("## Best Rows by Marginal Random-Effect BIC-Like Score")
    lines.append("")
    if summary.empty:
        lines.append("_No successful rows._")
    else:
        ranked = summary.sort_values(["BIC_random_marginal_base", "posterior_objective"]).copy()
        cols = [
            "frame",
            "random_mode",
            "prior_config",
            "model",
            "N",
            "k_base",
            "k_random",
            "dataset_sigma_mag",
            "survey_sigma_mag",
            "H0",
            "beta",
            "data_chi2_map",
            "prior_chi2",
            "marginal_score",
            "BIC_random_marginal_base",
        ]
        lines.append(markdown_table(ranked, cols, max_rows=24))
    lines.append("")
    lines.append("## H0=67.5 Beta-Free Focus")
    lines.append("")
    focus = summary[summary["model"].eq("FR_BETAfree_H0fixed675")].copy() if not summary.empty else pd.DataFrame()
    if focus.empty:
        lines.append("_No focus rows._")
    else:
        focus = focus.sort_values(["frame", "BIC_random_marginal_base", "posterior_objective"])
        cols = [
            "frame",
            "random_mode",
            "prior_config",
            "k_random",
            "dataset_sigma_mag",
            "survey_sigma_mag",
            "beta",
            "beta_minus_0p5",
            "data_chi2_map",
            "prior_chi2",
            "marginal_score",
            "BIC_random_marginal_base",
        ]
        lines.append(markdown_table(focus, cols, max_rows=40))
    lines.append("")
    lines.append("## Dataset Blocks")
    lines.append("")
    lines.append(markdown_table(blocks, list(blocks.columns) if not blocks.empty else []))
    lines.append("")
    lines.append("## Largest Offset Pulls")
    lines.append("")
    if offsets.empty:
        lines.append("_No offset rows._")
    else:
        pull = offsets.assign(abs_pull=offsets["pull_sigma"].abs()).sort_values("abs_pull", ascending=False)
        cols = [
            "frame",
            "random_mode",
            "prior_config",
            "model",
            "offset_group",
            "component_type",
            "prior_sigma_mag",
            "residual_offset_mag",
            "pull_sigma",
        ]
        lines.append(markdown_table(pull, cols, max_rows=40))
    lines.append("")
    lines.append("## Interpretation Notes")
    lines.append("")
    lines.append("- If beta near 0.5 survives only for loose calibration priors, the raw-MU signal remains calibration-dependent.")
    lines.append("- If dataset+survey priors improve data chi2 but are disfavored by the marginal random-effect score, the free-offset improvement is mostly Occam-penalized survey structure.")
    lines.append("- Large offset pulls identify which calibration blocks need direct audit before any model-promotion claim.")
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--frames", default="HD,HEL")
    parser.add_argument("--random-modes", default="none,dataset,dataset+survey")
    parser.add_argument(
        "--prior-configs",
        default="tight:0.02:0.01,moderate:0.05:0.02,loose:0.10:0.05,very_loose:0.20:0.10",
        help="Comma-separated name:dataset_sigma:survey_sigma configs in magnitudes.",
    )
    parser.add_argument("--z-min", type=float, default=0.01)
    parser.add_argument("--max-z", type=float, default=None)
    parser.add_argument("--min-n", type=int, default=10)
    parser.add_argument("--min-survey-n", type=int, default=10)
    parser.add_argument("--keep-calibrators", action="store_true")
    parser.add_argument("--allow-precision-submatrix", action="store_true")
    parser.add_argument("--h0-min", type=float, default=45.0)
    parser.add_argument("--h0-max", type=float, default=90.0)
    parser.add_argument("--beta-min", type=float, default=-0.25)
    parser.add_argument("--beta-max", type=float, default=1.25)
    parser.add_argument("--grid-steps", type=int, default=301)
    return parser


def main() -> None:
    cli = build_parser().parse_args()
    cli.outdir.mkdir(parents=True, exist_ok=True)

    frames = parse_list(cli.frames)
    random_modes = parse_list(cli.random_modes)
    prior_configs = parse_prior_configs(cli.prior_configs)

    rows: list[dict[str, object]] = []
    offset_table: list[dict[str, object]] = []
    block_rows: list[dict[str, object]] = []
    errors: list[dict[str, object]] = []

    for frame in frames:
        if frame not in FRAME_MAPS:
            errors.append({"frame": frame, "random_mode": "", "prior_config": "", "model": "", "error": f"Unknown frame {frame}"})
            continue
        try:
            blocks, these_block_rows = build_blocks(frame, cli)
            block_rows.extend(these_block_rows)
            z, mu, precision = combine_blocks(blocks)
        except Exception as exc:
            errors.append({"frame": frame, "random_mode": "", "prior_config": "", "model": "", "error": str(exc)})
            print(f"ERROR loading frame {frame}: {exc}", flush=True)
            continue

        for random_mode in random_modes:
            configs = [PriorConfig("none", 0.0, 0.0)] if random_mode == "none" else prior_configs
            for prior in configs:
                try:
                    design = build_random_design(blocks, random_mode, prior)
                    solver = prepare_solver(precision, design)
                except Exception as exc:
                    errors.append({"frame": frame, "random_mode": random_mode, "prior_config": prior.name, "model": "", "error": str(exc)})
                    print(f"ERROR design {frame} {random_mode} {prior.name}: {exc}", flush=True)
                    continue
                for model in MODELS:
                    run_label = f"joint_des_pantheon_union31_{frame}_{random_mode}_{prior.name}_{model.name}"
                    try:
                        fit, offsets, _profiled = fit_hierarchical_model(z, mu, solver, model, analysis_args(cli))
                        row = row_for_fit(run_label, frame, blocks, random_mode, prior, design, model, fit)
                        rows.append(row)
                        offset_table.extend(offset_rows(run_label, frame, blocks, random_mode, prior, model.name, design, offsets))
                        print(
                            f"{frame} {random_mode} {prior.name} {model.name}: "
                            f"post={float(row['posterior_objective']):.3f} "
                            f"marg={float(row['marginal_score']):.3f} "
                            f"beta={float(row['beta']):.5f}",
                            flush=True,
                        )
                    except Exception as exc:
                        errors.append({"frame": frame, "random_mode": random_mode, "prior_config": prior.name, "model": model.name, "error": str(exc)})
                        print(f"ERROR {frame} {random_mode} {prior.name} {model.name}: {exc}", flush=True)

    summary = pd.DataFrame(rows)
    offsets = pd.DataFrame(offset_table)
    blocks = pd.DataFrame(block_rows).drop_duplicates() if block_rows else pd.DataFrame()
    errors_df = pd.DataFrame(errors)

    summary_path = cli.outdir / "joint_rawmu_hierarchical_summary.csv"
    offsets_path = cli.outdir / "joint_rawmu_hierarchical_offsets.csv"
    blocks_path = cli.outdir / "joint_rawmu_hierarchical_blocks.csv"
    errors_path = cli.outdir / "joint_rawmu_hierarchical_errors.csv"
    report_path = cli.outdir / "joint_rawmu_hierarchical_report.md"
    config_path = cli.outdir / "joint_rawmu_hierarchical_config.json"

    summary.to_csv(summary_path, index=False)
    offsets.to_csv(offsets_path, index=False)
    blocks.to_csv(blocks_path, index=False)
    errors_df.to_csv(errors_path, index=False)
    config_path.write_text(json.dumps(vars(cli), indent=2, default=str), encoding="utf-8")
    write_report(report_path, summary, offsets, blocks, cli)

    print(f"Saved summary: {summary_path}")
    print(f"Saved offsets: {offsets_path}")
    print(f"Saved blocks: {blocks_path}")
    print(f"Saved report: {report_path}")
    if not errors_df.empty:
        print(f"Saved errors: {errors_path}")


if __name__ == "__main__":
    main()
