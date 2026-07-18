r"""FR alpha=0 calibration ladder and direct LCDM comparison.

This follow-up takes the clarified flat/no-expansion FR rule selected by the
previous dimming-sequence diagnostic,

    DL_FR(z) = (c/H0) z (1 + z/2)

and asks two promotion-gate questions:

1. Does the alpha=0 FR no-loss readout survive the 10/16/25/50 mmag survey
   calibration-budget ladder?
2. Under identical data blocks, redshift frames, covariance treatment, and
   calibration random effects, does FR no-loss beat LCDM?

Ranking columns match the July 14 raw-MU hierarchical offset audit:

    posterior_objective = data_chi2_at_MAP + prior_chi2
    marginal_score      = posterior_objective + random-effect log determinant
    BIC_base            = marginal_score + k_base log(N)

where k_base counts only cosmology shape parameters, because calibration
offsets are marginalized with explicit Gaussian priors.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

try:
    from scipy.optimize import minimize_scalar

    HAVE_SCIPY = True
except Exception:  # pragma: no cover
    minimize_scalar = None
    HAVE_SCIPY = False


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
SHARED_CODE = REPO_ROOT / "github_export" / "code" / "shared"
if str(SHARED_CODE) not in sys.path:
    sys.path.insert(0, str(SHARED_CODE))

from diagnose_pantheon_rawmu_fr import C_KMS  # noqa: E402
from fit_plamb_rawmu_nuisance import (  # noqa: E402
    PRESETS,
    DatasetBlock,
    combine_blocks,
    design_matrix,
    load_dataset,
    select_block,
)
from plamb_clock_distance import (  # noqa: E402
    PETER_LINEAR_REDSHIFT,
    clock_path_distance,
)


DEFAULT_OUTDIR = (
    REPO_ROOT
    / "plamb_runs"
    / "diagnostics"
    / "fr_c1pz_noloss_lcdm_calibration_ladder_20260716"
)

FRAME_MAPS: dict[str, dict[str, str]] = {
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


@dataclass(frozen=True)
class PriorConfig:
    name: str
    dataset_sigma: float
    survey_sigma: float


@dataclass(frozen=True)
class ModelSpec:
    name: str
    family: str
    param_name: str | None
    fixed_value: float | None
    default_value: float
    bounds: tuple[float, float]
    k_base: int


@dataclass(frozen=True)
class Solver:
    precision: np.ndarray
    x: np.ndarray | None
    labels: list[str]
    prior_sigmas: np.ndarray
    xtp: np.ndarray | None
    normal_post_inv: np.ndarray | None
    lambda_inv: np.ndarray | None
    logdet_random_effect: float


MODELS = [
    ModelSpec(
        "FR_C1PZ_ALPHA0_H0fixed675",
        "fr_c1pz_alpha0",
        None,
        None,
        0.0,
        (0.0, 0.0),
        0,
    ),
    ModelSpec(
        "LCDM_Om03_H0fixed675",
        "lcdm",
        "Om",
        0.3,
        0.3,
        (0.05, 0.6),
        0,
    ),
    ModelSpec(
        "LCDM_Omfree_H0fixed675",
        "lcdm",
        "Om",
        None,
        0.3,
        (0.05, 0.6),
        1,
    ),
]


def parse_list(raw: str) -> list[str]:
    return [item.strip() for item in str(raw).split(",") if item.strip()]


def parse_prior_configs(raw: str) -> list[PriorConfig]:
    configs: list[PriorConfig] = []
    for item in parse_list(raw):
        parts = item.split(":")
        if len(parts) != 3:
            raise ValueError(
                "Prior configs must be name:dataset_sigma:survey_sigma, "
                f"got {item!r}."
            )
        configs.append(PriorConfig(parts[0], float(parts[1]), float(parts[2])))
    if not configs:
        raise ValueError("At least one prior config is required.")
    return configs


def selected_presets(cli: argparse.Namespace):
    wanted = parse_list(cli.datasets)
    if wanted == ["all"]:
        wanted = ["pantheon", "des_raw", "union31"]
    out = []
    for key in wanted:
        if key not in PRESETS:
            raise ValueError(f"Unknown dataset preset {key!r}; valid keys include {sorted(PRESETS)}")
        out.append(PRESETS[key])
    return out


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
    )


def build_blocks(frame: str, cli: argparse.Namespace) -> tuple[list[DatasetBlock], list[dict[str, object]]]:
    z_map = FRAME_MAPS[frame]
    args = analysis_args(cli)
    blocks: list[DatasetBlock] = []
    rows: list[dict[str, object]] = []
    for spec in selected_presets(cli):
        loaded = load_dataset(spec, args)
        z_col = z_map[loaded.label]
        block = select_block(loaded, z_col, args)
        blocks.append(block)
        rows.append(
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
    return blocks, rows


def prior_sigmas_for_labels(labels: list[str], prior: PriorConfig) -> np.ndarray:
    sigmas = []
    for label in labels:
        if ":IDSURVEY_" in label:
            sigmas.append(prior.survey_sigma)
        else:
            sigmas.append(prior.dataset_sigma)
    return np.asarray(sigmas, dtype=float)


def prepare_solver(
    precision: np.ndarray,
    x: np.ndarray | None,
    labels: list[str],
    prior: PriorConfig | None,
) -> Solver:
    if x is None or x.shape[1] == 0 or prior is None:
        return Solver(
            precision=precision,
            x=x,
            labels=labels,
            prior_sigmas=np.full(len(labels), np.nan, dtype=float),
            xtp=None,
            normal_post_inv=None,
            lambda_inv=None,
            logdet_random_effect=0.0,
        )

    sigmas = prior_sigmas_for_labels(labels, prior)
    xtp = x.T @ precision
    normal_data = xtp @ x
    lambda_inv = np.diag(1.0 / np.square(sigmas))
    normal_post = normal_data + lambda_inv
    try:
        normal_post_inv = np.linalg.inv(normal_post)
    except np.linalg.LinAlgError:
        normal_post_inv = np.linalg.pinv(normal_post)

    sign_a, logdet_a = np.linalg.slogdet(normal_post)
    if sign_a <= 0:
        raise np.linalg.LinAlgError("Random-effect posterior normal matrix is not positive definite.")
    logdet_lambda = float(np.sum(np.log(np.square(sigmas))))
    logdet_random_effect = float(logdet_lambda + logdet_a)
    return Solver(
        precision=precision,
        x=x,
        labels=labels,
        prior_sigmas=sigmas,
        xtp=xtp,
        normal_post_inv=normal_post_inv,
        lambda_inv=lambda_inv,
        logdet_random_effect=logdet_random_effect,
    )


def lcdm_integral(z: np.ndarray, om: float, n_grid: int = 4096) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    if om <= 0.0 or om >= 1.5 or not np.all(np.isfinite(z)):
        return np.full_like(z, np.inf, dtype=float)
    zmax = float(np.nanmax(z)) if z.size else 0.0
    if zmax <= 0.0:
        return np.zeros_like(z)
    grid = np.linspace(0.0, zmax, max(256, int(n_grid)))
    ez = np.sqrt(om * np.power(1.0 + grid, 3.0) + (1.0 - om))
    if np.any(ez <= 0.0) or not np.all(np.isfinite(ez)):
        return np.full_like(z, np.inf, dtype=float)
    inv_e = 1.0 / ez
    step = grid[1] - grid[0]
    cumulative = np.empty_like(grid)
    cumulative[0] = 0.0
    cumulative[1:] = np.cumsum(0.5 * (inv_e[1:] + inv_e[:-1]) * step)
    return np.interp(z, grid, cumulative)


def mu_model(z: np.ndarray, h0: float, spec: ModelSpec, value: float | None) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    if spec.family == "fr_c1pz_alpha0":
        d_l = clock_path_distance(z, h0, C_KMS, PETER_LINEAR_REDSHIFT)
    elif spec.family == "lcdm":
        om = float(spec.fixed_value if spec.fixed_value is not None else value)
        d_l = (C_KMS / h0) * (1.0 + z) * lcdm_integral(z, om)
    else:
        raise ValueError(f"Unknown model family {spec.family!r}")
    if np.any(d_l <= 0.0) or not np.all(np.isfinite(d_l)):
        return np.full_like(z, np.inf, dtype=float)
    return 5.0 * np.log10(d_l) + 25.0


def score_for_value(
    z: np.ndarray,
    mu: np.ndarray,
    solver: Solver,
    h0: float,
    spec: ModelSpec,
    value: float | None,
) -> tuple[dict[str, float], np.ndarray, np.ndarray]:
    model = mu_model(z, h0, spec, value)
    if not np.all(np.isfinite(model)):
        return {
            "raw_chi2": float("inf"),
            "data_chi2_map": float("inf"),
            "prior_chi2": float("inf"),
            "posterior_objective": float("inf"),
            "marginal_score": float("inf"),
            "logdet_random_effect": solver.logdet_random_effect,
        }, np.zeros(0, dtype=float), np.full_like(mu, np.nan)

    residual = mu - model
    raw_chi2 = float(residual @ solver.precision @ residual)
    if (
        solver.x is None
        or solver.xtp is None
        or solver.normal_post_inv is None
        or solver.lambda_inv is None
    ):
        return {
            "raw_chi2": raw_chi2,
            "data_chi2_map": raw_chi2,
            "prior_chi2": 0.0,
            "posterior_objective": raw_chi2,
            "marginal_score": raw_chi2,
            "logdet_random_effect": 0.0,
        }, np.zeros(0, dtype=float), residual

    rhs = solver.xtp @ residual
    offsets = solver.normal_post_inv @ rhs
    profiled = residual - solver.x @ offsets
    data_chi2_map = float(profiled @ solver.precision @ profiled)
    prior_chi2 = float(offsets @ solver.lambda_inv @ offsets)
    posterior_objective = data_chi2_map + prior_chi2
    marginal_score = posterior_objective + solver.logdet_random_effect
    return {
        "raw_chi2": raw_chi2,
        "data_chi2_map": data_chi2_map,
        "prior_chi2": prior_chi2,
        "posterior_objective": posterior_objective,
        "marginal_score": marginal_score,
        "logdet_random_effect": solver.logdet_random_effect,
    }, offsets, profiled


def fit_model(
    z: np.ndarray,
    mu: np.ndarray,
    solver: Solver,
    h0: float,
    spec: ModelSpec,
    grid_steps: int,
) -> tuple[dict[str, object], np.ndarray, np.ndarray]:
    if spec.fixed_value is not None or spec.param_name is None:
        value = spec.fixed_value
        scores, offsets, profiled = score_for_value(z, mu, solver, h0, spec, value)
        return {
            "param_value": float(value) if value is not None else float("nan"),
            "method": "fixed",
            "success": np.isfinite(scores["posterior_objective"]),
            **scores,
        }, offsets, profiled

    def objective(v: float) -> float:
        scores, _offsets, _profiled = score_for_value(z, mu, solver, h0, spec, float(v))
        return scores["posterior_objective"]

    if HAVE_SCIPY:
        opt = minimize_scalar(objective, bounds=spec.bounds, method="bounded", options={"xatol": 1e-7})
        value = float(opt.x)
        success = bool(opt.success)
        method = f"minimize_{spec.param_name}"
    else:
        grid = np.linspace(spec.bounds[0], spec.bounds[1], grid_steps)
        vals = np.array([objective(v) for v in grid], dtype=float)
        value = float(grid[int(np.nanargmin(vals))])
        success = bool(np.isfinite(vals).any())
        method = f"grid_{spec.param_name}"

    scores, offsets, profiled = score_for_value(z, mu, solver, h0, spec, value)
    return {
        "param_value": value,
        "method": method,
        "success": bool(success and np.isfinite(scores["posterior_objective"])),
        **scores,
    }, offsets, profiled


def run_sequence(cli: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    frames = parse_list(cli.frames)
    offset_modes = parse_list(cli.offset_modes)
    priors = parse_prior_configs(cli.prior_configs)
    summary_rows: list[dict[str, object]] = []
    offset_rows: list[dict[str, object]] = []
    block_rows: list[dict[str, object]] = []

    for frame in frames:
        if frame not in FRAME_MAPS:
            raise ValueError(f"Unknown frame {frame!r}; valid frames are {sorted(FRAME_MAPS)}")
        blocks, rows = build_blocks(frame, cli)
        block_rows.extend(rows)
        z, mu, precision = combine_blocks(blocks)
        n = len(z)
        print(f"[frame] {frame}: N={n}", flush=True)

        for offset_mode in offset_modes:
            x, labels, design_note = design_matrix(blocks, offset_mode)
            mode_priors: list[PriorConfig | None] = [None] if offset_mode == "none" else priors
            for prior in mode_priors:
                prior_name = "none" if prior is None else prior.name
                solver = prepare_solver(precision, x, labels, prior)
                for spec in MODELS:
                    fit, offsets, profiled = fit_model(z, mu, solver, cli.h0, spec, cli.grid_steps)
                    k_base = int(spec.k_base)
                    bic_base = float(fit["marginal_score"]) + k_base * math.log(n)
                    aic_base = float(fit["marginal_score"]) + 2.0 * k_base
                    dof_base = max(n - k_base, 1)
                    max_pull = float("nan")
                    if prior is not None and offsets.size and solver.prior_sigmas.size:
                        max_pull = float(np.nanmax(np.abs(offsets / solver.prior_sigmas)))
                    rms = float(np.sqrt(np.nanmean(np.square(profiled)))) if profiled.size else float("nan")
                    summary_rows.append(
                        {
                            "run_date": datetime.now().isoformat(timespec="seconds"),
                            "frame": frame,
                            "datasets": "+".join(block.label for block in blocks),
                            "offset_mode": offset_mode,
                            "offset_design_note": design_note,
                            "prior_config": prior_name,
                            "dataset_sigma_mag": float(prior.dataset_sigma) if prior else 0.0,
                            "survey_sigma_mag": float(prior.survey_sigma) if prior else 0.0,
                            "model": spec.name,
                            "family": spec.family,
                            "H0": float(cli.h0),
                            "param_name": spec.param_name or "",
                            "param_value": float(fit["param_value"]),
                            "k_base": k_base,
                            "N": n,
                            "raw_chi2": float(fit["raw_chi2"]),
                            "data_chi2_map": float(fit["data_chi2_map"]),
                            "prior_chi2": float(fit["prior_chi2"]),
                            "posterior_objective": float(fit["posterior_objective"]),
                            "logdet_random_effect": float(fit["logdet_random_effect"]),
                            "marginal_score": float(fit["marginal_score"]),
                            "AIC_base": aic_base,
                            "BIC_base": bic_base,
                            "chi2_dof_base": float(fit["data_chi2_map"]) / dof_base,
                            "profiled_rms_mag": rms,
                            "max_abs_offset_pull_sigma": max_pull,
                            "method": str(fit["method"]),
                            "success": bool(fit["success"]),
                        }
                    )
                    for label, offset, sigma in zip(labels, offsets, solver.prior_sigmas):
                        offset_rows.append(
                            {
                                "frame": frame,
                                "offset_mode": offset_mode,
                                "prior_config": prior_name,
                                "model": spec.name,
                                "offset_label": label,
                                "offset_type": "survey" if ":IDSURVEY_" in label else "dataset",
                                "prior_sigma_mag": float(sigma),
                                "offset_mag": float(offset),
                                "pull_sigma": float(offset / sigma) if np.isfinite(sigma) and sigma > 0 else float("nan"),
                            }
                        )

    return pd.DataFrame(summary_rows), pd.DataFrame(offset_rows), pd.DataFrame(block_rows)


def markdown_table(df: pd.DataFrame, cols: list[str], max_rows: int | None = None) -> str:
    if df.empty:
        return "_No rows._"
    out = df[cols].copy()
    if max_rows is not None:
        out = out.head(max_rows)
    return out.to_markdown(index=False, floatfmt=".6g")


def comparison_rows(summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    keys = ["frame", "offset_mode", "prior_config"]
    for key, grp in summary.groupby(keys):
        fr = grp[grp["model"] == "FR_C1PZ_ALPHA0_H0fixed675"]
        lcdm_free = grp[grp["model"] == "LCDM_Omfree_H0fixed675"]
        lcdm_03 = grp[grp["model"] == "LCDM_Om03_H0fixed675"]
        if fr.empty:
            continue
        fr_row = fr.iloc[0]
        row: dict[str, object] = {
            "frame": key[0],
            "offset_mode": key[1],
            "prior_config": key[2],
            "N": int(fr_row["N"]),
            "FR_BIC_base": float(fr_row["BIC_base"]),
            "FR_posterior": float(fr_row["posterior_objective"]),
            "FR_max_pull": float(fr_row["max_abs_offset_pull_sigma"]),
        }
        if not lcdm_free.empty:
            lrow = lcdm_free.iloc[0]
            row.update(
                {
                    "LCDMfree_Om": float(lrow["param_value"]),
                    "LCDMfree_BIC_base": float(lrow["BIC_base"]),
                    "delta_BIC_FR_minus_LCDMfree": float(fr_row["BIC_base"] - lrow["BIC_base"]),
                    "delta_posterior_FR_minus_LCDMfree": float(fr_row["posterior_objective"] - lrow["posterior_objective"]),
                }
            )
        if not lcdm_03.empty:
            lrow = lcdm_03.iloc[0]
            row.update(
                {
                    "LCDM03_BIC_base": float(lrow["BIC_base"]),
                    "delta_BIC_FR_minus_LCDM03": float(fr_row["BIC_base"] - lrow["BIC_base"]),
                }
            )
        rows.append(row)
    return pd.DataFrame(rows)


def write_report(
    path: Path,
    summary: pd.DataFrame,
    offsets: pd.DataFrame,
    blocks: pd.DataFrame,
    comp: pd.DataFrame,
    cli: argparse.Namespace,
) -> None:
    lines: list[str] = []
    lines.append("# FR c(z)=c0(1+z) Alpha=0 Calibration Ladder vs LCDM")
    lines.append("")
    lines.append(f"Date: {datetime.now():%Y-%m-%d %H:%M:%S}")
    lines.append("")
    lines.append("## Test")
    lines.append("")
    lines.append("FR no-loss model:")
    lines.append("")
    lines.append("`DL_FR = (c/H0) z (1 + z/2)`")
    lines.append("")
    lines.append("LCDM comparison models:")
    lines.append("")
    lines.append("- `LCDM_Om03_H0fixed675`: flat LCDM with Om fixed at 0.3.")
    lines.append("- `LCDM_Omfree_H0fixed675`: flat LCDM with Om fitted.")
    lines.append("")
    lines.append(f"H0 fixed at `{cli.h0}` km/s/Mpc for all rows. Calibration offsets and priors are identical within each comparison cell.")
    lines.append("")

    if not blocks.empty:
        lines.append("## Loaded Blocks")
        lines.append("")
        lines.append(markdown_table(blocks, ["frame", "dataset", "z_col", "N_used", "cov_note", "subset_note"]))
        lines.append("")

    hd_ref = summary[
        (summary["frame"] == "HD")
        & (summary["offset_mode"] == "dataset+idsurvey")
        & (summary["model"] == "FR_C1PZ_ALPHA0_H0fixed675")
    ].sort_values("survey_sigma_mag")
    lines.append("## FR Alpha=0 Calibration-Budget Ladder")
    lines.append("")
    lines.append(markdown_table(
        hd_ref,
        [
            "prior_config",
            "dataset_sigma_mag",
            "survey_sigma_mag",
            "data_chi2_map",
            "prior_chi2",
            "posterior_objective",
            "marginal_score",
            "BIC_base",
            "profiled_rms_mag",
            "max_abs_offset_pull_sigma",
        ],
    ))
    lines.append("")

    lines.append("## Direct HD Model Comparison")
    lines.append("")
    hd_comp = comp[
        (comp["frame"] == "HD")
        & (comp["offset_mode"] == "dataset+idsurvey")
    ].sort_values("prior_config")
    lines.append(markdown_table(
        hd_comp,
        [
            "prior_config",
            "N",
            "FR_BIC_base",
            "LCDMfree_Om",
            "LCDMfree_BIC_base",
            "delta_BIC_FR_minus_LCDMfree",
            "delta_posterior_FR_minus_LCDMfree",
            "delta_BIC_FR_minus_LCDM03",
        ],
    ))
    lines.append("")
    lines.append("Negative `delta_BIC_FR_minus_LCDMfree` favors FR; positive favors LCDM.")
    lines.append("")

    lines.append("## Best Model By Cell")
    lines.append("")
    best = summary.sort_values("BIC_base").groupby(
        ["frame", "offset_mode", "prior_config"], as_index=False
    ).head(1)
    best = best.sort_values(["frame", "offset_mode", "prior_config"])
    lines.append(markdown_table(
        best,
        [
            "frame",
            "offset_mode",
            "prior_config",
            "model",
            "param_value",
            "BIC_base",
            "posterior_objective",
            "marginal_score",
        ],
    ))
    lines.append("")

    lines.append("## Largest HD Offset Pulls")
    lines.append("")
    pulls = offsets[
        (offsets["frame"] == "HD")
        & (offsets["offset_mode"] == "dataset+idsurvey")
        & (offsets["model"] == "FR_C1PZ_ALPHA0_H0fixed675")
        & np.isfinite(offsets["pull_sigma"])
    ].copy()
    if not pulls.empty:
        pulls["abs_pull"] = pulls["pull_sigma"].abs()
        pulls = pulls.sort_values("abs_pull", ascending=False)
        lines.append(markdown_table(
            pulls,
            [
                "prior_config",
                "offset_label",
                "offset_type",
                "prior_sigma_mag",
                "offset_mag",
                "pull_sigma",
            ],
            max_rows=20,
        ))
    else:
        lines.append("_No finite pull rows._")
    lines.append("")

    lines.append("## Readout Boundary")
    lines.append("")
    lines.append("This is still a promotion-gate diagnostic. A stable FR alpha=0 preference over added redshift-dimming factors does not by itself establish that FR beats LCDM. The decisive comparison here is the identical-offset FR-vs-LCDM BIC/base-score table, interpreted together with the calibration-budget dependence and offset pulls.")
    lines.append("")
    lines.append("## Output Files")
    lines.append("")
    lines.append("- `fr_c1pz_noloss_lcdm_ladder_summary.csv`")
    lines.append("- `fr_c1pz_noloss_lcdm_ladder_comparison.csv`")
    lines.append("- `fr_c1pz_noloss_lcdm_ladder_offsets.csv`")
    lines.append("- `fr_c1pz_noloss_lcdm_ladder_blocks.csv`")
    lines.append("- `fr_c1pz_noloss_lcdm_ladder_config.json`")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_readout(path: Path, comp: pd.DataFrame, summary: pd.DataFrame) -> None:
    lines: list[str] = []
    lines.append("# FR Alpha=0 Calibration Ladder vs LCDM Readout")
    lines.append("")
    lines.append("Date: 2026-07-16")
    lines.append("")
    lines.append("## Executive readout")
    lines.append("")
    hd = comp[(comp["frame"] == "HD") & (comp["offset_mode"] == "dataset+idsurvey")].copy()
    if hd.empty:
        lines.append("No HD dataset+survey comparison rows were produced.")
    else:
        bic_wins = (hd["delta_BIC_FR_minus_LCDMfree"] < 0).sum()
        posterior_wins = (hd["delta_posterior_FR_minus_LCDMfree"] < 0).sum()
        lines.append(f"FR alpha=0 beats LCDM Om-free by base BIC in `{bic_wins}` of `{len(hd)}` HD calibration-budget rows and by raw posterior objective in `{posterior_wins}` of `{len(hd)}` rows.")
        stress = hd[hd["prior_config"].astype(str).str.contains("050mmag", na=False)]
        if not stress.empty:
            row = stress.iloc[0]
            lines.append(
                "At the loosest 50/100 mmag stress budget, the distinction matters: "
                f"FR is better by base BIC (`{float(row['delta_BIC_FR_minus_LCDMfree']):.3f}`), "
                f"but LCDM has the slightly better raw posterior objective (`FR-LCDM = {float(row['delta_posterior_FR_minus_LCDMfree']):.3f}`)."
            )
        lines.append("")
        lines.append(markdown_table(
            hd.sort_values("prior_config"),
            [
                "prior_config",
                "FR_BIC_base",
                "LCDMfree_Om",
                "LCDMfree_BIC_base",
                "delta_BIC_FR_minus_LCDMfree",
                "delta_posterior_FR_minus_LCDMfree",
            ],
        ))
    lines.append("")
    lines.append("## FR calibration ladder")
    lines.append("")
    fr = summary[
        (summary["frame"] == "HD")
        & (summary["offset_mode"] == "dataset+idsurvey")
        & (summary["model"] == "FR_C1PZ_ALPHA0_H0fixed675")
    ].copy()
    lines.append(markdown_table(
        fr.sort_values("survey_sigma_mag"),
        [
            "prior_config",
            "data_chi2_map",
            "prior_chi2",
            "BIC_base",
            "profiled_rms_mag",
            "max_abs_offset_pull_sigma",
        ],
    ))
    lines.append("")
    lines.append("## Claim boundary")
    lines.append("")
    lines.append("A negative FR-minus-LCDM delta favors the clarified FR no-loss distance law under identical offsets. The result should still be treated as calibration-gated, because the tightest budgets expose large offset pulls.")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--datasets", default="pantheon,des_raw,union31")
    parser.add_argument("--frames", default="HD,HEL,CMB_PANTHEON_ONLY")
    parser.add_argument("--offset-modes", default="none,dataset+idsurvey")
    parser.add_argument(
        "--prior-configs",
        default=(
            "budget_010mmag_ds020:0.020:0.010,"
            "budget_016mmag_ds032:0.032:0.016,"
            "budget_025mmag_ds050:0.050:0.025,"
            "budget_050mmag_ds100_stress:0.100:0.050"
        ),
    )
    parser.add_argument("--h0", type=float, default=67.5)
    parser.add_argument("--z-min", type=float, default=0.01)
    parser.add_argument("--max-z", type=float, default=None)
    parser.add_argument("--min-n", type=int, default=10)
    parser.add_argument("--min-survey-n", type=int, default=20)
    parser.add_argument("--keep-calibrators", action="store_true")
    parser.add_argument("--allow-precision-submatrix", action="store_true")
    parser.add_argument("--grid-steps", type=int, default=401)
    return parser


def main() -> None:
    parser = build_parser()
    cli = parser.parse_args()
    cli.outdir.mkdir(parents=True, exist_ok=True)

    summary, offsets, blocks = run_sequence(cli)
    comp = comparison_rows(summary)

    summary_path = cli.outdir / "fr_c1pz_noloss_lcdm_ladder_summary.csv"
    comp_path = cli.outdir / "fr_c1pz_noloss_lcdm_ladder_comparison.csv"
    offsets_path = cli.outdir / "fr_c1pz_noloss_lcdm_ladder_offsets.csv"
    blocks_path = cli.outdir / "fr_c1pz_noloss_lcdm_ladder_blocks.csv"
    report_path = cli.outdir / "fr_c1pz_noloss_lcdm_ladder_report.md"
    readout_path = cli.outdir / "fr_c1pz_noloss_lcdm_ladder_readout.md"
    config_path = cli.outdir / "fr_c1pz_noloss_lcdm_ladder_config.json"

    summary.to_csv(summary_path, index=False)
    comp.to_csv(comp_path, index=False)
    offsets.to_csv(offsets_path, index=False)
    blocks.to_csv(blocks_path, index=False)
    write_report(report_path, summary, offsets, blocks, comp, cli)
    write_readout(readout_path, comp, summary)
    config_path.write_text(
        json.dumps(
            {
                "script": str(SCRIPT_PATH),
                "repo_root": str(REPO_ROOT),
                "outdir": str(cli.outdir),
                "datasets": cli.datasets,
                "frames": cli.frames,
                "offset_modes": cli.offset_modes,
                "prior_configs": cli.prior_configs,
                "h0": cli.h0,
                "z_min": cli.z_min,
                "max_z": cli.max_z,
                "min_n": cli.min_n,
                "min_survey_n": cli.min_survey_n,
                "models": [m.name for m in MODELS],
                "claim_boundary": "promotion-gate diagnostic only",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Saved summary: {summary_path}")
    print(f"Saved comparison: {comp_path}")
    print(f"Saved report: {report_path}")
    if not comp.empty:
        hd = comp[(comp["frame"] == "HD") & (comp["offset_mode"] == "dataset+idsurvey")]
        for row in hd.sort_values("prior_config").itertuples(index=False):
            print(
                f"HD {row.prior_config}: delta_BIC_FR_minus_LCDMfree="
                f"{float(row.delta_BIC_FR_minus_LCDMfree):.3f} "
                f"Om={float(row.LCDMfree_Om):.4f}"
            )


if __name__ == "__main__":
    main()
