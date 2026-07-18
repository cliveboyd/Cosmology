r"""FR c(z)=c0(1+z) no-loss dimming sequence for raw supernova MU.

This diagnostic tests Clive/Peter's clarified flat, no-expansion FR rule:

    c(z) = c0 (1 + z)
    r(z) = (c0 / H0) integral_0^z (1 + z') dz'
         = (c0 / H0) z (1 + z/2)

The test keeps the path term fixed and compares possible luminosity-distance
exponents:

    DL_alpha(z) = r(z) (1 + z)^alpha

where alpha=0 is pure no post-emission photon-energy loss/no extra dimming,
alpha=0.5 is one redshift-rate flux factor, and alpha=1 is two redshift flux
factors on the same flat path.

The loader and covariance handling are intentionally matched to the July 14
raw-MU audit via the live `fit_plamb_rawmu_nuisance.py` helpers.
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
    / "fr_c1pz_noloss_dimming_sequence_20260716"
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

ALPHA_LABELS = {
    0.0: "alpha0_no_extra_dimming",
    0.5: "alpha0p5_one_flux_factor",
    1.0: "alpha1_two_flux_factors",
}


@dataclass(frozen=True)
class PriorConfig:
    name: str
    dataset_sigma: float
    survey_sigma: float


def parse_list(raw: str) -> list[str]:
    return [item.strip() for item in str(raw).split(",") if item.strip()]


def parse_float_list(raw: str) -> list[float]:
    values = [float(item) for item in parse_list(raw)]
    if not values:
        raise ValueError("At least one alpha value is required.")
    return values


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


def mu_model_c1pz_noloss(z: np.ndarray, h0: float, alpha: float) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    path = clock_path_distance(z, h0, C_KMS, PETER_LINEAR_REDSHIFT)
    d_l = path * np.power(np.maximum(1.0 + z, 1e-12), float(alpha))
    if np.any(d_l <= 0.0) or not np.all(np.isfinite(d_l)):
        return np.full_like(z, np.inf, dtype=float)
    return 5.0 * np.log10(d_l) + 25.0


def prior_sigmas_for_labels(labels: list[str], prior: PriorConfig) -> np.ndarray:
    sigmas = []
    for label in labels:
        if ":IDSURVEY_" in label:
            sigmas.append(prior.survey_sigma)
        else:
            sigmas.append(prior.dataset_sigma)
    return np.asarray(sigmas, dtype=float)


def solve_offsets_with_optional_priors(
    residual: np.ndarray,
    precision: np.ndarray,
    x: np.ndarray | None,
    labels: list[str],
    prior: PriorConfig | None,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    if x is None or x.shape[1] == 0:
        data_chi2 = float(residual @ precision @ residual)
        return np.zeros(0, dtype=float), residual, {
            "data_chi2": data_chi2,
            "prior_chi2": 0.0,
            "posterior_objective": data_chi2,
            "logdet_random_effect": 0.0,
        }

    xtp = x.T @ precision
    normal = xtp @ x
    rhs = xtp @ residual

    if prior is not None:
        sigmas = prior_sigmas_for_labels(labels, prior)
        lambda_inv = np.diag(1.0 / np.square(sigmas))
        normal = normal + lambda_inv
    else:
        sigmas = np.full(len(labels), np.inf, dtype=float)
        lambda_inv = np.zeros_like(normal)

    try:
        offsets = np.linalg.solve(normal, rhs)
    except np.linalg.LinAlgError:
        offsets = np.linalg.pinv(normal) @ rhs

    profiled = residual - x @ offsets
    data_chi2 = float(profiled @ precision @ profiled)
    if prior is None:
        prior_chi2 = 0.0
        logdet_random_effect = 0.0
    else:
        prior_chi2 = float(offsets @ lambda_inv @ offsets)
        sign, logdet = np.linalg.slogdet(normal)
        logdet_random_effect = float(logdet) if sign > 0 else float("nan")
    return offsets, profiled, {
        "data_chi2": data_chi2,
        "prior_chi2": prior_chi2,
        "posterior_objective": data_chi2 + prior_chi2,
        "logdet_random_effect": logdet_random_effect,
    }


def score_model(
    z: np.ndarray,
    mu: np.ndarray,
    precision: np.ndarray,
    x: np.ndarray | None,
    labels: list[str],
    h0: float,
    alpha: float,
    prior: PriorConfig | None,
) -> tuple[dict[str, float], np.ndarray, np.ndarray]:
    model = mu_model_c1pz_noloss(z, h0, alpha)
    if not np.all(np.isfinite(model)):
        scores = {
            "data_chi2": float("inf"),
            "prior_chi2": float("nan"),
            "posterior_objective": float("inf"),
            "logdet_random_effect": float("nan"),
        }
        return scores, np.zeros(0, dtype=float), np.full_like(mu, np.nan)
    residual = mu - model
    offsets, profiled, scores = solve_offsets_with_optional_priors(
        residual, precision, x, labels, prior
    )
    return scores, offsets, profiled


def offset_rows(
    frame: str,
    offset_mode: str,
    prior_name: str,
    alpha: float,
    labels: list[str],
    offsets: np.ndarray,
    prior: PriorConfig | None,
) -> list[dict[str, object]]:
    sigmas = (
        prior_sigmas_for_labels(labels, prior)
        if prior is not None and labels
        else np.full(len(labels), np.nan, dtype=float)
    )
    rows: list[dict[str, object]] = []
    for label, offset, sigma in zip(labels, offsets, sigmas):
        rows.append(
            {
                "frame": frame,
                "offset_mode": offset_mode,
                "prior_config": prior_name,
                "alpha": alpha,
                "offset_label": label,
                "offset_type": "survey" if ":IDSURVEY_" in label else "dataset",
                "prior_sigma_mag": float(sigma),
                "offset_mag": float(offset),
                "pull_sigma": float(offset / sigma) if np.isfinite(sigma) and sigma > 0 else float("nan"),
            }
        )
    return rows


def run_sequence(cli: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    frames = parse_list(cli.frames)
    alphas = parse_float_list(cli.alphas)
    prior_configs = parse_prior_configs(cli.prior_configs)
    summary_rows: list[dict[str, object]] = []
    offset_out: list[dict[str, object]] = []
    block_rows: list[dict[str, object]] = []

    for frame in frames:
        if frame not in FRAME_MAPS:
            raise ValueError(f"Unknown frame {frame!r}; valid frames are {sorted(FRAME_MAPS)}")
        blocks, rows = build_blocks(frame, cli)
        block_rows.extend(rows)
        z, mu, precision = combine_blocks(blocks)
        n = len(z)
        print(f"[frame] {frame}: N={n}", flush=True)

        offset_modes = parse_list(cli.offset_modes)
        for offset_mode in offset_modes:
            x, labels, design_note = design_matrix(blocks, offset_mode)
            mode_priors: list[PriorConfig | None]
            if offset_mode == "none":
                mode_priors = [None]
            else:
                mode_priors = prior_configs

            for prior in mode_priors:
                prior_name = "none" if prior is None else prior.name
                for alpha in alphas:
                    scores, offsets, profiled = score_model(
                        z,
                        mu,
                        precision,
                        x,
                        labels,
                        cli.h0,
                        alpha,
                        prior,
                    )
                    k_base = 0
                    k_offsets = 0 if x is None else int(x.shape[1])
                    dof_data = max(n - k_base - (0 if prior is not None else k_offsets), 1)
                    bic_base = scores["posterior_objective"] + k_base * math.log(n)
                    rms_profiled = float(np.sqrt(np.nanmean(np.square(profiled)))) if profiled.size else float("nan")
                    max_abs_pull = float("nan")
                    if prior is not None and len(offsets):
                        sigmas = prior_sigmas_for_labels(labels, prior)
                        max_abs_pull = float(np.nanmax(np.abs(offsets / sigmas)))
                    summary_rows.append(
                        {
                            "run_date": datetime.now().isoformat(timespec="seconds"),
                            "model_family": "FR_C1PZ_NOLOSS",
                            "frame": frame,
                            "datasets": "+".join(block.label for block in blocks),
                            "offset_mode": offset_mode,
                            "offset_design_note": design_note,
                            "prior_config": prior_name,
                            "H0": float(cli.h0),
                            "alpha": float(alpha),
                            "alpha_label": ALPHA_LABELS.get(float(alpha), f"alpha{alpha:g}"),
                            "N": n,
                            "k_offsets": k_offsets,
                            "data_chi2": scores["data_chi2"],
                            "prior_chi2": scores["prior_chi2"],
                            "posterior_objective": scores["posterior_objective"],
                            "chi2_dof_data": scores["data_chi2"] / dof_data,
                            "bic_base": bic_base,
                            "profiled_rms_mag": rms_profiled,
                            "max_abs_offset_pull_sigma": max_abs_pull,
                            "path_distance": "(c/H0) * z * (1 + z/2)",
                            "luminosity_distance": "(c/H0) * z * (1 + z/2) * (1+z)^alpha",
                            "physics_note": (
                                "No post-emission photon energy loss; alpha controls any residual "
                                "flux/time-rate dimming factor applied to the flat c(z)=c0(1+z) path."
                            ),
                        }
                    )
                    offset_out.extend(
                        offset_rows(frame, offset_mode, prior_name, alpha, labels, offsets, prior)
                    )

    return pd.DataFrame(summary_rows), pd.DataFrame(offset_out), pd.DataFrame(block_rows)


def markdown_table(df: pd.DataFrame, cols: list[str], max_rows: int | None = None) -> str:
    if df.empty:
        return "_No rows._"
    out = df[cols].copy()
    if max_rows is not None:
        out = out.head(max_rows)
    return out.to_markdown(index=False, floatfmt=".6g")


def write_report(
    path: Path,
    summary: pd.DataFrame,
    offsets: pd.DataFrame,
    blocks: pd.DataFrame,
    cli: argparse.Namespace,
) -> None:
    lines: list[str] = []
    lines.append("# FR c(z)=c0(1+z) No-Loss Dimming Sequence")
    lines.append("")
    lines.append(f"Date: {datetime.now():%Y-%m-%d %H:%M:%S}")
    lines.append("")
    lines.append("## Model Under Test")
    lines.append("")
    lines.append("Flat/no-expansion path with no post-emission photon-energy loss:")
    lines.append("")
    lines.append("`r(z) = (c/H0) * integral_0^z (1+z') dz' = (c/H0) * z * (1+z/2)`")
    lines.append("")
    lines.append("The sequence tests:")
    lines.append("")
    lines.append("`DL_alpha(z) = r(z) * (1+z)^alpha`")
    lines.append("")
    lines.append("- `alpha=0`: no extra redshift dimming beyond geometric distance.")
    lines.append("- `alpha=0.5`: one flux/rate redshift factor.")
    lines.append("- `alpha=1`: two redshift flux factors on the same flat path.")
    lines.append("")
    lines.append(f"H0 fixed at `{cli.h0}` km/s/Mpc. Datasets and redshift frames follow the July 14 raw-MU audit policy.")
    lines.append("")

    if not blocks.empty:
        lines.append("## Loaded Blocks")
        lines.append("")
        lines.append(markdown_table(blocks, ["frame", "dataset", "z_col", "N_used", "cov_note", "subset_note"]))
        lines.append("")

    if not summary.empty:
        lines.append("## Primary HD Reference")
        lines.append("")
        hd = summary[(summary["frame"] == "HD") & (summary["offset_mode"] == "dataset+idsurvey")]
        if not hd.empty:
            hd = hd.sort_values(["prior_config", "posterior_objective"])
            lines.append(markdown_table(
                hd,
                [
                    "prior_config",
                    "alpha",
                    "N",
                    "data_chi2",
                    "prior_chi2",
                    "posterior_objective",
                    "chi2_dof_data",
                    "profiled_rms_mag",
                    "max_abs_offset_pull_sigma",
                ],
            ))
        else:
            lines.append("_No HD dataset+idsurvey rows._")
        lines.append("")

        lines.append("## No-Offset Control")
        lines.append("")
        noff = summary[summary["offset_mode"] == "none"].sort_values(["frame", "posterior_objective"])
        lines.append(markdown_table(
            noff,
            [
                "frame",
                "alpha",
                "N",
                "data_chi2",
                "posterior_objective",
                "chi2_dof_data",
                "profiled_rms_mag",
            ],
        ))
        lines.append("")

        lines.append("## Best Alpha By Frame And Offset Mode")
        lines.append("")
        best = summary.sort_values("posterior_objective").groupby(
            ["frame", "offset_mode", "prior_config"], as_index=False
        ).head(1)
        best = best.sort_values(["frame", "offset_mode", "prior_config"])
        lines.append(markdown_table(
            best,
            [
                "frame",
                "offset_mode",
                "prior_config",
                "alpha",
                "posterior_objective",
                "data_chi2",
                "prior_chi2",
                "profiled_rms_mag",
            ],
        ))
        lines.append("")

    if not offsets.empty:
        lines.append("## Largest HD Offset Pulls")
        lines.append("")
        pulls = offsets[
            (offsets["frame"] == "HD")
            & np.isfinite(offsets["pull_sigma"])
        ].copy()
        if not pulls.empty:
            pulls["abs_pull"] = pulls["pull_sigma"].abs()
            pulls = pulls.sort_values("abs_pull", ascending=False)
            lines.append(markdown_table(
                pulls,
                [
                    "offset_mode",
                    "prior_config",
                    "alpha",
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
    lines.append("This is a model-consistency diagnostic. It does not promote FR/PLAMB as a cosmological detection by itself. The relevant question is whether the clarified no-loss path prefers alpha=0, 0.5, or 1 under the same calibration and redshift-frame controls used in the existing raw-MU audit.")
    lines.append("")
    lines.append("## Output Files")
    lines.append("")
    lines.append("- `fr_c1pz_noloss_dimming_summary.csv`")
    lines.append("- `fr_c1pz_noloss_dimming_offsets.csv`")
    lines.append("- `fr_c1pz_noloss_dimming_blocks.csv`")
    lines.append("- `fr_c1pz_noloss_dimming_config.json`")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--datasets", default="pantheon,des_raw,union31")
    parser.add_argument("--frames", default="HD,HEL,CMB_PANTHEON_ONLY")
    parser.add_argument("--alphas", default="0,0.5,1")
    parser.add_argument("--offset-modes", default="none,dataset+idsurvey")
    parser.add_argument(
        "--prior-configs",
        default="budget_025mmag_ds050:0.05:0.025",
        help="Comma-separated name:dataset_sigma:survey_sigma configs.",
    )
    parser.add_argument("--h0", type=float, default=67.5)
    parser.add_argument("--z-min", type=float, default=0.01)
    parser.add_argument("--max-z", type=float, default=None)
    parser.add_argument("--min-n", type=int, default=10)
    parser.add_argument("--min-survey-n", type=int, default=20)
    parser.add_argument("--keep-calibrators", action="store_true")
    parser.add_argument("--allow-precision-submatrix", action="store_true")
    return parser


def main() -> None:
    parser = build_parser()
    cli = parser.parse_args()
    cli.outdir.mkdir(parents=True, exist_ok=True)

    summary, offsets, blocks = run_sequence(cli)

    summary_path = cli.outdir / "fr_c1pz_noloss_dimming_summary.csv"
    offsets_path = cli.outdir / "fr_c1pz_noloss_dimming_offsets.csv"
    blocks_path = cli.outdir / "fr_c1pz_noloss_dimming_blocks.csv"
    report_path = cli.outdir / "fr_c1pz_noloss_dimming_report.md"
    config_path = cli.outdir / "fr_c1pz_noloss_dimming_config.json"

    summary.to_csv(summary_path, index=False)
    offsets.to_csv(offsets_path, index=False)
    blocks.to_csv(blocks_path, index=False)
    write_report(report_path, summary, offsets, blocks, cli)
    config_path.write_text(
        json.dumps(
            {
                "script": str(SCRIPT_PATH),
                "repo_root": str(REPO_ROOT),
                "outdir": str(cli.outdir),
                "datasets": cli.datasets,
                "frames": cli.frames,
                "alphas": cli.alphas,
                "offset_modes": cli.offset_modes,
                "prior_configs": cli.prior_configs,
                "h0": cli.h0,
                "z_min": cli.z_min,
                "max_z": cli.max_z,
                "min_n": cli.min_n,
                "min_survey_n": cli.min_survey_n,
                "keep_calibrators": cli.keep_calibrators,
                "allow_precision_submatrix": cli.allow_precision_submatrix,
                "model": "DL=(c/H0)*z*(1+z/2)*(1+z)^alpha",
                "claim_boundary": "diagnostic only; not a FR/PLAMB promotion by itself",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Saved summary: {summary_path}")
    print(f"Saved offsets: {offsets_path}")
    print(f"Saved report: {report_path}")
    if not summary.empty:
        best = summary.sort_values("posterior_objective").groupby(
            ["frame", "offset_mode", "prior_config"], as_index=False
        ).head(1)
        for row in best.sort_values(["frame", "offset_mode", "prior_config"]).itertuples(index=False):
            print(
                f"{row.frame} {row.offset_mode} {row.prior_config}: "
                f"alpha={float(row.alpha):.3f} score={float(row.posterior_objective):.3f} "
                f"data_chi2={float(row.data_chi2):.3f}"
            )


if __name__ == "__main__":
    main()
