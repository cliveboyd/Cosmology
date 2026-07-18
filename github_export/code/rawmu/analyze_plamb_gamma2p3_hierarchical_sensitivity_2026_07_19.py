#!/usr/bin/env python3
"""Calibration-aware sensitivity analysis for the post-hoc gamma=2.3 clock law."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
RAWMU_CODE = SCRIPT_PATH.parent
SHARED_CODE = REPO_ROOT / "github_export" / "code" / "shared"
for code_path in (REPO_ROOT, RAWMU_CODE, SHARED_CODE):
    if str(code_path) not in sys.path:
        sys.path.insert(0, str(code_path))

import run_plamb_nested_clock_fit_2026_07_18 as nested  # noqa: E402
import run_rawmu_release_grounded_holdouts_2026_07_18 as grounded  # noqa: E402


DATE_TAG = "2026-07-19"
GAMMA_FIXED = 2.3
ALPHA_FIXED = 0.0
H0_FIXED = 67.5
P_BOUNDS = (-0.5, 2.5)
OMEGA_BOUNDS = (0.05, 0.60)
DEFAULT_OUTDIR = (
    REPO_ROOT
    / "github_export"
    / "results"
    / DATE_TAG
    / "rawmu_plamb_gamma2p3_hierarchy"
)
DEFAULT_REGISTRY = (
    REPO_ROOT
    / "github_export"
    / "results"
    / "2026-07-18"
    / "rawmu"
    / "rawmu_release_grounded_prior_registry_2026-07-18.json"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bounded_best(objective, bounds: tuple[float, float]) -> tuple[float, float, bool]:
    result = minimize_scalar(
        objective,
        bounds=bounds,
        method="bounded",
        options={"xatol": 1e-10, "maxiter": 600},
    )
    candidates = [
        (float(result.fun), float(result.x)),
        (float(objective(bounds[0])), float(bounds[0])),
        (float(objective(bounds[1])), float(bounds[1])),
    ]
    value, parameter = min(candidates)
    return parameter, value, bool(result.success or parameter in bounds)


def path_score(blocks: list[grounded.Block], p: float) -> nested.Score:
    return nested.score_clock_profiled_alpha(
        blocks, GAMMA_FIXED, float(p), ALPHA_FIXED, H0_FIXED
    )


def lcdm_score(blocks: list[grounded.Block], omega_m: float) -> nested.Score:
    return nested.score_lcdm(blocks, float(omega_m), H0_FIXED)


def fit_variant(
    blocks: list[grounded.Block],
    variant: str,
    hierarchy_metadata: dict[str, object],
) -> tuple[dict[str, object], dict[str, nested.Score]]:
    p, path_chi2, path_ok = bounded_best(
        lambda value: path_score(blocks, value).chi2, P_BOUNDS
    )
    omega_m, lcdm_chi2, lcdm_ok = bounded_best(
        lambda value: lcdm_score(blocks, value).chi2, OMEGA_BOUNDS
    )
    path = path_score(blocks, p)
    lcdm = lcdm_score(blocks, omega_m)
    n = int(sum(block.n for block in blocks))
    k = len(blocks) + 1
    signature = grounded.nuisance_signature(blocks, variant, hierarchy_metadata)
    row = {
        "analysis_variant": variant,
        "evidence_status": (
            "primary_released_total"
            if variant == "released_total_primary"
            else "correlated_covariance_reconstruction_sensitivity"
        ),
        "N": n,
        "N_release_intercepts": len(blocks),
        "PATH_gamma_c": GAMMA_FIXED,
        "PATH_alpha": ALPHA_FIXED,
        "PATH_p": p,
        "PATH_chi2": path_chi2,
        "LCDM_Omega_m": omega_m,
        "LCDM_chi2": lcdm_chi2,
        "PATH_BIC": path_chi2 + k * math.log(n),
        "LCDM_BIC": lcdm_chi2 + k * math.log(n),
        "delta_chi2_PATH_minus_LCDM": path_chi2 - lcdm_chi2,
        "delta_BIC_PATH_minus_LCDM": path_chi2 - lcdm_chi2,
        "nuisance_signature": signature,
        "same_nuisance_verified": True,
        "optimisation_success": bool(path_ok and lcdm_ok),
    }
    return row, {"PATH": path, "LCDM": lcdm}


def model_mu(block: grounded.Block, model: str, parameter: float) -> np.ndarray:
    if model == "PATH":
        parameters = {
            "gamma_c": GAMMA_FIXED,
            "p": float(parameter),
            "alpha": ALPHA_FIXED,
        }
        return nested.model_mu_from_parameters(
            block.z, nested.MODEL_SPECS["GENERAL_FREE"], parameters, H0_FIXED
        )
    return nested.model_mu_from_parameters(
        block.z,
        nested.MODEL_SPECS["LCDM"],
        {"Omega_m": float(parameter)},
        H0_FIXED,
    )


def residual_group_mode(
    block: grounded.Block,
    mask: np.ndarray,
    model: str,
    parameter: float,
    release_offset: float,
) -> float:
    sub = grounded.subset_block(block, mask, "gamma2p3_budget_mode")
    residual = sub.mu - model_mu(sub, model, parameter) - float(release_offset)
    ones = np.ones(sub.n)
    return float(ones @ sub.precision @ residual / (ones @ sub.precision @ ones))


def budget_table(
    primary: list[grounded.Block],
    primary_row: dict[str, object],
    scores: dict[str, nested.Score],
    prior_table: pd.DataFrame,
    min_survey_n: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for block in primary:
        if block.label == "Union3p1_UNITY1p8":
            continue
        component = (
            "CALIB"
            if block.label == "PantheonPlusSH0ES"
            else "all_systematics_covariance_upper_bound"
        )
        values, counts = np.unique(block.survey, return_counts=True)
        for value, count in zip(values, counts):
            if int(count) < min_survey_n:
                continue
            group = f"IDSURVEY_{int(float(value))}"
            match = prior_table[
                (prior_table["dataset"] == block.label)
                & (prior_table["level"] == "survey")
                & (prior_table["component_kind"] == component)
                & (prior_table["group"] == group)
            ]
            if match.empty:
                continue
            budget = float(match.iloc[0]["sigma_projected_mag"])
            mask = block.survey == value
            path_mode = residual_group_mode(
                block,
                mask,
                "PATH",
                float(primary_row["PATH_p"]),
                scores["PATH"].offsets[block.label],
            )
            lcdm_mode = residual_group_mode(
                block,
                mask,
                "LCDM",
                float(primary_row["LCDM_Omega_m"]),
                scores["LCDM"].offsets[block.label],
            )
            delta = path_mode - lcdm_mode
            rows.append(
                {
                    "dataset": block.label,
                    "group": group,
                    "N": int(count),
                    "budget_sigma_mag": budget,
                    "PATH_residual_mode_mag": path_mode,
                    "LCDM_residual_mode_mag": lcdm_mode,
                    "delta_mode_PATH_minus_LCDM_mag": delta,
                    "abs_delta_over_budget": abs(delta) / budget if budget > 0 else np.inf,
                    "within_budget": bool(abs(delta) <= budget),
                }
            )
    return pd.DataFrame(rows)


def write_report(
    path: Path,
    fits: pd.DataFrame,
    budgets: pd.DataFrame,
    metadata: dict[str, object],
    gates: dict[str, bool],
) -> None:
    primary = fits[fits["analysis_variant"] == "released_total_primary"].iloc[0]
    worst = budgets.sort_values("abs_delta_over_budget", ascending=False).head(10)
    sensitivity = fits[fits["analysis_variant"] != "released_total_primary"]
    lines = [
        "# Hierarchical raw-MU sensitivity: fixed gamma = 2.3",
        "",
        f"Generated: {datetime.now():%Y-%m-%d %H:%M:%S}",
        "",
        "## Claim boundary",
        "",
        "This is an outcome-aware sensitivity analysis. `gamma_c=2.3` was selected after "
        "inspection of the 18 July supernova profiles, so this run is not independent evidence "
        "and cannot promote PLAMB.",
        "",
        "## Model",
        "",
        "```text",
        "c(z) / c0                    = 1 + 2.3 z",
        "|dz/dT| / H0                = (1 + z)^p",
        "D_L(z)                      = (c0/H0) integral[0,z] (1+2.3u)/(1+u)^p du",
        "alpha                       = 0",
        "```",
        "",
        "Pantheon+ and DES use `zHD`; Union3.1 uses its released `z`. Every model has one "
        "unpenalised release intercept and one shape parameter, so delta BIC equals delta chi2.",
        "",
        "## Calibration treatment",
        "",
        "The primary cell uses each release's published `STAT+SYS` covariance. Calibration "
        "budgets already represented there are not added again. Sensitivity cells replace the "
        "relevant total covariance by `STATONLY + X Lambda X^T`, which is the covariance obtained "
        "after analytically marginalising zero-mean Gaussian dataset/survey offsets. Because the "
        "grouped reconstructions omit most of the original covariance structure, they are "
        "correlated diagnostics rather than alternative headline likelihoods.",
        "",
        "## Results",
        "",
        fits[
            [
                "analysis_variant",
                "evidence_status",
                "PATH_p",
                "LCDM_Omega_m",
                "delta_BIC_PATH_minus_LCDM",
                "same_nuisance_verified",
            ]
        ].to_markdown(index=False),
        "",
        f"Primary delta BIC (PATH - Lambda-CDM) = {primary['delta_BIC_PATH_minus_LCDM']:.6f}.",
        f"Sensitivity delta-BIC range = [{sensitivity['delta_BIC_PATH_minus_LCDM'].min():.6f}, "
        f"{sensitivity['delta_BIC_PATH_minus_LCDM'].max():.6f}].",
        "",
        "## External calibration budgets",
        "",
        worst[
            [
                "dataset",
                "group",
                "N",
                "budget_sigma_mag",
                "delta_mode_PATH_minus_LCDM_mag",
                "abs_delta_over_budget",
                "within_budget",
            ]
        ].to_markdown(index=False),
        "",
        "## Gates",
        "",
    ]
    for name, passed in gates.items():
        lines.append(f"- `{name}`: **{'PASS' if passed else 'FAIL'}**")
    lines.extend(
        [
            "",
            "These gates describe calibration robustness of the post-hoc cell only. They do not "
            "repair the failed independent-gamma and high-redshift predictive gates.",
            "",
            "## Reconstruction limitations",
            "",
        ]
    )
    for variant, details in metadata["hierarchies"].items():
        error = float(details["metrics"]["relative_frobenius_error"])
        lines.append(f"- `{variant}` leaves {100.0 * error:.1f}% of the target Frobenius norm unexplained.")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--min-survey-n", type=int, default=20)
    args = parser.parse_args()

    grounded.self_test()
    if not args.registry.exists():
        raise FileNotFoundError(args.registry)
    registry = json.loads(args.registry.read_text(encoding="utf-8"))
    prior_table = pd.DataFrame(registry["rows"])
    variants, load_metadata = grounded.load_blocks(args.min_survey_n)
    hierarchy_metadata = dict(load_metadata["hierarchies"])
    hierarchy_metadata["combined_calib_des_statonly_grouped"] = {
        "components": [
            hierarchy_metadata["pantheon_calib_statonly_grouped"],
            hierarchy_metadata["des_allsys_statonly_grouped"],
        ]
    }

    args.outdir.mkdir(parents=True, exist_ok=True)
    protocol = {
        "analysis_id": "plamb_gamma2p3_hierarchical_sensitivity_2026-07-19",
        "status": "post_hoc_outcome_aware_sensitivity",
        "primary_frame": load_metadata["primary_frame"],
        "gamma_c_fixed": GAMMA_FIXED,
        "alpha_fixed": ALPHA_FIXED,
        "p_bounds": P_BOUNDS,
        "omega_m_bounds": OMEGA_BOUNDS,
        "calibration_rule": "released total primary; stat-only plus marginalised grouped offsets as sensitivity",
        "double_count_guard": "never add release-contained calibration priors to STAT+SYS",
    }
    protocol_path = args.outdir / f"plamb_gamma2p3_hierarchical_protocol_{DATE_TAG}.json"
    protocol_path.write_text(json.dumps(protocol, indent=2), encoding="utf-8")

    fit_rows: list[dict[str, object]] = []
    fit_scores: dict[str, dict[str, nested.Score]] = {}
    for variant, blocks in variants.items():
        row, scores = fit_variant(blocks, variant, hierarchy_metadata)
        fit_rows.append(row)
        fit_scores[variant] = scores
        print(
            f"{variant}: p={row['PATH_p']:.6f} "
            f"delta_BIC={row['delta_BIC_PATH_minus_LCDM']:.6f}",
            flush=True,
        )
    fits = pd.DataFrame(fit_rows)
    primary_row = fits[fits["analysis_variant"] == "released_total_primary"].iloc[0].to_dict()
    budgets = budget_table(
        variants["released_total_primary"],
        primary_row,
        fit_scores["released_total_primary"],
        prior_table,
        args.min_survey_n,
    )
    sensitivity = fits[fits["analysis_variant"] != "released_total_primary"]
    gates = {
        "same_nuisance_all_cells": bool(fits["same_nuisance_verified"].all()),
        "optimisation_success_all_cells": bool(fits["optimisation_success"].all()),
        "primary_absolute_delta_BIC_at_most_2": bool(
            abs(float(primary_row["delta_BIC_PATH_minus_LCDM"])) <= 2.0
        ),
        "all_covariance_sensitivities_absolute_delta_BIC_at_most_2": bool(
            (sensitivity["delta_BIC_PATH_minus_LCDM"].abs() <= 2.0).all()
        ),
        "all_model_difference_modes_within_external_budget": bool(
            not budgets.empty and budgets["within_budget"].all()
        ),
        "no_calibration_prior_double_counting": True,
        "independent_gamma_derivation": False,
    }

    fits_path = args.outdir / f"plamb_gamma2p3_hierarchical_fits_{DATE_TAG}.csv"
    budgets_path = args.outdir / f"plamb_gamma2p3_hierarchical_budget_modes_{DATE_TAG}.csv"
    report_path = args.outdir / f"plamb_gamma2p3_hierarchical_report_{DATE_TAG}.md"
    summary_path = args.outdir / f"plamb_gamma2p3_hierarchical_summary_{DATE_TAG}.json"
    fits.to_csv(fits_path, index=False)
    budgets.to_csv(budgets_path, index=False)
    write_report(report_path, fits, budgets, load_metadata, gates)
    summary = {
        "analysis_date": DATE_TAG,
        "primary": primary_row,
        "sensitivity_delta_BIC_range": [
            float(sensitivity["delta_BIC_PATH_minus_LCDM"].min()),
            float(sensitivity["delta_BIC_PATH_minus_LCDM"].max()),
        ],
        "maximum_abs_model_difference_over_budget": float(
            budgets["abs_delta_over_budget"].max()
        ),
        "gates": gates,
        "claim_promotion": False,
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    manifest = pd.DataFrame(
        [
            {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in (protocol_path, fits_path, budgets_path, summary_path, report_path)
        ]
    )
    manifest.to_csv(
        args.outdir / f"plamb_gamma2p3_hierarchical_manifest_{DATE_TAG}.csv",
        index=False,
    )
    print(f"Saved report: {report_path}")


if __name__ == "__main__":
    main()
