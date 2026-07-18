#!/usr/bin/env python3
r"""Post hoc robustness audit of the locked SPARC FR environment study.

These diagnostics cannot alter the preregistered development gate. They inspect
component dependence, leave-one-out stability and host-covariate correlations
to determine whether the suggestive primary association merits a redesigned,
independently preregistered conventional-environment study.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[3]
LOCKED_PROGRAMME = ROOT / "github_export" / "code" / "sparc" / "run_sparc_fr_environment_asymmetry_2026-07-18.py"
RUN_DIR = ROOT / "plamb_runs" / "diagnostics" / "sparc_fr_environment_asymmetry_20260718"
EXPORT_DIR = ROOT / "github_export" / "results" / "2026-07-18" / "fr_environment_asymmetry"
PREREG = RUN_DIR / "sparc_fr_environment_preregistration.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_locked() -> Any:
    spec = importlib.util.spec_from_file_location("locked_fr_environment", LOCKED_PROGRAMME)
    if spec is None or spec.loader is None:
        raise ImportError(LOCKED_PROGRAMME)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def holm_adjust(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    adjusted = np.empty_like(p_values, dtype=float)
    running = 0.0
    count = len(p_values)
    for rank, index in enumerate(order):
        value = min(1.0, (count - rank) * float(p_values[index]))
        running = max(running, value)
        adjusted[index] = running
    return adjusted


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in frame.iterrows():
        values: list[str] = []
        for column in columns:
            value = row[column]
            if isinstance(value, (float, np.floating)):
                values.append(f"{float(value):.6g}" if math.isfinite(float(value)) else "")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def main() -> None:
    locked = load_locked()
    config = json.loads(PREREG.read_text(encoding="utf-8"))
    if sha256(LOCKED_PROGRAMME) != config["programme_sha256"]:
        raise RuntimeError("Locked programme no longer matches the preregistered hash")

    environment = pd.read_csv(RUN_DIR / "sparc_fr_environment_covariates_nonreserved.csv")
    outcomes = pd.read_csv(RUN_DIR / "sparc_fr_environment_development_outcomes.csv")
    catalogue = locked.host_catalogue()
    primary, _ = locked.build_test_frame(
        environment,
        catalogue,
        outcomes,
        locked.PRIMARY_DEFINITION,
        "combined_conventional",
    )

    component_rows: list[dict[str, Any]] = []
    fold_frames: list[pd.DataFrame] = []
    for component in locked.ENVIRONMENT_COMPONENTS:
        result, folds, _ = locked.cross_fitted_test(
            primary,
            component,
            f"posthoc_component__{component}",
        )
        component_rows.append(
            {
                "component": component,
                "partial_spearman_rho": result["partial_spearman_rho"],
                "permutation_p_two_sided_unadjusted": result["permutation_p_two_sided"],
                "cv_rmse_fractional_improvement": result["cv_rmse_fractional_improvement"],
                "coefficient_sign_fraction": result["environment_coefficient_sign_fraction"],
            }
        )
        fold_frames.append(folds)
    components = pd.DataFrame(component_rows)
    components["holm_p_four_components"] = holm_adjust(
        components["permutation_p_two_sided_unadjusted"].to_numpy(dtype=float)
    )

    primary_test = pd.read_csv(RUN_DIR / "sparc_fr_environment_primary_partial_residuals.csv")
    x = primary_test["environment_crossfit_residual"].to_numpy(dtype=float)
    y = primary_test["outcome_crossfit_residual"].to_numpy(dtype=float)
    original_rho = float(spearmanr(x, y).statistic)
    jackknife_rows: list[dict[str, Any]] = []
    for index, galaxy in enumerate(primary_test["Galaxy"].astype(str)):
        keep = np.arange(len(primary_test)) != index
        rho = float(spearmanr(x[keep], y[keep]).statistic)
        jackknife_rows.append(
            {
                "removed_galaxy": galaxy,
                "leave_one_out_rho": rho,
                "change_from_full_rho": rho - original_rho,
                "same_sign_as_full": bool(np.sign(rho) == np.sign(original_rho)),
            }
        )
    jackknife = pd.DataFrame(jackknife_rows).sort_values(
        "change_from_full_rho", key=lambda values: np.abs(values), ascending=False
    )

    quartile_frame = primary_test.copy()
    quartile_frame["environment_residual_quartile"] = pd.qcut(
        quartile_frame["environment_crossfit_residual"],
        4,
        labels=["Q1_low", "Q2", "Q3", "Q4_high"],
    )
    quartiles = (
        quartile_frame.groupby("environment_residual_quartile", observed=True)
        .agg(
            N=("Galaxy", "size"),
            median_environment_residual=("environment_crossfit_residual", "median"),
            mean_outcome_residual=("outcome_crossfit_residual", "mean"),
            median_outcome_residual=("outcome_crossfit_residual", "median"),
        )
        .reset_index()
    )

    host_rows: list[dict[str, Any]] = []
    for feature in locked.HOST_FEATURES:
        values = pd.to_numeric(primary[feature], errors="coerce")
        rho = float(spearmanr(primary["environment_score"], values, nan_policy="omit").statistic)
        host_rows.append({"host_feature": feature, "spearman_rho_with_environment_score": rho})
    host_correlations = pd.DataFrame(host_rows).sort_values(
        "spearman_rho_with_environment_score", key=lambda values: np.abs(values), ascending=False
    )

    components.to_csv(RUN_DIR / "sparc_fr_environment_posthoc_component_tests.csv", index=False)
    pd.concat(fold_frames, ignore_index=True).to_csv(
        RUN_DIR / "sparc_fr_environment_posthoc_component_fold_coefficients.csv", index=False
    )
    jackknife.to_csv(RUN_DIR / "sparc_fr_environment_posthoc_jackknife.csv", index=False)
    quartiles.to_csv(RUN_DIR / "sparc_fr_environment_posthoc_quartiles.csv", index=False)
    host_correlations.to_csv(RUN_DIR / "sparc_fr_environment_posthoc_host_correlations.csv", index=False)

    min_rho = float(jackknife["leave_one_out_rho"].min())
    max_rho = float(jackknife["leave_one_out_rho"].max())
    sign_fraction = float(jackknife["same_sign_as_full"].mean())
    component_columns = [
        "component",
        "partial_spearman_rho",
        "permutation_p_two_sided_unadjusted",
        "holm_p_four_components",
        "cv_rmse_fractional_improvement",
        "coefficient_sign_fraction",
    ]
    report = [
        "# SPARC FR Environment Result: Post Hoc Robustness Audit",
        "",
        "Date: 2026-07-18",
        "",
        "## Status",
        "",
        "This audit is post hoc. It cannot alter the failed preregistered development gate and cannot authorise access to the reserved replication outcomes.",
        "",
        "## Component Readout",
        "",
        markdown_table(components, component_columns),
        "",
        "The four component p-values are Holm-adjusted as one exploratory family. A component result is descriptive even if its adjusted p-value is small because the decomposition was examined after the composite result was known.",
        "",
        "## Leave-One-Out Stability",
        "",
        f"The full partial Spearman correlation is `{original_rho:.6g}`. Across all `{len(jackknife)}` single-galaxy removals, rho ranges from `{min_rho:.6g}` to `{max_rho:.6g}`, with the original sign retained in `{100.0 * sign_fraction:.1f}` per cent of removals.",
        "",
        "Largest absolute changes:",
        "",
        markdown_table(
            jackknife.head(10),
            ["removed_galaxy", "leave_one_out_rho", "change_from_full_rho", "same_sign_as_full"],
        ),
        "",
        "## Residual Quartiles",
        "",
        markdown_table(
            quartiles,
            [
                "environment_residual_quartile",
                "N",
                "median_environment_residual",
                "mean_outcome_residual",
                "median_outcome_residual",
            ],
        ),
        "",
        "## Host-Covariate Correlations",
        "",
        markdown_table(host_correlations, ["host_feature", "spearman_rho_with_environment_score"]),
        "",
        "## Interpretation",
        "",
        "The primary pattern should be carried forward only as a possible conventional environment term. The next defensible analysis would use an external group catalogue or a deeper volume-limited tracer, freeze a single dominant component suggested here, and test genuinely held-out galaxies. No antimatter interpretation follows from these diagnostics.",
        "",
    ]
    report_path = RUN_DIR / "sparc_fr_environment_posthoc_audit.md"
    report_path.write_text("\n".join(report), encoding="utf-8")

    EXPORT_DIR.mkdir(parents=True, exist_ok=True)
    for name in [
        "sparc_fr_environment_posthoc_component_tests.csv",
        "sparc_fr_environment_posthoc_component_fold_coefficients.csv",
        "sparc_fr_environment_posthoc_jackknife.csv",
        "sparc_fr_environment_posthoc_quartiles.csv",
        "sparc_fr_environment_posthoc_host_correlations.csv",
        "sparc_fr_environment_posthoc_audit.md",
    ]:
        source = RUN_DIR / name
        (EXPORT_DIR / name).write_bytes(source.read_bytes())

    manifest_rows: list[dict[str, Any]] = []
    for path in sorted(EXPORT_DIR.iterdir()):
        if path.name == "manifest.csv" or not path.is_file():
            continue
        manifest_rows.append(
            {
                "path": path.name,
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
                "role": "committed analysis product",
                "tracked_in_git": True,
            }
        )
    for path, role, tracked in [
        (LOCKED_PROGRAMME, "locked programme", True),
        (Path(__file__).resolve(), "post hoc audit programme", True),
        (locked.TWOMRS, "external 2MRS catalogue; URL in preregistration", False),
        (
            locked.POSITION_CACHE,
            "external NED resolution cache excluding reserved galaxies",
            False,
        ),
    ]:
        manifest_rows.append(
            {
                "path": (
                    str(path.relative_to(ROOT)).replace("\\", "/")
                    if tracked
                    else str(path)
                ),
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
                "role": role,
                "tracked_in_git": tracked,
            }
        )
    pd.DataFrame(manifest_rows).to_csv(EXPORT_DIR / "manifest.csv", index=False)
    print(f"Saved post hoc audit: {report_path}")
    print(f"Leave-one-out rho range: {min_rho:.6g} to {max_rho:.6g}")


if __name__ == "__main__":
    main()
