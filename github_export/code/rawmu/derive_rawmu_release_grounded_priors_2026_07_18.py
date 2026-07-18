r"""Derive release-grounded raw-MU calibration-mode budgets.

The output is a prior registry, not an instruction to add every projected
width to the likelihood. Release-internal calibration already represented in
STAT+SYS is marked ``covariance_contained`` to prevent double-counting.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from diagnose_pantheon_rawmu_fr import read_table  # noqa: E402
from fit_rawmu_fr_cov_likelihood import load_covariance_or_precision  # noqa: E402


DATE_TAG = "2026-07-18"
DEFAULT_OUTDIR = REPO_ROOT / "github_export" / "results" / DATE_TAG / "rawmu"

PANTHEON_DATA = REPO_ROOT / "Pantheon+SH0ES.dat"
PANTHEON_TOTAL = REPO_ROOT / "Pantheon+SH0ES_STAT+SYS.cov"
PANTHEON_COMPONENT_DIR = (
    REPO_ROOT / "external_datasets" / "current_cosmology_sources" / "PantheonPlusSH0ES"
)
PANTHEON_STAT = PANTHEON_COMPONENT_DIR / "Pantheon+SH0ES_STATONLY.cov"
PANTHEON_CALIB = PANTHEON_COMPONENT_DIR / "Pantheon+SH0ES_122221_CALIB.cov"
PANTHEON_CSP_RECAL = PANTHEON_COMPONENT_DIR / "Pantheon+SH0ES_122221_CSP_RECAL.cov"
PANTHEON_HSTCALSPEC = PANTHEON_COMPONENT_DIR / "Pantheon+SH0ES_122221_HSTCALSPEC.cov"
PANTHEON_FRAGILISTIC = PANTHEON_COMPONENT_DIR / "FRAGILISTIC_COVARIANCE.npz"

DES_DIR = REPO_ROOT / "external_datasets" / "current_cosmology_sources" / "DES_SN5YR_Dovekie"
DES_DATA = DES_DIR / "DES-Dovekie_HD.csv"
DES_TOTAL = DES_DIR / "STAT+SYS.npz"
DES_STAT = DES_DIR / "STATONLY.npz"

UNION_PATH = (
    REPO_ROOT
    / "external_datasets"
    / "current_cosmology_sources"
    / "Union3_release"
    / "mu_mat_union3.1_UNITY1.8_template_cosmo=2_0_mu.fits"
)

SOURCE_URLS = {
    "PantheonPlusSH0ES_release": (
        "https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/"
        "Pantheon%2B_Data/4_DISTANCES_AND_COVAR"
    ),
    "PantheonPlusSH0ES_systematic_groupings": (
        "https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/"
        "Pantheon%2B_Data/4_DISTANCES_AND_COVAR/sytematic_groupings"
    ),
    "PantheonPlus_calibration_paper": "https://arxiv.org/abs/2110.03486",
    "PantheonPlusSH0ES_STATONLY_raw": (
        "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
        "Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES_STATONLY.cov"
    ),
    "PantheonPlusSH0ES_CALIB_raw": (
        "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
        "Pantheon%2B_Data/4_DISTANCES_AND_COVAR/sytematic_groupings/"
        "Pantheon%2BSH0ES_122221_CALIB.cov"
    ),
    "PantheonPlusSH0ES_CSP_RECAL_raw": (
        "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
        "Pantheon%2B_Data/4_DISTANCES_AND_COVAR/sytematic_groupings/"
        "Pantheon%2BSH0ES_122221_CSP_RECAL.cov"
    ),
    "PantheonPlusSH0ES_HSTCALSPEC_raw": (
        "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
        "Pantheon%2B_Data/4_DISTANCES_AND_COVAR/sytematic_groupings/"
        "Pantheon%2BSH0ES_122221_HSTCALSPEC.cov"
    ),
    "PantheonPlusSH0ES_FRAGILISTIC_raw": (
        "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
        "Pantheon%2B_Data/2_CALIBRATION/FRAGILISTIC_COVARIANCE.npz"
    ),
    "DES_SN5YR_release": "https://github.com/des-science/DES-SN5YR",
    "Union3_paper": "https://doi.org/10.3847/1538-4357/adc0a5",
}


@dataclass(frozen=True)
class ProjectionSpec:
    dataset: str
    component_kind: str
    source_path: Path
    covariance_role: str
    additive_prior: bool
    budget_kind: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_record(path: Path, role: str, url_key: str | None = None) -> dict[str, object]:
    if not path.exists():
        raise FileNotFoundError(path)
    return {
        "role": role,
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "source_url": SOURCE_URLS.get(url_key or "", ""),
    }


def load_pantheon_covariance(path: Path, n: int) -> np.ndarray:
    flat = np.loadtxt(path, dtype=float).reshape(-1)
    if flat.size != n * n + 1 or int(round(flat[0])) != n:
        raise ValueError(f"Unexpected Pantheon+ covariance shape in {path}: {flat.size}")
    cov = flat[1:].reshape(n, n)
    return 0.5 * (cov + cov.T)


def invert_precision(path: Path, n: int) -> tuple[np.ndarray, np.ndarray, str]:
    cov, precision, note = load_covariance_or_precision(path, n)
    if precision is None:
        if cov is None:
            raise ValueError(f"No covariance or precision in {path}")
        precision = np.linalg.pinv(cov)
    try:
        covariance = np.linalg.inv(precision)
    except np.linalg.LinAlgError:
        covariance = np.linalg.pinv(precision)
    covariance = 0.5 * (covariance + covariance.T)
    return covariance, precision, note


def projected_sigma(covariance: np.ndarray, indices: np.ndarray) -> tuple[float, float]:
    if indices.size == 0:
        return float("nan"), float("nan")
    sub = covariance[np.ix_(indices, indices)]
    variance = float(np.sum(sub) / (indices.size * indices.size))
    clipped = max(variance, 0.0)
    return math.sqrt(clipped), variance


def grouped_hierarchy(
    component: np.ndarray,
    selected_idx: np.ndarray,
    raw_groups: np.ndarray,
    min_survey_n: int,
) -> tuple[float, dict[str, dict[str, object]], np.ndarray, dict[str, float]]:
    """Moment-match dataset and survey common modes to a covariance component."""
    target = component[np.ix_(selected_idx, selected_idx)]
    selected_groups = np.asarray(raw_groups[selected_idx], dtype=object)
    counts = pd.Series(selected_groups).value_counts()
    keep = set(counts[counts >= min_survey_n].index)
    pooled = np.asarray([value if value in keep else "OTHER" for value in selected_groups], dtype=object)
    labels = sorted(set(pooled), key=lambda value: (str(value) == "OTHER", str(value)))
    masks = {str(label): np.flatnonzero(pooled == label) for label in labels}

    n = len(selected_idx)
    within_sum = 0.0
    within_pairs = 0
    for idx in masks.values():
        within_sum += float(np.sum(target[np.ix_(idx, idx)]))
        within_pairs += int(idx.size * idx.size)
    cross_pairs = n * n - within_pairs
    cross_sum = float(np.sum(target)) - within_sum
    dataset_variance_raw = cross_sum / cross_pairs if cross_pairs > 0 else float(np.sum(target)) / (n * n)
    dataset_variance = max(dataset_variance_raw, 0.0)

    approximation = np.full_like(target, dataset_variance)
    survey_info: dict[str, dict[str, object]] = {}
    for label, idx in masks.items():
        projected_variance = float(np.sum(target[np.ix_(idx, idx)]) / (idx.size * idx.size))
        residual_variance = max(projected_variance - dataset_variance, 0.0)
        approximation[np.ix_(idx, idx)] += residual_variance
        survey_info[label] = {
            "indices_local": idx,
            "N": int(idx.size),
            "projected_variance_raw_mag2": projected_variance,
            "sigma_projected_mag": math.sqrt(max(projected_variance, 0.0)),
            "hierarchy_variance_mag2": residual_variance,
            "sigma_hierarchy_mag": math.sqrt(residual_variance),
        }

    residual = target - approximation
    target_norm = float(np.linalg.norm(target, ord="fro"))
    residual_norm = float(np.linalg.norm(residual, ord="fro"))
    relative_error = residual_norm / target_norm if target_norm > 0.0 else float("nan")
    metrics = {
        "target_frobenius_norm_mag2": target_norm,
        "residual_frobenius_norm_mag2": residual_norm,
        "relative_frobenius_error": relative_error,
        "captured_frobenius_power_fraction": 1.0 - relative_error**2 if np.isfinite(relative_error) else float("nan"),
        "diagonal_rmse_mag2": float(np.sqrt(np.mean(np.square(np.diag(residual))))),
        "dataset_variance_raw_mag2": dataset_variance_raw,
        "dataset_variance_mag2": dataset_variance,
        "dataset_sigma_hierarchy_mag": math.sqrt(dataset_variance),
        "N_groups": len(masks),
    }
    return math.sqrt(dataset_variance), survey_info, approximation, metrics


def make_row(
    spec: ProjectionSpec,
    level: str,
    group: str,
    n: int,
    sigma_projected: float,
    variance_raw: float,
    sigma_documented: float | None = None,
    sigma_hierarchy: float | None = None,
    explicit_base: str = "released_STAT+SYS",
    note: str = "",
) -> dict[str, object]:
    sigma_included = sigma_projected if np.isfinite(sigma_projected) else 0.0
    if sigma_documented is None or not spec.additive_prior:
        sigma_add = 0.0
    else:
        sigma_add = math.sqrt(max(sigma_documented**2 - sigma_included**2, 0.0))
    return {
        "dataset": spec.dataset,
        "level": level,
        "group": group,
        "N": int(n),
        "component_kind": spec.component_kind,
        "covariance_role": spec.covariance_role,
        "budget_kind": spec.budget_kind,
        "sigma_projected_mag": sigma_projected,
        "projected_variance_raw_mag2": variance_raw,
        "sigma_documented_mag": sigma_documented,
        "sigma_included_mag": sigma_included,
        "sigma_add_mag": sigma_add,
        "sigma_hierarchy_mag": sigma_hierarchy,
        "additive_prior": bool(spec.additive_prior and sigma_add > 0.0),
        "prior_action": "residual_sensitivity" if sigma_add > 0.0 else "covariance_contained",
        "explicit_base": explicit_base,
        "explicit_prior_role": "component_replacement_sensitivity" if sigma_hierarchy is not None else "none",
        "source_path": str(spec.source_path),
        "note": note,
    }


def derive_pantheon(
    min_survey_n: int,
) -> tuple[list[dict[str, object]], dict[str, object], list[dict[str, object]]]:
    df = read_table(PANTHEON_DATA)
    n_full = len(df)
    selected = (
        np.isfinite(pd.to_numeric(df["zHD"], errors="coerce"))
        & np.isfinite(pd.to_numeric(df["MU_SH0ES"], errors="coerce"))
        & (pd.to_numeric(df["zHD"], errors="coerce") > 0.01)
        & (pd.to_numeric(df["IS_CALIBRATOR"], errors="coerce").fillna(0.0) == 0.0)
    ).to_numpy()
    selected_idx = np.flatnonzero(selected)
    calibration_group = load_pantheon_covariance(PANTHEON_CALIB, n_full)
    csp_recal_group = load_pantheon_covariance(PANTHEON_CSP_RECAL, n_full)
    hstcalspec_group = load_pantheon_covariance(PANTHEON_HSTCALSPEC, n_full)
    stat = load_pantheon_covariance(PANTHEON_STAT, n_full)
    total = load_pantheon_covariance(PANTHEON_TOTAL, n_full)
    calibration = calibration_group - stat
    csp_recal = csp_recal_group - stat
    hstcalspec = hstcalspec_group - stat
    rows: list[dict[str, object]] = []
    surveys = pd.to_numeric(df["IDSURVEY"], errors="coerce").to_numpy()
    reconstructions: list[dict[str, object]] = []
    components = [
        ("CALIB", calibration, PANTHEON_CALIB, "primary_calibration_component"),
        ("CSP_RECAL", csp_recal, PANTHEON_CSP_RECAL, "calibration_component_control"),
        ("HSTCALSPEC", hstcalspec, PANTHEON_HSTCALSPEC, "calibration_component_control"),
        (
            "CALIB+CSP_RECAL+HSTCALSPEC",
            calibration + csp_recal + hstcalspec,
            PANTHEON_CALIB,
            "aggregate_calibration_component_control",
        ),
    ]
    for component_name, component, source_path, budget_kind in components:
        spec = ProjectionSpec(
            dataset="PantheonPlusSH0ES",
            component_kind=component_name,
            source_path=source_path,
            covariance_role="contained_in_STAT+SYS",
            additive_prior=False,
            budget_kind=budget_kind,
        )
        dataset_sigma, survey_info, approximation, metrics = grouped_hierarchy(
            component, selected_idx, surveys, min_survey_n
        )
        sigma, variance = projected_sigma(component, selected_idx)
        rows.append(
            make_row(
                spec,
                "dataset",
                "ALL",
                selected_idx.size,
                sigma,
                variance,
                sigma_hierarchy=dataset_sigma,
                explicit_base="STATONLY",
                note="Grouping increment is release grouping covariance minus STATONLY.",
            )
        )
        for label, info in survey_info.items():
            group = f"IDSURVEY_{int(float(label))}" if label != "OTHER" else "IDSURVEY_OTHER"
            rows.append(
                make_row(
                    spec,
                    "survey",
                    group,
                    int(info["N"]),
                    float(info["sigma_projected_mag"]),
                    float(info["projected_variance_raw_mag2"]),
                    sigma_hierarchy=float(info["sigma_hierarchy_mag"]),
                    explicit_base="STATONLY",
                    note="Grouping increment is release grouping covariance minus STATONLY.",
                )
            )
        stat_sub = stat[np.ix_(selected_idx, selected_idx)]
        total_sub = total[np.ix_(selected_idx, selected_idx)]
        grouped_reconstruction = stat_sub + approximation
        component_residual = component[np.ix_(selected_idx, selected_idx)] - approximation
        grouping_target = stat_sub + component[np.ix_(selected_idx, selected_idx)]
        reconstructions.append(
            {
                "dataset": "PantheonPlusSH0ES",
                "component": component_name,
                "reconstruction": "dataset+survey grouped offsets",
                "base": "STATONLY",
                "status": "stat_only_grouping_reconstruction_sensitivity",
                "grouping_increment_definition": "grouping_covariance-STATONLY",
                "grouping_covariance_relative_frobenius_error": float(
                    np.linalg.norm(component_residual, ord="fro")
                    / np.linalg.norm(grouping_target, ord="fro")
                ),
                "released_total_relative_frobenius_error": float(
                    np.linalg.norm(total_sub - grouped_reconstruction, ord="fro")
                    / np.linalg.norm(total_sub, ord="fro")
                ),
                "increment_diagonal_min_mag2": float(
                    np.min(np.diag(component[np.ix_(selected_idx, selected_idx)]))
                ),
                **metrics,
            }
        )

    sys_cov = total - stat
    denom = float(np.linalg.norm(sys_cov, ord="fro"))
    calib_fraction = float(np.linalg.norm(calibration, ord="fro") / denom) if denom > 0 else float("nan")
    audit = {
        "N_source": n_full,
        "N_selected": int(selected_idx.size),
        "total_symmetry_max_abs": float(np.max(np.abs(total - total.T))),
        "stat_symmetry_max_abs": float(np.max(np.abs(stat - stat.T))),
        "calibration_symmetry_max_abs": float(np.max(np.abs(calibration - calibration.T))),
        "calibration_to_all_systematics_frobenius_ratio": calib_fraction,
        "grouping_files_include_STATONLY": True,
        "calibration_increment_diagonal_min_mag2_selected": float(
            np.min(np.diag(calibration[np.ix_(selected_idx, selected_idx)]))
        ),
        "calibration_increment_diagonal_max_mag2_selected": float(
            np.max(np.diag(calibration[np.ix_(selected_idx, selected_idx)]))
        ),
        "CSP_RECAL_to_all_systematics_frobenius_ratio": float(np.linalg.norm(csp_recal, ord="fro") / denom),
        "HSTCALSPEC_to_all_systematics_frobenius_ratio": float(np.linalg.norm(hstcalspec, ord="fro") / denom),
    }
    fragilistic = np.load(PANTHEON_FRAGILISTIC)
    frag_cov = np.asarray(fragilistic["cov"], dtype=float)
    frag_sigma = np.sqrt(np.maximum(np.diag(frag_cov), 0.0))
    audit["FRAGILISTIC"] = {
        "shape": list(frag_cov.shape),
        "N_labels": int(len(fragilistic["labels"])),
        "sigma_min_mag": float(np.min(frag_sigma)),
        "sigma_median_mag": float(np.median(frag_sigma)),
        "sigma_max_mag": float(np.max(frag_sigma)),
        "eigenvalue_min_mag2": float(np.min(np.linalg.eigvalsh(frag_cov))),
        "likelihood_use": (
            "source calibration-parameter covariance only; not mapped directly without the "
            "release SN-distance Jacobian"
        ),
    }
    return rows, audit, reconstructions


def derive_des(
    min_survey_n: int,
) -> tuple[list[dict[str, object]], dict[str, object], list[dict[str, object]]]:
    df = read_table(DES_DATA)
    n = len(df)
    total, _p_total, total_note = invert_precision(DES_TOTAL, n)
    stat, _p_stat, stat_note = invert_precision(DES_STAT, n)
    systematics = 0.5 * ((total - stat) + (total - stat).T)

    z = pd.to_numeric(df["zHD"], errors="coerce").to_numpy(dtype=float)
    mu = pd.to_numeric(df["MU"], errors="coerce").to_numpy(dtype=float)
    selected = np.isfinite(z) & np.isfinite(mu) & (z > 0.01)
    selected_idx = np.flatnonzero(selected)
    spec = ProjectionSpec(
        dataset="DES_Dovekie_raw",
        component_kind="all_systematics_covariance_upper_bound",
        source_path=DES_TOTAL,
        covariance_role="contained_in_STAT+SYS",
        additive_prior=False,
        budget_kind="release_all_systematics_upper_bound",
    )
    rows: list[dict[str, object]] = []
    surveys = pd.to_numeric(df["IDSURVEY"], errors="coerce").to_numpy()
    dataset_sigma, survey_info, approximation, metrics = grouped_hierarchy(
        systematics, selected_idx, surveys, min_survey_n
    )
    sigma, variance = projected_sigma(systematics, selected_idx)
    rows.append(
        make_row(
            spec,
            "dataset",
            "ALL",
            selected_idx.size,
            sigma,
            variance,
            sigma_hierarchy=dataset_sigma,
            explicit_base="STATONLY",
            note="Upper bound only: C_total-C_stat is not calibration-only.",
        )
    )
    for label, info in survey_info.items():
        group = f"IDSURVEY_{int(float(label))}" if label != "OTHER" else "IDSURVEY_OTHER"
        rows.append(
            make_row(
                spec,
                "survey",
                group,
                int(info["N"]),
                float(info["sigma_projected_mag"]),
                float(info["projected_variance_raw_mag2"]),
                sigma_hierarchy=float(info["sigma_hierarchy_mag"]),
                explicit_base="STATONLY",
                note="Upper bound only: includes every released systematic component.",
            )
        )
    audit = {
        "N_source": n,
        "N_selected": int(selected_idx.size),
        "total_loader_note": total_note,
        "stat_loader_note": stat_note,
        "systematics_diagonal_min_mag2": float(np.min(np.diag(systematics))),
        "systematics_diagonal_max_mag2": float(np.max(np.diag(systematics))),
        "negative_systematics_diagonal_count": int(np.count_nonzero(np.diag(systematics) < 0.0)),
    }
    reconstruction = {
        "dataset": "DES_Dovekie_raw",
        "component": "ALL_SYSTEMATICS",
        "reconstruction": "dataset+survey grouped offsets",
        "base": "STATONLY",
        "status": "all_systematics_common_mode_sensitivity",
        "grouping_increment_definition": "inv(P_STAT+SYS)-inv(P_STATONLY)",
        "grouping_covariance_relative_frobenius_error": float(
            np.linalg.norm(systematics[np.ix_(selected_idx, selected_idx)] - approximation, ord="fro")
            / np.linalg.norm(total[np.ix_(selected_idx, selected_idx)], ord="fro")
        ),
        "released_total_relative_frobenius_error": float(
            np.linalg.norm(
                total[np.ix_(selected_idx, selected_idx)]
                - (stat[np.ix_(selected_idx, selected_idx)] + approximation),
                ord="fro",
            )
            / np.linalg.norm(total[np.ix_(selected_idx, selected_idx)], ord="fro")
        ),
        "increment_diagonal_min_mag2": float(np.min(np.diag(systematics))),
        **metrics,
    }
    return rows, audit, [reconstruction]


def derive_union() -> tuple[list[dict[str, object]], dict[str, object]]:
    from astropy.io import fits

    matrix = np.asarray(fits.getdata(UNION_PATH), dtype=float)
    precision = matrix[1:, 1:]
    try:
        covariance = np.linalg.inv(precision)
    except np.linalg.LinAlgError:
        covariance = np.linalg.pinv(precision)
    n = precision.shape[0]
    idx = np.arange(n)
    sigma_total, variance_total = projected_sigma(covariance, idx)
    spec = ProjectionSpec(
        dataset="Union3p1_UNITY1p8",
        component_kind="documented_zero_point_budget_in_compressed_covariance",
        source_path=UNION_PATH,
        covariance_role="contained_in_release_precision",
        additive_prior=False,
        budget_kind="documented_important_dataset_zero_point",
    )
    row = make_row(
        spec,
        "dataset",
        "ALL",
        n,
        sigma_projected=0.01,
        variance_raw=0.01**2,
        sigma_documented=0.01,
        note=(
            "The 0.01 mag documented scale is already propagated through UNITY calibration "
            f"derivatives; total compressed common-mode sigma is {sigma_total:.6g} mag."
        ),
    )
    audit = {
        "N_compressed_nodes": n,
        "matrix_shape": list(matrix.shape),
        "total_compressed_common_mode_sigma_mag": sigma_total,
        "total_compressed_common_mode_variance_mag2": variance_total,
    }
    return [row], audit


def write_report(
    path: Path,
    table: pd.DataFrame,
    reconstructions: pd.DataFrame,
    source_inventory: pd.DataFrame,
    audit: dict[str, object],
) -> None:
    dataset_rows = table[table["level"] == "dataset"].copy()
    survey_rows = table[table["level"] == "survey"].copy()
    lines = [
        "# Release-Grounded Raw-MU Calibration Prior Registry",
        "",
        f"Generated: {datetime.now():%Y-%m-%d %H:%M:%S}",
        "",
        "## Interpretation",
        "",
        "Projected widths are calibration-mode budgets. They are not automatically additive priors. "
        "Every row in this registry is already represented in the corresponding release `STAT+SYS` "
        "product, so `sigma_add_mag = 0` and the primary likelihood does not add it again.",
        "",
        "Pantheon+ uses its calibration-only covariance. DES uses `C_total-C_stat` only as an "
        "all-systematics upper bound. Union3.1 records its documented 0.01 mag important-dataset "
        "zero-point scale, already propagated in the compressed release covariance.",
        "",
        "## Dataset Modes",
        "",
        dataset_rows[
            [
                "dataset",
                "N",
                "component_kind",
                "sigma_projected_mag",
                "sigma_documented_mag",
                "sigma_hierarchy_mag",
                "sigma_add_mag",
                "prior_action",
            ]
        ].to_markdown(index=False),
        "",
        "## Survey Modes",
        "",
        survey_rows[
            [
                "dataset",
                "component_kind",
                "group",
                "N",
                "sigma_projected_mag",
                "sigma_hierarchy_mag",
                "explicit_base",
            ]
        ].to_markdown(index=False),
        "",
        "## Explicit-Hierarchy Reconstruction",
        "",
        reconstructions[
            [
                "dataset",
                "component",
                "base",
                "status",
                "relative_frobenius_error",
                "captured_frobenius_power_fraction",
                "grouping_covariance_relative_frobenius_error",
                "released_total_relative_frobenius_error",
                "diagonal_rmse_mag2",
            ]
        ].to_markdown(index=False),
        "",
        "Pantheon+ grouping files include `STATONLY`; the program first forms each grouping "
        "increment as `C_grouping - C_stat`, then reconstructs `C_stat + C_grouped_offsets`. "
        "DES starts from `STATONLY` and approximates all released systematics with common modes. "
        "The released total-covariance likelihood remains primary; reconstruction rows are sensitivities.",
        "",
        "## Covariance Audit",
        "",
        "```json",
        json.dumps(audit, indent=2, sort_keys=True),
        "```",
        "",
        "## Source Metadata",
        "",
        source_inventory[["name", "role", "bytes", "sha256", "source_url"]].to_markdown(index=False),
        "",
        "The Fragilistic file is a 102 x 102 calibration-parameter covariance. It is archived and "
        "hashed here, but is not mapped directly into SN-distance offsets without the corresponding "
        "release Jacobian; the released distance-covariance grouping remains the applicable product.",
        "",
    ]
    for label, url in SOURCE_URLS.items():
        lines.append(f"- [{label}]({url})")
    lines.extend(
        [
            "",
            "## Double-Counting Guard",
            "",
            "A residual-offset sensitivity may be run only for a demonstrably external component "
            "absent from the release covariance. It must use the identical design and prior for FR "
            "and Lambda-CDM and must not be presented as independent evidence.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--min-survey-n", type=int, default=20)
    return parser


def main() -> None:
    cli = build_parser().parse_args()
    cli.outdir.mkdir(parents=True, exist_ok=True)

    pantheon_rows, pantheon_audit, pantheon_reconstruction = derive_pantheon(cli.min_survey_n)
    des_rows, des_audit, des_reconstruction = derive_des(cli.min_survey_n)
    union_rows, union_audit = derive_union()
    table = pd.DataFrame(pantheon_rows + des_rows + union_rows)
    table = table.sort_values(["dataset", "level", "group"]).reset_index(drop=True)
    reconstructions = pd.DataFrame(pantheon_reconstruction + des_reconstruction)

    sources = {
        "PantheonPlusSH0ES_data": source_record(PANTHEON_DATA, "distance table", "PantheonPlusSH0ES_release"),
        "PantheonPlusSH0ES_total": source_record(PANTHEON_TOTAL, "STAT+SYS covariance", "PantheonPlusSH0ES_release"),
        "PantheonPlusSH0ES_stat": source_record(
            PANTHEON_STAT, "STATONLY covariance", "PantheonPlusSH0ES_STATONLY_raw"
        ),
        "PantheonPlusSH0ES_calibration": source_record(
            PANTHEON_CALIB, "calibration-only covariance", "PantheonPlusSH0ES_CALIB_raw"
        ),
        "PantheonPlusSH0ES_CSP_RECAL": source_record(
            PANTHEON_CSP_RECAL, "CSP recalibration covariance", "PantheonPlusSH0ES_CSP_RECAL_raw"
        ),
        "PantheonPlusSH0ES_HSTCALSPEC": source_record(
            PANTHEON_HSTCALSPEC, "HST CALSPEC covariance", "PantheonPlusSH0ES_HSTCALSPEC_raw"
        ),
        "PantheonPlusSH0ES_FRAGILISTIC": source_record(
            PANTHEON_FRAGILISTIC,
            "102-parameter Fragilistic calibration covariance",
            "PantheonPlusSH0ES_FRAGILISTIC_raw",
        ),
        "DES_data": source_record(DES_DATA, "distance table", "DES_SN5YR_release"),
        "DES_total": source_record(DES_TOTAL, "STAT+SYS precision", "DES_SN5YR_release"),
        "DES_stat": source_record(DES_STAT, "STATONLY precision", "DES_SN5YR_release"),
        "Union3p1": source_record(UNION_PATH, "compressed distance/precision", "Union3_paper"),
    }
    audit = {
        "PantheonPlusSH0ES": pantheon_audit,
        "DES_Dovekie_raw": des_audit,
        "Union3p1_UNITY1p8": union_audit,
    }
    source_inventory = pd.DataFrame(
        [{"name": name, **record} for name, record in sources.items()]
    ).sort_values("name")
    registry = {
        "analysis_id": "rawmu_release_grounded_calibration_2026-07-18",
        "generated": datetime.now().isoformat(timespec="seconds"),
        "primary_frame": {"PantheonPlusSH0ES": "zHD", "DES_Dovekie_raw": "zHD", "Union3p1_UNITY1p8": "z"},
        "double_count_guard": "release-internal projected modes are not added to STAT+SYS",
        "residual_sensitivity_is_independent_evidence": False,
        "source_urls": SOURCE_URLS,
        "sources": sources,
        "audit": audit,
        "reconstructions": reconstructions.replace({np.nan: None}).to_dict(orient="records"),
        "rows": table.replace({np.nan: None}).to_dict(orient="records"),
    }

    csv_path = cli.outdir / f"rawmu_release_grounded_prior_registry_{DATE_TAG}.csv"
    reconstruction_path = cli.outdir / f"rawmu_release_grounded_covariance_reconstruction_{DATE_TAG}.csv"
    source_inventory_path = cli.outdir / f"rawmu_release_grounded_source_inventory_{DATE_TAG}.csv"
    json_path = cli.outdir / f"rawmu_release_grounded_prior_registry_{DATE_TAG}.json"
    report_path = cli.outdir / f"rawmu_release_grounded_prior_derivation_report_{DATE_TAG}.md"
    table.to_csv(csv_path, index=False)
    reconstructions.to_csv(reconstruction_path, index=False)
    source_inventory.to_csv(source_inventory_path, index=False)
    json_path.write_text(json.dumps(registry, indent=2, sort_keys=True), encoding="utf-8")
    write_report(report_path, table, reconstructions, source_inventory, audit)
    print(f"Saved prior table: {csv_path}")
    print(f"Saved covariance reconstruction: {reconstruction_path}")
    print(f"Saved source inventory: {source_inventory_path}")
    print(f"Saved prior registry: {json_path}")
    print(f"Saved derivation report: {report_path}")


if __name__ == "__main__":
    main()
