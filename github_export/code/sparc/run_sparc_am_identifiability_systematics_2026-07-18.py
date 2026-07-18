#!/usr/bin/env python3
r"""Preregister and run the SPARC antimatter-identifiability systematics screen.

The six previously excluded persistent-negative galaxies are a fixed development
set. Controls are matched without using RAR or PLAMB residuals. The programme
first records whether the implemented rotation-curve likelihood contains a
matter/antimatter discriminator, then profiles conventional galaxy nuisances
under a locked sensitivity matrix.

This is a fixed-global-parameter MAP profile screen. It is not a posterior model
probability, a discovery test, or evidence that any galaxy contains antimatter.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment, minimize
from scipy.stats import binomtest


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import fit_sparc_hierarchical_map as fit
from diagnose_plamb_sparc_h0_nuisance import acceleration_cH0_over_2pi
from diagnose_plamb_sparc_rotation import KPC_M, load_dataset


RUN_DATE = "2026-07-18"
TASK_ROOT = Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu")
OUTPUTS = TASK_ROOT / "outputs"
SAMPLE = ROOT / "external_datasets" / "current_cosmology_sources" / "SPARC" / "sparc_galaxy_sample.csv"
POINTS = ROOT / "external_datasets" / "current_cosmology_sources" / "SPARC" / "sparc_rotation_points.csv"
MAP_SUMMARY = (
    ROOT
    / "plamb_runs"
    / "diagnostics"
    / "sparc_hierarchical_map"
    / "optical_depth_hierarchical_20260714"
    / "sparc_hierarchical_map_summary.csv"
)
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "sparc_am_identifiability_systematics_20260718"
DEFAULT_EXPORT = ROOT / "github_export" / "results" / "2026-07-18" / "sparc_am"

CANDIDATES = (
    "UGC06787",
    "UGC03580",
    "UGC02487",
    "NGC5985",
    "UGC09133",
    "NGC2903",
)
RESERVED_REPLICATION = (
    "NGC2841",
    "UGC02953",
    "UGC07399",
    "UGC00128",
    "NGC2403",
)
MODELS = ("RAR", "PLAMB_OPTICAL_DEPTH_KAPPA_P")
H0 = 65.0
SPARC_H0_REF = 73.0
MATCHES_PER_CANDIDATE = 2


@dataclass(frozen=True)
class Scenario:
    name: str
    description: str
    sigma_ln_ml: float = 0.25
    ml_mode: str = "coupled"
    distance_floor_frac: float = 0.03
    distance_prior: str = "gaussian"
    hubble_prior_mode: str = "model_h0_rescale"
    inclination_nuisance: bool = False
    velocity_jitter_kms: float = 0.0
    radial_correlation: float = 0.0
    gas_scale_sigma: float = 0.0


SCENARIOS = (
    Scenario("baseline", "Published coupled stellar scaling, Gaussian distance prior and 3 km/s error floor."),
    Scenario("stellar_ml_tight", "Tighter 15 per cent log stellar mass-to-light prior.", sigma_ln_ml=0.15),
    Scenario("stellar_ml_wide", "Wider 50 per cent log stellar mass-to-light prior.", sigma_ln_ml=0.50),
    Scenario(
        "disk_bulge_decoupled",
        "Independent galaxy-level disc and bulge mass-to-light multipliers.",
        ml_mode="separate",
    ),
    Scenario(
        "distance_published_frame",
        "No model-H0 rescaling of Hubble-flow distance-prior centres.",
        hubble_prior_mode="published",
    ),
    Scenario(
        "distance_floor_10pct",
        "Minimum fractional distance uncertainty increased to 10 per cent.",
        distance_floor_frac=0.10,
    ),
    Scenario(
        "distance_student_t3",
        "Robust Student-t(3) distance-pull penalty.",
        distance_prior="student_t3",
    ),
    Scenario(
        "inclination_nuisance",
        "First-order line-of-sight inclination nuisance with the SPARC inclination error.",
        inclination_nuisance=True,
    ),
    Scenario(
        "velocity_jitter_5kms",
        "Five km/s non-circular-motion jitter added in quadrature.",
        velocity_jitter_kms=5.0,
    ),
    Scenario(
        "radial_correlation_rho035",
        "Within-galaxy AR(1) radial residual correlation fixed at rho=0.35.",
        radial_correlation=0.35,
    ),
    Scenario(
        "gas_scale_20pct",
        "Galaxy-level gas contribution scaling with a 20 per cent log prior.",
        gas_scale_sigma=0.20,
    ),
    Scenario(
        "combined_conventional",
        "Combined wide separate stellar, robust distance, inclination, jitter, radial-correlation and gas terms.",
        sigma_ln_ml=0.50,
        ml_mode="separate",
        distance_floor_frac=0.10,
        distance_prior="student_t3",
        inclination_nuisance=True,
        velocity_jitter_kms=5.0,
        radial_correlation=0.35,
        gas_scale_sigma=0.20,
    ),
)

MATCH_FEATURES = (
    "T",
    "log_L36",
    "log_SBeff",
    "log_MHI",
    "log_Vflat",
    "Inc_deg",
    "frac_distance_error",
    "bulge_proxy",
    "log_n_points",
    "Q",
)
MATCH_WEIGHTS = np.asarray([0.75, 1.25, 1.00, 1.00, 1.25, 0.75, 0.75, 1.25, 0.75, 0.50])


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, Any]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        values: list[str] = []
        for column in columns:
            value = row.get(column, "")
            if isinstance(value, float):
                values.append(f"{value:.6g}" if math.isfinite(value) else "")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def transformed_catalogue(sample: pd.DataFrame, points: pd.DataFrame) -> pd.DataFrame:
    sample = sample.copy()
    points = points.copy()
    for column in ["Vgas_kms", "Vdisk_kms", "Vbul_kms"]:
        points[column] = pd.to_numeric(points[column], errors="coerce").fillna(0.0)
    vgas2 = np.abs(points["Vgas_kms"] * points["Vgas_kms"])
    vdisk2 = 0.5 * points["Vdisk_kms"] ** 2
    vbul2 = 0.7 * points["Vbul_kms"] ** 2
    denom = np.maximum(vgas2 + vdisk2 + vbul2, 1e-12)
    points["bulge_fraction_proxy"] = vbul2 / denom
    point_summary = points.groupby("Galaxy").agg(
        bulge_proxy=("bulge_fraction_proxy", "median"),
        n_points=("Galaxy", "size"),
    )
    sample = sample.merge(point_summary, left_on="Galaxy", right_index=True, how="left")
    sample["has_bulge"] = sample["bulge_proxy"].fillna(0.0) > 1e-12
    sample["frac_distance_error"] = pd.to_numeric(sample["e_D_Mpc"], errors="coerce") / np.maximum(
        pd.to_numeric(sample["D_Mpc"], errors="coerce"), 1e-12
    )
    for source, target in [
        ("L36_1e9_Lsun", "log_L36"),
        ("SBeff_Lsun_pc2", "log_SBeff"),
        ("MHI_1e9_Msun", "log_MHI"),
        ("Vflat_kms", "log_Vflat"),
    ]:
        sample[target] = np.log10(np.maximum(pd.to_numeric(sample[source], errors="coerce"), 1e-12))
    sample["log_n_points"] = np.log10(np.maximum(pd.to_numeric(sample["n_points"], errors="coerce"), 1.0))
    for column in ["T", "Inc_deg", "Q", "f_D"]:
        sample[column] = pd.to_numeric(sample[column], errors="coerce")
    return sample


def robust_standardise(frame: pd.DataFrame, columns: tuple[str, ...]) -> tuple[pd.DataFrame, dict[str, dict[str, float]]]:
    result = pd.DataFrame(index=frame.index)
    scaling: dict[str, dict[str, float]] = {}
    for column in columns:
        values = pd.to_numeric(frame[column], errors="coerce")
        median = float(values.median())
        q25 = float(values.quantile(0.25))
        q75 = float(values.quantile(0.75))
        scale = q75 - q25
        if not math.isfinite(scale) or scale <= 1e-12:
            scale = float(values.std(ddof=0))
        if not math.isfinite(scale) or scale <= 1e-12:
            scale = 1.0
        result[column] = (values.fillna(median) - median) / scale
        scaling[column] = {"median": median, "scale_iqr": scale}
    return result, scaling


def select_controls(catalogue: pd.DataFrame) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    candidate_rows = catalogue[catalogue["Galaxy"].isin(CANDIDATES)].copy()
    if set(candidate_rows["Galaxy"]) != set(CANDIDATES):
        missing = sorted(set(CANDIDATES) - set(candidate_rows["Galaxy"]))
        raise ValueError(f"Candidate galaxies missing from SPARC catalogue: {missing}")
    candidate_distance_methods = set(candidate_rows["f_D"].astype(int))
    if candidate_distance_methods != {1}:
        raise ValueError(f"Locked candidate set no longer has a common f_D=1 distance method: {candidate_distance_methods}")
    eligible = catalogue[
        (catalogue["Q"] <= 2)
        & (catalogue["f_D"] == 1)
        & (~catalogue["Galaxy"].isin(CANDIDATES))
        & (~catalogue["Galaxy"].isin(RESERVED_REPLICATION))
    ].copy()
    combined = pd.concat([candidate_rows, eligible], ignore_index=True)
    standardised, scaling = robust_standardise(combined, MATCH_FEATURES)
    candidate_z = standardised.iloc[: len(candidate_rows)].to_numpy(dtype=float)
    eligible_z = standardised.iloc[len(candidate_rows) :].to_numpy(dtype=float)
    slots = np.repeat(np.arange(len(candidate_rows)), MATCHES_PER_CANDIDATE)
    weighted_candidate = candidate_z[slots] * np.sqrt(MATCH_WEIGHTS)[None, :]
    weighted_eligible = eligible_z * np.sqrt(MATCH_WEIGHTS)[None, :]
    cost = np.sqrt(np.sum((weighted_candidate[:, None, :] - weighted_eligible[None, :, :]) ** 2, axis=2))
    candidate_has_bulge = candidate_rows["has_bulge"].to_numpy(dtype=bool)[slots]
    eligible_has_bulge = eligible["has_bulge"].to_numpy(dtype=bool)
    cost += 1.0e6 * (candidate_has_bulge[:, None] != eligible_has_bulge[None, :])
    row_index, column_index = linear_sum_assignment(cost)
    if len(row_index) != len(slots):
        raise RuntimeError("Could not allocate the preregistered number of unique controls")
    assignments: list[dict[str, Any]] = []
    for row, column in sorted(zip(row_index, column_index), key=lambda item: item[0]):
        candidate_position = int(slots[row])
        candidate = str(candidate_rows.iloc[candidate_position]["Galaxy"])
        control = str(eligible.iloc[column]["Galaxy"])
        assignments.append(
            {
                "candidate": candidate,
                "match_slot": 1 + sum(item["candidate"] == candidate for item in assignments),
                "control": control,
                "weighted_match_distance": float(cost[row, column]),
            }
        )
    selected_controls = [item["control"] for item in assignments]
    candidate_table = candidate_rows.set_index("Galaxy")
    control_table = eligible[eligible["Galaxy"].isin(selected_controls)].set_index("Galaxy")
    balance_rows: list[dict[str, Any]] = []
    all_z = standardised.copy()
    all_z["Galaxy"] = combined["Galaxy"].values
    for feature in MATCH_FEATURES:
        candidate_values = all_z[all_z["Galaxy"].isin(CANDIDATES)][feature]
        control_values = all_z[all_z["Galaxy"].isin(selected_controls)][feature]
        balance_rows.append(
            {
                "feature": feature,
                "candidate_mean_robust_z": float(candidate_values.mean()),
                "control_mean_robust_z": float(control_values.mean()),
                "standardised_mean_difference": float(candidate_values.mean() - control_values.mean()),
                "candidate_median_original": float(pd.to_numeric(candidate_table[feature], errors="coerce").median()),
                "control_median_original": float(pd.to_numeric(control_table[feature], errors="coerce").median()),
            }
        )
    metadata = {
        "eligible_control_population": int(len(eligible)),
        "matches_per_candidate": MATCHES_PER_CANDIDATE,
        "unique_controls": len(set(selected_controls)),
        "features": list(MATCH_FEATURES),
        "weights": {name: float(weight) for name, weight in zip(MATCH_FEATURES, MATCH_WEIGHTS)},
        "robust_scaling": scaling,
        "outcome_blind": True,
        "exact_match": "presence or absence of a non-zero SPARC bulge component",
        "excluded_replication_galaxies": list(RESERVED_REPLICATION),
    }
    return assignments, balance_rows, metadata


def candidate_inventory(catalogue: pd.DataFrame) -> list[dict[str, Any]]:
    columns = [
        "Galaxy",
        "T",
        "D_Mpc",
        "e_D_Mpc",
        "f_D",
        "Inc_deg",
        "e_Inc_deg",
        "L36_1e9_Lsun",
        "SBeff_Lsun_pc2",
        "MHI_1e9_Msun",
        "Vflat_kms",
        "Q",
        "bulge_proxy",
        "n_points",
    ]
    return catalogue[catalogue["Galaxy"].isin(CANDIDATES)][columns].sort_values("Galaxy").to_dict("records")


def identifiability_matrix() -> list[dict[str, Any]]:
    return [
        {
            "observable_or_term": "Hydrogen/antihydrogen atomic spectrum",
            "matter_to_antimatter_transformation": "Invariant under CPT to current experimental precision",
            "implemented_in_sparc_likelihood": "No particle-sign variable",
            "identifies_antimatter": "No",
        },
        {
            "observable_or_term": "Attractive gravitational rotation curve",
            "matter_to_antimatter_transformation": "No sign reversal in standard gravity; ALPHA-g is consistent with downward attraction",
            "implemented_in_sparc_likelihood": "Positive gbar and positive predicted speed",
            "identifies_antimatter": "No",
        },
        {
            "observable_or_term": "PLAMB optical depth tau=(gbar/g0)^p",
            "matter_to_antimatter_transformation": "Unspecified by the supplied FR material",
            "implemented_in_sparc_likelihood": "Same positive tau for every galaxy",
            "identifies_antimatter": "No",
        },
        {
            "observable_or_term": "PLAMB predicted acceleration gbar/(1-exp(-tau))",
            "matter_to_antimatter_transformation": "Unspecified by the supplied FR material",
            "implemented_in_sparc_likelihood": "Single branch with no charge-conjugation state",
            "identifies_antimatter": "No",
        },
        {
            "observable_or_term": "Matter-antimatter boundary annihilation",
            "matter_to_antimatter_transformation": "Produces an independent gamma-ray and secondary-particle signature at contact",
            "implemented_in_sparc_likelihood": "Absent",
            "identifies_antimatter": "Potentially, with external data",
        },
    ]


def write_identifiability_document(path: Path, matrix: list[dict[str, Any]]) -> None:
    lines = [
        "# FR Matter-to-Antimatter Identifiability Preregistration",
        "",
        f"Date: {RUN_DATE}",
        "",
        "## Question",
        "",
        "Can the currently implemented SPARC RAR/PLAMB rotation-curve likelihood identify whether a galaxy is constructed from matter or antimatter?",
        "",
        "## Implemented Equations",
        "",
        r"\[",
        r"\begin{aligned}",
        r"\tau_g       &= \left(\frac{g_{{\rm bar},g}}{g_0}\right)^p, \\",
        r"g_{{\rm pred},g} &= \frac{g_{{\rm bar},g}}{1-\exp(-\tau_g)}.",
        r"\end{aligned}",
        r"\]",
        "",
        "Every input is positive and the code contains no matter/antimatter state. Therefore, under the implemented likelihood,",
        "",
        r"\[",
        r"\begin{aligned}",
        r"\mathcal{L}_{\rm RC}(D_g\mid M,\theta)       &= \mathcal{L}_{\rm RC}(D_g\mid \bar M,\theta), \\",
        r"B^{\rm RC}_{\bar M,M}                         &= 1, \\",
        r"P(\bar M\mid D_{\rm RC})                     &= P(\bar M).",
        r"\end{aligned}",
        r"\]",
        "",
        "The final equality holds only for the present rotation-curve likelihood. It does not rule out an FR matter-antimatter theory; it shows that an additional, independently derived transformation law is required before the likelihood can test one.",
        "",
        "## Locked Transformation Matrix",
        "",
        markdown_table(matrix, ["observable_or_term", "matter_to_antimatter_transformation", "implemented_in_sparc_likelihood", "identifies_antimatter"]),
        "",
        "## Identifiability Decision",
        "",
        "**Gate 1 result: NOT IDENTIFIABLE in the current SPARC likelihood.**",
        "",
        "No arbitrary sign flip will be fitted. A future FR antimatter branch must specify before data fitting:",
        "",
        "1. which field or source variable changes under charge conjugation;",
        "2. whether the change affects gbar, g0, kappa, p or a separate field term;",
        "3. the predicted radial shape and sign;",
        "4. its matter-matter, antimatter-antimatter and matter-antimatter limits; and",
        "5. a compatible annihilation/contact model for external gamma-ray testing.",
        "",
        "## External Physics Boundary",
        "",
        "Antihydrogen spectroscopy is consistent with hydrogen, and ALPHA-g observes behaviour consistent with attractive terrestrial gravity. The independent astronomical discriminator is therefore expected at matter-antimatter contact, not in an otherwise isolated optical or 21-cm rotation curve.",
        "",
        "References:",
        "",
        "- https://www.nature.com/articles/nature21040",
        "- https://www.nature.com/articles/s41586-023-06527-1",
        "- https://arxiv.org/abs/2103.10073",
        "- https://arxiv.org/abs/0808.1122",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def write_preregistration(path: Path, config: dict[str, Any], assignments: list[dict[str, Any]]) -> None:
    lines = [
        "# SPARC Antimatter-Candidate Conventional-Systematics Preregistration",
        "",
        f"- Written before outcome profiling: `{config['written_before_outcomes_utc']}`",
        f"- Candidate galaxies: `{len(CANDIDATES)}`",
        f"- Unique matched controls: `{len({row['control'] for row in assignments})}`",
        f"- Reserved replication galaxies: `{len(RESERVED_REPLICATION)}`",
        "",
        "## Claim Boundary",
        "",
        "The candidate set was selected previously from model residuals. Its profile scores are therefore descriptive and cannot supply an unbiased discovery p-value. Controls are selected here without using any RAR or PLAMB outcome. The reserved replication set is not examined by this programme.",
        "",
        "Gate 1 has already determined that the current rotation-curve likelihood is matter/antimatter non-identifiable. Consequently, persistence in this matrix means only that a conventional-systematics challenge did not erase a model-relative anomaly. It does not identify antimatter.",
        "",
        "## Locked Candidate Set",
        "",
        ", ".join(f"`{name}`" for name in CANDIDATES),
        "",
        "## Locked Control Matching",
        "",
        f"{MATCHES_PER_CANDIDATE} unique f_D=1, Q<=2 controls per candidate are assigned by minimum total weighted robust-standardised Euclidean distance, with exact matching on the presence or absence of a SPARC bulge component. Candidate residuals and model scores are excluded from matching. The matching variables are morphology, luminosity, effective surface brightness, H I mass, flat speed, inclination, fractional distance error, bulge proxy, point count and quality.",
        "",
        markdown_table(assignments, ["candidate", "match_slot", "control", "weighted_match_distance"]),
        "",
        "## Locked Sensitivity Matrix",
        "",
        markdown_table(config["scenarios"], ["name", "description"]),
        "",
        "## Profile Method",
        "",
        "The all-Q2 July 14 MAP global RAR and PLAMB parameters are held fixed. For every candidate and control, local nuisances are re-profiled independently in every scenario. Positive Delta objective = objective_PLAMB - objective_RAR means PLAMB fits worse than RAR. Both the data term and the data-plus-prior objective are retained.",
        "",
        "A candidate is labelled persistent conventional tension only when PLAMB is worse in at least 80 per cent of locked scenarios and remains worse in the combined-conventional scenario. This label is a robustness flag, not an antimatter classification.",
        "",
        "## Reproducibility Locks",
        "",
        "```json",
        json.dumps(config, indent=2),
        "```",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def load_global_parameters() -> dict[str, dict[str, float]]:
    summary = pd.read_csv(MAP_SUMMARY)
    selected = summary[(summary["split"] == "all_Q2") & summary["model"].isin(MODELS)]
    if set(selected["model"]) != set(MODELS):
        raise ValueError("Could not find both locked all-Q2 MAP model rows")
    output: dict[str, dict[str, float]] = {}
    for _, row in selected.iterrows():
        output[str(row["model"])] = {
            "ydisk": float(row["ydisk"]),
            "ybul": float(row["ybul"]),
            "log10_gdagger": float(row["log10_gdagger"]) if pd.notna(row["log10_gdagger"]) else float("nan"),
            "log10_kappa": float(row["log10_kappa"]) if pd.notna(row["log10_kappa"]) else float("nan"),
            "bridge_exponent": float(row["bridge_exponent"]) if pd.notna(row["bridge_exponent"]) else float("nan"),
            "H0": float(row["H0_trial"]),
        }
    return output


def correlation_precision(size: int, rho: float) -> np.ndarray | None:
    if rho <= 0.0:
        return None
    indices = np.arange(size)
    correlation = rho ** np.abs(indices[:, None] - indices[None, :])
    return np.linalg.inv(correlation)


def deterministic_optimisation_starts(
    centre: np.ndarray,
    bounds: list[tuple[float, float]],
    key: str,
) -> list[np.ndarray]:
    starts = [np.asarray(centre, dtype=float)]
    for index, (lower, upper) in enumerate(bounds):
        for fraction in (0.25, 0.75):
            trial = np.asarray(centre, dtype=float).copy()
            trial[index] = lower + fraction * (upper - lower)
            starts.append(trial)
    seed = int(hashlib.sha256(key.encode("utf-8")).hexdigest()[:16], 16) % (2**32)
    rng = np.random.default_rng(seed)
    lower = np.asarray([item[0] for item in bounds], dtype=float)
    upper = np.asarray([item[1] for item in bounds], dtype=float)
    for _ in range(6):
        starts.append(lower + rng.uniform(0.10, 0.90, size=len(bounds)) * (upper - lower))
    unique: list[np.ndarray] = []
    seen: set[tuple[float, ...]] = set()
    for start in starts:
        signature = tuple(np.round(start, 12))
        if signature not in seen:
            seen.add(signature)
            unique.append(start)
    return unique


def profile_galaxy(
    galaxy: str,
    model: str,
    scenario: Scenario,
    data: Any,
    sample_info: dict[str, dict[str, float]],
    globals_: dict[str, float],
) -> dict[str, Any]:
    mask = np.asarray(data.galaxy == galaxy, dtype=bool)
    order = np.argsort(data.rad_kpc[mask])
    rad0 = np.asarray(data.rad_kpc[mask][order], dtype=float)
    vobs = np.asarray(data.vobs[mask][order], dtype=float)
    sigma0 = np.asarray(data.sigma_v[mask][order], dtype=float)
    vgas0 = np.asarray(data.vgas[mask][order], dtype=float)
    vdisk0 = np.asarray(data.vdisk[mask][order], dtype=float)
    vbul0 = np.asarray(data.vbul[mask][order], dtype=float)
    item = sample_info[galaxy]
    d = max(item["D_Mpc"], 1e-12)
    sigma_logd = max(item["e_D_Mpc"] / d, scenario.distance_floor_frac)
    mu_logd = 0.0
    if scenario.hubble_prior_mode == "model_h0_rescale" and int(round(item["f_D"])) == 1:
        mu_logd = math.log(SPARC_H0_REF / globals_["H0"])
    inc0_deg = float(item["Inc_deg"])
    inc_sigma_deg = max(float(item["e_Inc_deg"]), 1.0)
    has_bulge = bool(np.any(np.abs(vbul0) > 1e-12))
    precision = correlation_precision(len(rad0), scenario.radial_correlation)

    parameter_names = ["logd", "logeta_disk"]
    start = [mu_logd, 0.0]
    distance_half_width = max((6.0 if scenario.distance_prior == "student_t3" else 4.0) * sigma_logd, 0.08)
    bounds: list[tuple[float, float]] = [
        (mu_logd - distance_half_width, mu_logd + distance_half_width),
        (-4.0 * scenario.sigma_ln_ml, 4.0 * scenario.sigma_ln_ml),
    ]
    if scenario.ml_mode == "separate" and has_bulge:
        parameter_names.append("logeta_bulge")
        start.append(0.0)
        bounds.append((-4.0 * scenario.sigma_ln_ml, 4.0 * scenario.sigma_ln_ml))
    if scenario.inclination_nuisance:
        parameter_names.append("delta_inc_deg")
        start.append(0.0)
        lower = max(-3.0 * inc_sigma_deg, 15.0 - inc0_deg)
        upper = min(3.0 * inc_sigma_deg, 90.0 - inc0_deg)
        if lower > upper:
            raise ValueError(
                f"Invalid inclination bounds for {galaxy}: published={inc0_deg}, "
                f"sigma={inc_sigma_deg}, bounds=({lower}, {upper})"
            )
        bounds.append((lower, upper))
    if scenario.gas_scale_sigma > 0.0:
        parameter_names.append("logeta_gas")
        start.append(0.0)
        bounds.append((-4.0 * scenario.gas_scale_sigma, 4.0 * scenario.gas_scale_sigma))
    index = {name: position for position, name in enumerate(parameter_names)}

    def evaluate(theta: np.ndarray) -> tuple[float, float, float, dict[str, float]]:
        logd = float(theta[index["logd"]])
        logeta_disk = float(theta[index["logeta_disk"]])
        logeta_bulge = logeta_disk
        if "logeta_bulge" in index:
            logeta_bulge = float(theta[index["logeta_bulge"]])
        delta_inc_deg = float(theta[index["delta_inc_deg"]]) if "delta_inc_deg" in index else 0.0
        logeta_gas = float(theta[index["logeta_gas"]]) if "logeta_gas" in index else 0.0
        scale = math.exp(logd)
        radius = rad0 * scale
        vgas = vgas0 * math.sqrt(scale * math.exp(logeta_gas))
        vdisk = vdisk0 * math.sqrt(scale * math.exp(logeta_disk))
        vbul = vbul0 * math.sqrt(scale * math.exp(logeta_bulge))
        v2_bar = vgas * np.abs(vgas) + globals_["ydisk"] * vdisk**2 + globals_["ybul"] * vbul**2
        gbar = np.maximum(v2_bar, 1e-18) * 1.0e6 / (radius * KPC_M)
        if model == "RAR":
            tau = np.sqrt(np.maximum(gbar / (10.0 ** globals_["log10_gdagger"]), 1e-30))
        elif model == "PLAMB_OPTICAL_DEPTH_KAPPA_P":
            g0 = (10.0 ** globals_["log10_kappa"]) * acceleration_cH0_over_2pi(globals_["H0"])
            tau = np.maximum(gbar / g0, 1e-30) ** globals_["bridge_exponent"]
        else:
            raise ValueError(model)
        gpred = gbar / np.maximum(1.0 - np.exp(-tau), 1e-12)
        prediction = np.sqrt(np.maximum(gpred, 0.0) * radius * KPC_M) / 1000.0
        sigma = np.sqrt(sigma0**2 + scenario.velocity_jitter_kms**2)
        if scenario.inclination_nuisance:
            sin_published = math.sin(math.radians(inc0_deg))
            sin_trial = math.sin(math.radians(inc0_deg + delta_inc_deg))
            residual = (vobs * sin_published - prediction * sin_trial) / (sigma * sin_published)
        else:
            residual = (vobs - prediction) / sigma
        data_objective = float(residual @ residual if precision is None else residual @ precision @ residual)
        distance_pull = (logd - mu_logd) / sigma_logd
        if scenario.distance_prior == "student_t3":
            distance_penalty = 4.0 * math.log1p(distance_pull**2 / 3.0)
        else:
            distance_penalty = distance_pull**2
        prior = distance_penalty + (logeta_disk / scenario.sigma_ln_ml) ** 2
        if "logeta_bulge" in index:
            prior += (logeta_bulge / scenario.sigma_ln_ml) ** 2
        if "delta_inc_deg" in index:
            prior += (delta_inc_deg / inc_sigma_deg) ** 2
        if "logeta_gas" in index:
            prior += (logeta_gas / scenario.gas_scale_sigma) ** 2
        details = {
            "logd": logd,
            "distance_pull": distance_pull,
            "distance_multiplier": math.exp(logd),
            "logeta_disk": logeta_disk,
            "disk_ml_multiplier": math.exp(logeta_disk),
            "logeta_bulge": logeta_bulge,
            "bulge_ml_multiplier": math.exp(logeta_bulge),
            "delta_inc_deg": delta_inc_deg,
            "inclination_pull": delta_inc_deg / inc_sigma_deg,
            "trial_inclination_deg": inc0_deg + delta_inc_deg,
            "logeta_gas": logeta_gas,
            "gas_multiplier": math.exp(logeta_gas),
        }
        return data_objective + prior, data_objective, prior, details

    starts = deterministic_optimisation_starts(
        np.asarray(start, dtype=float),
        bounds,
        f"{RUN_DATE}|{scenario.name}|{model}|{galaxy}",
    )
    results = [
        minimize(
            lambda theta: evaluate(np.asarray(theta, dtype=float))[0],
            trial,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 700, "maxfun": 30000, "ftol": 1e-10, "maxls": 50},
        )
        for trial in starts
    ]
    objectives = [evaluate(np.asarray(item.x, dtype=float))[0] for item in results]
    best_index = int(np.argmin(objectives))
    result = results[best_index]
    total, data_objective, prior, details = evaluate(np.asarray(result.x, dtype=float))
    return {
        "scenario": scenario.name,
        "model": model,
        "galaxy": galaxy,
        "N_points": len(rad0),
        "success": bool(result.success),
        "message": str(result.message),
        "optimisation_starts": len(starts),
        "successful_starts": sum(bool(item.success) for item in results),
        "best_start_index": best_index,
        "profile_objective_start_spread": float(max(objectives) - min(objectives)),
        "profile_objective": total,
        "data_objective": data_objective,
        "prior_objective": prior,
        **details,
    }


def profile_matrix(
    scenarios: tuple[Scenario, ...],
    assignments: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    data = load_dataset(SAMPLE, POINTS, quality_max=2, distance_method="all", err_floor=3.0)
    sample = pd.read_csv(SAMPLE)
    sample_info: dict[str, dict[str, float]] = {}
    for _, row in sample.iterrows():
        sample_info[str(row["Galaxy"])] = {
            key: float(row[key]) for key in ["D_Mpc", "e_D_Mpc", "f_D", "Inc_deg", "e_Inc_deg"]
        }
    globals_by_model = load_global_parameters()
    controls = sorted({row["control"] for row in assignments})
    galaxies = list(CANDIDATES) + controls
    available = set(str(value) for value in data.galaxy)
    missing = sorted(set(galaxies) - available)
    if missing:
        raise ValueError(f"Profile galaxies missing from Q<=2 rotation data: {missing}")
    rows: list[dict[str, Any]] = []
    total = len(scenarios) * len(galaxies) * len(MODELS)
    completed = 0
    for scenario in scenarios:
        print(f"=== {scenario.name} ===", flush=True)
        for galaxy in galaxies:
            for model in MODELS:
                rows.append(profile_galaxy(galaxy, model, scenario, data, sample_info, globals_by_model[model]))
                completed += 1
            if completed % 24 == 0 or completed == total:
                print(f"profiled {completed}/{total}", flush=True)
    return rows


def paired_deltas(profile_rows: list[dict[str, Any]], assignments: list[dict[str, Any]]) -> list[dict[str, Any]]:
    lookup = {(row["scenario"], row["galaxy"], row["model"]): row for row in profile_rows}
    control_to_candidate = {row["control"]: row["candidate"] for row in assignments}
    galaxies = sorted(set(CANDIDATES) | set(control_to_candidate))
    rows: list[dict[str, Any]] = []
    for scenario in sorted({row["scenario"] for row in profile_rows}):
        for galaxy in galaxies:
            rar = lookup[(scenario, galaxy, "RAR")]
            plamb = lookup[(scenario, galaxy, "PLAMB_OPTICAL_DEPTH_KAPPA_P")]
            rows.append(
                {
                    "scenario": scenario,
                    "galaxy": galaxy,
                    "group": "candidate" if galaxy in CANDIDATES else "matched_control",
                    "matched_candidate": galaxy if galaxy in CANDIDATES else control_to_candidate[galaxy],
                    "N_points": rar["N_points"],
                    "delta_profile_objective_plamb_minus_rar": plamb["profile_objective"] - rar["profile_objective"],
                    "delta_data_objective_plamb_minus_rar": plamb["data_objective"] - rar["data_objective"],
                    "rar_profile_objective": rar["profile_objective"],
                    "plamb_profile_objective": plamb["profile_objective"],
                    "rar_distance_pull": rar["distance_pull"],
                    "plamb_distance_pull": plamb["distance_pull"],
                    "rar_disk_ml_multiplier": rar["disk_ml_multiplier"],
                    "plamb_disk_ml_multiplier": plamb["disk_ml_multiplier"],
                    "rar_bulge_ml_multiplier": rar["bulge_ml_multiplier"],
                    "plamb_bulge_ml_multiplier": plamb["bulge_ml_multiplier"],
                    "rar_delta_inc_deg": rar["delta_inc_deg"],
                    "plamb_delta_inc_deg": plamb["delta_inc_deg"],
                    "rar_gas_multiplier": rar["gas_multiplier"],
                    "plamb_gas_multiplier": plamb["gas_multiplier"],
                    "both_optimisers_succeeded": bool(rar["success"] and plamb["success"]),
                }
            )
    return rows


def summarise_deltas(delta_rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    frame = pd.DataFrame(delta_rows)
    summary_rows: list[dict[str, Any]] = []
    for scenario, selected in frame.groupby("scenario"):
        candidates = selected[selected["group"] == "candidate"]
        controls = selected[selected["group"] == "matched_control"]
        contrasts: list[float] = []
        for candidate in CANDIDATES:
            candidate_delta = float(candidates[candidates["galaxy"] == candidate]["delta_profile_objective_plamb_minus_rar"].iloc[0])
            block_controls = controls[controls["matched_candidate"] == candidate]
            contrasts.append(candidate_delta - float(block_controls["delta_profile_objective_plamb_minus_rar"].mean()))
        positive = sum(value > 0.0 for value in contrasts)
        sign_p = float(binomtest(positive, len(contrasts), 0.5, alternative="greater").pvalue)
        summary_rows.append(
            {
                "scenario": scenario,
                "candidate_mean_delta": float(candidates["delta_profile_objective_plamb_minus_rar"].mean()),
                "candidate_median_delta": float(candidates["delta_profile_objective_plamb_minus_rar"].median()),
                "control_mean_delta": float(controls["delta_profile_objective_plamb_minus_rar"].mean()),
                "control_median_delta": float(controls["delta_profile_objective_plamb_minus_rar"].median()),
                "candidate_minus_control_mean": float(candidates["delta_profile_objective_plamb_minus_rar"].mean() - controls["delta_profile_objective_plamb_minus_rar"].mean()),
                "mean_matched_block_contrast": float(np.mean(contrasts)),
                "median_matched_block_contrast": float(np.median(contrasts)),
                "positive_candidate_blocks": positive,
                "descriptive_one_sided_sign_p": sign_p,
                "candidates_plamb_worse": int(np.sum(candidates["delta_profile_objective_plamb_minus_rar"] > 0.0)),
                "controls_plamb_worse": int(np.sum(controls["delta_profile_objective_plamb_minus_rar"] > 0.0)),
                "optimiser_failures": int(np.sum(~selected["both_optimisers_succeeded"])),
            }
        )
    stability_rows: list[dict[str, Any]] = []
    candidate_frame = frame[frame["group"] == "candidate"]
    for candidate in CANDIDATES:
        selected = candidate_frame[candidate_frame["galaxy"] == candidate]
        combined = selected[selected["scenario"] == "combined_conventional"]
        combined_delta = float(combined["delta_profile_objective_plamb_minus_rar"].iloc[0]) if len(combined) else float("nan")
        fraction = float(np.mean(selected["delta_profile_objective_plamb_minus_rar"] > 0.0))
        if fraction >= 0.80 and combined_delta > 0.0:
            status = "persistent_conventional_tension"
        elif fraction <= 0.50 or combined_delta <= 0.0:
            status = "erased_or_reversed"
        else:
            status = "mixed_sensitivity"
        stability_rows.append(
            {
                "galaxy": candidate,
                "scenarios": len(selected),
                "fraction_scenarios_plamb_worse": fraction,
                "minimum_delta_profile_objective": float(selected["delta_profile_objective_plamb_minus_rar"].min()),
                "maximum_delta_profile_objective": float(selected["delta_profile_objective_plamb_minus_rar"].max()),
                "baseline_delta_profile_objective": float(selected[selected["scenario"] == "baseline"]["delta_profile_objective_plamb_minus_rar"].iloc[0]),
                "combined_delta_profile_objective": combined_delta,
                "robustness_status": status,
            }
        )
    return summary_rows, stability_rows


def make_plots(delta_rows: list[dict[str, Any]], summary_rows: list[dict[str, Any]], out_dir: Path) -> list[Path]:
    frame = pd.DataFrame(delta_rows)
    scenario_order = [scenario.name for scenario in SCENARIOS if scenario.name in set(frame["scenario"])]
    candidate_frame = frame[frame["group"] == "candidate"]
    matrix = candidate_frame.pivot(index="galaxy", columns="scenario", values="delta_profile_objective_plamb_minus_rar")
    matrix = matrix.reindex(index=list(CANDIDATES), columns=scenario_order)
    limit = max(float(np.nanpercentile(np.abs(matrix.to_numpy()), 95)), 1.0)
    fig, ax = plt.subplots(figsize=(13.0, 5.4))
    image = ax.imshow(matrix.to_numpy(), aspect="auto", cmap="coolwarm", vmin=-limit, vmax=limit)
    ax.set_xticks(np.arange(len(scenario_order)), [name.replace("_", "\n") for name in scenario_order], fontsize=7)
    ax.set_yticks(np.arange(len(CANDIDATES)), list(CANDIDATES))
    ax.set_title("Candidate profile objective: PLAMB minus RAR")
    ax.set_xlabel("Locked conventional-systematics scenario")
    ax.set_ylabel("Development-set galaxy")
    colourbar = fig.colorbar(image, ax=ax, pad=0.02)
    colourbar.set_label("Positive means PLAMB worse")
    fig.tight_layout()
    heatmap = out_dir / "sparc_am_candidate_systematics_heatmap.png"
    fig.savefig(heatmap, dpi=180)
    plt.close(fig)

    summary = pd.DataFrame(summary_rows).set_index("scenario").reindex(scenario_order)
    fig, ax = plt.subplots(figsize=(12.0, 5.2))
    x = np.arange(len(summary))
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.bar(x, summary["mean_matched_block_contrast"], color="#35618f")
    ax.set_xticks(x, [name.replace("_", "\n") for name in scenario_order], fontsize=7)
    ax.set_ylabel("Mean candidate minus matched-control Delta objective")
    ax.set_title("Matched-block conventional-systematics contrast")
    fig.tight_layout()
    contrast_plot = out_dir / "sparc_am_matched_block_contrast.png"
    fig.savefig(contrast_plot, dpi=180)
    plt.close(fig)
    return [heatmap, contrast_plot]


def write_report(
    path: Path,
    config: dict[str, Any],
    balance_rows: list[dict[str, Any]],
    summary_rows: list[dict[str, Any]],
    stability_rows: list[dict[str, Any]],
) -> None:
    persistent = [row["galaxy"] for row in stability_rows if row["robustness_status"] == "persistent_conventional_tension"]
    erased = [row["galaxy"] for row in stability_rows if row["robustness_status"] == "erased_or_reversed"]
    max_smd = max(abs(float(row["standardised_mean_difference"])) for row in balance_rows)
    combined = next(row for row in summary_rows if row["scenario"] == "combined_conventional")
    lines = [
        "# SPARC Antimatter Identifiability and Conventional-Systematics Report",
        "",
        f"Date: {RUN_DATE}",
        f"Completed: `{utc_now()}`",
        "",
        "## Bottom Line",
        "",
        "The current RAR/PLAMB SPARC likelihood is matter/antimatter non-identifiable: it has no particle-sign or charge-conjugation branch. Rotation-curve persistence can therefore nominate an anomaly for further study, but cannot identify antimatter.",
        "",
        f"After the locked conventional-systematics profile, `{len(persistent)}` of six development galaxies meet the preregistered persistent-tension rule: {', '.join(persistent) if persistent else 'none'}. `{len(erased)}` are erased or reversed: {', '.join(erased) if erased else 'none'}.",
        "",
        f"In the combined-conventional scenario, the mean candidate-minus-matched-control Delta objective is `{float(combined['mean_matched_block_contrast']):.6g}`. Positive values mean that the candidate set remains relatively more difficult for PLAMB than its matched controls.",
        "",
        "## Gate 1: FR Antimatter Identifiability",
        "",
        "**NOT IDENTIFIABLE.** The supplied FR/SPARC implementation does not specify a matter-to-antimatter transformation for gbar, g0, kappa, p or an additional field. An arbitrary fitted sign is prohibited because it would be an empirical residual absorber rather than an FR prediction.",
        "",
        "## Outcome-Blind Matched Controls",
        "",
        f"The maximum absolute robust-standardised mean difference across the ten locked matching variables is `{max_smd:.4g}`. Matching used no RAR or PLAMB residual.",
        "",
        markdown_table(balance_rows, ["feature", "candidate_median_original", "control_median_original", "standardised_mean_difference"]),
        "",
        "## Conventional-Systematics Matrix",
        "",
        "Positive Delta objective means the fixed-global PLAMB branch profiles worse than fixed-global RAR. The sign-test column is descriptive only because the six candidates were selected previously from their residuals.",
        "",
        markdown_table(
            summary_rows,
            [
                "scenario",
                "candidate_mean_delta",
                "control_mean_delta",
                "mean_matched_block_contrast",
                "positive_candidate_blocks",
                "descriptive_one_sided_sign_p",
                "optimiser_failures",
            ],
        ),
        "",
        "## Candidate Stability",
        "",
        markdown_table(
            stability_rows,
            [
                "galaxy",
                "fraction_scenarios_plamb_worse",
                "baseline_delta_profile_objective",
                "combined_delta_profile_objective",
                "robustness_status",
            ],
        ),
        "",
        "## Interpretation",
        "",
        "A persistent flag means only that this particular nuisance matrix did not remove the candidate's relative PLAMB tension. It is not evidence for antimatter. A reversed or unstable flag supports an ordinary calibration, geometry or baryonic-decomposition explanation.",
        "",
        "The candidate profiles hold the July 14 all-Q2 global parameters fixed. This provides a transparent and symmetric local challenge but is not a refitted hierarchical posterior. Publication-grade follow-up would require a newly derived FR antimatter transformation, a converged hierarchical nuisance model, and validation on the untouched replication galaxies.",
        "",
        "## Locked Claim Boundary",
        "",
        "- Do not label any SPARC galaxy as matter or antimatter from this analysis.",
        "- Do not treat the descriptive sign p-value as a discovery probability.",
        "- Do not inspect the reserved replication set until an FR transformation and decision rule are frozen.",
        "- External gamma-ray analysis is an independent contact-annihilation test, not confirmation of a rotation-curve label by itself.",
        "",
        "## Configuration",
        "",
        "```json",
        json.dumps(config, indent=2),
        "```",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def build_config(scenarios: tuple[Scenario, ...], matching: dict[str, Any]) -> dict[str, Any]:
    return {
        "date": RUN_DATE,
        "written_before_outcomes_utc": utc_now(),
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "sample": str(SAMPLE),
        "sample_sha256": sha256(SAMPLE),
        "points": str(POINTS),
        "points_sha256": sha256(POINTS),
        "map_summary": str(MAP_SUMMARY),
        "map_summary_sha256": sha256(MAP_SUMMARY),
        "candidates": list(CANDIDATES),
        "reserved_replication": list(RESERVED_REPLICATION),
        "models": list(MODELS),
        "global_parameter_treatment": "fixed at July 14 all-Q2 hierarchical MAP values",
        "matching": matching,
        "scenarios": [asdict(scenario) for scenario in scenarios],
        "persistence_gate": "PLAMB profile objective worse in >=80% scenarios and in combined_conventional",
        "identifiability_gate": "current likelihood must contain an independently derived FR charge-conjugation transformation; current result is NOT IDENTIFIABLE",
        "selection_inference": "development-set descriptive only; no discovery p-value",
        "optimisation": "deterministic multi-start L-BFGS-B: centre, per-parameter quartile starts and six hash-seeded joint starts; minimum objective retained",
    }


def preregister(out_dir: Path, scenarios: tuple[Scenario, ...]) -> None:
    out_dir.mkdir(parents=True, exist_ok=False)
    sample = pd.read_csv(SAMPLE)
    points = pd.read_csv(POINTS)
    catalogue = transformed_catalogue(sample, points)
    assignments, balance_rows, matching = select_controls(catalogue)
    matrix = identifiability_matrix()
    config = build_config(scenarios, matching)
    write_csv(out_dir / "sparc_am_candidate_inventory.csv", candidate_inventory(catalogue))
    write_csv(out_dir / "sparc_am_matched_controls.csv", assignments)
    write_csv(out_dir / "sparc_am_matching_balance.csv", balance_rows)
    write_csv(out_dir / "sparc_am_identifiability_matrix.csv", matrix)
    (out_dir / "sparc_am_systematics_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    write_identifiability_document(out_dir / "sparc_am_identifiability_preregistration.md", matrix)
    write_preregistration(out_dir / "sparc_am_systematics_preregistration.md", config, assignments)
    print(f"Preregistered before outcomes: {out_dir}", flush=True)


def execute(out_dir: Path, export_dir: Path, scenarios: tuple[Scenario, ...], copy_to_outputs: bool) -> None:
    required = [
        out_dir / "sparc_am_matched_controls.csv",
        out_dir / "sparc_am_matching_balance.csv",
        out_dir / "sparc_am_systematics_config.json",
        out_dir / "sparc_am_identifiability_preregistration.md",
        out_dir / "sparc_am_systematics_preregistration.md",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Run --preregister-only first; missing: {missing}")
    config = json.loads((out_dir / "sparc_am_systematics_config.json").read_text(encoding="utf-8"))
    locked_scenarios = tuple(Scenario(**item) for item in config["scenarios"])
    if [asdict(item) for item in locked_scenarios] != [asdict(item) for item in scenarios]:
        raise ValueError("Runtime scenario matrix differs from the preregistered matrix")
    assignments = pd.read_csv(out_dir / "sparc_am_matched_controls.csv").to_dict("records")
    balance_rows = pd.read_csv(out_dir / "sparc_am_matching_balance.csv").to_dict("records")
    profile_rows = profile_matrix(scenarios, assignments)
    delta_rows = paired_deltas(profile_rows, assignments)
    summary_rows, stability_rows = summarise_deltas(delta_rows)
    profile_path = out_dir / "sparc_am_systematics_profile_fits.csv"
    delta_path = out_dir / "sparc_am_systematics_galaxy_deltas.csv"
    summary_path = out_dir / "sparc_am_systematics_summary.csv"
    stability_path = out_dir / "sparc_am_candidate_stability.csv"
    report_path = out_dir / "sparc_am_identifiability_systematics_report.md"
    write_csv(profile_path, profile_rows)
    write_csv(delta_path, delta_rows)
    write_csv(summary_path, summary_rows)
    write_csv(stability_path, stability_rows)
    plot_paths = make_plots(delta_rows, summary_rows, out_dir)
    write_report(report_path, config, balance_rows, summary_rows, stability_rows)
    result_files = sorted(path for path in out_dir.iterdir() if path.is_file() and path.name != "manifest.csv")
    write_csv(
        out_dir / "manifest.csv",
        [{"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)} for path in result_files],
    )
    result_files.append(out_dir / "manifest.csv")
    export_dir.mkdir(parents=True, exist_ok=True)
    for path in result_files:
        shutil.copy2(path, export_dir / path.name)
    if copy_to_outputs:
        OUTPUTS.mkdir(parents=True, exist_ok=True)
        for path in [report_path, summary_path, stability_path, *plot_paths]:
            shutil.copy2(path, OUTPUTS / path.name)
    print(f"Saved report: {report_path}", flush=True)
    print(f"Exported reproducibility bundle: {export_dir}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="SPARC antimatter identifiability and conventional-systematics screen.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--export-dir", type=Path, default=DEFAULT_EXPORT)
    parser.add_argument("--preregister-only", action="store_true")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--quick", action="store_true", help="Development smoke run using baseline and combined scenarios only.")
    parser.add_argument("--copy-to-outputs", action=argparse.BooleanOptionalAction, default=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.preregister_only == args.execute:
        raise ValueError("Choose exactly one of --preregister-only or --execute")
    scenarios = (SCENARIOS[0], SCENARIOS[-1]) if args.quick else SCENARIOS
    if args.preregister_only:
        preregister(args.out_dir, scenarios)
    else:
        execute(args.out_dir, args.export_dir, scenarios, args.copy_to_outputs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
