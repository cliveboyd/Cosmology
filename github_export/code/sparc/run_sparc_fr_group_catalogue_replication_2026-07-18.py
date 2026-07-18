#!/usr/bin/env python3
r"""Preregistered SPARC FR replication using a published galaxy-group catalogue.

The acquire and preregister phases do not read outcome values. The analyse phase
uses the previously frozen non-reserved PLAMB-minus-RAR development outcomes.
The five reserved galaxies are never requested, matched, profiled or opened.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[3]
RUN_DATE = "2026-07-18"
SEED = 2026071802

SOURCE_DIR = (
    ROOT
    / "external_datasets"
    / "current_cosmology_sources"
    / "Kourkchi_Tully_2017_Groups"
)
SOURCE_README = SOURCE_DIR / "ReadMe"
GROUP_TABLE = SOURCE_DIR / "table2.dat.gz"
MEMBER_TABLE = SOURCE_DIR / "table3.dat.gz"
SAMPLE = ROOT / "external_datasets" / "current_cosmology_sources" / "SPARC" / "sparc_galaxy_sample.csv"
POINTS = ROOT / "external_datasets" / "current_cosmology_sources" / "SPARC" / "sparc_rotation_points.csv"
POSITION_CACHE = (
    ROOT
    / "external_datasets"
    / "current_cosmology_sources"
    / "SPARC"
    / "sparc_ned_positions_environment_nonreserved_20260718.csv"
)
PARENT_PROGRAMME = (
    ROOT
    / "github_export"
    / "code"
    / "sparc"
    / "run_sparc_fr_environment_asymmetry_2026-07-18.py"
)
OUTCOME_TABLE = (
    ROOT
    / "github_export"
    / "results"
    / "2026-07-18"
    / "fr_environment_asymmetry"
    / "sparc_fr_environment_development_outcomes.csv"
)
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "sparc_fr_group_catalogue_replication_20260718"
DEFAULT_EXPORT = (
    ROOT
    / "github_export"
    / "results"
    / "2026-07-18"
    / "fr_group_catalogue_replication"
)

CATALOGUE_URL = "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/843/16/"
CATALOGUE_QUERY = "https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJ/843/16"
PAPER_URL = "https://doi.org/10.3847/1538-4357/aa76db"

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

MATCH_COORDINATE_ARCMIN = 1.0
MATCH_NAME_ARCMIN = 5.0
MATCH_DV_KMS = 400.0
NAME_MATCH_DV_KMS = 500.0
SHIFTED_DV_KMS = 500.0
PERMUTATIONS = 20000
OUTER_FOLDS = 10
RIDGE_ALPHA = 1.0

HOST_FEATURES = (
    "T",
    "log_L36",
    "log_SBeff",
    "log_MHI",
    "log_Vflat",
    "Inc_deg",
    "frac_distance_error",
    "bulge_proxy",
    "log_n_points",
    "D_Mpc",
    "abs_gal_b_deg",
)

SAMPLE_DEFINITIONS = (
    {
        "name": "primary_10_40_mpc",
        "distance_min_mpc": 10.0,
        "distance_max_mpc": 40.0,
        "abs_galactic_latitude_min_deg": 10.0,
        "min_log_L36": None,
    },
    {
        "name": "near_5_40_mpc",
        "distance_min_mpc": 5.0,
        "distance_max_mpc": 40.0,
        "abs_galactic_latitude_min_deg": 10.0,
        "min_log_L36": None,
    },
    {
        "name": "core_10_35_mpc",
        "distance_min_mpc": 10.0,
        "distance_max_mpc": 35.0,
        "abs_galactic_latitude_min_deg": 10.0,
        "min_log_L36": None,
    },
    {
        "name": "latitude_strict_10_40_mpc",
        "distance_min_mpc": 10.0,
        "distance_max_mpc": 40.0,
        "abs_galactic_latitude_min_deg": 20.0,
        "min_log_L36": None,
    },
    {
        "name": "luminous_10_40_mpc",
        "distance_min_mpc": 10.0,
        "distance_max_mpc": 40.0,
        "abs_galactic_latitude_min_deg": 10.0,
        "min_log_L36": 0.5,
    },
)
PRIMARY_SAMPLE = SAMPLE_DEFINITIONS[0]["name"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_parent() -> Any:
    spec = importlib.util.spec_from_file_location("locked_fr_environment_parent", PARENT_PROGRAMME)
    if spec is None or spec.loader is None:
        raise ImportError(PARENT_PROGRAMME)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def normalise_name(value: Any) -> str:
    name = re.sub(r"[^A-Z0-9]", "", str(value).upper())
    match = re.match(r"^(UGC|NGC|IC|DDO|ESO)0+([0-9].*)$", name)
    return f"{match.group(1)}{match.group(2)}" if match else name


def read_groups() -> pd.DataFrame:
    columns = [
        "PGC1",
        "PGC1_association",
        "Nm",
        "GLON",
        "GLAT",
        "SGLON",
        "SGLAT",
        "Ksmag",
        "logK",
        "HRV",
        "VLS",
        "Nd",
        "Dist",
        "e_Dist",
        "sigmaL",
        "sigmaV",
        "R2t",
        "Rg",
        "logMK",
        "logMd",
    ]
    colspecs = [
        (0, 7),
        (8, 15),
        (16, 19),
        (20, 28),
        (29, 37),
        (38, 46),
        (47, 55),
        (56, 61),
        (62, 67),
        (68, 72),
        (73, 77),
        (78, 81),
        (82, 87),
        (88, 90),
        (91, 94),
        (95, 98),
        (99, 104),
        (105, 110),
        (111, 117),
        (118, 124),
    ]
    frame = pd.read_fwf(
        GROUP_TABLE,
        compression="gzip",
        colspecs=colspecs,
        names=columns,
    )
    for column in columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def read_members() -> pd.DataFrame:
    columns = [
        "PGC",
        "Name",
        "RAdeg",
        "DEdeg",
        "GLON",
        "GLAT",
        "SGLON",
        "SGLAT",
        "TT",
        "Bmag",
        "Ksmag",
        "logK",
        "HRV",
        "VLS",
        "Dist",
        "e_Dist",
        "PGC1",
    ]
    colspecs = [
        (0, 7),
        (8, 36),
        (37, 45),
        (46, 54),
        (55, 63),
        (64, 72),
        (73, 81),
        (82, 90),
        (91, 95),
        (96, 101),
        (102, 107),
        (108, 113),
        (114, 118),
        (119, 123),
        (124, 130),
        (131, 133),
        (134, 141),
    ]
    frame = pd.read_fwf(
        MEMBER_TABLE,
        compression="gzip",
        colspecs=colspecs,
        names=columns,
    )
    for column in columns:
        if column != "Name":
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["normalised_name"] = frame["Name"].map(normalise_name)
    return frame


def shifted_group_proxy(
    ra_deg: float,
    dec_deg: float,
    cz_kms: float,
    groups: pd.DataFrame,
    group_coordinates: SkyCoord,
) -> tuple[float, float, float, float]:
    valid = np.isfinite(groups["HRV"]) & (np.abs(groups["HRV"] - cz_kms) <= SHIFTED_DV_KMS)
    if not np.any(valid):
        return float("nan"), float("nan"), float("nan"), float("nan")
    target = SkyCoord(ra=((ra_deg + 90.0) % 360.0) * u.deg, dec=dec_deg * u.deg)
    valid_indices = np.flatnonzero(valid.to_numpy())
    separations = target.separation(group_coordinates[valid_indices]).deg
    selected_index = valid_indices[int(np.argmin(separations))]
    row = groups.iloc[selected_index]
    richness = float(row["Nm"])
    return (
        math.log(max(richness, 1.0)),
        float(row["logMK"]),
        float(separations[np.argmin(separations)]),
        float(row["PGC1"]),
    )


def match_catalogue() -> tuple[pd.DataFrame, dict[str, Any]]:
    sample = pd.read_csv(SAMPLE)
    positions = pd.read_csv(POSITION_CACHE)
    groups = read_groups()
    members = read_members()
    targets = sample.loc[
        (pd.to_numeric(sample["Q"], errors="coerce") <= 2)
        & (~sample["Galaxy"].astype(str).isin(RESERVED_REPLICATION))
    ].copy()
    if set(targets["Galaxy"].astype(str)) & set(RESERVED_REPLICATION):
        raise AssertionError("Reserved galaxy entered target list")
    positions = positions[positions["Galaxy"].isin(targets["Galaxy"])].copy()
    merged_targets = targets.merge(positions, on="Galaxy", how="left", suffixes=("", "_ned"))
    member_coordinates = SkyCoord(
        ra=members["RAdeg"].to_numpy(dtype=float) * u.deg,
        dec=members["DEdeg"].to_numpy(dtype=float) * u.deg,
    )
    group_coordinates = SkyCoord(
        l=groups["GLON"].to_numpy(dtype=float) * u.deg,
        b=groups["GLAT"].to_numpy(dtype=float) * u.deg,
        frame="galactic",
    ).icrs
    group_lookup = groups.set_index("PGC1")
    rows: list[dict[str, Any]] = []
    for _, target_row in merged_targets.sort_values("Galaxy").iterrows():
        galaxy = str(target_row["Galaxy"])
        resolved = str(target_row.get("resolved", "")).lower() in {"true", "1"}
        ra = float(target_row.get("ra_deg", np.nan))
        dec = float(target_row.get("dec_deg", np.nan))
        cz = float(target_row.get("cz_kms", np.nan))
        galactic_latitude = float("nan")
        if resolved and math.isfinite(ra) and math.isfinite(dec):
            galactic_latitude = float(
                SkyCoord(ra=ra * u.deg, dec=dec * u.deg).galactic.b.deg
            )
        base: dict[str, Any] = {
            "Galaxy": galaxy,
            "Q": int(target_row["Q"]),
            "D_Mpc": float(target_row["D_Mpc"]),
            "resolved": resolved,
            "ra_deg": ra,
            "dec_deg": dec,
            "central_cz_kms": cz,
            "abs_gal_b_deg": abs(galactic_latitude),
            "matched": False,
            "match_rule": "unmatched",
            "match_ambiguous": False,
        }
        if not resolved or not all(math.isfinite(value) for value in [ra, dec, cz]):
            rows.append(base)
            continue
        target_coordinate = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
        separations = target_coordinate.separation(member_coordinates).arcmin
        velocity_offsets = np.abs(members["HRV"].to_numpy(dtype=float) - cz)
        nearest = int(np.nanargmin(separations))
        normalised_targets = {
            normalise_name(galaxy),
            normalise_name(target_row.get("preferred_name", "")),
        }
        name_match = members["normalised_name"].isin(normalised_targets).to_numpy()
        coordinate_match = (separations <= MATCH_COORDINATE_ARCMIN) & (
            velocity_offsets <= MATCH_DV_KMS
        )
        named_match = name_match & (separations <= MATCH_NAME_ARCMIN) & (
            velocity_offsets <= NAME_MATCH_DV_KMS
        )
        accepted_indices = np.flatnonzero(coordinate_match | named_match)
        selected_index: int | None = None
        rule = "unmatched"
        if np.any(named_match):
            named_indices = np.flatnonzero(named_match)
            selected_index = int(named_indices[np.argmin(separations[named_indices])])
            rule = "name_coordinate_velocity"
        elif np.any(coordinate_match):
            coordinate_indices = np.flatnonzero(coordinate_match)
            selected_index = int(coordinate_indices[np.argmin(separations[coordinate_indices])])
            rule = "coordinate_velocity"
        candidate_index = selected_index if selected_index is not None else nearest
        member = members.iloc[candidate_index]
        base.update(
            {
                "matched": selected_index is not None,
                "match_rule": rule,
                "match_ambiguous": bool(len(accepted_indices) > 1 and not np.any(named_match)),
                "catalogue_member_name": str(member["Name"]),
                "member_pgc": int(member["PGC"]) if math.isfinite(float(member["PGC"])) else np.nan,
                "group_pgc1": int(member["PGC1"]) if math.isfinite(float(member["PGC1"])) else np.nan,
                "match_separation_arcmin": float(separations[candidate_index]),
                "match_velocity_offset_kms": float(velocity_offsets[candidate_index]),
                "catalogue_member_hrv_kms": float(member["HRV"]),
                "catalogue_member_vls_kms": float(member["VLS"]),
            }
        )
        if selected_index is not None and float(member["PGC1"]) in group_lookup.index:
            group = group_lookup.loc[float(member["PGC1"])]
            if isinstance(group, pd.DataFrame):
                group = group.iloc[0]
            richness = float(group["Nm"])
            base.update(
                {
                    "group_richness": richness,
                    "group_companions": richness - 1.0,
                    "group_richness_score": math.log(max(richness, 1.0)),
                    "grouped_indicator": float(richness >= 2.0),
                    "group_log_mass_luminosity": float(group["logMK"]),
                    "group_log_mass_dynamical": float(group["logMd"]),
                    "group_distance_mpc": float(group["Dist"]),
                    "group_r2t_mpc": float(group["R2t"]),
                }
            )
        shifted_richness, shifted_mass, shifted_sep, shifted_pgc = shifted_group_proxy(
            ra,
            dec,
            cz,
            groups,
            group_coordinates,
        )
        base.update(
            {
                "shifted_group_richness_score": shifted_richness,
                "shifted_group_log_mass_luminosity": shifted_mass,
                "shifted_nearest_group_separation_deg": shifted_sep,
                "shifted_group_pgc1": shifted_pgc,
            }
        )
        rows.append(base)
    frame = pd.DataFrame(rows).sort_values("Galaxy").reset_index(drop=True)
    if set(frame["Galaxy"]) & set(RESERVED_REPLICATION):
        raise AssertionError("Reserved galaxy entered group covariates")
    accepted = frame[frame["matched"].astype(bool)]
    audit = {
        "date": RUN_DATE,
        "created_utc": utc_now(),
        "outcomes_read": False,
        "catalogue": "Kourkchi and Tully 2017, J/ApJ/843/16",
        "source_url": CATALOGUE_URL,
        "member_rows": int(len(members)),
        "group_rows": int(len(groups)),
        "target_rows_nonreserved_q2": int(len(frame)),
        "accepted_matches": int(len(accepted)),
        "unmatched": int(np.sum(~frame["matched"].astype(bool))),
        "ambiguous_matches": int(np.sum(frame["match_ambiguous"].astype(bool))),
        "median_match_separation_arcmin": float(accepted["match_separation_arcmin"].median()),
        "max_match_separation_arcmin": float(accepted["match_separation_arcmin"].max()),
        "median_velocity_offset_kms": float(accepted["match_velocity_offset_kms"].median()),
        "reserved_targets": int(np.sum(frame["Galaxy"].isin(RESERVED_REPLICATION))),
        "file_hashes": {
            "ReadMe": sha256(SOURCE_README),
            "table2.dat.gz": sha256(GROUP_TABLE),
            "table3.dat.gz": sha256(MEMBER_TABLE),
            "NED_position_cache": sha256(POSITION_CACHE),
        },
    }
    return frame, audit


def sample_mask(frame: pd.DataFrame, definition: dict[str, Any]) -> pd.Series:
    mask = (
        frame["matched"].astype(bool)
        & (~frame["match_ambiguous"].astype(bool))
        & (pd.to_numeric(frame["Q"], errors="coerce") <= 2)
        & (pd.to_numeric(frame["D_Mpc"], errors="coerce") >= definition["distance_min_mpc"])
        & (pd.to_numeric(frame["D_Mpc"], errors="coerce") <= definition["distance_max_mpc"])
        & (
            pd.to_numeric(frame["abs_gal_b_deg"], errors="coerce")
            >= definition["abs_galactic_latitude_min_deg"]
        )
        & (~frame["Galaxy"].isin(CANDIDATES))
        & (~frame["Galaxy"].isin(RESERVED_REPLICATION))
    )
    if definition["min_log_L36"] is not None:
        mask &= pd.to_numeric(frame["log_L36"], errors="coerce") >= float(
            definition["min_log_L36"]
        )
    return mask


def preregistration_config(covariates: pd.DataFrame, programme: Path) -> dict[str, Any]:
    parent = load_parent()
    host = parent.host_catalogue()
    coverage = covariates.merge(host[["Galaxy", "log_L36"]], on="Galaxy", how="left")
    counts = {
        definition["name"]: int(np.sum(sample_mask(coverage, definition)))
        for definition in SAMPLE_DEFINITIONS
    }
    return {
        "date": RUN_DATE,
        "written_utc": utc_now(),
        "programme": str(programme),
        "programme_sha256": sha256(programme),
        "outcome_values_read_in_this_programme_before_lock": False,
        "study_status": (
            "Independent-predictor replication on previously generated development outcomes; "
            "not a held-out-outcome replication"
        ),
        "known_prior_result": (
            "The earlier 2MRS composite test had rho=-0.307254, p=0.0117994 and "
            "2.398 per cent cross-validated improvement, and failed its locked gate."
        ),
        "catalogue": {
            "name": "Kourkchi and Tully 2017 Galaxy Groups within 3500 km/s",
            "vizier": CATALOGUE_QUERY,
            "paper": PAPER_URL,
            "source_directory": CATALOGUE_URL,
            "construction_independence": (
                "Published group-finding construction independent of this project's 2MRS neighbour-count algorithm"
            ),
            "non_independence_limit": (
                "The catalogue uses overlapping nearby-galaxy and Ks-band information, including 2MRS completeness; "
                "it is not an independent raw survey."
            ),
        },
        "input_hashes": {
            "sample": sha256(SAMPLE),
            "points": sha256(POINTS),
            "positions": sha256(POSITION_CACHE),
            "catalogue_readme": sha256(SOURCE_README),
            "catalogue_groups": sha256(GROUP_TABLE),
            "catalogue_members": sha256(MEMBER_TABLE),
            "covariates": sha256(DEFAULT_OUT / "sparc_fr_group_catalogue_covariates_nonreserved.csv"),
            "parent_programme": sha256(PARENT_PROGRAMME),
            "frozen_outcomes": sha256(OUTCOME_TABLE),
        },
        "exclusions": {
            "previously_outcome_selected": list(CANDIDATES),
            "reserved_replication": list(RESERVED_REPLICATION),
        },
        "matching": {
            "coordinate_arcmin_max": MATCH_COORDINATE_ARCMIN,
            "coordinate_velocity_kms_max": MATCH_DV_KMS,
            "name_coordinate_arcmin_max": MATCH_NAME_ARCMIN,
            "name_velocity_kms_max": NAME_MATCH_DV_KMS,
            "ambiguous_matches_excluded": True,
        },
        "samples": {
            "definitions": list(SAMPLE_DEFINITIONS),
            "outcome_blind_counts": counts,
        },
        "predictors": {
            "primary": "group_richness_score = ln(N_group), including N_group=1 singles",
            "secondary_not_in_gate": [
                "group_log_mass_luminosity",
                "grouped_indicator",
            ],
            "negative_control": (
                "richness of the nearest published group after RA is shifted by +90 degrees, "
                "restricted to |Delta v_hel|<=500 km/s"
            ),
        },
        "outcomes": {
            "primary": "frozen combined_conventional signed_log1p PLAMB-minus-RAR profile contrast per point",
            "secondary": "same frozen contrast under the baseline nuisance scenario",
            "no_new_rotation_curve_profiling": True,
        },
        "host_controls": list(HOST_FEATURES) + ["f_D one-hot indicators"],
        "statistics": {
            "partial_test": "Spearman correlation of cross-fitted host-control residuals",
            "ridge_alpha": RIDGE_ALPHA,
            "outer_folds": OUTER_FOLDS,
            "fold_seed": SEED,
            "permutations": PERMUTATIONS,
            "predictive_metric": "cross-fitted RMSE change after adding group richness",
        },
        "development_gate": {
            "all_required": True,
            "primary_permutation_p_max": 0.01,
            "primary_cv_rmse_fractional_improvement_min": 0.05,
            "coefficient_sign_fraction_min": 0.80,
            "same_partial_correlation_sign_in": [
                "near_5_40_mpc",
                "core_10_35_mpc",
                "latitude_strict_10_40_mpc",
                "luminous_10_40_mpc",
                "primary_baseline_outcome",
            ],
            "negative_control": "absolute shifted-sky rho must be smaller than absolute actual rho",
        },
        "replication_policy": (
            "Reserved outcomes remain sealed regardless of this result. A separate held-out replication "
            "preregistration and explicit unsealing decision are required."
        ),
        "claim_boundary": [
            "A pass supports a conventional group-environment predictor, not antimatter.",
            "No galaxy matter/antimatter labels may be fitted or reported.",
            "A fail ends this catalogue branch without unsealing the reserved galaxies.",
        ],
    }


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in frame.iterrows():
        rendered: list[str] = []
        for column in columns:
            value = row[column]
            if isinstance(value, (float, np.floating)):
                rendered.append(f"{float(value):.6g}" if math.isfinite(float(value)) else "")
            else:
                rendered.append(str(value))
        lines.append("| " + " | ".join(rendered) + " |")
    return "\n".join(lines)


def write_preregistration(path: Path, config: dict[str, Any]) -> None:
    counts = pd.DataFrame(
        [
            {"sample": name, "N_outcome_blind": count}
            for name, count in config["samples"]["outcome_blind_counts"].items()
        ]
    )
    lines = [
        "# SPARC FR Published-Group-Catalogue Replication Preregistration",
        "",
        f"Date: {RUN_DATE}",
        f"Locked: `{config['written_utc']}`",
        "",
        "## Status",
        "",
        config["study_status"],
        "",
        config["known_prior_result"],
        "",
        "The group covariates and sample coverage were constructed without reading outcome values in this programme. The development outcomes themselves were generated previously, so this is not a new held-out-outcome experiment.",
        "",
        "## Catalogue",
        "",
        f"The predictor comes from [Kourkchi and Tully 2017]({PAPER_URL}); fixed tables are archived by [VizieR]({CATALOGUE_QUERY}). The catalogue lists 15,004 galaxies and group assignments within 3500 km/s.",
        "",
        config["catalogue"]["non_independence_limit"],
        "",
        "## Locked Samples",
        "",
        markdown_table(counts, ["sample", "N_outcome_blind"]),
        "",
        "The primary sample is Q<=2, 10-40 Mpc and |b|>=10 degrees. Six outcome-selected galaxies and all five reserved galaxies are excluded. Accepted coordinates must satisfy the locked coordinate/velocity match; ambiguous matches are excluded.",
        "",
        "## Primary Predictor",
        "",
        "The single primary predictor is `group_richness_score = ln(N_group)`, where catalogue singles have N_group=1 and score zero. It avoids selecting a predictor using the development outcome and is less directly tied to Ks luminosity than catalogue halo mass.",
        "",
        "Luminosity-derived group mass and grouped/single status are secondary diagnostics and cannot pass the gate. A +90-degree RA shifted group-richness assignment is the negative control.",
        "",
        "## Outcome And Controls",
        "",
        "The outcome is the frozen signed-log transformed PLAMB-minus-RAR profile-objective contrast per rotation-curve point under the combined-conventional nuisance scenario. The baseline scenario is a sensitivity control. Host morphology, luminosity, surface brightness, gas mass, velocity, inclination, distance uncertainty, bulge proxy, point count, distance, Galactic latitude and distance-method indicators are controlled by cross-fitting.",
        "",
        "## Locked Gate",
        "",
        "Every condition is required: two-sided permutation p<=0.01, cross-fitted RMSE improvement>=5 per cent, one coefficient sign in at least 8/10 folds, the same partial-correlation sign in all listed distance/latitude/luminosity and baseline-outcome sensitivities, and a smaller absolute shifted-sky rho.",
        "",
        "Reserved outcomes remain sealed even after a pass. A separate held-out replication preregistration is mandatory.",
        "",
        "## Claim Boundary",
        "",
        "A pass would support a conventional group-environment predictor. It would not identify antimatter, estimate a signed matter/antimatter field or authorise matter/antimatter galaxy labels.",
        "",
        f"Programme SHA-256: `{config['programme_sha256']}`",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def host_matrix(frame: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    matrix = pd.DataFrame(index=frame.index)
    for column in HOST_FEATURES:
        values = pd.to_numeric(frame[column], errors="coerce")
        matrix[column] = values.fillna(float(values.median()))
    distance_dummies = pd.get_dummies(
        pd.to_numeric(frame["f_D"], errors="coerce").fillna(-1).astype(int),
        prefix="f_D",
        dtype=float,
    )
    matrix = pd.concat([matrix, distance_dummies], axis=1)
    return matrix.to_numpy(dtype=float), list(matrix.columns)


def permutation_p(x: np.ndarray, y: np.ndarray, label: str) -> tuple[float, float]:
    rho = float(spearmanr(x, y).statistic)
    if not math.isfinite(rho):
        return float("nan"), float("nan")
    xr = pd.Series(x).rank(method="average").to_numpy(dtype=float)
    yr = pd.Series(y).rank(method="average").to_numpy(dtype=float)
    xr -= np.mean(xr)
    yr -= np.mean(yr)
    denominator = math.sqrt(float(np.sum(xr**2) * np.sum(yr**2)))
    observed = float(np.dot(xr, yr) / denominator)
    offset = int(hashlib.sha256(label.encode("utf-8")).hexdigest()[:8], 16)
    rng = np.random.default_rng(SEED + offset)
    exceed = 0
    for _ in range(PERMUTATIONS):
        permuted = rng.permutation(yr)
        if abs(float(np.dot(xr, permuted) / denominator)) >= abs(observed) - 1.0e-15:
            exceed += 1
    return rho, (1.0 + exceed) / (1.0 + PERMUTATIONS)


def cross_fitted_test(
    frame: pd.DataFrame,
    predictor: str,
    label: str,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    y = frame["signed_log1p_delta_per_point"].to_numpy(dtype=float)
    environment = frame[predictor].to_numpy(dtype=float)
    x, feature_names = host_matrix(frame)
    n_splits = min(OUTER_FOLDS, max(3, len(frame) // 5))
    splitter = KFold(n_splits=n_splits, shuffle=True, random_state=SEED)
    baseline_prediction = np.full(len(frame), np.nan)
    augmented_prediction = np.full(len(frame), np.nan)
    environment_prediction = np.full(len(frame), np.nan)
    fold_rows: list[dict[str, Any]] = []
    for fold_index, (train, test) in enumerate(splitter.split(x), start=1):
        scaler = StandardScaler().fit(x[train])
        x_train = scaler.transform(x[train])
        x_test = scaler.transform(x[test])
        baseline = Ridge(alpha=RIDGE_ALPHA).fit(x_train, y[train])
        baseline_prediction[test] = baseline.predict(x_test)
        environment_model = Ridge(alpha=RIDGE_ALPHA).fit(x_train, environment[train])
        environment_prediction[test] = environment_model.predict(x_test)
        augmented_scaler = StandardScaler().fit(
            np.column_stack([x[train], environment[train]])
        )
        augmented = Ridge(alpha=RIDGE_ALPHA).fit(
            augmented_scaler.transform(np.column_stack([x[train], environment[train]])),
            y[train],
        )
        augmented_prediction[test] = augmented.predict(
            augmented_scaler.transform(np.column_stack([x[test], environment[test]]))
        )
        fold_rows.append(
            {
                "test": label,
                "fold": fold_index,
                "N_train": len(train),
                "N_test": len(test),
                "environment_coefficient_standardised": float(augmented.coef_[-1]),
            }
        )
    outcome_residual = y - baseline_prediction
    environment_residual = environment - environment_prediction
    rho, p_value = permutation_p(environment_residual, outcome_residual, label)
    baseline_rmse = math.sqrt(mean_squared_error(y, baseline_prediction))
    augmented_rmse = math.sqrt(mean_squared_error(y, augmented_prediction))
    improvement = (baseline_rmse - augmented_rmse) / baseline_rmse
    coefficients = np.asarray(
        [row["environment_coefficient_standardised"] for row in fold_rows], dtype=float
    )
    sign_fraction = max(float(np.mean(coefficients > 0)), float(np.mean(coefficients < 0)))
    result = {
        "test": label,
        "predictor": predictor,
        "N": len(frame),
        "partial_spearman_rho": rho,
        "permutation_p_two_sided": p_value,
        "baseline_cv_rmse": baseline_rmse,
        "environment_cv_rmse": augmented_rmse,
        "cv_rmse_fractional_improvement": improvement,
        "environment_coefficient_median": float(np.median(coefficients)),
        "environment_coefficient_sign_fraction": sign_fraction,
        "outer_folds_used": n_splits,
        "host_features": ";".join(feature_names),
    }
    residuals = pd.DataFrame(
        {
            "Galaxy": frame["Galaxy"].to_numpy(),
            "test": label,
            "outcome": y,
            "environment_score": environment,
            "outcome_crossfit_residual": outcome_residual,
            "environment_crossfit_residual": environment_residual,
            "baseline_prediction": baseline_prediction,
            "augmented_prediction": augmented_prediction,
        }
    )
    return result, pd.DataFrame(fold_rows), residuals


def build_frame(
    covariates: pd.DataFrame,
    host: pd.DataFrame,
    outcomes: pd.DataFrame,
    definition: dict[str, Any],
    scenario: str,
) -> pd.DataFrame:
    frame = covariates.merge(host, on="Galaxy", how="inner", suffixes=("", "_host"))
    frame = frame.merge(
        outcomes[outcomes["scenario"] == scenario],
        on="Galaxy",
        how="inner",
    )
    frame = frame[sample_mask(frame, definition)].copy()
    frame = frame[frame["rar_success"].astype(bool) & frame["plamb_success"].astype(bool)]
    frame = frame.sort_values("Galaxy").reset_index(drop=True)
    excluded = set(frame["Galaxy"]) & (set(CANDIDATES) | set(RESERVED_REPLICATION))
    if excluded:
        raise AssertionError(f"Excluded galaxies entered analysis: {sorted(excluded)}")
    return frame


def evaluate(
    covariates: pd.DataFrame,
    host: pd.DataFrame,
    outcomes: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    results: list[dict[str, Any]] = []
    folds: list[pd.DataFrame] = []
    residuals: list[pd.DataFrame] = []
    primary_frame = pd.DataFrame()
    for definition in SAMPLE_DEFINITIONS:
        frame = build_frame(covariates, host, outcomes, definition, "combined_conventional")
        label = f"{definition['name']}__combined__group_richness"
        result, fold, residual = cross_fitted_test(frame, "group_richness_score", label)
        result.update(
            {
                "sample": definition["name"],
                "scenario": "combined_conventional",
                "control": "actual",
            }
        )
        results.append(result)
        folds.append(fold)
        residuals.append(residual)
        if definition["name"] == PRIMARY_SAMPLE:
            primary_frame = frame
    baseline_frame = build_frame(
        covariates,
        host,
        outcomes,
        SAMPLE_DEFINITIONS[0],
        "baseline",
    )
    baseline_result, baseline_folds, baseline_residuals = cross_fitted_test(
        baseline_frame,
        "group_richness_score",
        "primary_10_40_mpc__baseline__group_richness",
    )
    baseline_result.update(
        {"sample": PRIMARY_SAMPLE, "scenario": "baseline", "control": "actual"}
    )
    results.append(baseline_result)
    folds.append(baseline_folds)
    residuals.append(baseline_residuals)
    for predictor, control in [
        ("shifted_group_richness_score", "shifted_ra90"),
        ("group_log_mass_luminosity", "secondary_group_mass"),
        ("grouped_indicator", "secondary_grouped_indicator"),
    ]:
        result, fold, residual = cross_fitted_test(
            primary_frame,
            predictor,
            f"primary_10_40_mpc__combined__{control}",
        )
        result.update(
            {"sample": PRIMARY_SAMPLE, "scenario": "combined_conventional", "control": control}
        )
        results.append(result)
        folds.append(fold)
        residuals.append(residual)
    summary = pd.DataFrame(results)
    primary = summary[
        (summary["sample"] == PRIMARY_SAMPLE)
        & (summary["scenario"] == "combined_conventional")
        & (summary["control"] == "actual")
    ].iloc[0]
    shifted = summary[summary["control"] == "shifted_ra90"].iloc[0]
    primary_sign = int(np.sign(float(primary["partial_spearman_rho"])))
    sensitivity_rows = summary[
        (
            summary["sample"].isin(
                [
                    "near_5_40_mpc",
                    "core_10_35_mpc",
                    "latitude_strict_10_40_mpc",
                    "luminous_10_40_mpc",
                ]
            )
            & (summary["control"] == "actual")
        )
        | ((summary["sample"] == PRIMARY_SAMPLE) & (summary["scenario"] == "baseline"))
    ]
    same_sign = bool(
        primary_sign != 0
        and all(
            int(np.sign(float(value))) == primary_sign
            for value in sensitivity_rows["partial_spearman_rho"]
        )
    )
    conditions = {
        "primary_p": bool(float(primary["permutation_p_two_sided"]) <= 0.01),
        "cv_improvement": bool(float(primary["cv_rmse_fractional_improvement"]) >= 0.05),
        "coefficient_sign": bool(
            float(primary["environment_coefficient_sign_fraction"]) >= 0.80
        ),
        "sensitivity_signs": same_sign,
        "negative_control": bool(
            abs(float(shifted["partial_spearman_rho"]))
            < abs(float(primary["partial_spearman_rho"]))
        ),
    }
    gate = {
        "date": RUN_DATE,
        "development_gate_passed": bool(all(conditions.values())),
        "conditions": conditions,
        "replication_unsealed": False,
        "replication_policy": "A separate held-out preregistration is required even after a pass.",
    }
    return (
        summary,
        pd.concat(folds, ignore_index=True),
        pd.concat(residuals, ignore_index=True),
        primary_frame,
        gate,
    )


def make_plot(residuals: pd.DataFrame, path: Path) -> None:
    x = residuals["environment_crossfit_residual"].to_numpy(dtype=float)
    y = residuals["outcome_crossfit_residual"].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.scatter(x, y, s=34, color="#2f6f6d", alpha=0.8, edgecolor="white", linewidth=0.4)
    if len(x) >= 2 and np.std(x) > 0:
        slope, intercept = np.polyfit(x, y, 1)
        grid = np.linspace(float(np.min(x)), float(np.max(x)), 100)
        ax.plot(grid, intercept + slope * grid, color="#b54435", linewidth=2.0)
    ax.axhline(0, color="#555555", linewidth=0.8)
    ax.axvline(0, color="#555555", linewidth=0.8)
    ax.set_xlabel("Published group richness after host controls")
    ax.set_ylabel("PLAMB-RAR outcome after host controls")
    ax.set_title("SPARC FR Published-Group Replication")
    ax.grid(alpha=0.18)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_report(path: Path, summary: pd.DataFrame, gate: dict[str, Any]) -> None:
    primary = summary[
        (summary["sample"] == PRIMARY_SAMPLE)
        & (summary["scenario"] == "combined_conventional")
        & (summary["control"] == "actual")
    ].iloc[0]
    decision = "PASS" if gate["development_gate_passed"] else "FAIL"
    columns = [
        "sample",
        "scenario",
        "control",
        "predictor",
        "N",
        "partial_spearman_rho",
        "permutation_p_two_sided",
        "cv_rmse_fractional_improvement",
        "environment_coefficient_sign_fraction",
    ]
    gate_frame = pd.DataFrame(
        [{"condition": key, "passed": value} for key, value in gate["conditions"].items()]
    )
    lines = [
        "# SPARC FR Published-Group-Catalogue Replication",
        "",
        f"Date: {RUN_DATE}",
        f"Completed: `{utc_now()}`",
        "",
        "## Bottom Line",
        "",
        f"The preregistered development gate is **{decision}**. For the primary published group-richness predictor, partial Spearman rho=`{float(primary['partial_spearman_rho']):.6g}`, two-sided permutation p=`{float(primary['permutation_p_two_sided']):.6g}`, and cross-fitted RMSE change=`{100.0 * float(primary['cv_rmse_fractional_improvement']):.3f}` per cent.",
        "",
        "This is an independent-predictor replication on previously generated development outcomes. It is not an independent held-out-outcome replication. Kourkchi-Tully group construction is independent of this project's neighbour-count algorithm, but its nearby-galaxy and Ks-band inputs overlap 2MRS.",
        "",
        "## Locked Test Matrix",
        "",
        markdown_table(summary, columns),
        "",
        "## Gate",
        "",
        markdown_table(gate_frame, ["condition", "passed"]),
        "",
        "The five reserved galaxies remain sealed. No result in this programme identifies antimatter or authorises matter/antimatter labels.",
        "",
        "## Sources",
        "",
        f"[Kourkchi and Tully 2017]({PAPER_URL}); [VizieR fixed catalogue tables]({CATALOGUE_QUERY}).",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def export_bundle(out_dir: Path, export_dir: Path, programme: Path) -> None:
    export_dir.mkdir(parents=True, exist_ok=True)
    for path in out_dir.iterdir():
        if path.is_file():
            shutil.copy2(path, export_dir / path.name)
    rows: list[dict[str, Any]] = []
    for path in sorted(export_dir.iterdir()):
        if not path.is_file() or path.name == "manifest.csv":
            continue
        rows.append(
            {
                "path": path.name,
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
                "role": "committed analysis product",
                "tracked_in_git": True,
            }
        )
    for path, role, tracked in [
        (programme, "locked replication programme", True),
        (PARENT_PROGRAMME, "locked parent environment programme", True),
        (SOURCE_README, "external VizieR catalogue documentation", False),
        (GROUP_TABLE, "external VizieR group table", False),
        (MEMBER_TABLE, "external VizieR member table", False),
        (POSITION_CACHE, "external NED position cache excluding reserved galaxies", False),
        (OUTCOME_TABLE, "frozen non-reserved development outcomes", True),
    ]:
        rows.append(
            {
                "path": str(path.relative_to(ROOT)).replace("\\", "/") if tracked else str(path),
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
                "role": role,
                "tracked_in_git": tracked,
            }
        )
    pd.DataFrame(rows).to_csv(export_dir / "manifest.csv", index=False)


def acquire_phase(out_dir: Path, export_dir: Path) -> None:
    for path in [SOURCE_README, GROUP_TABLE, MEMBER_TABLE, SAMPLE, POSITION_CACHE]:
        if not path.exists():
            raise FileNotFoundError(path)
    out_dir.mkdir(parents=True, exist_ok=True)
    covariates, audit = match_catalogue()
    covariate_path = out_dir / "sparc_fr_group_catalogue_covariates_nonreserved.csv"
    covariates.to_csv(covariate_path, index=False)
    audit_path = out_dir / "sparc_fr_group_catalogue_acquisition_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    export_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(covariate_path, export_dir / covariate_path.name)
    shutil.copy2(audit_path, export_dir / audit_path.name)
    print(f"Saved group covariates: {covariate_path}")
    print(f"Accepted matches: {audit['accepted_matches']}/{audit['target_rows_nonreserved_q2']}")
    print(f"Ambiguous matches: {audit['ambiguous_matches']}")
    print("Reserved replication targets: 0")


def preregister_phase(out_dir: Path, export_dir: Path, programme: Path) -> None:
    covariate_path = out_dir / "sparc_fr_group_catalogue_covariates_nonreserved.csv"
    if not covariate_path.exists():
        raise FileNotFoundError("Run --phase acquire first")
    covariates = pd.read_csv(covariate_path)
    config = preregistration_config(covariates, programme)
    json_path = out_dir / "sparc_fr_group_catalogue_preregistration.json"
    markdown_path = out_dir / "sparc_fr_group_catalogue_preregistration.md"
    json_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    write_preregistration(markdown_path, config)
    export_dir.mkdir(parents=True, exist_ok=True)
    for path in [
        covariate_path,
        out_dir / "sparc_fr_group_catalogue_acquisition_audit.json",
        json_path,
        markdown_path,
    ]:
        shutil.copy2(path, export_dir / path.name)
    print(f"Saved preregistration: {markdown_path}")
    print(f"Programme SHA-256: {config['programme_sha256']}")
    for name, count in config["samples"]["outcome_blind_counts"].items():
        print(f"{name}: N={count}")
    print("Outcome values read before lock: no")


def analyse_phase(out_dir: Path, export_dir: Path, programme: Path) -> None:
    prereg_path = out_dir / "sparc_fr_group_catalogue_preregistration.json"
    covariate_path = out_dir / "sparc_fr_group_catalogue_covariates_nonreserved.csv"
    if not prereg_path.exists():
        raise FileNotFoundError("Run and commit --phase preregister first")
    config = json.loads(prereg_path.read_text(encoding="utf-8"))
    if sha256(programme) != config["programme_sha256"]:
        raise RuntimeError("Programme hash no longer matches preregistration")
    locked_paths = {
        "sample": SAMPLE,
        "points": POINTS,
        "positions": POSITION_CACHE,
        "catalogue_readme": SOURCE_README,
        "catalogue_groups": GROUP_TABLE,
        "catalogue_members": MEMBER_TABLE,
        "covariates": covariate_path,
        "parent_programme": PARENT_PROGRAMME,
        "frozen_outcomes": OUTCOME_TABLE,
    }
    for key, path in locked_paths.items():
        if sha256(path) != config["input_hashes"][key]:
            raise RuntimeError(f"Input hash changed after preregistration: {key}")
    covariates = pd.read_csv(covariate_path)
    outcomes = pd.read_csv(OUTCOME_TABLE)
    excluded_outcomes = set(outcomes["Galaxy"]) & (
        set(CANDIDATES) | set(RESERVED_REPLICATION)
    )
    if excluded_outcomes:
        raise AssertionError(f"Excluded outcomes present: {sorted(excluded_outcomes)}")
    parent = load_parent()
    host = parent.host_catalogue()
    summary, folds, residuals, primary_frame, gate = evaluate(covariates, host, outcomes)
    summary.to_csv(out_dir / "sparc_fr_group_catalogue_test_summary.csv", index=False)
    folds.to_csv(out_dir / "sparc_fr_group_catalogue_fold_coefficients.csv", index=False)
    residuals.to_csv(out_dir / "sparc_fr_group_catalogue_test_residuals.csv", index=False)
    primary_label = "primary_10_40_mpc__combined__group_richness"
    primary_residuals = residuals[residuals["test"] == primary_label].copy()
    primary_residuals.to_csv(
        out_dir / "sparc_fr_group_catalogue_primary_residuals.csv", index=False
    )
    primary_frame.to_csv(out_dir / "sparc_fr_group_catalogue_primary_frame.csv", index=False)
    (out_dir / "sparc_fr_group_catalogue_gate.json").write_text(
        json.dumps(gate, indent=2) + "\n", encoding="utf-8"
    )
    make_plot(primary_residuals, out_dir / "sparc_fr_group_catalogue_primary.png")
    write_report(out_dir / "sparc_fr_group_catalogue_report.md", summary, gate)
    export_bundle(out_dir, export_dir, programme)
    print(f"Saved report: {out_dir / 'sparc_fr_group_catalogue_report.md'}")
    print(f"Development gate passed: {gate['development_gate_passed']}")
    print("Reserved replication outcomes opened: no")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=["acquire", "preregister", "analyse"], required=True)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--export-dir", type=Path, default=DEFAULT_EXPORT)
    args = parser.parse_args()
    programme = Path(__file__).resolve()
    if args.phase == "acquire":
        acquire_phase(args.out_dir, args.export_dir)
    elif args.phase == "preregister":
        preregister_phase(args.out_dir, args.export_dir, programme)
    else:
        analyse_phase(args.out_dir, args.export_dir, programme)


if __name__ == "__main__":
    main()
