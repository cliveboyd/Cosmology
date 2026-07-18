#!/usr/bin/env python3
r"""Run the preregistered SPARC FR environmental-asymmetry study.

The programme has three deliberately separated phases:

1. ``acquire`` resolves non-reserved SPARC coordinates/redshifts through NED
   and constructs charge-blind 2MRS environment covariates without opening
   rotation-curve outcomes.
2. ``preregister`` freezes samples, outcomes, tests, gates and file hashes.
3. ``analyse`` verifies every lock before profiling model outcomes.

The six previously selected galaxies are excluded from the primary association
sample. The five reserved replication galaxies are excluded from every phase.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import shutil
import sys
import time
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

RUN_DATE = "2026-07-18"
SEED = 20260718
C_KMS = 299792.458

TASK_ROOT = Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu")
SAMPLE = ROOT / "external_datasets" / "current_cosmology_sources" / "SPARC" / "sparc_galaxy_sample.csv"
POINTS = ROOT / "external_datasets" / "current_cosmology_sources" / "SPARC" / "sparc_rotation_points.csv"
TWOMRS = (
    ROOT
    / "external_datasets"
    / "current_cosmology_sources"
    / "2MRS"
    / "2mrs_table3_environment_20260718.tsv"
)
POSITION_CACHE = (
    ROOT
    / "external_datasets"
    / "current_cosmology_sources"
    / "SPARC"
    / "sparc_ned_positions_environment_nonreserved_20260718.csv"
)
MAP_SUMMARY = (
    ROOT
    / "plamb_runs"
    / "diagnostics"
    / "sparc_hierarchical_map"
    / "optical_depth_hierarchical_20260714"
    / "sparc_hierarchical_map_summary.csv"
)
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "sparc_fr_environment_asymmetry_20260718"
DEFAULT_EXPORT = ROOT / "github_export" / "results" / "2026-07-18" / "fr_environment_asymmetry"

AM_PROGRAMME = ROOT / "github_export" / "code" / "sparc" / "run_sparc_am_identifiability_systematics_2026-07-18.py"
FR_RULE = (
    ROOT
    / "github_export"
    / "results"
    / "2026-07-18"
    / "fr_charge_conjugation"
    / "fr_charge_conjugation_rule.json"
)

TWOMRS_URL = (
    "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?"
    "-source=J%2FApJS%2F199%2F26%2Ftable3&"
    "-out=ID,RAJ2000,DEJ2000,GLON,GLAT,Ktmag,cz,SimbadName&-out.max=unlimited"
)
TWOMRS_LANDING = "https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJS/199/26"
TWOMRS_PAPER = "https://doi.org/10.1088/0067-0049/199/2/26"
NED_LOOKUP = "https://ned.ipac.caltech.edu/srs/ObjectLookup"
NED_GUIDE = "https://ned.ipac.caltech.edu/Documents/Guides/Interface/ObjectLookup"

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

PRIMARY_DEFINITION = {
    "name": "primary_10_40_mpc_mk22",
    "distance_min_mpc": 10.0,
    "distance_max_mpc": 40.0,
    "abs_galactic_latitude_min_deg": 10.0,
    "absolute_k_limit": -22.0,
}
SENSITIVITY_DEFINITIONS = (
    {
        "name": "near_5_40_mpc_mk22",
        "distance_min_mpc": 5.0,
        "distance_max_mpc": 40.0,
        "abs_galactic_latitude_min_deg": 10.0,
        "absolute_k_limit": -22.0,
    },
    {
        "name": "wide_10_80_mpc_mk23",
        "distance_min_mpc": 10.0,
        "distance_max_mpc": 80.0,
        "abs_galactic_latitude_min_deg": 10.0,
        "absolute_k_limit": -23.0,
    },
    {
        "name": "latitude_strict_10_40_mpc_mk22",
        "distance_min_mpc": 10.0,
        "distance_max_mpc": 40.0,
        "abs_galactic_latitude_min_deg": 20.0,
        "absolute_k_limit": -22.0,
    },
)

ENVIRONMENT_DV_KMS = 500.0
SIDEBAND_DV_MIN_KMS = 1000.0
SIDEBAND_DV_MAX_KMS = 2000.0
ENVIRONMENT_RADIUS_MPC = 5.0
INNER_RADIUS_MPC = 2.0
TIDAL_SOFTENING_MPC = 0.05
SELF_ANGLE_DEG = 0.5 / 60.0
SELF_DV_KMS = 300.0
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
    "log1p_sideband_count_5mpc",
)

ENVIRONMENT_COMPONENTS = (
    "log1p_count_2mpc",
    "log1p_count_5mpc",
    "minus_log10_nearest_mpc",
    "log10_tidal_strength",
)


def load_am_module() -> Any:
    spec = importlib.util.spec_from_file_location("sparc_am_systematics_locked", AM_PROGRAMME)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load {AM_PROGRAMME}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_csv(path: Path, rows: Iterable[dict[str, Any]] | pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(rows, pd.DataFrame):
        rows.to_csv(path, index=False)
        return
    materialised = list(rows)
    if not materialised:
        path.write_text("", encoding="utf-8")
        return
    fields: list[str] = []
    for row in materialised:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(materialised)


def markdown_table(rows: Iterable[dict[str, Any]], columns: list[str]) -> str:
    rows = list(rows)
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
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


def ensure_twomrs(path: Path) -> None:
    if path.exists() and path.stat().st_size > 1_000_000:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    response = requests.get(TWOMRS_URL, timeout=180)
    response.raise_for_status()
    path.write_bytes(response.content)
    if path.stat().st_size <= 1_000_000:
        raise RuntimeError(f"2MRS download is unexpectedly small: {path.stat().st_size} bytes")


def parse_twomrs(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    required = ["ID", "RAJ2000", "DEJ2000", "GLON", "GLAT", "Ktmag", "cz", "SimbadName"]
    missing = [column for column in required if column not in frame]
    if missing:
        raise ValueError(f"2MRS columns missing: {missing}")
    for column in ["RAJ2000", "DEJ2000", "GLON", "GLAT", "Ktmag", "cz"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame.dropna(subset=["RAJ2000", "DEJ2000", "Ktmag", "cz"]).copy()
    frame = frame[np.isfinite(frame["cz"])].reset_index(drop=True)
    return frame


def ned_lookup(session: requests.Session, name: str) -> dict[str, Any]:
    payload = {"json": json.dumps({"name": {"v": name}})}
    error = ""
    for attempt in range(4):
        try:
            response = session.post(NED_LOOKUP, data=payload, timeout=45)
            response.raise_for_status()
            result = response.json()
            preferred = result.get("Preferred") or {}
            position = preferred.get("Position") or {}
            redshift = preferred.get("Redshift") or {}
            ra = float(position["RA"])
            dec = float(position["Dec"])
            z_value = redshift.get("Value")
            z = float(z_value) if z_value is not None else float("nan")
            return {
                "Galaxy": name,
                "resolved": True,
                "resolver": "NASA_IPAC_NED_ObjectLookup",
                "preferred_name": preferred.get("Name", ""),
                "ra_deg": ra,
                "dec_deg": dec,
                "redshift": z,
                "cz_kms": z * C_KMS if math.isfinite(z) else float("nan"),
                "object_type": (preferred.get("ObjType") or {}).get("Value", ""),
                "position_reference": position.get("RefCode", ""),
                "redshift_reference": redshift.get("RefCode", ""),
                "resolution_error": "",
            }
        except Exception as exc:  # Network and catalogue failures are retained in the audit table.
            error = f"{type(exc).__name__}: {exc}"
            time.sleep(0.5 * (attempt + 1))
    return {
        "Galaxy": name,
        "resolved": False,
        "resolver": "NASA_IPAC_NED_ObjectLookup",
        "preferred_name": "",
        "ra_deg": float("nan"),
        "dec_deg": float("nan"),
        "redshift": float("nan"),
        "cz_kms": float("nan"),
        "object_type": "",
        "position_reference": "",
        "redshift_reference": "",
        "resolution_error": error,
    }


def acquire_positions(sample: pd.DataFrame, path: Path) -> pd.DataFrame:
    target_names = sorted(
        set(sample.loc[sample["Q"] <= 2, "Galaxy"].astype(str)) - set(RESERVED_REPLICATION)
    )
    if set(target_names) & set(RESERVED_REPLICATION):
        raise AssertionError("Reserved replication name entered the resolver target list")
    existing = pd.DataFrame()
    if path.exists() and path.stat().st_size > 0:
        existing = pd.read_csv(path)
    completed = set(existing.loc[existing.get("resolved", False).astype(str).str.lower() == "true", "Galaxy"].astype(str)) if len(existing) else set()
    rows = existing.to_dict("records") if len(existing) else []
    pending = [name for name in target_names if name not in completed]
    if pending:
        session = requests.Session()
        session.headers.update({"User-Agent": "SPARC-FR-environment-audit/1.0 research reproducibility"})
        for index, name in enumerate(pending, start=1):
            rows = [row for row in rows if str(row.get("Galaxy")) != name]
            rows.append(ned_lookup(session, name))
            if index % 20 == 0 or index == len(pending):
                print(f"NED resolved {index}/{len(pending)} pending names", flush=True)
                write_csv(path, pd.DataFrame(rows).sort_values("Galaxy"))
    output = pd.DataFrame(rows)
    output = output[output["Galaxy"].isin(target_names)].sort_values("Galaxy").reset_index(drop=True)
    if set(output["Galaxy"]) & set(RESERVED_REPLICATION):
        raise AssertionError("Reserved replication galaxy was written to the position cache")
    write_csv(path, output)
    return output


def angular_geometry(ra_deg: float, dec_deg: float, catalogue: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    ra = np.deg2rad(catalogue["RAJ2000"].to_numpy(dtype=float))
    dec = np.deg2rad(catalogue["DEJ2000"].to_numpy(dtype=float))
    ra0 = math.radians(ra_deg)
    dec0 = math.radians(dec_deg)
    cos_angle = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    angle = np.arccos(np.clip(cos_angle, -1.0, 1.0))
    return angle, np.rad2deg(angle)


def environment_metrics_at_position(
    ra_deg: float,
    dec_deg: float,
    distance_mpc: float,
    central_cz_kms: float,
    catalogue: pd.DataFrame,
    absolute_k_limit: float,
) -> dict[str, float]:
    angle, angle_deg = angular_geometry(ra_deg, dec_deg, catalogue)
    projected_mpc = distance_mpc * np.sin(angle)
    delta_v = np.abs(catalogue["cz"].to_numpy(dtype=float) - central_cz_kms)
    absolute_k = catalogue["Ktmag"].to_numpy(dtype=float) - 5.0 * math.log10(distance_mpc) - 25.0
    bright = absolute_k <= absolute_k_limit
    is_self = (angle_deg <= SELF_ANGLE_DEG) & (delta_v <= SELF_DV_KMS)
    physical = bright & (~is_self) & (delta_v <= ENVIRONMENT_DV_KMS)
    within_2 = physical & (projected_mpc <= INNER_RADIUS_MPC)
    within_5 = physical & (projected_mpc <= ENVIRONMENT_RADIUS_MPC)
    sideband = (
        bright
        & (~is_self)
        & (delta_v >= SIDEBAND_DV_MIN_KMS)
        & (delta_v < SIDEBAND_DV_MAX_KMS)
        & (projected_mpc <= ENVIRONMENT_RADIUS_MPC)
    )
    distances = projected_mpc[within_5]
    nearest = float(np.min(distances)) if distances.size else ENVIRONMENT_RADIUS_MPC
    luminosity_relative = 10.0 ** (-0.4 * (absolute_k + 23.0))
    tidal = float(
        np.sum(
            luminosity_relative[within_5]
            / np.maximum(projected_mpc[within_5] ** 2 + TIDAL_SOFTENING_MPC**2, 1e-12) ** 1.5
        )
    )
    return {
        "count_2mpc": int(np.sum(within_2)),
        "count_5mpc": int(np.sum(within_5)),
        "nearest_mpc": nearest,
        "tidal_strength": tidal,
        "sideband_count_5mpc": int(np.sum(sideband)),
    }


def build_environment_covariates(
    sample: pd.DataFrame,
    positions: pd.DataFrame,
    catalogue: pd.DataFrame,
) -> pd.DataFrame:
    position_lookup = positions.set_index("Galaxy")
    rows: list[dict[str, Any]] = []
    thresholds = sorted({float(PRIMARY_DEFINITION["absolute_k_limit"])} | {float(item["absolute_k_limit"]) for item in SENSITIVITY_DEFINITIONS})
    for _, galaxy_row in sample[sample["Q"] <= 2].sort_values("Galaxy").iterrows():
        galaxy = str(galaxy_row["Galaxy"])
        if galaxy in RESERVED_REPLICATION:
            continue
        position = position_lookup.loc[galaxy] if galaxy in position_lookup.index else None
        resolved = position is not None and str(position.get("resolved", "")).lower() in {"true", "1"}
        ra = float(position["ra_deg"]) if resolved else float("nan")
        dec = float(position["dec_deg"]) if resolved else float("nan")
        cz = float(position["cz_kms"]) if resolved else float("nan")
        gal_b = float("nan")
        if resolved and math.isfinite(ra) and math.isfinite(dec):
            gal_b = float(SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs").galactic.b.deg)
        output: dict[str, Any] = {
            "Galaxy": galaxy,
            "Q": int(galaxy_row["Q"]),
            "D_Mpc": float(galaxy_row["D_Mpc"]),
            "resolved": resolved,
            "preferred_name": position["preferred_name"] if resolved else "",
            "ra_deg": ra,
            "dec_deg": dec,
            "gal_b_deg": gal_b,
            "abs_gal_b_deg": abs(gal_b) if math.isfinite(gal_b) else float("nan"),
            "central_cz_kms": cz,
        }
        for threshold in thresholds:
            label = f"mk{abs(int(round(threshold)))}"
            if resolved and all(math.isfinite(value) for value in [ra, dec, cz]) and float(galaxy_row["D_Mpc"]) > 0:
                metrics = environment_metrics_at_position(
                    ra,
                    dec,
                    float(galaxy_row["D_Mpc"]),
                    cz,
                    catalogue,
                    threshold,
                )
                shifted = environment_metrics_at_position(
                    (ra + 90.0) % 360.0,
                    dec,
                    float(galaxy_row["D_Mpc"]),
                    cz,
                    catalogue,
                    threshold,
                )
            else:
                metrics = {key: float("nan") for key in ["count_2mpc", "count_5mpc", "nearest_mpc", "tidal_strength", "sideband_count_5mpc"]}
                shifted = dict(metrics)
            for key, value in metrics.items():
                output[f"{key}_{label}"] = value
            for key, value in shifted.items():
                output[f"shifted_{key}_{label}"] = value
        rows.append(output)
    frame = pd.DataFrame(rows).sort_values("Galaxy").reset_index(drop=True)
    if set(frame["Galaxy"]) & set(RESERVED_REPLICATION):
        raise AssertionError("Reserved replication galaxy entered environment covariates")
    return frame


def sample_mask(frame: pd.DataFrame, definition: dict[str, Any], exclude_candidates: bool = True) -> pd.Series:
    mask = (
        (frame["Q"] <= 2)
        & frame["resolved"].astype(bool)
        & np.isfinite(pd.to_numeric(frame["central_cz_kms"], errors="coerce"))
        & (frame["D_Mpc"] >= float(definition["distance_min_mpc"]))
        & (frame["D_Mpc"] <= float(definition["distance_max_mpc"]))
        & (frame["abs_gal_b_deg"] >= float(definition["abs_galactic_latitude_min_deg"]))
        & (~frame["Galaxy"].isin(RESERVED_REPLICATION))
    )
    if exclude_candidates:
        mask &= ~frame["Galaxy"].isin(CANDIDATES)
    return mask


def preregistration_config(environment: pd.DataFrame, programme: Path) -> dict[str, Any]:
    definitions = [PRIMARY_DEFINITION, *SENSITIVITY_DEFINITIONS]
    sample_counts = {
        str(item["name"]): int(np.sum(sample_mask(environment, item, exclude_candidates=True)))
        for item in definitions
    }
    return {
        "date": RUN_DATE,
        "written_before_outcomes_utc": utc_now(),
        "programme": str(programme),
        "programme_sha256": sha256(programme),
        "outcome_opened_before_lock": False,
        "theory_boundary": (
            "FR complete charge conjugation is charge-even. This study tests whether an independently observed "
            "charge-blind environment predicts model residuals; it does not classify antimatter."
        ),
        "development_exclusions": {
            "previously_outcome_selected": list(CANDIDATES),
            "reserved_replication": list(RESERVED_REPLICATION),
        },
        "catalogues": {
            "central_positions_and_redshifts": {
                "source": "NASA/IPAC NED ObjectLookup",
                "url": NED_GUIDE,
                "cache": str(POSITION_CACHE),
            },
            "environment_tracer": {
                "source": "2MASS Redshift Survey table 3",
                "url": TWOMRS_LANDING,
                "paper": TWOMRS_PAPER,
                "catalogue": str(TWOMRS),
                "catalogue_limit": "K_s <= 11.75 mag; 97.6 per cent redshift complete over 91 per cent of sky",
            },
        },
        "input_hashes": {
            "sample": sha256(SAMPLE),
            "points": sha256(POINTS),
            "map_summary": sha256(MAP_SUMMARY),
            "2mrs": sha256(TWOMRS),
            "positions": sha256(POSITION_CACHE),
            "fr_rule": sha256(FR_RULE),
            "am_profile_programme": sha256(AM_PROGRAMME),
        },
        "samples": {
            "primary": PRIMARY_DEFINITION,
            "sensitivities": list(SENSITIVITY_DEFINITIONS),
            "outcome_blind_counts": sample_counts,
            "quality": "SPARC Q <= 2; finite NED position/redshift; stated distance and latitude cuts",
        },
        "environment": {
            "line_of_sight_window_kms": ENVIRONMENT_DV_KMS,
            "projected_radii_mpc": [INNER_RADIUS_MPC, ENVIRONMENT_RADIUS_MPC],
            "sideband_kms": [SIDEBAND_DV_MIN_KMS, SIDEBAND_DV_MAX_KMS],
            "self_exclusion": {"angle_deg": SELF_ANGLE_DEG, "delta_v_kms": SELF_DV_KMS},
            "tidal_softening_mpc": TIDAL_SOFTENING_MPC,
            "composite_components": list(ENVIRONMENT_COMPONENTS),
            "composite_rule": (
                "mean of robust-standardised log1p N(<2 Mpc), log1p N(<5 Mpc), "
                "-log10 nearest-neighbour distance and log10 softened K-luminosity tidal sum"
            ),
            "negative_control": "same proxy at RA shifted by +90 degrees with Dec and central redshift fixed",
        },
        "outcomes": {
            "primary": (
                "signed_log1p((profile_objective_PLAMB - profile_objective_RAR) / N_points) "
                "under the locked combined_conventional nuisance scenario"
            ),
            "secondary": "same transformed contrast under the baseline nuisance scenario",
            "global_parameters": "fixed July 14 all-Q2 MAP values",
            "optimisation": "deterministic multi-start galaxy profile inherited unchanged from the July 18 systematics audit",
        },
        "host_controls": list(HOST_FEATURES) + ["f_D one-hot indicators"],
        "statistics": {
            "ridge_alpha": RIDGE_ALPHA,
            "outer_folds": OUTER_FOLDS,
            "fold_seed": SEED,
            "partial_test": "Spearman correlation of cross-fitted host-control residuals",
            "permutations": PERMUTATIONS,
            "p_value": "two-sided deterministic permutation p=(1+exceedances)/(1+permutations)",
            "predictive_metric": "cross-fitted RMSE baseline versus baseline plus environment composite",
        },
        "development_gate_for_future_replication": {
            "all_required": True,
            "primary_partial_permutation_p_max": 0.01,
            "primary_cv_rmse_fractional_improvement_min": 0.05,
            "environment_coefficient_sign_fraction_min": 0.80,
            "same_partial_correlation_sign_in": [
                "wide_10_80_mpc_mk23",
                "latitude_strict_10_40_mpc_mk22",
                "baseline_outcome_primary_sample",
            ],
            "negative_control": "absolute shifted-position partial rho must be smaller than actual-environment rho",
            "replication_policy": "never unseal automatically; write a separate replication preregistration first",
        },
        "claim_boundary": [
            "A positive result would show ordinary environment dependence of model residuals, not antimatter.",
            "No fitted matter/antimatter galaxy sign is permitted.",
            "Candidate-galaxy environment comparisons are descriptive because those galaxies were outcome-selected.",
            "Reserved replication galaxies remain unopened in this programme.",
        ],
    }


def write_preregistration(path: Path, config: dict[str, Any]) -> None:
    counts = [
        {"sample": key, "N_outcome_blind": value}
        for key, value in config["samples"]["outcome_blind_counts"].items()
    ]
    lines = [
        "# SPARC FR Environmental-Asymmetry Preregistration",
        "",
        f"Date: {RUN_DATE}",
        f"Written before outcome profiling: `{config['written_before_outcomes_utc']}`",
        "",
        "## Scientific Question",
        "",
        "Does an independently measured, charge-blind galaxy environment predict where the fixed-global PLAMB bridge differs from RAR after conventional galaxy properties and nuisance structure are controlled?",
        "",
        "This is the only rotation-curve consequence left open by the FR charge-conjugation derivation: a host-only change in a fixed asymmetric environment can alter the magnitude of the FR background asymmetry. The observable 2MRS environment is not a signed matter/antimatter field and cannot identify antimatter.",
        "",
        "## Locked Samples",
        "",
        markdown_table(counts, ["sample", "N_outcome_blind"]),
        "",
        "The primary sample uses SPARC Q<=2 galaxies at 10-40 Mpc, |b|>=10 degrees and a 2MRS neighbour threshold M_K<=-22. The six previously outcome-selected development galaxies and all five reserved replication galaxies are excluded from primary inference.",
        "",
        "## Environment Composite",
        "",
        r"$$",
        r"\begin{aligned}",
        r"E_1 &= \log(1+N_{2\,{\rm Mpc}}), \\",
        r"E_2 &= \log(1+N_{5\,{\rm Mpc}}), \\",
        r"E_3 &= -\log_{10} R_{\rm nearest}, \\",
        r"E_4 &= \log_{10}\!\left[10^{-6}+\sum_j\frac{L_{K,j}/L_{K,-23}}{(R_j^2+0.05^2)^{3/2}}\right], \\",
        r"E   &= \frac{1}{4}\sum_{k=1}^{4}\frac{E_k-\operatorname{median}(E_k)}{\operatorname{IQR}(E_k)}.",
        r"\end{aligned}",
        r"$$",
        "",
        "Neighbours must lie within |Delta v|<=500 km/s. A 1000<=|Delta v|<2000 km/s sideband enters the host-control model as a projection/selection control. The negative-control environment repeats the calculation at RA+90 degrees.",
        "",
        "## Locked Outcome",
        "",
        r"$$",
        r"\begin{aligned}",
        r"\Delta_g &= \frac{Q_{{\rm PLAMB},g}-Q_{{\rm RAR},g}}{N_g}, \\",
        r"Y_g      &= \operatorname{sign}(\Delta_g)\log(1+|\Delta_g|).",
        r"\end{aligned}",
        r"$$",
        "",
        "The primary outcome uses the combined-conventional nuisance profile. The baseline profile is secondary. Positive values mean that PLAMB fits worse than RAR.",
        "",
        "## Locked Test",
        "",
        "Ridge models with fixed alpha=1 use ten deterministic folds. The primary association is a two-sided 20,000-permutation Spearman test of cross-fitted outcome and environment residuals after the locked host controls. Predictive value is the fractional reduction in cross-fitted RMSE when E is added.",
        "",
        "## Development Gate",
        "",
        "Every condition below is required before a separate replication preregistration may be written:",
        "",
        "1. primary permutation p<=0.01;",
        "2. cross-fitted RMSE improvement>=5 per cent;",
        "3. environment coefficient has one sign in at least 8/10 folds;",
        "4. the partial-correlation sign agrees in the wide-distance, strict-latitude and baseline-outcome controls; and",
        "5. the shifted-position negative-control |rho| is smaller than the actual-environment |rho|.",
        "",
        "Passing this gate would establish an environment residual worth replication. It would not establish a signed antimatter field.",
        "",
        "## Reproducibility Lock",
        "",
        "```json",
        json.dumps(config, indent=2),
        "```",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def copy_preregistration(out_dir: Path, export_dir: Path) -> None:
    export_dir.mkdir(parents=True, exist_ok=True)
    for name in [
        "sparc_fr_environment_preregistration.json",
        "sparc_fr_environment_preregistration.md",
        "sparc_fr_environment_covariates_nonreserved.csv",
        "sparc_fr_environment_acquisition_audit.json",
    ]:
        source = out_dir / name
        if source.exists():
            shutil.copy2(source, export_dir / name)


def acquire_phase(out_dir: Path, export_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    sample = pd.read_csv(SAMPLE)
    for column in ["Q", "D_Mpc"]:
        sample[column] = pd.to_numeric(sample[column], errors="coerce")
    ensure_twomrs(TWOMRS)
    twomrs = parse_twomrs(TWOMRS)
    positions = acquire_positions(sample, POSITION_CACHE)
    environment = build_environment_covariates(sample, positions, twomrs)
    environment_path = out_dir / "sparc_fr_environment_covariates_nonreserved.csv"
    write_csv(environment_path, environment)
    resolved = positions[positions["resolved"].astype(str).str.lower().isin({"true", "1"})]
    finite_cz = resolved[np.isfinite(pd.to_numeric(resolved["cz_kms"], errors="coerce"))]
    audit = {
        "date": RUN_DATE,
        "completed_utc": utc_now(),
        "outcomes_opened": False,
        "2mrs": {
            "path": str(TWOMRS),
            "url": TWOMRS_URL,
            "sha256": sha256(TWOMRS),
            "bytes": TWOMRS.stat().st_size,
            "rows_with_finite_position_magnitude_redshift": int(len(twomrs)),
        },
        "ned_positions": {
            "path": str(POSITION_CACHE),
            "url": NED_GUIDE,
            "sha256": sha256(POSITION_CACHE),
            "requested_nonreserved_q2": int(len(positions)),
            "resolved_positions": int(len(resolved)),
            "finite_redshifts": int(len(finite_cz)),
            "unresolved": sorted(set(positions.loc[~positions.index.isin(resolved.index), "Galaxy"].astype(str))),
        },
        "reserved_replication_requested": False,
        "reserved_replication_present_in_covariates": bool(
            set(environment["Galaxy"].astype(str)) & set(RESERVED_REPLICATION)
        ),
    }
    (out_dir / "sparc_fr_environment_acquisition_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    copy_preregistration(out_dir, export_dir)
    print(f"Saved outcome-blind environment covariates: {environment_path}")
    print(f"NED finite redshifts: {len(finite_cz)}/{len(positions)}")
    print("Reserved replication galaxies requested: no")


def preregister_phase(out_dir: Path, export_dir: Path, programme: Path) -> None:
    environment_path = out_dir / "sparc_fr_environment_covariates_nonreserved.csv"
    if not environment_path.exists():
        raise FileNotFoundError("Run --phase acquire before preregistration")
    required = [SAMPLE, POINTS, TWOMRS, POSITION_CACHE, MAP_SUMMARY, FR_RULE, AM_PROGRAMME]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Required inputs missing: {missing}")
    environment = pd.read_csv(environment_path)
    config = preregistration_config(environment, programme)
    config["input_hashes"]["environment_covariates"] = sha256(environment_path)
    json_path = out_dir / "sparc_fr_environment_preregistration.json"
    md_path = out_dir / "sparc_fr_environment_preregistration.md"
    json_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    write_preregistration(md_path, config)
    copy_preregistration(out_dir, export_dir)
    print(f"Saved preregistration: {md_path}")
    print(f"Locked programme SHA-256: {config['programme_sha256']}")
    print("Rotation-curve outcomes opened: no")


def verify_locks(config: dict[str, Any], programme: Path, environment_path: Path) -> None:
    current_programme = sha256(programme)
    if current_programme != config["programme_sha256"]:
        raise RuntimeError(
            "Programme hash changed after preregistration: "
            f"locked={config['programme_sha256']} current={current_programme}"
        )
    paths = {
        "sample": SAMPLE,
        "points": POINTS,
        "map_summary": MAP_SUMMARY,
        "2mrs": TWOMRS,
        "positions": POSITION_CACHE,
        "fr_rule": FR_RULE,
        "am_profile_programme": AM_PROGRAMME,
        "environment_covariates": environment_path,
    }
    for key, path in paths.items():
        actual = sha256(path)
        expected = config["input_hashes"][key]
        if actual != expected:
            raise RuntimeError(f"Input hash changed for {key}: locked={expected} current={actual}")


def transformed_environment_components(
    frame: pd.DataFrame,
    definition: dict[str, Any],
    scaling: dict[str, dict[str, float]] | None = None,
) -> tuple[pd.DataFrame, dict[str, dict[str, float]]]:
    output = frame.copy()
    label = f"mk{abs(int(round(float(definition['absolute_k_limit']))))}"
    actual_raw = {
        "log1p_count_2mpc": np.log1p(pd.to_numeric(output[f"count_2mpc_{label}"], errors="coerce")),
        "log1p_count_5mpc": np.log1p(pd.to_numeric(output[f"count_5mpc_{label}"], errors="coerce")),
        "minus_log10_nearest_mpc": -np.log10(
            np.maximum(pd.to_numeric(output[f"nearest_mpc_{label}"], errors="coerce"), 0.01)
        ),
        "log10_tidal_strength": np.log10(
            np.maximum(pd.to_numeric(output[f"tidal_strength_{label}"], errors="coerce"), 0.0) + 1.0e-6
        ),
    }
    shifted_raw = {
        "log1p_count_2mpc": np.log1p(
            pd.to_numeric(output[f"shifted_count_2mpc_{label}"], errors="coerce")
        ),
        "log1p_count_5mpc": np.log1p(
            pd.to_numeric(output[f"shifted_count_5mpc_{label}"], errors="coerce")
        ),
        "minus_log10_nearest_mpc": -np.log10(
            np.maximum(pd.to_numeric(output[f"shifted_nearest_mpc_{label}"], errors="coerce"), 0.01)
        ),
        "log10_tidal_strength": np.log10(
            np.maximum(pd.to_numeric(output[f"shifted_tidal_strength_{label}"], errors="coerce"), 0.0)
            + 1.0e-6
        ),
    }
    if scaling is None:
        scaling = {}
        for name, values in actual_raw.items():
            finite = values[np.isfinite(values)]
            median = float(np.median(finite))
            q25, q75 = np.quantile(finite, [0.25, 0.75])
            scale = float(q75 - q25)
            if not math.isfinite(scale) or scale <= 1.0e-12:
                scale = float(np.std(finite, ddof=0))
            if not math.isfinite(scale) or scale <= 1.0e-12:
                scale = 1.0
            scaling[name] = {"median": median, "iqr": scale}
    actual_z: list[np.ndarray] = []
    shifted_z: list[np.ndarray] = []
    for name in ENVIRONMENT_COMPONENTS:
        centre = scaling[name]["median"]
        scale = scaling[name]["iqr"]
        output[name] = actual_raw[name]
        output[f"shifted_{name}"] = shifted_raw[name]
        actual_z.append((actual_raw[name].to_numpy(dtype=float) - centre) / scale)
        shifted_z.append((shifted_raw[name].to_numpy(dtype=float) - centre) / scale)
    output["environment_score"] = np.mean(np.column_stack(actual_z), axis=1)
    output["shifted_environment_score"] = np.mean(np.column_stack(shifted_z), axis=1)
    output["log1p_sideband_count_5mpc"] = np.log1p(
        pd.to_numeric(output[f"sideband_count_5mpc_{label}"], errors="coerce")
    )
    return output, scaling


def development_names(environment: pd.DataFrame) -> list[str]:
    definitions = [PRIMARY_DEFINITION, *SENSITIVITY_DEFINITIONS]
    names: set[str] = set()
    for definition in definitions:
        names.update(environment.loc[sample_mask(environment, definition, True), "Galaxy"].astype(str))
    if names & set(CANDIDATES):
        raise AssertionError("Previously outcome-selected candidate entered primary/sensitivity outcomes")
    if names & set(RESERVED_REPLICATION):
        raise AssertionError("Reserved replication galaxy entered development outcomes")
    return sorted(names)


def profile_development_outcomes(names: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    am = load_am_module()
    data = am.load_dataset(SAMPLE, POINTS, quality_max=2, distance_method="all", err_floor=3.0)
    sample = pd.read_csv(SAMPLE)
    sample_info: dict[str, dict[str, float]] = {}
    for _, row in sample.iterrows():
        sample_info[str(row["Galaxy"])] = {
            key: float(row[key]) for key in ["D_Mpc", "e_D_Mpc", "f_D", "Inc_deg", "e_Inc_deg"]
        }
    available = set(str(value) for value in data.galaxy)
    missing = sorted(set(names) - available)
    if missing:
        raise ValueError(f"Development galaxies missing from Q<=2 rotation data: {missing}")
    globals_by_model = am.load_global_parameters()
    scenarios = {
        scenario.name: scenario
        for scenario in am.SCENARIOS
        if scenario.name in {"baseline", "combined_conventional"}
    }
    if set(scenarios) != {"baseline", "combined_conventional"}:
        raise ValueError("Could not load both locked nuisance scenarios")
    profile_rows: list[dict[str, Any]] = []
    total = len(names) * len(scenarios) * len(am.MODELS)
    complete = 0
    for scenario_name in ["baseline", "combined_conventional"]:
        scenario = scenarios[scenario_name]
        print(f"=== profiling {scenario_name} ===", flush=True)
        for galaxy in names:
            for model in am.MODELS:
                row = am.profile_galaxy(
                    galaxy,
                    model,
                    scenario,
                    data,
                    sample_info,
                    globals_by_model[model],
                )
                profile_rows.append(row)
                complete += 1
            if complete % 20 == 0 or complete == total:
                print(f"profiled {complete}/{total}", flush=True)
    profile = pd.DataFrame(profile_rows)
    paired_rows: list[dict[str, Any]] = []
    for scenario_name in ["baseline", "combined_conventional"]:
        selected = profile[profile["scenario"] == scenario_name]
        for galaxy in names:
            rows = selected[selected["galaxy"] == galaxy].set_index("model")
            rar = rows.loc["RAR"]
            plamb = rows.loc["PLAMB_OPTICAL_DEPTH_KAPPA_P"]
            n_points = int(rar["N_points"])
            delta = float(plamb["profile_objective"] - rar["profile_objective"])
            per_point = delta / max(n_points, 1)
            transformed = math.copysign(math.log1p(abs(per_point)), per_point) if per_point != 0 else 0.0
            paired_rows.append(
                {
                    "Galaxy": galaxy,
                    "scenario": scenario_name,
                    "N_points": n_points,
                    "rar_success": bool(rar["success"]),
                    "plamb_success": bool(plamb["success"]),
                    "rar_profile_objective": float(rar["profile_objective"]),
                    "plamb_profile_objective": float(plamb["profile_objective"]),
                    "delta_objective_plamb_minus_rar": delta,
                    "delta_objective_per_point": per_point,
                    "signed_log1p_delta_per_point": transformed,
                }
            )
    paired = pd.DataFrame(paired_rows)
    if set(paired["Galaxy"]) & (set(CANDIDATES) | set(RESERVED_REPLICATION)):
        raise AssertionError("Excluded galaxy entered outcome table")
    return profile, paired


def host_catalogue() -> pd.DataFrame:
    am = load_am_module()
    sample = pd.read_csv(SAMPLE)
    points = pd.read_csv(POINTS)
    catalogue = am.transformed_catalogue(sample, points)
    catalogue["D_Mpc"] = pd.to_numeric(catalogue["D_Mpc"], errors="coerce")
    return catalogue


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


def deterministic_permutation_p(x: np.ndarray, y: np.ndarray, seed_offset: int = 0) -> tuple[float, float]:
    rho = float(spearmanr(x, y).statistic)
    if not math.isfinite(rho):
        return float("nan"), float("nan")
    x_rank = pd.Series(x).rank(method="average").to_numpy(dtype=float)
    y_rank = pd.Series(y).rank(method="average").to_numpy(dtype=float)
    x_rank -= np.mean(x_rank)
    y_rank -= np.mean(y_rank)
    denominator = math.sqrt(float(np.sum(x_rank**2) * np.sum(y_rank**2)))
    observed = float(np.dot(x_rank, y_rank) / denominator)
    rng = np.random.default_rng(SEED + seed_offset)
    exceed = 0
    for _ in range(PERMUTATIONS):
        permuted = rng.permutation(y_rank)
        value = float(np.dot(x_rank, permuted) / denominator)
        if abs(value) >= abs(observed) - 1.0e-15:
            exceed += 1
    return rho, (1.0 + exceed) / (1.0 + PERMUTATIONS)


def cross_fitted_test(
    frame: pd.DataFrame,
    environment_column: str,
    test_label: str,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    y = frame["signed_log1p_delta_per_point"].to_numpy(dtype=float)
    environment = frame[environment_column].to_numpy(dtype=float)
    x, feature_names = host_matrix(frame)
    n_splits = min(OUTER_FOLDS, max(3, len(frame) // 5))
    folds = KFold(n_splits=n_splits, shuffle=True, random_state=SEED)
    baseline_prediction = np.full(len(frame), np.nan)
    augmented_prediction = np.full(len(frame), np.nan)
    environment_prediction = np.full(len(frame), np.nan)
    fold_rows: list[dict[str, Any]] = []
    for fold_index, (train, test) in enumerate(folds.split(x), start=1):
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
        augmented_train = augmented_scaler.transform(
            np.column_stack([x[train], environment[train]])
        )
        augmented_test = augmented_scaler.transform(
            np.column_stack([x[test], environment[test]])
        )
        augmented = Ridge(alpha=RIDGE_ALPHA).fit(augmented_train, y[train])
        augmented_prediction[test] = augmented.predict(augmented_test)
        fold_rows.append(
            {
                "test": test_label,
                "fold": fold_index,
                "N_train": len(train),
                "N_test": len(test),
                "environment_coefficient_standardised": float(augmented.coef_[-1]),
            }
        )
    y_residual = y - baseline_prediction
    environment_residual = environment - environment_prediction
    rho, p_value = deterministic_permutation_p(
        environment_residual,
        y_residual,
        seed_offset=int(hashlib.sha256(test_label.encode("utf-8")).hexdigest()[:8], 16),
    )
    baseline_rmse = math.sqrt(mean_squared_error(y, baseline_prediction))
    augmented_rmse = math.sqrt(mean_squared_error(y, augmented_prediction))
    improvement = (baseline_rmse - augmented_rmse) / baseline_rmse if baseline_rmse > 0 else float("nan")
    coefficients = np.asarray([row["environment_coefficient_standardised"] for row in fold_rows], dtype=float)
    positive = float(np.mean(coefficients > 0.0))
    negative = float(np.mean(coefficients < 0.0))
    sign_fraction = max(positive, negative)
    result = {
        "test": test_label,
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
            "test": test_label,
            "outcome": y,
            "environment_score": environment,
            "outcome_crossfit_residual": y_residual,
            "environment_crossfit_residual": environment_residual,
            "baseline_prediction": baseline_prediction,
            "environment_prediction": augmented_prediction,
        }
    )
    return result, pd.DataFrame(fold_rows), residuals


def build_test_frame(
    environment: pd.DataFrame,
    catalogue: pd.DataFrame,
    outcomes: pd.DataFrame,
    definition: dict[str, Any],
    scenario: str,
) -> tuple[pd.DataFrame, dict[str, dict[str, float]]]:
    selected = environment.loc[sample_mask(environment, definition, True)].copy()
    selected, scaling = transformed_environment_components(selected, definition)
    selected = selected.merge(catalogue, on="Galaxy", how="inner", suffixes=("", "_sample"))
    selected = selected.merge(
        outcomes[outcomes["scenario"] == scenario],
        on="Galaxy",
        how="inner",
    )
    selected = selected[selected["rar_success"] & selected["plamb_success"]].copy()
    selected = selected.sort_values("Galaxy").reset_index(drop=True)
    if set(selected["Galaxy"]) & (set(CANDIDATES) | set(RESERVED_REPLICATION)):
        raise AssertionError("Excluded galaxy entered test frame")
    return selected, scaling


def evaluate_tests(
    environment: pd.DataFrame,
    catalogue: pd.DataFrame,
    outcomes: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any], dict[str, dict[str, float]]]:
    definitions = [PRIMARY_DEFINITION, *SENSITIVITY_DEFINITIONS]
    summary_rows: list[dict[str, Any]] = []
    fold_frames: list[pd.DataFrame] = []
    residual_frames: list[pd.DataFrame] = []
    primary_frame = pd.DataFrame()
    primary_scaling: dict[str, dict[str, float]] = {}
    for definition in definitions:
        frame, scaling = build_test_frame(
            environment,
            catalogue,
            outcomes,
            definition,
            "combined_conventional",
        )
        label = f"{definition['name']}__combined_conventional"
        result, folds, residuals = cross_fitted_test(frame, "environment_score", label)
        result.update({"sample": definition["name"], "scenario": "combined_conventional", "control": "actual"})
        summary_rows.append(result)
        fold_frames.append(folds)
        residual_frames.append(residuals)
        if definition["name"] == PRIMARY_DEFINITION["name"]:
            primary_frame = frame
            primary_scaling = scaling
            shifted_result, shifted_folds, shifted_residuals = cross_fitted_test(
                frame,
                "shifted_environment_score",
                f"{definition['name']}__combined_conventional__shifted_ra90",
            )
            shifted_result.update(
                {"sample": definition["name"], "scenario": "combined_conventional", "control": "shifted_ra90"}
            )
            summary_rows.append(shifted_result)
            fold_frames.append(shifted_folds)
            residual_frames.append(shifted_residuals)
    baseline_frame, _ = build_test_frame(
        environment,
        catalogue,
        outcomes,
        PRIMARY_DEFINITION,
        "baseline",
    )
    baseline_result, baseline_folds, baseline_residuals = cross_fitted_test(
        baseline_frame,
        "environment_score",
        f"{PRIMARY_DEFINITION['name']}__baseline",
    )
    baseline_result.update(
        {"sample": PRIMARY_DEFINITION["name"], "scenario": "baseline", "control": "actual"}
    )
    summary_rows.append(baseline_result)
    fold_frames.append(baseline_folds)
    residual_frames.append(baseline_residuals)

    summary = pd.DataFrame(summary_rows)
    primary = summary[
        (summary["sample"] == PRIMARY_DEFINITION["name"])
        & (summary["scenario"] == "combined_conventional")
        & (summary["control"] == "actual")
    ].iloc[0]
    shifted = summary[
        (summary["sample"] == PRIMARY_DEFINITION["name"])
        & (summary["scenario"] == "combined_conventional")
        & (summary["control"] == "shifted_ra90")
    ].iloc[0]
    sign = int(np.sign(float(primary["partial_spearman_rho"])))
    agreement_labels = [
        ("wide_10_80_mpc_mk23", "combined_conventional"),
        ("latitude_strict_10_40_mpc_mk22", "combined_conventional"),
        (PRIMARY_DEFINITION["name"], "baseline"),
    ]
    agreement = []
    for sample_name, scenario in agreement_labels:
        row = summary[
            (summary["sample"] == sample_name)
            & (summary["scenario"] == scenario)
            & (summary["control"] == "actual")
        ].iloc[0]
        agreement.append(int(np.sign(float(row["partial_spearman_rho"]))) == sign and sign != 0)
    conditions = {
        "primary_p": bool(float(primary["permutation_p_two_sided"]) <= 0.01),
        "cv_improvement": bool(float(primary["cv_rmse_fractional_improvement"]) >= 0.05),
        "coefficient_sign": bool(float(primary["environment_coefficient_sign_fraction"]) >= 0.80),
        "sensitivity_signs": bool(all(agreement)),
        "negative_control": bool(
            abs(float(shifted["partial_spearman_rho"])) < abs(float(primary["partial_spearman_rho"]))
        ),
    }
    gate = {
        "date": RUN_DATE,
        "development_gate_passed": bool(all(conditions.values())),
        "conditions": conditions,
        "replication_unsealed": False,
        "replication_policy": "A separate preregistration is required even if every development condition passes.",
    }
    return (
        summary,
        pd.concat(fold_frames, ignore_index=True),
        pd.concat(residual_frames, ignore_index=True),
        primary_frame,
        gate,
        primary_scaling,
    )


def candidate_environment_description(
    environment: pd.DataFrame,
    primary_scaling: dict[str, dict[str, float]],
) -> pd.DataFrame:
    matched_path = (
        ROOT
        / "github_export"
        / "results"
        / "2026-07-18"
        / "sparc_am"
        / "sparc_am_matched_controls.csv"
    )
    if not matched_path.exists():
        return pd.DataFrame()
    assignments = pd.read_csv(matched_path)
    names = set(CANDIDATES) | set(assignments["control"].astype(str))
    selected = environment[environment["Galaxy"].isin(names)].copy()
    selected, _ = transformed_environment_components(selected, PRIMARY_DEFINITION, scaling=primary_scaling)
    lookup = selected.set_index("Galaxy")
    rows: list[dict[str, Any]] = []
    for candidate in CANDIDATES:
        controls = assignments.loc[assignments["candidate"] == candidate, "control"].astype(str).tolist()
        candidate_score = float(lookup.loc[candidate, "environment_score"]) if candidate in lookup.index else float("nan")
        control_scores = [float(lookup.loc[name, "environment_score"]) for name in controls if name in lookup.index]
        rows.append(
            {
                "candidate": candidate,
                "candidate_environment_score": candidate_score,
                "matched_controls": ";".join(controls),
                "matched_control_mean_environment_score": float(np.mean(control_scores)) if control_scores else float("nan"),
                "candidate_minus_control_environment_score": candidate_score - float(np.mean(control_scores))
                if control_scores and math.isfinite(candidate_score)
                else float("nan"),
                "interpretation": "descriptive_only_outcome_selected_candidate",
            }
        )
    return pd.DataFrame(rows)


def make_primary_plot(primary_residuals: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    x = primary_residuals["environment_crossfit_residual"].to_numpy(dtype=float)
    y = primary_residuals["outcome_crossfit_residual"].to_numpy(dtype=float)
    ax.scatter(x, y, s=32, color="#16697a", alpha=0.78, edgecolor="white", linewidth=0.4)
    if len(x) >= 2 and np.std(x) > 0:
        slope, intercept = np.polyfit(x, y, 1)
        grid = np.linspace(float(np.min(x)), float(np.max(x)), 100)
        ax.plot(grid, intercept + slope * grid, color="#c44536", linewidth=2.0)
    ax.axhline(0.0, color="#555555", linewidth=0.8)
    ax.axvline(0.0, color="#555555", linewidth=0.8)
    ax.set_xlabel("Environment composite after host controls")
    ax.set_ylabel("PLAMB-RAR outcome after host controls")
    ax.set_title("SPARC FR Environment Test: Primary Development Sample")
    ax.grid(alpha=0.18)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_report(
    path: Path,
    summary: pd.DataFrame,
    gate: dict[str, Any],
    candidate_description: pd.DataFrame,
    profile: pd.DataFrame,
) -> None:
    primary = summary[
        (summary["sample"] == PRIMARY_DEFINITION["name"])
        & (summary["scenario"] == "combined_conventional")
        & (summary["control"] == "actual")
    ].iloc[0]
    shifted = summary[
        (summary["sample"] == PRIMARY_DEFINITION["name"])
        & (summary["scenario"] == "combined_conventional")
        & (summary["control"] == "shifted_ra90")
    ].iloc[0]
    optimiser_failures = int(np.sum(~profile["success"].astype(bool)))
    decision = "PASS" if gate["development_gate_passed"] else "FAIL"
    result_rows = summary[
        [
            "sample",
            "scenario",
            "control",
            "N",
            "partial_spearman_rho",
            "permutation_p_two_sided",
            "cv_rmse_fractional_improvement",
            "environment_coefficient_sign_fraction",
        ]
    ].to_dict("records")
    candidate_rows = candidate_description.to_dict("records") if len(candidate_description) else []
    lines = [
        "# SPARC FR Environmental-Asymmetry Study",
        "",
        f"Date: {RUN_DATE}",
        f"Completed: `{utc_now()}`",
        "",
        "## Bottom Line",
        "",
        f"The preregistered development gate is **{decision}**. The primary charge-blind 2MRS environment association has partial Spearman rho=`{float(primary['partial_spearman_rho']):.6g}`, two-sided permutation p=`{float(primary['permutation_p_two_sided']):.6g}`, and cross-fitted RMSE change=`{100.0 * float(primary['cv_rmse_fractional_improvement']):.3f}` per cent.",
        "",
        "This test cannot identify antimatter. It asks whether ordinary observed environment predicts the residual difference between the current PLAMB bridge and RAR after conventional host controls. A positive association would motivate an environment-conditioned FR response; it would not measure a signed matter/antimatter background.",
        "",
        "## Primary Test",
        "",
        f"The primary sample contains `{int(primary['N'])}` previously unselected, non-reserved Q<=2 galaxies at 10-40 Mpc and |b|>=10 degrees. Neighbours are 2MRS galaxies with M_K<=-22 and |Delta v|<=500 km/s.",
        "",
        f"The RA+90-degree negative control has rho=`{float(shifted['partial_spearman_rho']):.6g}` and p=`{float(shifted['permutation_p_two_sided']):.6g}`.",
        "",
        "## Locked Test Matrix",
        "",
        markdown_table(
            result_rows,
            [
                "sample",
                "scenario",
                "control",
                "N",
                "partial_spearman_rho",
                "permutation_p_two_sided",
                "cv_rmse_fractional_improvement",
                "environment_coefficient_sign_fraction",
            ],
        ),
        "",
        "## Gate Conditions",
        "",
        markdown_table(
            [{"condition": key, "passed": value} for key, value in gate["conditions"].items()],
            ["condition", "passed"],
        ),
        "",
        "Replication was not unsealed. Even a passing development gate requires a separate frozen replication programme.",
        "",
        "## Previously Selected Galaxies",
        "",
        "The following environment comparison is descriptive only because the six galaxies were selected previously from model outcomes. No candidate outcome enters the primary test.",
        "",
        markdown_table(
            candidate_rows,
            [
                "candidate",
                "candidate_environment_score",
                "matched_control_mean_environment_score",
                "candidate_minus_control_environment_score",
            ],
        ) if candidate_rows else "No candidate comparison was available.",
        "",
        "## Data and Model Boundary",
        "",
        f"2MRS is a K_s<=11.75 near-infrared redshift catalogue with 97.6 per cent redshift completeness over 91 per cent of the sky. [2MRS catalogue]({TWOMRS_LANDING}); [Huchra et al. 2012]({TWOMRS_PAPER}). Central positions and redshifts were resolved independently through [NASA/IPAC NED]({NED_GUIDE}).",
        "",
        "The proxy is charge-blind: it traces luminous galaxy environment, not matter-minus-antimatter sign. Its principal limitations are redshift-space projection, 2MRS luminosity selection, peculiar velocities in the nearby volume and the use of K-band luminosity as a mass tracer.",
        "",
        f"Deterministic profile optimiser failures: `{optimiser_failures}` out of `{len(profile)}` model-galaxy-scenario fits.",
        "",
        "## Locked Claim Boundary",
        "",
        "- Do not label any galaxy as antimatter from this study.",
        "- Do not reinterpret a charge-blind density association as a signed FR background.",
        "- Do not inspect the five reserved replication outcomes under this programme.",
        "- If the gate fails, retain conventional environment only as a possible nuisance/control path.",
        "- If the gate passes, preregister an independent replication before opening any reserved outcome.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def export_analysis(out_dir: Path, export_dir: Path, programme: Path) -> None:
    export_dir.mkdir(parents=True, exist_ok=True)
    names = [
        "sparc_fr_environment_acquisition_audit.json",
        "sparc_fr_environment_covariates_nonreserved.csv",
        "sparc_fr_environment_preregistration.json",
        "sparc_fr_environment_preregistration.md",
        "sparc_fr_environment_profile_fits.csv",
        "sparc_fr_environment_development_outcomes.csv",
        "sparc_fr_environment_test_summary.csv",
        "sparc_fr_environment_fold_coefficients.csv",
        "sparc_fr_environment_test_residuals.csv",
        "sparc_fr_environment_primary_partial_residuals.csv",
        "sparc_fr_environment_candidate_descriptive.csv",
        "sparc_fr_environment_gate.json",
        "sparc_fr_environment_report.md",
        "sparc_fr_environment_primary_partial.png",
    ]
    for name in names:
        source = out_dir / name
        if source.exists():
            shutil.copy2(source, export_dir / name)
    rows: list[dict[str, Any]] = []
    for path in sorted(export_dir.iterdir()):
        if path.name == "manifest.csv" or not path.is_file():
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
    rows.extend(
        [
            {
                "path": str(programme.relative_to(ROOT)).replace("\\", "/"),
                "sha256": sha256(programme),
                "bytes": programme.stat().st_size,
                "role": "locked programme",
                "tracked_in_git": True,
            },
            {
                "path": str(TWOMRS),
                "sha256": sha256(TWOMRS),
                "bytes": TWOMRS.stat().st_size,
                "role": "external 2MRS catalogue; URL in preregistration",
                "tracked_in_git": False,
            },
            {
                "path": str(POSITION_CACHE),
                "sha256": sha256(POSITION_CACHE),
                "bytes": POSITION_CACHE.stat().st_size,
                "role": "external NED resolution cache excluding reserved galaxies",
                "tracked_in_git": False,
            },
        ]
    )
    write_csv(export_dir / "manifest.csv", rows)


def analyse_phase(out_dir: Path, export_dir: Path, programme: Path) -> None:
    prereg_path = out_dir / "sparc_fr_environment_preregistration.json"
    environment_path = out_dir / "sparc_fr_environment_covariates_nonreserved.csv"
    if not prereg_path.exists():
        raise FileNotFoundError("Run and commit --phase preregister before analysis")
    config = json.loads(prereg_path.read_text(encoding="utf-8"))
    verify_locks(config, programme, environment_path)
    environment = pd.read_csv(environment_path)
    names = development_names(environment)
    print(f"Locked development outcomes to profile: {len(names)} galaxies", flush=True)
    profile, outcomes = profile_development_outcomes(names)
    catalogue = host_catalogue()
    summary, folds, residuals, primary_frame, gate, primary_scaling = evaluate_tests(
        environment,
        catalogue,
        outcomes,
    )
    candidate_description = candidate_environment_description(environment, primary_scaling)

    write_csv(out_dir / "sparc_fr_environment_profile_fits.csv", profile)
    write_csv(out_dir / "sparc_fr_environment_development_outcomes.csv", outcomes)
    write_csv(out_dir / "sparc_fr_environment_test_summary.csv", summary)
    write_csv(out_dir / "sparc_fr_environment_fold_coefficients.csv", folds)
    write_csv(out_dir / "sparc_fr_environment_test_residuals.csv", residuals)
    primary_label = f"{PRIMARY_DEFINITION['name']}__combined_conventional"
    primary_residuals = residuals[residuals["test"] == primary_label].copy()
    write_csv(out_dir / "sparc_fr_environment_primary_partial_residuals.csv", primary_residuals)
    write_csv(out_dir / "sparc_fr_environment_candidate_descriptive.csv", candidate_description)
    (out_dir / "sparc_fr_environment_gate.json").write_text(
        json.dumps(gate, indent=2) + "\n", encoding="utf-8"
    )
    make_primary_plot(primary_residuals, out_dir / "sparc_fr_environment_primary_partial.png")
    write_report(
        out_dir / "sparc_fr_environment_report.md",
        summary,
        gate,
        candidate_description,
        profile,
    )
    export_analysis(out_dir, export_dir, programme)
    print(f"Saved report: {out_dir / 'sparc_fr_environment_report.md'}")
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
