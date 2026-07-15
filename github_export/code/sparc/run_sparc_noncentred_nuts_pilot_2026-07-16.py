#!/usr/bin/env python3
r"""Run a preregistered non-centred NUTS pilot for filtered all-Q2 SPARC.

This is a sampler-geometry pilot, not a model-comparison result. It keeps the
July 14 likelihood, priors, H0 choice, data filter, and physical model formulas
fixed while replacing the blocked random-walk sampler with PyMC plus nutpie
NUTS. Galaxy distance and stellar-M/L nuisances are sampled as standardised
latent variables and transformed into the original coordinates.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

os.environ.setdefault("PYTENSOR_FLAGS", "cxx=")

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm
import pytensor
import pytensor.tensor as pt
import nutpie

import fit_sparc_hierarchical_map as fit
import sample_sparc_hierarchical_posterior as legacy
from diagnose_plamb_sparc_h0_nuisance import acceleration_cH0_over_2pi
from diagnose_plamb_sparc_rotation import KPC_M, Dataset


ROOT = Path(__file__).resolve().parents[3]
TASK_ROOT = Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu")
OUTPUTS = TASK_ROOT / "outputs"
FILTERED_INPUT_DIR = TASK_ROOT / "work" / "sparc_filtered_persistent_negative"
DEFAULT_SAMPLE = FILTERED_INPUT_DIR / "sparc_galaxy_sample_without_persistent_negative_high_leverage.csv"
DEFAULT_POINTS = FILTERED_INPUT_DIR / "sparc_rotation_points_without_persistent_negative_high_leverage.csv"
DEFAULT_MAP_DIR = (
    ROOT
    / "plamb_runs"
    / "diagnostics"
    / "sparc_hierarchical_map"
    / "optical_depth_minus_persistent_negative_20260714"
)
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "sparc_noncentred_nuts_pilot_20260716"
SUPPORTED_MODELS = ("RAR", "PLAMB_OPTICAL_DEPTH_KAPPA_P")
RUN_DATE = "2026-07-16"


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


def subset_dataset(data: Dataset, mask: np.ndarray) -> Dataset:
    return Dataset(
        galaxy=data.galaxy[mask],
        rad_kpc=data.rad_kpc[mask],
        vobs=data.vobs[mask],
        sigma_v=data.sigma_v[mask],
        vgas=data.vgas[mask],
        vdisk=data.vdisk[mask],
        vbul=data.vbul[mask],
        quality=data.quality[mask],
        distance_method=data.distance_method[mask],
        inc_deg=data.inc_deg[mask],
        d_mpc=data.d_mpc[mask],
    )


def filtered_prepared(prep: fit.Prepared, selected: list[str]) -> fit.Prepared:
    selected = sorted(selected)
    keep = set(selected)
    mask = np.asarray([str(galaxy) in keep for galaxy in prep.data.galaxy], dtype=bool)
    data = subset_dataset(prep.data, mask)
    old_index = {galaxy: index for index, galaxy in enumerate(prep.galaxies)}
    selected_index = np.asarray([old_index[galaxy] for galaxy in selected], dtype=int)
    new_index = {galaxy: index for index, galaxy in enumerate(selected)}
    gal_idx = np.asarray([new_index[str(galaxy)] for galaxy in data.galaxy], dtype=int)
    return fit.Prepared(
        data=data,
        galaxies=selected,
        gal_idx=gal_idx,
        sigma_logd=prep.sigma_logd[selected_index],
        mu_logd_by_h0={h0: values[selected_index] for h0, values in prep.mu_logd_by_h0.items()},
        sigma_ln_ml=prep.sigma_ln_ml,
    )


def stratum_allocation(counts: dict[str, int], total: int) -> dict[str, int]:
    if total < len(counts):
        raise ValueError(f"pilot-galaxies={total} is smaller than the {len(counts)} populated strata")
    grand_total = sum(counts.values())
    ideal = {key: total * count / grand_total for key, count in counts.items()}
    allocation = {key: min(counts[key], max(1, int(math.floor(ideal[key])))) for key in counts}
    while sum(allocation.values()) > total:
        candidates = [key for key in counts if allocation[key] > 1]
        key = max(candidates, key=lambda item: (allocation[item] - ideal[item], item))
        allocation[key] -= 1
    while sum(allocation.values()) < total:
        candidates = [key for key in counts if allocation[key] < counts[key]]
        if not candidates:
            raise ValueError("Could not allocate the requested pilot subset")
        key = max(candidates, key=lambda item: (ideal[item] - allocation[item], counts[item], item))
        allocation[key] += 1
    return allocation


def select_pilot_galaxies(prep: fit.Prepared, total: int, seed: int) -> tuple[list[str], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    for galaxy in prep.galaxies:
        mask = prep.data.galaxy == galaxy
        distance_method = int(round(float(np.median(prep.data.distance_method[mask]))))
        inclination = float(np.median(prep.data.inc_deg[mask]))
        distance_group = "hubble" if distance_method == 1 else "non_hubble"
        inclination_group = "high_inc" if inclination >= 60.0 else "low_inc"
        stratum = f"{distance_group}__{inclination_group}"
        selection_hash = hashlib.sha256(f"{seed}|{galaxy}".encode("utf-8")).hexdigest()
        rows.append(
            {
                "galaxy": galaxy,
                "stratum": stratum,
                "distance_method": distance_method,
                "median_inclination_deg": inclination,
                "n_points": int(np.sum(mask)),
                "selection_hash": selection_hash,
            }
        )
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["stratum"]] = counts.get(row["stratum"], 0) + 1
    allocation = stratum_allocation(counts, total)
    selected: set[str] = set()
    for stratum, number in sorted(allocation.items()):
        candidates = sorted(
            (row for row in rows if row["stratum"] == stratum),
            key=lambda row: (row["selection_hash"], row["galaxy"]),
        )
        selected.update(row["galaxy"] for row in candidates[:number])
    for row in rows:
        row["stratum_population"] = counts[row["stratum"]]
        row["stratum_pilot_allocation"] = allocation[row["stratum"]]
        row["selected"] = row["galaxy"] in selected
    return sorted(selected), rows


def parse_models(text: str) -> list[str]:
    models = [part.strip() for part in text.split(",") if part.strip()]
    unknown = [model for model in models if model not in SUPPORTED_MODELS]
    if unknown:
        raise ValueError(f"Unsupported models: {unknown}")
    return models


def pytensor_velocity_prediction(
    radius_kpc: Any,
    vgas_data: Any,
    vdisk_data: Any,
    vbul_data: Any,
    gal_idx: Any,
    logd: Any,
    logeta: Any,
    ydisk: Any,
    ybul: Any,
    model_name: str,
    h0: float,
    log10_gdagger: Any | None = None,
    log10_kappa: Any | None = None,
    bridge_exponent: Any | None = None,
) -> Any:
    distance_scale = pt.exp(logd[gal_idx])
    stellar_scale = pt.exp(logeta[gal_idx])
    root_distance = pt.sqrt(distance_scale)
    root_stellar = pt.sqrt(stellar_scale)
    radius_scaled = radius_kpc * distance_scale
    vgas = vgas_data * root_distance
    vdisk = vdisk_data * root_distance * root_stellar
    vbul = vbul_data * root_distance * root_stellar
    v2_bar = vgas * pt.abs(vgas) + ydisk * vdisk**2 + ybul * vbul**2
    gbar = pt.maximum(v2_bar, 1e-18) * 1.0e6 / (radius_scaled * KPC_M)

    if model_name == "RAR":
        if log10_gdagger is None:
            raise ValueError("RAR requires log10_gdagger")
        tau = pt.sqrt(pt.maximum(gbar / (10.0**log10_gdagger), 1e-30))
    elif model_name == "PLAMB_OPTICAL_DEPTH_KAPPA_P":
        if log10_kappa is None or bridge_exponent is None:
            raise ValueError("PLAMB optical-depth kappa,p requires both global parameters")
        g0 = (10.0**log10_kappa) * acceleration_cH0_over_2pi(h0)
        tau = pt.maximum(gbar / g0, 1e-30) ** bridge_exponent
    else:  # pragma: no cover - guarded by parse_models
        raise ValueError(model_name)

    gpred = gbar / pt.maximum(1.0 - pt.exp(-tau), 1e-12)
    return pt.sqrt(pt.maximum(gpred, 0.0) * radius_scaled * KPC_M) / 1000.0


def build_model(
    prep: fit.Prepared,
    model_name: str,
    h0: float,
    sigma_global_ml: float,
) -> pm.Model:
    n_galaxies = len(prep.galaxies)
    n_observations = len(prep.data.vobs)
    mu_logd = prep.mu_logd_by_h0[h0]
    log_ydisk_mu = math.log(0.5)
    log_ybul_mu = math.log(0.7)
    z_ydisk_lower = (math.log(0.05) - log_ydisk_mu) / sigma_global_ml
    z_ydisk_upper = (math.log(1.50) - log_ydisk_mu) / sigma_global_ml
    z_ybul_lower = (math.log(0.05) - log_ybul_mu) / sigma_global_ml
    z_ybul_upper = (math.log(1.80) - log_ybul_mu) / sigma_global_ml

    coords = {
        # Numeric coordinates avoid pandas Arrow-string metadata incompatibility
        # in nutpie; names remain locked in the subset and preregistration files.
        "galaxy": np.arange(n_galaxies, dtype=int),
        "observation": np.arange(n_observations, dtype=int),
    }
    with pm.Model(coords=coords) as model:
        gal_idx = pm.Data("galaxy_index", prep.gal_idx, dims="observation")
        radius_kpc = pm.Data("radius_kpc", prep.data.rad_kpc, dims="observation")
        vgas_data = pm.Data("vgas_data", prep.data.vgas, dims="observation")
        vdisk_data = pm.Data("vdisk_data", prep.data.vdisk, dims="observation")
        vbul_data = pm.Data("vbul_data", prep.data.vbul, dims="observation")
        sigma_v = pm.Data("sigma_v", prep.data.sigma_v, dims="observation")
        mu_logd_data = pm.Data("mu_logd", mu_logd, dims="galaxy")
        sigma_logd_data = pm.Data("sigma_logd", prep.sigma_logd, dims="galaxy")

        z_ydisk = pm.TruncatedNormal(
            "z_ydisk",
            mu=0.0,
            sigma=1.0,
            lower=z_ydisk_lower,
            upper=z_ydisk_upper,
        )
        z_ybul = pm.TruncatedNormal(
            "z_ybul",
            mu=0.0,
            sigma=1.0,
            lower=z_ybul_lower,
            upper=z_ybul_upper,
        )
        log_ydisk = pm.Deterministic("log_ydisk", log_ydisk_mu + sigma_global_ml * z_ydisk)
        log_ybul = pm.Deterministic("log_ybul", log_ybul_mu + sigma_global_ml * z_ybul)
        ydisk = pm.Deterministic("ydisk", pt.exp(log_ydisk))
        ybul = pm.Deterministic("ybul", pt.exp(log_ybul))

        z_logd = pm.TruncatedNormal(
            "z_logd",
            mu=0.0,
            sigma=1.0,
            lower=-4.0,
            upper=4.0,
            dims="galaxy",
        )
        z_logeta = pm.TruncatedNormal(
            "z_logeta",
            mu=0.0,
            sigma=1.0,
            lower=-4.0,
            upper=4.0,
            dims="galaxy",
        )
        logd = pm.Deterministic("logd", mu_logd_data + sigma_logd_data * z_logd, dims="galaxy")
        logeta = pm.Deterministic("logeta", prep.sigma_ln_ml * z_logeta, dims="galaxy")

        if model_name == "RAR":
            log10_gdagger = pm.Uniform("log10_gdagger", lower=-11.6, upper=-9.2)
            velocity_pred = pytensor_velocity_prediction(
                radius_kpc,
                vgas_data,
                vdisk_data,
                vbul_data,
                gal_idx,
                logd,
                logeta,
                ydisk,
                ybul,
                model_name,
                h0,
                log10_gdagger=log10_gdagger,
            )
        elif model_name == "PLAMB_OPTICAL_DEPTH_KAPPA_P":
            log10_kappa = pm.Uniform("log10_kappa", lower=-0.50, upper=0.50)
            bridge_exponent = pm.Uniform("bridge_exponent", lower=0.15, upper=1.20)
            velocity_pred = pytensor_velocity_prediction(
                radius_kpc,
                vgas_data,
                vdisk_data,
                vbul_data,
                gal_idx,
                logd,
                logeta,
                ydisk,
                ybul,
                model_name,
                h0,
                log10_kappa=log10_kappa,
                bridge_exponent=bridge_exponent,
            )
        else:  # pragma: no cover - guarded by parse_models
            raise ValueError(model_name)
        pm.Normal("vobs", mu=velocity_pred, sigma=sigma_v, observed=prep.data.vobs, dims="observation")

    return model


def initial_legacy_state(prep: fit.Prepared, model_name: str, h0: float) -> np.ndarray:
    params = [math.log(0.5), math.log(0.7)]
    if model_name == "RAR":
        params.append((-11.6 + -9.2) / 2.0)
    elif model_name == "PLAMB_OPTICAL_DEPTH_KAPPA_P":
        params.extend([0.0, (0.15 + 1.20) / 2.0])
    params.extend(prep.mu_logd_by_h0[h0].tolist())
    params.extend([0.0] * len(prep.galaxies))
    return np.asarray(params, dtype=float)


def equivalence_check(
    prep: fit.Prepared,
    model_name: str,
    h0: float,
    sigma_global_ml: float,
) -> dict[str, Any]:
    state = initial_legacy_state(prep, model_name, h0)
    ydisk = pt.dscalar("equivalence_ydisk")
    ybul = pt.dscalar("equivalence_ybul")
    logd = pt.dvector("equivalence_logd")
    logeta = pt.dvector("equivalence_logeta")
    fixed = [
        pt.as_tensor_variable(prep.data.rad_kpc),
        pt.as_tensor_variable(prep.data.vgas),
        pt.as_tensor_variable(prep.data.vdisk),
        pt.as_tensor_variable(prep.data.vbul),
        pt.as_tensor_variable(prep.gal_idx),
        logd,
        logeta,
        ydisk,
        ybul,
        model_name,
        h0,
    ]
    if model_name == "RAR":
        global_parameter = pt.dscalar("equivalence_log10_gdagger")
        graph_expression = pytensor_velocity_prediction(*fixed, log10_gdagger=global_parameter)
        graph_predict = pytensor.function([ydisk, ybul, logd, logeta, global_parameter], graph_expression)
        graph_velocity = np.asarray(
            graph_predict(0.5, 0.7, prep.mu_logd_by_h0[h0], np.zeros(len(prep.galaxies)), state[2]),
            dtype=float,
        )
    else:
        log10_kappa = pt.dscalar("equivalence_log10_kappa")
        bridge_exponent = pt.dscalar("equivalence_bridge_exponent")
        graph_expression = pytensor_velocity_prediction(
            *fixed,
            log10_kappa=log10_kappa,
            bridge_exponent=bridge_exponent,
        )
        graph_predict = pytensor.function(
            [ydisk, ybul, logd, logeta, log10_kappa, bridge_exponent],
            graph_expression,
        )
        graph_velocity = np.asarray(
            graph_predict(
                0.5,
                0.7,
                prep.mu_logd_by_h0[h0],
                np.zeros(len(prep.galaxies)),
                state[2],
                state[3],
            ),
            dtype=float,
        )
    legacy_velocity, _meta = fit.predict(state, model_name, prep, h0, 0.0, 0.0)
    legacy_total, legacy_chi2, legacy_prior, _meta = fit.objective(
        state,
        model_name,
        prep,
        h0,
        0.0,
        0.0,
        sigma_global_ml,
    )
    graph_chi2 = float(np.sum(((prep.data.vobs - graph_velocity) / prep.data.sigma_v) ** 2))
    max_abs = float(np.max(np.abs(graph_velocity - legacy_velocity)))
    max_rel = float(np.max(np.abs(graph_velocity - legacy_velocity) / np.maximum(np.abs(legacy_velocity), 1e-12)))
    chi2_abs = abs(graph_chi2 - legacy_chi2)
    passed = max_abs <= 1e-8 and max_rel <= 1e-10 and chi2_abs <= 1e-7
    return {
        "model": model_name,
        "max_abs_prediction_difference_kms": max_abs,
        "max_relative_prediction_difference": max_rel,
        "graph_chi2_data": graph_chi2,
        "legacy_chi2_data": legacy_chi2,
        "chi2_data_absolute_difference": chi2_abs,
        "legacy_chi2_prior": legacy_prior,
        "legacy_total_objective": legacy_total,
        "pass": passed,
    }


def gate_variable_names(model_name: str) -> list[str]:
    names = ["log_ydisk", "log_ybul", "logd", "logeta"]
    if model_name == "RAR":
        names.append("log10_gdagger")
    else:
        names.extend(["log10_kappa", "bridge_exponent"])
    return names


def data_variable(idata: az.InferenceData, names: list[str]) -> Any | None:
    for name in names:
        if name in idata.sample_stats:
            return idata.sample_stats[name]
    return None


def serialisable_attribute(value: Any) -> Any:
    allowed = (str, bytes, int, float, np.number, np.ndarray)
    if isinstance(value, allowed):
        return value
    if isinstance(value, (list, tuple)) and all(isinstance(item, allowed) for item in value):
        return value
    return json.dumps(value, sort_keys=True, default=str)


def sanitise_netcdf_attributes(idata: az.InferenceData) -> None:
    containers: list[Any] = [idata]
    containers.extend(getattr(idata, group) for group in idata.groups())
    for container in containers:
        attributes = container.attrs
        for key, value in list(attributes.items()):
            attributes[key] = serialisable_attribute(value)


def flattened_dataset_values(dataset: Any) -> np.ndarray:
    arrays = [np.asarray(variable.values, dtype=float).ravel() for variable in dataset.data_vars.values()]
    return np.concatenate(arrays) if arrays else np.asarray([], dtype=float)


def diagnose(
    idata: az.InferenceData,
    model_name: str,
    args: argparse.Namespace,
    elapsed_seconds: float,
    equivalence_pass: bool,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    gate_variables = gate_variable_names(model_name)
    summary = az.summary(idata, var_names=gate_variables, kind="all", round_to="none")
    max_rhat = float(np.max(flattened_dataset_values(az.rhat(idata, var_names=gate_variables, method="rank"))))
    min_ess_bulk = float(np.min(flattened_dataset_values(az.ess(idata, var_names=gate_variables, method="bulk"))))
    min_ess_tail = float(np.min(flattened_dataset_values(az.ess(idata, var_names=gate_variables, method="tail"))))
    bfmi_values = np.asarray(az.bfmi(idata), dtype=float)
    minimum_bfmi = float(np.min(bfmi_values))
    divergence_data = data_variable(idata, ["diverging", "divergence"])
    divergences = int(np.asarray(divergence_data).sum()) if divergence_data is not None else 0
    tree_depth_data = data_variable(idata, ["tree_depth", "depth"])
    maximum_tree_depth = int(np.asarray(tree_depth_data).max()) if tree_depth_data is not None else -1
    max_depth_hits = (
        int(np.sum(np.asarray(tree_depth_data) >= args.max_tree_depth)) if tree_depth_data is not None else -1
    )

    strict_rhat_pass = max_rhat <= args.rhat_target
    ceiling_rhat_pass = max_rhat <= args.rhat_ceiling
    ess_pass = min_ess_bulk >= args.min_ess and min_ess_tail >= args.min_ess
    divergence_pass = divergences == 0
    energy_pass = minimum_bfmi >= args.min_bfmi
    strict_pass = equivalence_pass and strict_rhat_pass and ess_pass and divergence_pass and energy_pass
    conditional_pass = equivalence_pass and ceiling_rhat_pass and ess_pass and divergence_pass and energy_pass
    if strict_pass:
        status = "PASS"
    elif conditional_pass:
        status = "CONDITIONAL_EXTEND_PILOT"
    else:
        status = "FAIL"
    gate = {
        "model": model_name,
        "status": status,
        "ready_for_full_run": strict_pass,
        "likelihood_equivalence_pass": equivalence_pass,
        "max_rhat": max_rhat,
        "rhat_target": args.rhat_target,
        "rhat_absolute_ceiling": args.rhat_ceiling,
        "strict_rhat_pass": strict_rhat_pass,
        "ceiling_rhat_pass": ceiling_rhat_pass,
        "min_bulk_ess": min_ess_bulk,
        "min_tail_ess": min_ess_tail,
        "minimum_ess_required": args.min_ess,
        "ess_pass": ess_pass,
        "divergences": divergences,
        "divergence_pass": divergence_pass,
        "minimum_bfmi": minimum_bfmi,
        "minimum_bfmi_required": args.min_bfmi,
        "energy_pass": energy_pass,
        "maximum_tree_depth_observed": maximum_tree_depth,
        "configured_max_tree_depth": args.max_tree_depth,
        "max_tree_depth_hits": max_depth_hits,
        "elapsed_seconds": elapsed_seconds,
        "chains": args.chains,
        "tune_per_chain": args.tune,
        "draws_per_chain": args.draws,
    }
    summary = summary.reset_index().rename(columns={"index": "parameter"})
    return summary, gate


def preregistration_markdown(config: dict[str, Any]) -> str:
    gates = config["gates"]
    selected = config["selected_galaxies"]
    return "\n".join(
        [
            "# SPARC Non-centred NUTS Pilot Preregistration",
            "",
            f"**Written before sampling:** `{config['written_before_sampling_utc']}`  ",
            f"**Selected galaxies:** `{len(selected)}`  ",
            f"**Models:** `{', '.join(config['models'])}`",
            "",
            "## Fixed Scope",
            "",
            "This is a sampler-geometry pilot, not a PLAMB-versus-RAR evidence calculation. The filtered all-Q2 inputs, physical likelihood, fixed H0, and July prior scales are unchanged.",
            "",
            "The nuisance coordinates are transformed as",
            "",
            r"\[",
            r"\begin{aligned}",
            r"z_{D,g}       &\sim \mathcal{N}(0,1), & \log D_g       &= \mu_{D,g} + \sigma_{D,g} z_{D,g}, \\",
            r"z_{\eta,g}    &\sim \mathcal{N}(0,1), & \log \eta_g   &= \sigma_{\log\eta} z_{\eta,g}, \\",
            r"z_{\Upsilon d}&\sim \mathcal{N}(0,1), & \log\Upsilon_d&= \log(0.5)+\sigma_{\Upsilon}z_{\Upsilon d}, \\",
            r"z_{\Upsilon b}&\sim \mathcal{N}(0,1), & \log\Upsilon_b&= \log(0.7)+\sigma_{\Upsilon}z_{\Upsilon b}.",
            r"\end{aligned}",
            r"\]",
            "",
            "The same four-standard-deviation nuisance truncation and original global physical bounds are retained.",
            "",
            "## Outcome-blind Subset",
            "",
            f"The pilot subset is selected by SHA-256 ranking with seed `{config['selection_seed']}`, stratified only by published SPARC distance method and inclination below/above 60 degrees. No PLAMB-RAR residual or fit score enters selection.",
            "",
            "Selected galaxies:",
            "",
            ", ".join(f"`{galaxy}`" for galaxy in selected),
            "",
            "## Locked Gates",
            "",
            f"- likelihood-equivalence check must pass at absolute velocity tolerance `{gates['prediction_abs_tolerance_kms']}` km/s;",
            "- zero post-warmup divergences;",
            f"- rank-normalised split R-hat target `<= {gates['rhat_target']}` and absolute ceiling `<= {gates['rhat_absolute_ceiling']}`;",
            f"- minimum bulk and tail ESS `>= {gates['minimum_bulk_and_tail_ess']}`; and",
            f"- minimum chain E-BFMI `>= {gates['minimum_bfmi']}`.",
            "",
            "A model between the R-hat target and ceiling is labelled `CONDITIONAL_EXTEND_PILOT`, not ready for the full run. Any failure of equivalence, divergences, ESS, energy, or the R-hat ceiling is a pilot failure.",
            "",
        ]
    )


def write_report(
    path: Path,
    config: dict[str, Any],
    equivalence_rows: list[dict[str, Any]],
    gates: list[dict[str, Any]],
    global_rows: list[dict[str, Any]],
) -> None:
    lines = [
        "# SPARC Non-centred NUTS Pilot Report",
        "",
        f"**Completed:** `{utc_now()}`  ",
        f"**Pilot galaxies:** `{len(config['selected_galaxies'])}`  ",
        f"**Pilot points:** `{config['pilot_points']}`  ",
        f"**Chains:** `{config['sampling']['chains']}`",
        "",
        "## Decision",
        "",
    ]
    if all(bool(gate["ready_for_full_run"]) for gate in gates):
        lines.append("Both models pass every strict preregistered geometry gate. The implementation is ready for a separately registered full filtered all-Q2 run.")
    elif all(gate["status"] != "FAIL" for gate in gates):
        lines.append("No model fails the absolute gates, but at least one misses the strict R-hat target. Extend the pilot before authorising a full run.")
    else:
        lines.append("At least one model fails a preregistered pilot gate. Do not launch the full run until the failure is understood and the remedy is registered.")
    lines.extend(
        [
            "",
            "This pilot assesses sampler geometry only. It is not held-out prediction and must not be used as a PLAMB-versus-RAR evidence statement.",
            "",
            "## Gate Results",
            "",
            "| Model | Status | Max R-hat | Min bulk ESS | Min tail ESS | Divergences | Min E-BFMI |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for gate in gates:
        lines.append(
            f"| {gate['model']} | {gate['status']} | {gate['max_rhat']:.5f} | "
            f"{gate['min_bulk_ess']:.1f} | {gate['min_tail_ess']:.1f} | {gate['divergences']} | "
            f"{gate['minimum_bfmi']:.4f} |"
        )
    lines.extend(
        [
            "",
            "## Likelihood Equivalence",
            "",
            "| Model | Max velocity difference (km/s) | Chi-squared difference | Pass |",
            "| --- | ---: | ---: | --- |",
        ]
    )
    for row in equivalence_rows:
        lines.append(
            f"| {row['model']} | {float(row['max_abs_prediction_difference_kms']):.3e} | "
            f"{float(row['chi2_data_absolute_difference']):.3e} | {row['pass']} |"
        )
    lines.extend(
        [
            "",
            "## Global Posterior Readout",
            "",
            "These values are subset-level sampler checks, not full-sample estimates.",
            "",
            "| Model | Parameter | Mean | SD | 3% HDI | 97% HDI | R-hat | Bulk ESS | Tail ESS |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in global_rows:
        lines.append(
            f"| {row['model']} | {row['parameter']} | {row['mean']:.6g} | {row['sd']:.6g} | "
            f"{row['hdi_3%']:.6g} | {row['hdi_97%']:.6g} | {row['r_hat']:.5f} | "
            f"{row['ess_bulk']:.1f} | {row['ess_tail']:.1f} |"
        )
    lines.extend(
        [
            "",
            "## Reproducibility",
            "",
            f"- PyMC: `{config['environment']['pymc']}`",
            f"- nutpie: `{config['environment']['nutpie']}`",
            f"- ArviZ: `{config['environment']['arviz']}`",
            f"- Python: `{config['environment']['python']}`",
            f"- selection seed: `{config['selection_seed']}`",
            f"- sampling seed: `{config['sampling']['seed']}`",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def global_summary_rows(model_name: str, idata: az.InferenceData) -> list[dict[str, Any]]:
    variables = ["ydisk", "ybul"]
    if model_name == "RAR":
        variables.append("log10_gdagger")
    else:
        variables.extend(["log10_kappa", "bridge_exponent"])
    summary = az.summary(idata, var_names=variables, kind="all", round_to="none")
    rows: list[dict[str, Any]] = []
    for parameter, row in summary.iterrows():
        rows.append(
            {
                "model": model_name,
                "parameter": parameter,
                "mean": float(row["mean"]),
                "sd": float(row["sd"]),
                "hdi_3%": float(row["hdi_3%"]),
                "hdi_97%": float(row["hdi_97%"]),
                "r_hat": float(row["r_hat"]),
                "ess_bulk": float(row["ess_bulk"]),
                "ess_tail": float(row["ess_tail"]),
            }
        )
    return rows


def analyse_existing_results(args: argparse.Namespace, models: list[str]) -> int:
    prereg_json = args.out_dir / "sparc_noncentred_nuts_pilot_preregistration.json"
    equivalence_path = args.out_dir / "sparc_noncentred_nuts_likelihood_equivalence.csv"
    if not prereg_json.exists() or not equivalence_path.exists():
        raise FileNotFoundError("Analysis-only mode requires the existing preregistration and equivalence files")
    prereg = json.loads(prereg_json.read_text(encoding="utf-8"))
    registered_sampling = prereg["sampling"]
    args.chains = int(registered_sampling["chains"])
    args.tune = int(registered_sampling["tune_per_chain"])
    args.draws = int(registered_sampling["draws_per_chain"])
    args.max_tree_depth = int(registered_sampling["max_tree_depth"])
    registered_gates = prereg["gates"]
    args.rhat_target = float(registered_gates["rhat_target"])
    args.rhat_ceiling = float(registered_gates["rhat_absolute_ceiling"])
    args.min_ess = float(registered_gates["minimum_bulk_and_tail_ess"])
    args.min_bfmi = float(registered_gates["minimum_bfmi"])
    equivalence_rows = pd.read_csv(equivalence_path).to_dict(orient="records")
    equivalence_by_model = {row["model"]: row for row in equivalence_rows}
    gates: list[dict[str, Any]] = []
    global_rows: list[dict[str, Any]] = []
    for model_name in models:
        trace_path = args.out_dir / f"sparc_noncentred_nuts_{model_name}_trace.nc"
        if not trace_path.exists():
            raise FileNotFoundError(trace_path)
        idata = az.from_netcdf(trace_path)
        old_gate_path = args.out_dir / f"sparc_noncentred_nuts_{model_name}_gate.json"
        old_gate = json.loads(old_gate_path.read_text(encoding="utf-8")) if old_gate_path.exists() else {}
        summary, gate = diagnose(
            idata,
            model_name,
            args,
            float(old_gate.get("elapsed_seconds", 0.0)),
            str(equivalence_by_model[model_name]["pass"]).strip().lower() in {"true", "1", "yes"},
        )
        gate["compile_seconds"] = float(old_gate.get("compile_seconds", 0.0))
        summary.to_csv(args.out_dir / f"sparc_noncentred_nuts_{model_name}_diagnostics.csv", index=False)
        old_gate_path.write_text(json.dumps(gate, indent=2), encoding="utf-8")
        gates.append(gate)
        global_rows.extend(global_summary_rows(model_name, idata))
    gate_path = args.out_dir / "sparc_noncentred_nuts_pilot_gates.csv"
    globals_path = args.out_dir / "sparc_noncentred_nuts_pilot_global_summary.csv"
    report_path = args.out_dir / "sparc_noncentred_nuts_pilot_report.md"
    write_csv(gate_path, gates)
    write_csv(globals_path, global_rows)
    write_report(report_path, prereg, equivalence_rows, gates, global_rows)
    create_manifest(args.out_dir)
    if args.copy_to_outputs:
        OUTPUTS.mkdir(parents=True, exist_ok=True)
        prereg_md = args.out_dir / "sparc_noncentred_nuts_pilot_preregistration.md"
        subset_path = args.out_dir / "sparc_nuts_pilot_locked_subset.csv"
        for path in [prereg_md, subset_path, equivalence_path, gate_path, globals_path, report_path]:
            shutil.copy2(path, OUTPUTS / f"{path.stem}_{RUN_DATE}{path.suffix}")
    print(f"Reanalysed existing traces: {report_path}", flush=True)
    return 0


def create_manifest(out_dir: Path) -> Path:
    files = sorted(path for path in out_dir.iterdir() if path.is_file() and path.name != "manifest.csv")
    rows = [
        {
            "file": path.name,
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in files
    ]
    manifest = out_dir / "manifest.csv"
    write_csv(manifest, rows)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Non-centred PyMC/nutpie NUTS pilot for filtered all-Q2 SPARC.")
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--points", type=Path, default=DEFAULT_POINTS)
    parser.add_argument("--map-dir", type=Path, default=DEFAULT_MAP_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--models", default=",".join(SUPPORTED_MODELS))
    parser.add_argument("--pilot-galaxies", type=int, default=24)
    parser.add_argument("--selection-seed", type=int, default=160726)
    parser.add_argument("--seed", type=int, default=16072601)
    parser.add_argument("--chains", type=int, default=4)
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--tune", type=int, default=1000)
    parser.add_argument("--draws", type=int, default=1000)
    parser.add_argument("--target-accept", type=float, default=0.90)
    parser.add_argument("--max-tree-depth", type=int, default=12)
    parser.add_argument("--sigma-ln-ml", type=float, default=0.25)
    parser.add_argument("--sigma-global-ml", type=float, default=0.55)
    parser.add_argument("--distance-floor-frac", type=float, default=0.03)
    parser.add_argument("--hubble-prior-mode", choices=["model_h0_rescale", "published"], default="model_h0_rescale")
    parser.add_argument("--err-floor-kms", type=float, default=3.0)
    parser.add_argument("--rhat-target", type=float, default=1.01)
    parser.add_argument("--rhat-ceiling", type=float, default=1.05)
    parser.add_argument("--min-ess", type=float, default=400.0)
    parser.add_argument("--min-bfmi", type=float, default=0.30)
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--analyse-existing", action="store_true")
    parser.add_argument("--copy-to-outputs", action=argparse.BooleanOptionalAction, default=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    models = parse_models(args.models)
    if args.analyse_existing:
        return analyse_existing_results(args, models)
    if args.chains != 4 and not args.preflight_only:
        raise ValueError("The registered pilot requires exactly four chains")
    required = [args.sample, args.points, args.map_dir / "sparc_hierarchical_map_summary.csv"]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required inputs: {missing}")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    map_rows = legacy.load_map_rows(args.map_dir)
    h0_by_model = {model: legacy.parse_float(map_rows[("all_Q2", model)]["H0_trial"]) for model in models}
    if len(set(h0_by_model.values())) != 1:
        raise ValueError(f"The pilot requires a common fixed H0; found {h0_by_model}")
    h0 = next(iter(h0_by_model.values()))
    split_config = next(item for item in fit.split_configs() if item["split"] == "all_Q2")
    full_prep = fit.prepare(
        args.sample,
        args.points,
        int(split_config["quality_max"]),
        str(split_config["distance_method"]),
        args.err_floor_kms,
        [h0],
        args.sigma_ln_ml,
        args.distance_floor_frac,
        args.hubble_prior_mode,
    )
    selected, selection_rows = select_pilot_galaxies(full_prep, args.pilot_galaxies, args.selection_seed)
    prep = filtered_prepared(full_prep, selected)
    subset_path = args.out_dir / "sparc_nuts_pilot_locked_subset.csv"
    write_csv(subset_path, selection_rows)

    environment = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "pymc": pm.__version__,
        "pytensor": __import__("pytensor").__version__,
        "nutpie": nutpie.__version__,
        "arviz": az.__version__,
    }
    prereg = {
        "date": RUN_DATE,
        "written_before_sampling_utc": utc_now(),
        "analysis": "filtered_all_Q2_noncentred_nuts_geometry_pilot",
        "scientific_use": "sampler geometry only; no model promotion or held-out evidence",
        "models": models,
        "sample": str(args.sample),
        "points": str(args.points),
        "map_dir": str(args.map_dir),
        "input_sha256": {
            "sample": sha256(args.sample),
            "points": sha256(args.points),
            "map_summary": sha256(args.map_dir / "sparc_hierarchical_map_summary.csv"),
            "programme": sha256(Path(__file__)),
        },
        "full_filtered_galaxies": len(full_prep.galaxies),
        "full_filtered_points": len(full_prep.data.vobs),
        "selected_galaxies": selected,
        "pilot_points": len(prep.data.vobs),
        "selection_seed": args.selection_seed,
        "selection_rule": "SHA-256 rank within distance-method x inclination(<60/>=60 deg) strata; proportional allocation; no fit residuals",
        "fixed_h0": h0,
        "priors": {
            "log_ydisk": "Normal(log(0.5), 0.55), truncated to log(0.05)..log(1.50)",
            "log_ybul": "Normal(log(0.7), 0.55), truncated to log(0.05)..log(1.80)",
            "distance_nuisance": "Normal(mu_logd(H0), sigma_logd), truncated at +/-4 sigma",
            "stellar_ml_nuisance": "Normal(0, 0.25), truncated at +/-4 sigma",
            "RAR_log10_gdagger": "Uniform(-11.6, -9.2)",
            "PLAMB_log10_kappa": "Uniform(-0.50, 0.50)",
            "PLAMB_bridge_exponent": "Uniform(0.15, 1.20)",
        },
        "sampling": {
            "sampler": "nutpie NUTS compiled from PyMC model",
            "chains": args.chains,
            "cores": args.cores,
            "tune_per_chain": args.tune,
            "draws_per_chain": args.draws,
            "target_accept": args.target_accept,
            "max_tree_depth": args.max_tree_depth,
            "seed": args.seed,
        },
        "gates": {
            "prediction_abs_tolerance_kms": 1e-8,
            "rhat_target": args.rhat_target,
            "rhat_absolute_ceiling": args.rhat_ceiling,
            "minimum_bulk_and_tail_ess": args.min_ess,
            "maximum_divergences": 0,
            "minimum_bfmi": args.min_bfmi,
        },
        "environment": environment,
    }
    prereg_json = args.out_dir / "sparc_noncentred_nuts_pilot_preregistration.json"
    prereg_md = args.out_dir / "sparc_noncentred_nuts_pilot_preregistration.md"
    prereg_json.write_text(json.dumps(prereg, indent=2), encoding="utf-8")
    prereg_md.write_text(preregistration_markdown(prereg), encoding="utf-8")
    print(f"Preregistered {len(selected)} galaxies and {len(prep.data.vobs)} points", flush=True)

    equivalence_rows: list[dict[str, Any]] = []
    built_models: dict[str, pm.Model] = {}
    for model_index, model_name in enumerate(models):
        pymc_model = build_model(prep, model_name, h0, args.sigma_global_ml)
        equivalence = equivalence_check(
            prep,
            model_name,
            h0,
            args.sigma_global_ml,
        )
        equivalence_rows.append(equivalence)
        built_models[model_name] = pymc_model
        print(
            f"{model_name} equivalence max_abs={equivalence['max_abs_prediction_difference_kms']:.3e} "
            f"chi2_diff={equivalence['chi2_data_absolute_difference']:.3e} pass={equivalence['pass']}",
            flush=True,
        )
        if not equivalence["pass"]:
            raise RuntimeError(f"Likelihood-equivalence check failed for {model_name}")
    equivalence_path = args.out_dir / "sparc_noncentred_nuts_likelihood_equivalence.csv"
    write_csv(equivalence_path, equivalence_rows)

    if args.preflight_only:
        create_manifest(args.out_dir)
        print(f"Preflight complete: {args.out_dir}", flush=True)
        return 0

    gates: list[dict[str, Any]] = []
    global_rows: list[dict[str, Any]] = []
    for model_index, model_name in enumerate(models):
        print(
            f"=== Sampling {model_name}: chains={args.chains} tune={args.tune} draws={args.draws} "
            f"galaxies={len(prep.galaxies)} points={len(prep.data.vobs)} ===",
            flush=True,
        )
        compile_start = time.perf_counter()
        compiled = nutpie.compile_pymc_model(built_models[model_name], backend="numba")
        compile_seconds = time.perf_counter() - compile_start
        sample_start = time.perf_counter()
        idata = nutpie.sample(
            compiled,
            draws=args.draws,
            tune=args.tune,
            chains=args.chains,
            cores=args.cores,
            seed=args.seed + model_index * 100003,
            target_accept=args.target_accept,
            maxdepth=args.max_tree_depth,
            save_warmup=True,
            progress_bar=True,
        )
        sample_seconds = time.perf_counter() - sample_start
        trace_path = args.out_dir / f"sparc_noncentred_nuts_{model_name}_trace.nc"
        sanitise_netcdf_attributes(idata)
        idata.to_netcdf(trace_path)
        summary, gate = diagnose(
            idata,
            model_name,
            args,
            sample_seconds,
            equivalence_rows[model_index]["pass"],
        )
        gate["compile_seconds"] = compile_seconds
        summary.to_csv(args.out_dir / f"sparc_noncentred_nuts_{model_name}_diagnostics.csv", index=False)
        (args.out_dir / f"sparc_noncentred_nuts_{model_name}_gate.json").write_text(
            json.dumps(gate, indent=2),
            encoding="utf-8",
        )
        gates.append(gate)
        global_rows.extend(global_summary_rows(model_name, idata))
        print(
            f"{model_name} status={gate['status']} max_rhat={gate['max_rhat']:.5f} "
            f"min_bulk_ess={gate['min_bulk_ess']:.1f} min_tail_ess={gate['min_tail_ess']:.1f} "
            f"divergences={gate['divergences']} min_bfmi={gate['minimum_bfmi']:.4f}",
            flush=True,
        )

    gate_path = args.out_dir / "sparc_noncentred_nuts_pilot_gates.csv"
    globals_path = args.out_dir / "sparc_noncentred_nuts_pilot_global_summary.csv"
    report_path = args.out_dir / "sparc_noncentred_nuts_pilot_report.md"
    write_csv(gate_path, gates)
    write_csv(globals_path, global_rows)
    write_report(report_path, prereg, equivalence_rows, gates, global_rows)
    create_manifest(args.out_dir)

    if args.copy_to_outputs:
        OUTPUTS.mkdir(parents=True, exist_ok=True)
        for path in [prereg_md, subset_path, equivalence_path, gate_path, globals_path, report_path]:
            shutil.copy2(path, OUTPUTS / f"{path.stem}_{RUN_DATE}{path.suffix}")
    print(f"Saved report: {report_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
