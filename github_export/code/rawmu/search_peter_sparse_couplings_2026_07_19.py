#!/usr/bin/env python3
"""Enumerate sparse repairs to the Peter clock/redshift/source closure."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
DATE_TAG = "2026-07-19"
ANALYSIS_ID = "peter_sparse_coupling_search_2026-07-19"
EXPECTED_PROTOCOL_SHA256 = (
    "2a3287de46ea6cda9e4905bb2575ff32f70f20e74124d50767272db0c18b3530"
)
DEFAULT_OUTDIR = (
    REPO_ROOT
    / "github_export"
    / "results"
    / DATE_TAG
    / "peter_sparse_coupling_search"
)
DEFAULT_PROTOCOL = (
    DEFAULT_OUTDIR / f"peter_sparse_coupling_search_protocol_{DATE_TAG}.json"
)
PRIOR_SUMMARY = (
    REPO_ROOT
    / "github_export"
    / "results"
    / DATE_TAG
    / "peter_clock_law_v1"
    / f"peter_clock_law_v1_summary_{DATE_TAG}.json"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_json(value):
    if isinstance(value, dict):
        return {str(key): safe_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe_json(item) for item in value]
    if isinstance(value, np.ndarray):
        return [safe_json(item) for item in value.tolist()]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def endpoint_factor(exponent: float, one_plus_zmax: float) -> float:
    ratio = float(one_plus_zmax**float(exponent))
    return max(ratio, 1.0 / ratio)


def maximum_factor(values: dict[str, float], one_plus_zmax: float) -> float:
    if not values:
        return 1.0
    return max(endpoint_factor(value, one_plus_zmax) for value in values.values())


def verify_inputs(protocol_path: Path, prior_summary_path: Path) -> tuple[dict, dict]:
    if not protocol_path.exists():
        raise FileNotFoundError(protocol_path)
    if sha256_file(protocol_path) != EXPECTED_PROTOCOL_SHA256:
        raise RuntimeError("Frozen protocol hash mismatch")
    if not prior_summary_path.exists():
        raise FileNotFoundError(prior_summary_path)
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    prior = json.loads(prior_summary_path.read_text(encoding="utf-8"))
    known = protocol["known_inputs"]
    checks = {
        "exact_Peter_path_delta_BIC_vs_flat_LCDM": prior["primary"]
        ["PETER_NOLOSS_LCONST"]["delta_BIC_vs_LCDM"],
        "maximum_supernova_redshift": prior["protocol"]["known_before_registration"]
        .get("maximum_supernova_redshift", known["maximum_supernova_redshift"]),
        "minimum_registered_fixed_s_tension_sigma": prior["derived"]
        ["minimum_registered_fixed_s_tension_sigma"],
        "sample_N": prior["protocol"]["known_before_registration"]
        ["fixed_path_sample_N"],
        "supernova_composite_s": prior["primary"]["COMPOSITE_S_FREE"]
        ["parameter_value"],
        "supernova_composite_s_sigma": prior["primary"]["COMPOSITE_S_FREE"]
        ["parameter_sigma"],
        "supernova_free_s_delta_BIC_vs_flat_LCDM": prior["primary"]
        ["COMPOSITE_S_FREE"]["delta_BIC_vs_LCDM"],
    }
    for key, observed in checks.items():
        expected = known[key]
        if not np.isclose(float(observed), float(expected), rtol=0.0, atol=1e-11):
            raise RuntimeError(
                f"Registered prior input changed for {key}: {observed} != {expected}"
            )
    return protocol, prior


def stated_evaluator(
    values: dict[str, float], known: dict, one_plus_zmax: float
) -> dict[str, object]:
    delta_q = values.get("delta_Q_atomic", 0.0)
    sigma = values.get("sigma_spectral", 0.0)
    tau = values.get("tau_photon_distance", 0.0)
    ell = values.get("ell_source", 0.0)
    q = -1.0 + delta_q
    z_power = 1.0 + sigma - q
    s_power = tau - 0.5 * ell
    b_pred = 1.0
    b_sigma = math.hypot(
        known["DES_time_dilation_sigma_stat"],
        known["DES_time_dilation_sigma_sys"],
    )
    return {
        "atomic_Q_final_exponent": q,
        "spectroscopic_redshift_exponent": z_power,
        "supernova_composite_s": s_power,
        "duration_b_predicted": b_pred,
        "duration_tension_sigma": abs(b_pred - known["DES_time_dilation_b"])
        / b_sigma,
        "time_dilation_gate_pass": abs(b_pred - known["DES_time_dilation_b"])
        <= 2.0 * b_sigma,
        "Q_endpoint_ratio": one_plus_zmax**q,
        "spectral_endpoint_ratio": one_plus_zmax**sigma,
        "photon_distance_endpoint_ratio": one_plus_zmax**tau,
        "source_luminosity_endpoint_ratio": one_plus_zmax**ell,
        "maximum_dimensionless_endpoint_change_factor": maximum_factor(
            {"Q": q, "spectral": sigma}, one_plus_zmax
        ),
    }


def microscopic_quantities(values: dict[str, float]) -> dict[str, float]:
    delta_h = values.get("delta_h", 0.0)
    delta_g = values.get("delta_G", 0.0)
    delta_mp = values.get("delta_m_p", 0.0)
    delta_me = values.get("delta_m_e", 0.0)
    delta_alpha = values.get("delta_alpha", 0.0)
    sigma = values.get("sigma_spectral", 0.0)
    tau = values.get("tau_photon_distance", 0.0)
    eta = values.get("eta_Ni_luminosity", 0.0)
    duration = values.get("delta_duration_extra", 0.0)

    h = delta_h
    g = delta_g
    mp = -1.0 + delta_mp
    me = -1.0 + delta_me
    alpha = delta_alpha
    c = 1.0
    q = 2.0 * alpha + me + 2.0 * c - h
    redshift = c + sigma - q
    mch = 1.5 * (h + c - g) - 2.0 * mp
    ell = mch + eta
    s_power = tau - 0.5 * ell
    kappa = 2.0 * alpha + 2.0 * h - 2.0 * me - 2.0 * c - mp
    diffusion = 0.5 * (kappa + mch - 2.0 * c)
    duration_b = 1.0 + diffusion + duration
    return {
        "h_exponent": h,
        "G_exponent": g,
        "m_p_exponent": mp,
        "m_e_exponent": me,
        "alpha_exponent": alpha,
        "atomic_Q_exponent": q,
        "spectroscopic_redshift_exponent": redshift,
        "M_Ch_exponent": mch,
        "source_luminosity_exponent": ell,
        "supernova_composite_s": s_power,
        "opacity_exponent": kappa,
        "diffusion_time_exponent": diffusion,
        "duration_b_predicted": duration_b,
        "me_over_mp_exponent": me - mp,
        "proton_gravity_coupling_exponent": g + 2.0 * mp - h - c,
        "proton_gravity_coupling_delta_exponent": g + 2.0 * delta_mp - h,
    }


def microscopic_evaluator(
    values: dict[str, float], known: dict, one_plus_zmax: float, diffusion: bool
) -> dict[str, object]:
    quantities = microscopic_quantities(values)
    dimensionless_changes = {
        "alpha": quantities["alpha_exponent"],
        "me_over_mp": quantities["me_over_mp_exponent"],
        "proton_gravity": quantities["proton_gravity_coupling_delta_exponent"],
    }
    dim_factor = maximum_factor(dimensionless_changes, one_plus_zmax)
    b_sigma = math.hypot(
        known["DES_time_dilation_sigma_stat"],
        known["DES_time_dilation_sigma_sys"],
    )
    b_pred = quantities["duration_b_predicted"]
    result: dict[str, object] = dict(quantities)
    result.update(
        {
            "alpha_endpoint_ratio": one_plus_zmax
            ** quantities["alpha_exponent"],
            "me_over_mp_endpoint_ratio": one_plus_zmax
            ** quantities["me_over_mp_exponent"],
            "proton_gravity_endpoint_ratio_relative_to_baseline": one_plus_zmax
            ** quantities["proton_gravity_coupling_delta_exponent"],
            "maximum_dimensionless_endpoint_change_factor": dim_factor,
            "dimensionless_drift_warning": dim_factor
            > known.get("dimensionless_drift_warning_factor_at_zmax", 1.1),
            "time_dilation_evaluated": diffusion,
            "duration_tension_sigma": abs(b_pred - known["DES_time_dilation_b"])
            / b_sigma,
            "time_dilation_gate_pass": bool(
                diffusion
                and abs(b_pred - known["DES_time_dilation_b"]) <= 2.0 * b_sigma
            ),
        }
    )
    return result


def enumerate_sparse(
    system: str,
    variables: list[str],
    matrix: np.ndarray,
    target: np.ndarray,
    maximum_active: int,
    tolerance: float,
    one_plus_zmax: float,
    evaluator,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for active_count in range(1, maximum_active + 1):
        for indices in itertools.combinations(range(len(variables)), active_count):
            submatrix = matrix[:, indices]
            solution, _, rank, singular_values = np.linalg.lstsq(
                submatrix, target, rcond=None
            )
            predicted = submatrix @ solution
            residual = predicted - target
            maximum_residual = float(np.max(np.abs(residual)))
            if maximum_residual > tolerance:
                continue
            if np.any(np.abs(solution) <= tolerance):
                continue
            active_names = [variables[index] for index in indices]
            values = {name: float(value) for name, value in zip(active_names, solution)}
            full_column_rank = int(rank) == active_count
            candidate = {
                "system": system,
                "candidate_id": f"{system}:{'+'.join(active_names)}",
                "active_count": active_count,
                "active_couplings": ";".join(active_names),
                "active_values_json": json.dumps(values, sort_keys=True),
                "full_column_rank": full_column_rank,
                "matrix_rank": int(rank),
                "constraint_count": int(matrix.shape[0]),
                "maximum_constraint_residual": maximum_residual,
                "L1_exponent_norm": float(np.sum(np.abs(solution))),
                "L2_exponent_norm": float(np.linalg.norm(solution)),
                "maximum_absolute_exponent": float(np.max(np.abs(solution))),
                "minimum_singular_value": (
                    float(np.min(singular_values)) if len(singular_values) else np.nan
                ),
                "condition_number": (
                    float(np.max(singular_values) / np.min(singular_values))
                    if len(singular_values) and np.min(singular_values) > 0.0
                    else np.inf
                ),
                "maximum_active_endpoint_change_factor": maximum_factor(
                    values, one_plus_zmax
                ),
            }
            candidate.update(evaluator(values))
            candidate["supernova_reference_delta_BIC_vs_flat_LCDM"] = 83.76092105304724
            candidate["supernova_gate_pass"] = False
            candidate["action_level_model_supplied"] = False
            candidate["promotion_gate_pass"] = False
            rows.append(candidate)
    if not rows:
        raise RuntimeError(f"No exact sparse candidates found for {system}")
    frame = pd.DataFrame(rows)
    return frame.sort_values(
        [
            "active_count",
            "full_column_rank",
            "maximum_active_endpoint_change_factor",
            "L2_exponent_norm",
            "active_couplings",
        ],
        ascending=[True, False, True, True, True],
    ).reset_index(drop=True)


def constraint_table(
    systems: list[tuple[str, list[str], list[str], np.ndarray, np.ndarray]]
) -> pd.DataFrame:
    rows = []
    for system, variables, constraints, matrix, target in systems:
        for row_index, constraint in enumerate(constraints):
            for column_index, variable in enumerate(variables):
                rows.append(
                    {
                        "system": system,
                        "constraint": constraint,
                        "variable": variable,
                        "coefficient": float(matrix[row_index, column_index]),
                        "target": float(target[row_index]),
                    }
                )
    return pd.DataFrame(rows)


def identified_minimum(frame: pd.DataFrame) -> pd.DataFrame:
    identified = frame[frame["full_column_rank"]].copy()
    minimum = int(identified["active_count"].min())
    return identified[identified["active_count"] == minimum].copy()


def best_candidates(frames: list[pd.DataFrame], rows_per_system: int = 8) -> pd.DataFrame:
    selected = []
    for frame in frames:
        minimum = identified_minimum(frame)
        selected.append(minimum.head(rows_per_system))
    return pd.concat(selected, ignore_index=True, sort=False)


def format_cell(value) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, (bool, np.bool_)):
        return "PASS" if value else "FAIL"
    if isinstance(value, (float, np.floating)):
        if value == 0.0:
            return "0"
        if abs(value) >= 1e4 or abs(value) < 1e-4:
            return f"{value:.4e}"
        return f"{value:.6f}"
    return str(value).replace("|", "\\|")


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    view = frame.loc[:, [column for column in columns if column in frame.columns]]
    header = "| " + " | ".join(view.columns) + " |"
    separator = "| " + " | ".join("---" for _ in view.columns) + " |"
    body = [
        "| " + " | ".join(format_cell(value) for value in row) + " |"
        for row in view.itertuples(index=False, name=None)
    ]
    return "\n".join([header, separator, *body])


def make_plot(path: Path, frames: list[pd.DataFrame]) -> None:
    colours = {
        "stated_phenomenology": "#006d77",
        "microscopic_closure": "#bc4749",
        "diffusion_control": "#6a4c93",
    }
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.4), constrained_layout=True)
    for frame in frames:
        system = str(frame.iloc[0]["system"])
        axes[0].scatter(
            frame["active_count"],
            frame["maximum_active_endpoint_change_factor"],
            s=np.where(frame["full_column_rank"], 46, 18),
            alpha=np.where(frame["full_column_rank"], 0.82, 0.22),
            color=colours[system],
            label=system.replace("_", " "),
        )
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Active coupling count")
    axes[0].set_ylabel("Largest active-coupling change factor at zmax")
    axes[0].set_title("All exact algebraic families")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8)

    labels = []
    factors = []
    bar_colours = []
    for frame in frames:
        minimum = identified_minimum(frame).iloc[0]
        labels.append(str(minimum["system"]).replace("_", " "))
        factors.append(float(minimum["maximum_active_endpoint_change_factor"]))
        bar_colours.append(colours[str(minimum["system"])])
    axes[1].barh(labels, factors, color=bar_colours)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Smallest endpoint change factor among minimum-rank repairs")
    axes[1].set_title("Best-ranked identified candidate")
    axes[1].grid(axis="x", alpha=0.25)
    for index, value in enumerate(factors):
        axes[1].text(value * 1.03, index, f"{value:.3g}x", va="center", fontsize=9)
    fig.suptitle("Peter clock-law sparse-coupling closure search", fontsize=14)
    fig.savefig(path, dpi=190)
    plt.close(fig)


def write_report(
    path: Path,
    protocol_hash: str,
    known: dict,
    stated: pd.DataFrame,
    microscopic: pd.DataFrame,
    diffusion: pd.DataFrame,
) -> None:
    stated_min = identified_minimum(stated)
    micro_min = identified_minimum(microscopic)
    diffusion_min = identified_minimum(diffusion)
    stated_best = stated_min.iloc[0]
    micro_best = micro_min.iloc[0]
    diffusion_best = diffusion_min.iloc[0]
    lines = [
        "# Peter clock-law sparse-coupling closure search",
        "",
        f"Generated: {datetime.now():%Y-%m-%d %H:%M:%S}",
        "",
        "## Decision",
        "",
        "**THE STATED EQUATIONS REQUIRE AT LEAST TWO INDEPENDENT REPAIRS; THE "
        "ILLUSTRATIVE ATOMIC/CHANDRASEKHAR/DIFFUSION CLOSURE REQUIRES AT LEAST "
        "THREE. NO CANDIDATE PASSES THE EXISTING SUPERNOVA OR ACTION-LEVEL "
        "PROMOTION GATES.**",
        "",
        f"The smallest stated-phenomenology closure uses "
        f"`{int(stated_best['active_count'])}` active exponents. The microscopic "
        f"atomic plus luminosity closure also needs "
        f"`{int(micro_best['active_count'])}`, while adding the simplified "
        f"diffusion-time equation raises the minimum to "
        f"`{int(diffusion_best['active_count'])}`. This is evidence for missing "
        "structure in the explicit closure, not evidence that any enumerated "
        "power law is physical.",
        "",
        f"Every candidate inherits the best possible Peter-path composite "
        f"supernova result `Delta BIC={known['supernova_free_s_delta_BIC_vs_flat_LCDM']:+.6f}` "
        "relative to flat Lambda-CDM. The frozen promotion threshold is "
        "`abs(Delta BIC)<=2`, so no scalar-field integration or new supernova "
        "sampling is warranted from this screen.",
        "",
        "## Frozen status",
        "",
        f"Protocol SHA-256: `{protocol_hash}`.",
        "",
        "The earlier supernova results were disclosed before enumeration. The "
        "new output is the sparse rank, endpoint consequences and explicit "
        "separation of identified solutions from rank-deficient families.",
        "",
        "## Exponent convention",
        "",
        "$$",
        r"\begin{aligned}",
        r"Z                              &= 1+z, \\",
        r"X(z)/X(0)                      &= Z^{x_X}, \\",
        rf"Z_{{\rm max}}                  &= {1.0 + known['maximum_supernova_redshift']:.6f}.",
        r"\end{aligned}",
        "$$",
        "",
        "## Stated phenomenology",
        "",
        "The baseline assertion `x_Q=-1` is allowed either to change directly "
        "or to be accompanied by an extra spectral coupling. Photon transfer "
        "and source luminosity are kept separate even though supernova "
        "magnitudes identify only their combination.",
        "",
        "$$",
        r"\begin{aligned}",
        r"x_Q                            &= -1+\delta_Q, \\",
        r"-\delta_Q+\sigma              &= -1, \\",
        r"\tau-\frac{\ell}{2}           &= 0.0376049011.",
        r"\end{aligned}",
        "$$",
        "",
        f"There are `{len(stated_min)}` full-rank, minimum-count stated "
        "solutions. The best-ranked endpoint repair is "
        f"`{stated_best['active_couplings']}` with values "
        f"`{stated_best['active_values_json']}`.",
        "",
        markdown_table(
            stated_min,
            [
                "active_couplings",
                "active_values_json",
                "maximum_active_endpoint_change_factor",
                "atomic_Q_final_exponent",
                "supernova_composite_s",
                "duration_b_predicted",
            ],
        ),
        "",
        "The economical branches either replace the asserted atomic clock law "
        "with `Q=constant`, or retain it by inserting a new spectral factor "
        "proportional to `(1+z)^-1`. A second independent photon/source law is "
        "then still required to reproduce the measured supernova composite.",
        "",
        "## Microscopic closure",
        "",
        "$$",
        r"\begin{aligned}",
        r"Q                              &\propto \frac{\alpha^2m_ec^2}{h}, \\",
        r"M_{\rm Ch}                     &\propto \left(\frac{\hbar c}{G}\right)^{3/2}m_p^{-2}, \\",
        r"L_{\rm Ia}                     &\propto M_{\rm Ch}Z^{\eta_{\rm Ni}}, \\",
        r"\delta_h-\delta_{m_e}-2\delta_\alpha+\sigma &= 1, \\",
        r"\tau-\frac{3}{4}\delta_h+\frac{3}{4}\delta_G+\delta_{m_p}-\frac{1}{2}\eta_{\rm Ni} &= 1.7876049011.",
        r"\end{aligned}",
        "$$",
        "",
        f"There are `{len(micro_min)}` full-rank minimum-count microscopic "
        f"solutions. The best-ranked active-set endpoint factor is "
        f"`{float(micro_best['maximum_active_endpoint_change_factor']):.6g}` "
        f"for `{micro_best['active_couplings']}`. This ranking is descriptive: "
        "a small dimensional exponent norm is not a physical prior.",
        "",
        markdown_table(
            micro_min.head(12),
            [
                "active_couplings",
                "active_values_json",
                "maximum_active_endpoint_change_factor",
                "maximum_dimensionless_endpoint_change_factor",
                "atomic_Q_exponent",
                "M_Ch_exponent",
                "supernova_composite_s",
            ],
        ),
        "",
        "Order-unity changes in `alpha`, `m_e/m_p` or "
        "`G m_p^2/(hbar c)` are warnings, not accepted solutions. A local "
        "atomic-clock bound cannot be translated into a cosmological exponent "
        "without the scalar field's time history and environmental response.",
        "",
        "## Diffusion-time control",
        "",
        "$$",
        r"\begin{aligned}",
        r"\kappa                         &\propto \frac{\alpha^2h^2}{m_e^2c^2m_p}, \\",
        r"t_{\rm diff}                   &\propto \left(\frac{\kappa M_{\rm Ch}}{vc}\right)^{1/2}, \\",
        r"v                              &\propto c, \\",
        r"\frac{7}{4}\delta_h-\frac{3}{4}\delta_G-\frac{3}{2}\delta_{m_p}-\delta_{m_e}+\delta_\alpha+\delta_{\rm duration} &= -1.247.",
        r"\end{aligned}",
        "$$",
        "",
        f"The simplified light-curve system has `{len(diffusion_min)}` "
        f"full-rank solutions at the minimum active count of "
        f"`{int(diffusion_best['active_count'])}`. Its best-ranked member is "
        f"`{diffusion_best['active_couplings']}` with values "
        f"`{diffusion_best['active_values_json']}`.",
        "",
        markdown_table(
            diffusion_min.head(12),
            [
                "active_couplings",
                "active_values_json",
                "maximum_active_endpoint_change_factor",
                "maximum_dimensionless_endpoint_change_factor",
                "duration_b_predicted",
                "duration_tension_sigma",
                "time_dilation_gate_pass",
            ],
        ),
        "",
        "This branch is only an opacity/diffusion sensitivity. A credible Type "
        "Ia calculation must replace it with radiative transfer, nickel yield, "
        "ejecta structure, colour standardisation and selection effects.",
        "",
        "## Gate matrix",
        "",
        "| Gate | Result | Reason |",
        "| --- | --- | --- |",
        "| Algebraic stated closure | PASS | Exact two-coupling solutions exist. |",
        "| Algebraic microscopic closure | PASS | Exact two-coupling solutions exist. |",
        "| DES time dilation in diffusion control | PASS by construction | The third equation targets the measured central value. |",
        f"| Released-covariance supernova gate | FAIL | Best Peter-path composite has Delta BIC `{known['supernova_free_s_delta_BIC_vs_flat_LCDM']:+.6f}`. |",
        "| Dimensionless external constraints | NOT YET TESTED | Requires an action and scalar history; dimensional exponents alone are not observables. |",
        "| Action-level closure | FAIL | No scalar, matter and electromagnetic action was supplied. |",
        "| Promotion | FAIL | Statistical and theory gates remain open. |",
        "",
        "## What the search identifies",
        "",
        "1. A single missing scalar power cannot repair both the redshift/atomic "
        "identity and the supernova photon/source relation.",
        "2. At least two independent combinations are required before light-curve "
        "physics; the simplified diffusion equation requires a third.",
        "3. The least elaborate phenomenological branch is `Q=constant` plus a "
        "small composite flux correction. Retaining `Q=(1+z)^-1` instead needs "
        "an additional inverse spectral factor, which is observationally "
        "degenerate at the redshift-identity level.",
        "4. The supernova radial-shape failure remains. Algebraic closure does "
        "not make the existing distance law competitive.",
        "",
        "## Required next derivation",
        "",
        "Before further cosmological sampling, specify one covariant action with "
        "a dimensionless scalar and derive from it: (i) the observable "
        "spectroscopic redshift; (ii) atomic transition ratios; (iii) photon "
        "geodesics, arrival rates and opacity; (iv) dimensionless gravitational "
        "couplings; and (v) Type Ia source and diffusion laws. The resulting "
        "functions must be fixed before fitting and must reproduce local clocks, "
        "supernova time dilation and distance-duality controls.",
        "",
        "## Reproduction",
        "",
        "```powershell",
        "python github_export/code/rawmu/search_peter_sparse_couplings_2026_07_19.py",
        "```",
        "",
        "## Primary research links",
        "",
        "- [Bimetric varying-speed-of-light cosmology](https://arxiv.org/abs/gr-qc/0202012)",
        "- [Varying-alpha action and dimensionless observables](https://arxiv.org/abs/gr-qc/0208081)",
        "- [Dimensionless-constant interpretation caveat](https://arxiv.org/abs/hep-th/0208093)",
        "- [DES supernova time-dilation measurement](https://arxiv.org/abs/2406.05050)",
        "- [Type Ia radiative-transfer modelling](https://arxiv.org/abs/astro-ph/0609562)",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def write_manifest(outdir: Path, paths: list[Path]) -> Path:
    rows = []
    for path in sorted(set(paths), key=lambda item: item.as_posix()):
        relative = path.relative_to(REPO_ROOT).as_posix()
        rows.append(
            {"path": relative, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
        )
    manifest = outdir / f"peter_sparse_coupling_search_manifest_{DATE_TAG}.csv"
    pd.DataFrame(rows).to_csv(manifest, index=False)
    return manifest


def self_test() -> None:
    quantities = microscopic_quantities({})
    expected = {
        "atomic_Q_exponent": 1.0,
        "spectroscopic_redshift_exponent": 0.0,
        "M_Ch_exponent": 3.5,
        "supernova_composite_s": -1.75,
        "opacity_exponent": 1.0,
        "diffusion_time_exponent": 1.25,
        "duration_b_predicted": 2.25,
    }
    for key, value in expected.items():
        if not np.isclose(quantities[key], value, rtol=0.0, atol=1e-13):
            raise AssertionError(f"Microscopic baseline changed for {key}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--prior-summary", type=Path, default=PRIOR_SUMMARY)
    args = parser.parse_args()

    self_test()
    protocol, _ = verify_inputs(args.protocol, args.prior_summary)
    args.outdir.mkdir(parents=True, exist_ok=True)
    known = dict(protocol["known_inputs"])
    known["dimensionless_drift_warning_factor_at_zmax"] = protocol[
        "external_gate_policy"
    ]["dimensionless_drift_warning_factor_at_zmax"]
    tolerance = protocol["enumeration"][
        "acceptance_absolute_constraint_residual"
    ]
    one_plus_zmax = 1.0 + known["maximum_supernova_redshift"]

    stated_variables = [
        "delta_Q_atomic",
        "sigma_spectral",
        "tau_photon_distance",
        "ell_source",
    ]
    stated_constraints = ["spectroscopic_redshift", "supernova_composite"]
    stated_matrix = np.array([[-1.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, -0.5]])
    stated_target = np.array([-1.0, known["supernova_composite_s"]])

    micro_variables = [
        "delta_h",
        "delta_G",
        "delta_m_p",
        "delta_m_e",
        "delta_alpha",
        "sigma_spectral",
        "tau_photon_distance",
        "eta_Ni_luminosity",
    ]
    micro_constraints = ["spectroscopic_redshift", "supernova_composite"]
    micro_matrix = np.array(
        [
            [1.0, 0.0, 0.0, -1.0, -2.0, 1.0, 0.0, 0.0],
            [-0.75, 0.75, 1.0, 0.0, 0.0, 0.0, 1.0, -0.5],
        ]
    )
    micro_target = np.array(protocol["microscopic_constraints"]["target_vector_for_perturbations"])

    diffusion_variables = micro_variables + ["delta_duration_extra"]
    diffusion_constraints = micro_constraints + ["light_curve_duration"]
    diffusion_matrix = np.vstack(
        [
            np.pad(micro_matrix, ((0, 0), (0, 1))),
            np.array([1.75, -0.75, -1.5, -1.0, 1.0, 0.0, 0.0, 0.0, 1.0]),
        ]
    )
    diffusion_target = np.array(
        [
            *protocol["microscopic_constraints"]["target_vector_for_perturbations"],
            protocol["diffusion_control"]["target_perturbation"],
        ]
    )

    stated = enumerate_sparse(
        "stated_phenomenology",
        stated_variables,
        stated_matrix,
        stated_target,
        protocol["enumeration"]["stated_maximum_active_couplings"],
        tolerance,
        one_plus_zmax,
        lambda values: stated_evaluator(values, known, one_plus_zmax),
    )
    microscopic = enumerate_sparse(
        "microscopic_closure",
        micro_variables,
        micro_matrix,
        micro_target,
        protocol["enumeration"]["microscopic_maximum_active_couplings"],
        tolerance,
        one_plus_zmax,
        lambda values: microscopic_evaluator(
            values, known, one_plus_zmax, diffusion=False
        ),
    )
    diffusion = enumerate_sparse(
        "diffusion_control",
        diffusion_variables,
        diffusion_matrix,
        diffusion_target,
        protocol["enumeration"]["diffusion_maximum_active_couplings"],
        tolerance,
        one_plus_zmax,
        lambda values: microscopic_evaluator(
            values, known, one_plus_zmax, diffusion=True
        ),
    )

    systems = [
        (
            "stated_phenomenology",
            stated_variables,
            stated_constraints,
            stated_matrix,
            stated_target,
        ),
        (
            "microscopic_closure",
            micro_variables,
            micro_constraints,
            micro_matrix,
            micro_target,
        ),
        (
            "diffusion_control",
            diffusion_variables,
            diffusion_constraints,
            diffusion_matrix,
            diffusion_target,
        ),
    ]
    constraints = constraint_table(systems)
    best = best_candidates([stated, microscopic, diffusion])

    stated_path = args.outdir / f"peter_sparse_coupling_stated_candidates_{DATE_TAG}.csv"
    micro_path = args.outdir / f"peter_sparse_coupling_microscopic_candidates_{DATE_TAG}.csv"
    diffusion_path = args.outdir / f"peter_sparse_coupling_diffusion_candidates_{DATE_TAG}.csv"
    constraint_path = args.outdir / f"peter_sparse_coupling_constraint_matrix_{DATE_TAG}.csv"
    best_path = args.outdir / f"peter_sparse_coupling_best_candidates_{DATE_TAG}.csv"
    summary_path = args.outdir / f"peter_sparse_coupling_search_summary_{DATE_TAG}.json"
    report_path = args.outdir / f"peter_sparse_coupling_search_report_{DATE_TAG}.md"
    plot_path = args.outdir / f"peter_sparse_coupling_search_readout_{DATE_TAG}.png"

    stated.to_csv(stated_path, index=False)
    microscopic.to_csv(micro_path, index=False)
    diffusion.to_csv(diffusion_path, index=False)
    constraints.to_csv(constraint_path, index=False)
    best.to_csv(best_path, index=False)
    make_plot(plot_path, [stated, microscopic, diffusion])
    write_report(
        report_path,
        sha256_file(args.protocol),
        known,
        stated,
        microscopic,
        diffusion,
    )

    frame_summary = {}
    for frame in (stated, microscopic, diffusion):
        system = str(frame.iloc[0]["system"])
        identified = frame[frame["full_column_rank"]]
        minimum = identified_minimum(frame)
        frame_summary[system] = {
            "exact_family_count": len(frame),
            "identified_candidate_count": len(identified),
            "minimum_identified_active_count": int(minimum.iloc[0]["active_count"]),
            "minimum_identified_candidate_count": len(minimum),
            "best_ranked_candidate": safe_json(minimum.iloc[0].to_dict()),
        }
    summary = {
        "analysis_id": ANALYSIS_ID,
        "generated": datetime.now().isoformat(timespec="seconds"),
        "protocol_sha256": sha256_file(args.protocol),
        "prior_summary_sha256": sha256_file(args.prior_summary),
        "known_inputs": known,
        "systems": frame_summary,
        "gates": {
            "frozen_protocol_hash_verified": True,
            "prior_input_regression_verified": True,
            "stated_algebraic_closure_exists": True,
            "microscopic_algebraic_closure_exists": True,
            "diffusion_algebraic_closure_exists": True,
            "released_covariance_supernova_gate": False,
            "action_level_model_supplied": False,
            "external_dimensionless_constant_gates_applied": False,
            "promote_candidate": False,
        },
        "decision": "NO_PROMOTION_MINIMUM_TWO_COUPLINGS_THREE_WITH_DIFFUSION",
        "claim_boundary": protocol["claim_boundary"],
    }
    summary_path.write_text(
        json.dumps(safe_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    manifest = write_manifest(
        args.outdir,
        [
            SCRIPT_PATH,
            args.protocol,
            args.outdir / f"peter_sparse_coupling_search_protocol_{DATE_TAG}.md",
            args.prior_summary,
            stated_path,
            micro_path,
            diffusion_path,
            constraint_path,
            best_path,
            summary_path,
            report_path,
            plot_path,
        ],
    )
    print(
        "Minimum identified active counts: "
        f"stated={frame_summary['stated_phenomenology']['minimum_identified_active_count']} "
        f"microscopic={frame_summary['microscopic_closure']['minimum_identified_active_count']} "
        f"diffusion={frame_summary['diffusion_control']['minimum_identified_active_count']}",
        flush=True,
    )
    print(f"Decision: {summary['decision']}", flush=True)
    print(f"Saved report: {report_path}", flush=True)
    print(f"Saved manifest: {manifest}", flush=True)


if __name__ == "__main__":
    main()
