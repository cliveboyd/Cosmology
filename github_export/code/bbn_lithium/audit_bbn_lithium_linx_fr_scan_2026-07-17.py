#!/usr/bin/env python3
"""Audit the FR/LINX lithium scan catalogue without running LINX.

This programme reconstructs the registered proposal catalogue, checks point-ID
and score integrity, and quantifies proposal artefacts. It is read-only with
respect to the long-running LINX catalogue.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
SCAN_PROGRAM = SCRIPT_PATH.with_name("analyze_bbn_lithium_linx_fr_network_2026-07-16.py")
SCAN_DIR = (
    REPO_ROOT
    / "plamb_runs"
    / "diagnostics"
    / "bbn_lithium_linx_key_fr_network_overnight_20260716"
)
DEFAULT_OUTDIR = REPO_ROOT / "github_export" / "results" / "2026-07-17" / "bbn_lithium"
DEFAULT_LINX_PYTHON = REPO_ROOT / "plamb_runs" / "envs" / "linx_py311" / "python.exe"
TAU_REFERENCE_SECONDS = 879.4
ANCHOR_SETS = {
    "registered": {"D": (2.508, 0.029), "He": (0.245, 0.003), "tau": (879.4, 0.6)},
    "deuterium_2024": {"D": (2.533, 0.024), "He": (0.245, 0.003), "tau": (879.4, 0.6)},
    "helium_lbt_2026": {"D": (2.508, 0.029), "He": (0.2458, 0.0013), "tau": (879.4, 0.6)},
    "neutron_pdg_2025": {"D": (2.508, 0.029), "He": (0.245, 0.003), "tau": (878.4, 0.5)},
    "combined_current": {"D": (2.533, 0.024), "He": (0.2458, 0.0013), "tau": (878.4, 0.5)},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-program", type=Path, default=SCAN_PROGRAM)
    parser.add_argument("--scan-dir", type=Path, default=SCAN_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--linx-python", type=Path, default=DEFAULT_LINX_PYTHON)
    parser.add_argument("--n-random", type=int, default=30_000)
    parser.add_argument("--seed", type=int, default=20_260_716)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_scan_module(path: Path):
    spec = importlib.util.spec_from_file_location("fr_linx_scan_for_audit", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load scan programme: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def reconstruct_catalogue(module: Any, available_reactions: list[str], n_random: int, seed: int) -> dict[str, Any]:
    captured = []
    original_add_point = module.add_point

    def capture(points, family, eta_fac, clock_factor, tau_n_fac, q_values):
        original_add_point(points, family, eta_fac, clock_factor, tau_n_fac, q_values)
        captured.append(points[-1])

    module.add_point = capture
    try:
        built = module.build_points("overnight", n_random, seed)
    finally:
        module.add_point = original_add_point

    final_points = module.normalise_points_for_reactions(built, available_reactions)
    available = set(available_reactions)
    groups: dict[str, list[tuple[Any, ...]]] = defaultdict(list)
    group_families: dict[str, set[str]] = defaultdict(set)
    for point in captured:
        q_values = {name: value for name, value in point.q_values.items() if name in available}
        point_id = module.make_point_id(
            point.family,
            point.eta_fac,
            point.clock_factor,
            point.tau_n_fac,
            q_values,
        )
        semantic = (
            point.family,
            float(point.eta_fac),
            float(point.clock_factor),
            float(point.tau_n_fac),
            tuple(sorted((name, float(value)) for name, value in q_values.items())),
        )
        groups[point_id].append(semantic)
        group_families[point_id].add(point.family)

    duplicate_groups = {point_id: values for point_id, values in groups.items() if len(values) > 1}
    semantic_collision_groups = {
        point_id: values
        for point_id, values in duplicate_groups.items()
        if len(set(values)) > 1
    }
    random_collision_groups = {
        point_id: values
        for point_id, values in semantic_collision_groups.items()
        if group_families[point_id] == {"random_targeted_li_be"}
    }
    identical_collapse_rows = sum(len(values) - 1 for values in duplicate_groups.values() if len(set(values)) == 1)

    final_by_family = Counter(point.family for point in final_points)
    return {
        "raw_points_before_id_dedup": len(captured),
        "points_after_initial_id_dedup": len(built),
        "points_after_network_normalisation": len(final_points),
        "normalised_duplicate_id_groups": len(duplicate_groups),
        "identical_rows_removed_by_network_normalisation": identical_collapse_rows,
        "semantic_point_id_collision_groups": len(semantic_collision_groups),
        "random_rounding_collision_groups": len(random_collision_groups),
        "family_counts_expected": dict(sorted(final_by_family.items())),
        "expected_ids": {point.point_id for point in final_points},
        "final_points": final_points,
    }


def latest_successes(raw: pd.DataFrame) -> tuple[pd.DataFrame, set[str], set[str]]:
    frame = raw.copy()
    frame["_row_order"] = np.arange(len(frame))
    ok = frame.loc[frame["status"].eq("ok")].copy()
    ok = ok.sort_values("_row_order").drop_duplicates("point_id", keep="last")
    ok_ids = set(ok["point_id"].dropna().astype(str))
    failed_ids = set(frame.loc[frame["status"].eq("failed"), "point_id"].dropna().astype(str))
    return ok.drop(columns="_row_order"), ok_ids, failed_ids - ok_ids


def numeric(frame: pd.DataFrame, columns: list[str]) -> None:
    for column in columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")


def finite_max(values: pd.Series) -> float:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    return float(np.max(array)) if len(array) else math.nan


def verify_scores(frame: pd.DataFrame, reaction_columns: list[str]) -> dict[str, float]:
    q2 = (frame[reaction_columns].astype(float) ** 2).sum(axis=1)
    chi2_d = ((frame["D_H_1e5"] - 2.508) / 0.029) ** 2
    chi2_he = ((frame["Yp_mass"] - 0.245) / 0.003) ** 2
    chi2_li = ((frame["Li7_H_1e10"] - 1.45) / 0.25) ** 2
    chi2_li_lower = (np.maximum(0.0, 1.45 - frame["Li7_H_1e10"]) / 0.25) ** 2
    chi2_eta = ((frame["eta_fac"] - 1.0) / (0.04 / 6.12)) ** 2
    tau_sigma_fac = 0.6 / 879.4
    chi2_tau = ((frame["tau_n_fac"] - 1.0) / tau_sigma_fac) ** 2
    total = chi2_d + chi2_he + chi2_li + chi2_eta + chi2_tau + q2
    total_lower = chi2_d + chi2_he + chi2_li_lower + chi2_eta + chi2_tau + q2
    return {
        "max_abs_q_pull_chi2_error": finite_max((frame["q_pull_chi2"] - q2).abs()),
        "max_abs_chi2_D_error": finite_max((frame["chi2_D"] - chi2_d).abs()),
        "max_abs_chi2_He4_error": finite_max((frame["chi2_He4"] - chi2_he).abs()),
        "max_abs_chi2_Li_measurement_error": finite_max((frame["chi2_Li_measurement"] - chi2_li).abs()),
        "max_abs_chi2_eta_error": finite_max((frame["chi2_eta_cmb"] - chi2_eta).abs()),
        "max_abs_chi2_tau_error": finite_max((frame["chi2_tau"] - chi2_tau).abs()),
        "max_abs_total_measurement_error": finite_max((frame["chi2_total_li_measurement"] - total).abs()),
        "max_abs_total_lower_bound_error": finite_max((frame["chi2_total_li_lower_bound"] - total_lower).abs()),
    }


def proposal_metrics(frame: pd.DataFrame, reaction_columns: list[str]) -> dict[str, Any]:
    random_rows = frame.loc[frame["family"].eq("random_targeted_li_be")].copy()
    q = random_rows[reaction_columns].astype(float)
    all_zero = (q.abs() <= 1e-15).all(axis=1)
    eta_ok = (random_rows["eta_fac"] - 1.0).abs() <= 2.0 * (0.04 / 6.12)
    tau_ok = (random_rows["tau_n_fac"] - 1.0).abs() <= 2.0 * (0.6 / 879.4)
    su2_base = random_rows["clock_factor"].ge(1.0) & eta_ok & tau_ok

    by_reaction = {}
    for column in reaction_columns:
        values = q[column].to_numpy(dtype=float)
        nonzero = values[np.abs(values) > 1e-15]
        by_reaction[column.removeprefix("q_")] = {
            "fraction_exact_zero": float(np.mean(np.abs(values) <= 1e-15)),
            "sample_sd_including_zero": float(np.std(values, ddof=1)),
            "sample_sd_nonzero": float(np.std(nonzero, ddof=1)) if len(nonzero) > 1 else math.nan,
            "n_at_clip_abs_4": int(np.sum(np.abs(values) >= 4.0 - 1e-12)),
        }

    return {
        "n_random_rows": int(len(random_rows)),
        "n_random_all_available_q_exact_zero": int(all_zero.sum()),
        "fraction_random_all_available_q_exact_zero": float(all_zero.mean()),
        "n_random_q0_inside_registered_su2_base": int((all_zero & su2_base).sum()),
        "eta_fac_min": float(random_rows["eta_fac"].min()),
        "eta_fac_max": float(random_rows["eta_fac"].max()),
        "clock_factor_min": float(random_rows["clock_factor"].min()),
        "clock_factor_max": float(random_rows["clock_factor"].max()),
        "tau_n_fac_min": float(random_rows["tau_n_fac"].min()),
        "tau_n_fac_max": float(random_rows["tau_n_fac"].max()),
        "proposal_q_sigma_before_threshold": 1.35,
        "target_q_prior_sigma_used_in_score": 1.0,
        "q_zeroing_threshold_abs": 0.25,
        "reaction_metrics": by_reaction,
    }


def score_decomposition(frame: pd.DataFrame, reaction_columns: list[str]) -> pd.DataFrame:
    rows = []
    selectors = {
        "best_lithium_measurement_objective": frame["chi2_total_li_measurement"].idxmin(),
        "best_lithium_lower_bound_objective": frame["chi2_total_li_lower_bound"].idxmin(),
    }
    d_he = frame.loc[
        ((frame["D_H_1e5"] - 2.508).abs() <= 2 * 0.029)
        & ((frame["Yp_mass"] - 0.245).abs() <= 2 * 0.003)
    ]
    if not d_he.empty:
        selectors["minimum_lithium_given_D_He_2sigma"] = d_he["Li7_H_1e10"].idxmin()

    for label, index in selectors.items():
        row = frame.loc[index]
        item = {
            "selection": label,
            "point_id": row["point_id"],
            "family": row["family"],
            "eta_fac": row["eta_fac"],
            "clock_factor": row["clock_factor"],
            "delta_neff_initial": row["delta_neff_initial"],
            "tau_n_fac": row["tau_n_fac"],
            "D_H_1e5": row["D_H_1e5"],
            "Yp_mass": row["Yp_mass"],
            "Li7_H_1e10": row["Li7_H_1e10"],
            "chi2_D": row["chi2_D"],
            "chi2_He4": row["chi2_He4"],
            "chi2_Li_measurement": row["chi2_Li_measurement"],
            "chi2_eta_cmb": row["chi2_eta_cmb"],
            "chi2_tau": row["chi2_tau"],
            "q_pull_chi2": row["q_pull_chi2"],
            "chi2_total_li_measurement": row["chi2_total_li_measurement"],
            "chi2_total_li_lower_bound": row["chi2_total_li_lower_bound"],
        }
        for column in reaction_columns:
            item[column] = row[column]
        rows.append(item)
    return pd.DataFrame(rows)


def git_revision(path: Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def python_environment(python_path: Path) -> dict[str, Any]:
    packages = ["numpy", "jax", "jaxlib", "diffrax", "equinox", "lineax", "optimistix", "interpax"]
    code = (
        "import importlib.metadata as m,json,sys;"
        f"names={packages!r};"
        "print(json.dumps({'executable':sys.executable,'python':sys.version,"
        "'packages':{name:m.version(name) for name in names}}))"
    )
    try:
        output = subprocess.check_output(
            [str(python_path), "-c", code],
            text=True,
            stderr=subprocess.STDOUT,
        )
        return json.loads(output)
    except (OSError, subprocess.CalledProcessError, json.JSONDecodeError) as exc:
        return {"executable": str(python_path), "error": repr(exc)}


def write_report(path: Path, metrics: dict[str, Any], decomposition: pd.DataFrame) -> None:
    catalogue = metrics["catalogue"]
    proposal = metrics["proposal"]
    score = metrics["score_integrity"]
    best = decomposition.loc[decomposition["selection"].eq("best_lithium_measurement_objective")].iloc[0]
    lines = [
        "# FR LINX Lithium Scan: Empirical Audit",
        "",
        f"Generated: {metrics['generated_utc']}",
        "",
        "## Catalogue Integrity",
        "",
        f"- Successful unique rows read: `{catalogue['n_unique_ok']}`; unresolved failures: `{catalogue['n_unresolved']}`.",
        f"- Reconstructed registered target: `{catalogue['n_expected']}` points; missing IDs: `{catalogue['n_expected_missing']}`; unexpected IDs: `{catalogue['n_unexpected_ok']}`.",
        f"- Semantic point-ID collision groups: `{catalogue['semantic_point_id_collision_groups']}`; random rounding collision groups: `{catalogue['random_rounding_collision_groups']}`.",
        f"- Exact duplicate rows removed after key-network normalisation: `{catalogue['identical_rows_removed_by_network_normalisation']}`.",
        "",
        "## Score Integrity",
        "",
        f"- Maximum absolute reconstructed total-score discrepancy: `{score['max_abs_total_measurement_error']:.3g}`.",
        f"- Maximum absolute reconstructed lower-bound-score discrepancy: `{score['max_abs_total_lower_bound_error']:.3g}`.",
        "- These are floating-point reconstruction checks, not a validation of the adopted likelihood.",
        "",
        "## Proposal Diagnostics",
        "",
        f"- Nuclear pulls were proposed from `N(0, {proposal['proposal_q_sigma_before_threshold']}^2)` but scored against `N(0, {proposal['target_q_prior_sigma_used_in_score']}^2)`.",
        f"- Pulls with `|q| <= {proposal['q_zeroing_threshold_abs']}` were replaced by exactly zero.",
        f"- Random rows with all three available Li/Be pulls exactly zero: `{proposal['n_random_all_available_q_exact_zero']}` / `{proposal['n_random_rows']}`.",
        f"- Such random q=0 rows inside the registered SU2 base cuts: `{proposal['n_random_q0_inside_registered_su2_base']}`.",
        "- The q=0 expansion-only scenario therefore contains proposal-generated point mass as well as the deliberate fixed-rate grid.",
        "",
        "## Current Best Penalised Row",
        "",
        f"- Total measurement objective: `{best['chi2_total_li_measurement']:.6g}`.",
        f"- Predictions: D/H x1e5 `{best['D_H_1e5']:.6g}`, Yp `{best['Yp_mass']:.6g}`, Li/H x1e10 `{best['Li7_H_1e10']:.6g}`.",
        f"- Controls: eta factor `{best['eta_fac']:.6g}`, S `{best['clock_factor']:.6g}`, initial Delta Neff `{best['delta_neff_initial']:.6g}`, q-pull penalty `{best['q_pull_chi2']:.6g}`.",
        "",
        "## Interpretation Boundary",
        "",
        "The catalogue arithmetic and ID coverage pass the checks above when the final target is complete. The scan remains a proposal-driven, penalised sensitivity study. It is not posterior sampling, evidence calculation, a goodness-of-fit chi-squared with a fixed number of degrees of freedom, or a self-consistent FR/fundamental-constant BBN model.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def within(values: pd.Series, mean: float, sigma: float, n_sigma: float = 2.0) -> pd.Series:
    return (values - mean).abs().le(n_sigma * sigma)


def anchor_sensitivity(frame: pd.DataFrame) -> pd.DataFrame:
    eta_gate = within(frame["eta_fac"], 1.0, 0.04 / 6.12)
    li_gate = within(frame["Li7_H_1e10"], 1.45, 0.25)
    base_without_updated_terms = (
        frame["chi2_total_li_measurement"]
        - frame["chi2_D"]
        - frame["chi2_He4"]
        - frame["chi2_tau"]
    )
    rows: list[dict[str, Any]] = []
    for anchor_name, anchors in ANCHOR_SETS.items():
        d_mean, d_sigma = anchors["D"]
        he_mean, he_sigma = anchors["He"]
        tau_seconds, tau_sigma_seconds = anchors["tau"]
        tau_mean = tau_seconds / TAU_REFERENCE_SECONDS
        tau_sigma = tau_sigma_seconds / TAU_REFERENCE_SECONDS
        chi2_d = ((frame["D_H_1e5"] - d_mean) / d_sigma) ** 2
        chi2_he = ((frame["Yp_mass"] - he_mean) / he_sigma) ** 2
        chi2_tau = ((frame["tau_n_fac"] - tau_mean) / tau_sigma) ** 2
        reweighted = base_without_updated_terms + chi2_d + chi2_he + chi2_tau
        base = frame["clock_factor"].ge(1.0) & eta_gate & within(frame["tau_n_fac"], tau_mean, tau_sigma)
        abundance_gate = within(frame["D_H_1e5"], d_mean, d_sigma) & within(
            frame["Yp_mass"], he_mean, he_sigma
        )
        scenarios = {
            "su2_expansion_only": base & frame["q_pull_chi2"].le(1e-12),
            "su2_plus_modest_rate_controls": base & frame["q_pull_chi2"].le(9.0),
            "su2_plus_scanned_rate_controls": base,
            "all_scanned_controls": pd.Series(True, index=frame.index),
        }
        for scenario, mask in scenarios.items():
            joint = mask & abundance_gate
            triple = joint & li_gate
            best_index = reweighted.loc[mask].idxmin() if mask.any() else None
            min_li_index = frame.loc[joint, "Li7_H_1e10"].idxmin() if joint.any() else None
            rows.append(
                {
                    "anchor_set": anchor_name,
                    "scenario": scenario,
                    "D_mean": d_mean,
                    "D_sigma": d_sigma,
                    "He_mean": he_mean,
                    "He_sigma": he_sigma,
                    "tau_seconds_mean": tau_seconds,
                    "tau_seconds_sigma": tau_sigma_seconds,
                    "n_rows": int(mask.sum()),
                    "n_D_He_2sigma": int(joint.sum()),
                    "n_D_He_Li_2sigma": int(triple.sum()),
                    "best_reweighted_chi2": None if best_index is None else float(reweighted.at[best_index]),
                    "best_point_id": None if best_index is None else frame.at[best_index, "point_id"],
                    "minimum_Li7_H_1e10_given_D_He": (
                        None if min_li_index is None else float(frame.at[min_li_index, "Li7_H_1e10"])
                    ),
                    "minimum_li_point_id": None if min_li_index is None else frame.at[min_li_index, "point_id"],
                }
            )
    return pd.DataFrame(rows)


def write_anchor_report(path: Path, sensitivity: pd.DataFrame) -> None:
    combined = sensitivity.loc[sensitivity["anchor_set"].eq("combined_current")]

    def fmt(value: float) -> str:
        return "NA" if pd.isna(value) else f"{value:.4g}"

    lines = [
        "# FR LINX Anchor Sensitivity Reweight",
        "",
        "This is a likelihood-only reweight of the completed abundance catalogue. It does not rerun LINX, add nuclear theory covariance, or replace the registered gate.",
        "",
        "The combined-current case uses D/H x 1e5 = `2.533 +/- 0.024`, Yp = `0.2458 +/- 0.0013`, and neutron lifetime = `878.4 +/- 0.5 s`. Lithium remains `1.45 +/- 0.25` in units of Li/H x 1e10.",
        "",
        "## Combined-Current Readout",
        "",
        "| Scenario | Rows | D+He pass | D+He+Li pass | Best reweighted chi2 | Min Li given D+He |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in combined.itertuples(index=False):
        lines.append(
            f"| {row.scenario} | {row.n_rows} | {row.n_D_He_2sigma} | "
            f"{row.n_D_He_Li_2sigma} | {fmt(row.best_reweighted_chi2)} | "
            f"{fmt(row.minimum_Li7_H_1e10_given_D_He)} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The updated anchors are a sensitivity analysis, not a new prior default for an FR cosmology. They retain the standard-background baryon treatment and do not marginalise the fixed deuterium-sensitive reaction pulls. The table should therefore be read as the direction and scale of likelihood sensitivity only.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_manifest(path: Path, entries: list[tuple[str, Path]]) -> None:
    rows = []
    for role, entry in entries:
        rows.append(
            {
                "role": role,
                "path": str(entry),
                "bytes": entry.stat().st_size,
                "sha256": sha256(entry),
            }
        )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["role", "path", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    cli = parse_args()
    cli.outdir.mkdir(parents=True, exist_ok=True)
    points_path = cli.scan_dir / "bbn_lithium_linx_fr_network_points.csv"
    config_path = cli.scan_dir / "bbn_lithium_linx_fr_network_config.json"
    if not points_path.exists() or not config_path.exists():
        raise FileNotFoundError(f"Missing scan inputs under {cli.scan_dir}")

    config = json.loads(config_path.read_text(encoding="utf-8"))
    available = list(config["reaction_names_available"])
    module = load_scan_module(cli.scan_program)
    reconstructed = reconstruct_catalogue(module, available, cli.n_random, cli.seed)

    raw = pd.read_csv(points_path, low_memory=False)
    frame, ok_ids, unresolved = latest_successes(raw)
    reaction_columns = [f"q_{name}" for name in available]
    numeric_columns = [
        "eta_fac",
        "clock_factor",
        "delta_neff_initial",
        "tau_n_fac",
        "D_H_1e5",
        "Yp_mass",
        "Li7_H_1e10",
        "q_pull_chi2",
        "chi2_D",
        "chi2_He4",
        "chi2_Li_measurement",
        "chi2_Li_lower_bound",
        "chi2_eta_cmb",
        "chi2_tau",
        "chi2_total_li_measurement",
        "chi2_total_li_lower_bound",
        *reaction_columns,
    ]
    numeric(frame, numeric_columns)
    expected_ids = reconstructed.pop("expected_ids")
    reconstructed.pop("final_points")

    decomposition = score_decomposition(frame, reaction_columns)
    metrics = {
        "generated_utc": datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "inputs": {
            "scan_program": str(cli.scan_program),
            "scan_program_sha256": sha256(cli.scan_program),
            "points_csv": str(points_path),
            "points_csv_sha256_at_audit": sha256(points_path),
            "config_json": str(config_path),
            "config_json_sha256_at_audit": sha256(config_path),
            "linx_path": config.get("linx_path"),
            "linx_git_revision": git_revision(Path(config["linx_path"])),
            "linx_python_environment": python_environment(cli.linx_python),
            "network": config.get("network"),
            "rtol": 1e-5,
            "atol": 1e-9,
            "sampling_key": 150,
        },
        "catalogue": {
            "n_raw_rows": int(len(raw)),
            "n_unique_ok": int(len(ok_ids)),
            "n_unresolved": int(len(unresolved)),
            "n_expected": int(len(expected_ids)),
            "n_expected_missing": int(len(expected_ids - ok_ids)),
            "n_unexpected_ok": int(len(ok_ids - expected_ids)),
            **reconstructed,
        },
        "score_integrity": verify_scores(frame, reaction_columns),
        "proposal": proposal_metrics(frame, ["q_He3aBe7g", "q_Be7nLi7p", "q_Li7paa"]),
        "family_counts_observed": dict(sorted(Counter(frame["family"]).items())),
    }

    metrics_path = cli.outdir / "fr_linx_scan_audit_metrics_2026-07-17.json"
    decomposition_path = cli.outdir / "fr_linx_scan_audit_score_decomposition_2026-07-17.csv"
    report_path = cli.outdir / "fr_linx_scan_empirical_audit_2026-07-17.md"
    sensitivity_path = cli.outdir / "fr_linx_anchor_sensitivity_2026-07-17.csv"
    sensitivity_report_path = cli.outdir / "fr_linx_anchor_sensitivity_report_2026-07-17.md"
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    decomposition.to_csv(decomposition_path, index=False)
    write_report(report_path, metrics, decomposition)
    sensitivity = anchor_sensitivity(frame)
    sensitivity.to_csv(sensitivity_path, index=False)
    write_anchor_report(sensitivity_report_path, sensitivity)
    manifest_path = cli.outdir / "fr_linx_scan_audit_manifest_2026-07-17.csv"
    detailed_report = cli.outdir / "fr_linx_analysis_program_detailed_audit_2026-07-17.md"
    archived_config = cli.outdir / "fr_linx_scan_config_2026-07-17.json"
    entries = [
        ("scan_program", cli.scan_program),
        ("audit_program", SCRIPT_PATH),
        ("scan_config", config_path),
        ("scan_catalogue_snapshot", points_path),
        ("audit_metrics", metrics_path),
        ("score_decomposition", decomposition_path),
        ("empirical_report", report_path),
        ("anchor_sensitivity", sensitivity_path),
        ("anchor_sensitivity_report", sensitivity_report_path),
    ]
    if detailed_report.exists():
        entries.append(("detailed_audit_report", detailed_report))
    if archived_config.exists():
        entries.append(("scan_config_archived", archived_config))
    for result_path in sorted(cli.outdir.glob("su2_bbn_lithium_*_2026-07-17.*")):
        entries.append((f"final_gate_frontier:{result_path.name}", result_path))
    write_manifest(manifest_path, entries)
    print(f"Saved metrics: {metrics_path}")
    print(f"Saved score decomposition: {decomposition_path}")
    print(f"Saved report: {report_path}")
    print(f"Saved anchor sensitivity: {sensitivity_path}")
    print(f"Saved anchor sensitivity report: {sensitivity_report_path}")
    print(f"Saved manifest: {manifest_path}")


if __name__ == "__main__":
    main()
