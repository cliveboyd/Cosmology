#!/usr/bin/env python3
"""Build an SU2-compatible gate from the selected-network LINX BBN scan."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
TASK_ROOT = Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu")
DEFAULT_SCAN_DIR = (
    ROOT
    / "plamb_runs"
    / "diagnostics"
    / "bbn_lithium_linx_key_fr_network_overnight_20260716"
)
DEFAULT_OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "su2_bbn_lithium_gate_20260717"
DEFAULT_OUTPUTS = TASK_ROOT / "outputs"
RUN_DATE = "2026-07-17"

OBSERVED = {
    "D_H_1e5": (2.508, 0.029),
    "Yp_mass": (0.245, 0.003),
    "Li7_H_1e10": (1.45, 0.25),
    "eta_fac": (1.0, 0.04 / 6.12),
    "tau_n_fac": (1.0, 0.6 / 879.4),
}

REACTION_COLUMNS = ("q_He3aBe7g", "q_Be7nLi7p", "q_Li7paa")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-dir", type=Path, default=DEFAULT_SCAN_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--outputs", type=Path, default=DEFAULT_OUTPUTS)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def latest_status(points_path: Path) -> tuple[pd.DataFrame, int, int, int]:
    raw = pd.read_csv(points_path, low_memory=False)
    raw["_row_order"] = np.arange(len(raw))
    ok = raw.loc[raw["status"].eq("ok")].copy()
    ok = ok.sort_values("_row_order").drop_duplicates("point_id", keep="last")
    ok_ids = set(ok["point_id"].astype(str))
    attempted_ids = set(raw["point_id"].dropna().astype(str))
    failed_ids = set(raw.loc[raw["status"].eq("failed"), "point_id"].dropna().astype(str)) - ok_ids
    return ok.drop(columns="_row_order"), len(ok_ids), len(failed_ids), len(attempted_ids)


def within_sigma(values: pd.Series, key: str, n_sigma: float = 2.0) -> pd.Series:
    mean, sigma = OBSERVED[key]
    return (values.astype(float) - mean).abs() <= n_sigma * sigma


def finite_float(row: pd.Series | None, key: str) -> float:
    if row is None:
        return math.nan
    value = float(row[key])
    return value if math.isfinite(value) else math.nan


def scenario_row(name: str, description: str, frame: pd.DataFrame, abundance_gate: pd.Series) -> tuple[dict[str, Any], dict[str, Any]]:
    subset = frame.copy()
    joint = subset.loc[abundance_gate.reindex(subset.index, fill_value=False)].copy()
    li_gate = within_sigma(joint["Li7_H_1e10"], "Li7_H_1e10") if not joint.empty else pd.Series(dtype=bool)
    joint_li = joint.loc[li_gate] if not joint.empty else joint

    best = None if subset.empty else subset.loc[subset["chi2_total_li_measurement"].astype(float).idxmin()]
    best_lower = None if subset.empty else subset.loc[subset["chi2_total_li_lower_bound"].astype(float).idxmin()]
    min_li = None if joint.empty else joint.loc[joint["Li7_H_1e10"].astype(float).idxmin()]

    summary = {
        "scenario": name,
        "description": description,
        "n_rows": int(len(subset)),
        "n_D_He_2sigma": int(len(joint)),
        "n_D_He_Li_2sigma": int(len(joint_li)),
        "lithium_gate_pass": bool(len(joint_li) > 0),
        "best_chi2_li_measurement": finite_float(best, "chi2_total_li_measurement"),
        "best_chi2_li_lower_bound": finite_float(best_lower, "chi2_total_li_lower_bound"),
        "best_Li7_H_1e10": finite_float(best, "Li7_H_1e10"),
        "best_D_H_1e5": finite_float(best, "D_H_1e5"),
        "best_Yp": finite_float(best, "Yp_mass"),
        "minimum_Li7_H_1e10_given_D_He": finite_float(min_li, "Li7_H_1e10"),
        "minimum_li_eta_fac": finite_float(min_li, "eta_fac"),
        "minimum_li_clock_factor": finite_float(min_li, "clock_factor"),
        "minimum_li_delta_neff": finite_float(min_li, "delta_neff_initial"),
        "minimum_li_q_pull_chi2": finite_float(min_li, "q_pull_chi2"),
    }
    best_detail = {"scenario": name}
    for key in (
        "point_id",
        "family",
        "eta_fac",
        "clock_factor",
        "delta_neff_initial",
        "tau_n_fac",
        "D_H_1e5",
        "Yp_mass",
        "Li7_H_1e10",
        "stellar_depletion_needed",
        "chi2_total_li_measurement",
        "chi2_total_li_lower_bound",
        "q_pull_chi2",
        *REACTION_COLUMNS,
    ):
        best_detail[key] = None if best is None else best.get(key)
    return summary, best_detail


def influence_rows(frame: pd.DataFrame) -> list[dict[str, Any]]:
    random_rows = frame.loc[frame["family"].eq("random_targeted_li_be")].copy()
    if len(random_rows) < 20:
        return []
    predictors = ["eta_fac", "clock_factor", "tau_n_fac", *REACTION_COLUMNS]
    outcomes = ["D_H_1e5", "Yp_mass", "Li7_H_1e10"]
    x = random_rows[predictors].astype(float).to_numpy()
    x_std = np.std(x, axis=0, ddof=1)
    valid_predictor = x_std > 0
    x_z = (x[:, valid_predictor] - np.mean(x[:, valid_predictor], axis=0)) / x_std[valid_predictor]
    x_design = np.column_stack([np.ones(len(x_z)), x_z])
    kept = np.asarray(predictors)[valid_predictor]
    rows: list[dict[str, Any]] = []
    for outcome in outcomes:
        y = random_rows[outcome].astype(float).to_numpy()
        y_sd = float(np.std(y, ddof=1))
        if not math.isfinite(y_sd) or y_sd == 0:
            continue
        y_z = (y - np.mean(y)) / y_sd
        beta, *_ = np.linalg.lstsq(x_design, y_z, rcond=None)
        fitted = x_design @ beta
        r2 = 1.0 - float(np.sum((y_z - fitted) ** 2) / np.sum((y_z - np.mean(y_z)) ** 2))
        for predictor, coefficient in zip(kept, beta[1:], strict=True):
            rows.append(
                {
                    "outcome": outcome,
                    "predictor": predictor,
                    "standardised_linear_coefficient": float(coefficient),
                    "model_r2": r2,
                    "n_rows": int(len(random_rows)),
                }
            )
    return rows


def fmt(value: Any, digits: int = 4) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
    return f"{number:.{digits}g}" if math.isfinite(number) else "NA"


def write_report(
    path: Path,
    summaries: list[dict[str, Any]],
    n_ok: int,
    n_unresolved: int,
    n_attempted: int,
    requested: int,
) -> None:
    by_name = {row["scenario"]: row for row in summaries}
    expansion = by_name["su2_expansion_only"]
    modest = by_name["su2_plus_modest_rate_controls"]
    unrestricted = by_name["su2_plus_scanned_rate_controls"]
    complete = n_ok == requested and n_unresolved == 0
    status = "complete" if complete else "interim"
    lines = [
        "# SU2-Compatible BBN Lithium Gate",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        f"Status: **{status}**; successful unique LINX points `{n_ok} / {requested}`, attempted `{n_attempted}`, unresolved `{n_unresolved}`.",
        "",
        "## Scope",
        "",
        "This is an effective SU2 gate, not a first-principles SU2 nucleosynthesis calculation. A non-negative SU2-like early-radiation contribution is represented only through an expansion factor and its equivalent initial radiation shift:",
        "",
        r"\[",
        r"\begin{aligned}",
        r"S                    &= H_{\rm BBN}/H_{\rm SBBN}, \\",
        r"\Delta N_{\rm eff}  &= \frac{43}{7}\left(S^2-1\right), \\",
        r"\eta_{\rm fac}      &= \eta_{\rm BBN}/\eta_{\rm CMB}, \\",
        r"\tau_{n,{\rm fac}}  &= \tau_n/\tau_{n,0}.",
        r"\end{aligned}",
        r"\]",
        "",
        "The SU2-compatible rows require `S >= 1`, two-sigma CMB consistency in `eta_fac`, and two-sigma neutron-lifetime consistency. D/H and He-4 must each pass their two-sigma abundance gates before lithium is assessed. Nuclear-rate pulls are controls and are not attributed to SU2 without an explicit coupling model.",
        "",
        "## Scenario Gates",
        "",
        "| Scenario | Rows | D+He pass | D+He+Li pass | Best chi2, Li measurement | Min Li given D+He |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in summaries:
        lines.append(
            f"| {row['scenario']} | {row['n_rows']} | {row['n_D_He_2sigma']} | "
            f"{row['n_D_He_Li_2sigma']} | {fmt(row['best_chi2_li_measurement'])} | "
            f"{fmt(row['minimum_Li7_H_1e10_given_D_He'])} |"
        )
    lines.extend(
        [
            "",
            "## Readout",
            "",
            f"- Expansion-only proxy: best measurement chi-squared `{fmt(expansion['best_chi2_li_measurement'])}`; lithium gate pass `{expansion['lithium_gate_pass']}`.",
            f"- With modest joint rate controls (`sum(q_i^2) <= 9`): best chi-squared `{fmt(modest['best_chi2_li_measurement'])}`; lithium gate pass `{modest['lithium_gate_pass']}`.",
            f"- With the full scanned rate-control range: best chi-squared `{fmt(unrestricted['best_chi2_li_measurement'])}`; lithium gate pass `{unrestricted['lithium_gate_pass']}`.",
            "",
            "A lithium pass requires simultaneous two-sigma agreement with D/H, He-4 and the lithium plateau. A lower total chi-squared that still misses this joint gate is not a lithium solution.",
            "",
            "## Claim Boundary",
            "",
            "The scan can veto an expansion-only SU2 explanation within the tested effective-clock mapping. It cannot validate or exclude a fully specified SU2 model that changes binding energies, neutron-proton mass splitting, weak rates and nuclear Q-values self-consistently. The selected LINX reaction pulls remain survey-style nuisance controls for nuclear physics, not evidence for SU2.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    cli = parse_args()
    cli.outdir.mkdir(parents=True, exist_ok=True)
    cli.outputs.mkdir(parents=True, exist_ok=True)
    points_path = cli.scan_dir / "bbn_lithium_linx_fr_network_points.csv"
    config_path = cli.scan_dir / "bbn_lithium_linx_fr_network_config.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    requested = int(config.get("n_points_requested", 30503))
    # Chunk runs temporarily write smaller maxima; the registered full target remains fixed.
    requested = max(requested, 30503)

    frame, n_ok, n_unresolved, n_attempted = latest_status(points_path)
    for key in (
        "eta_fac",
        "clock_factor",
        "tau_n_fac",
        "D_H_1e5",
        "Yp_mass",
        "Li7_H_1e10",
        "q_pull_chi2",
        "chi2_total_li_measurement",
        "chi2_total_li_lower_bound",
    ):
        frame[key] = pd.to_numeric(frame[key], errors="coerce")

    base = (
        frame["clock_factor"].ge(1.0)
        & within_sigma(frame["eta_fac"], "eta_fac")
        & within_sigma(frame["tau_n_fac"], "tau_n_fac")
    )
    abundance_gate = within_sigma(frame["D_H_1e5"], "D_H_1e5") & within_sigma(frame["Yp_mass"], "Yp_mass")
    scenarios = [
        (
            "su2_expansion_only",
            "Non-negative early-radiation/clock proxy with no nuclear-rate pull.",
            frame.loc[base & frame["q_pull_chi2"].le(1e-12)],
        ),
        (
            "su2_plus_modest_rate_controls",
            "Same SU2 proxy with joint selected-rate pull sum(q_i^2) <= 9.",
            frame.loc[base & frame["q_pull_chi2"].le(9.0)],
        ),
        (
            "su2_plus_scanned_rate_controls",
            "Same SU2 proxy with the complete scanned selected-rate range.",
            frame.loc[base],
        ),
        (
            "all_scanned_controls",
            "All successful LINX rows, including negative clock shifts; diagnostic control only.",
            frame,
        ),
    ]

    summaries: list[dict[str, Any]] = []
    best_rows: list[dict[str, Any]] = []
    for name, description, subset in scenarios:
        summary, best = scenario_row(name, description, subset, abundance_gate)
        summaries.append(summary)
        best_rows.append(best)

    summary_path = cli.outdir / "su2_bbn_lithium_gate_summary.csv"
    best_path = cli.outdir / "su2_bbn_lithium_gate_best_rows.csv"
    influence_path = cli.outdir / "su2_bbn_lithium_gate_linear_influence.csv"
    report_path = cli.outdir / "su2_bbn_lithium_gate_report.md"
    pd.DataFrame(summaries).to_csv(summary_path, index=False)
    pd.DataFrame(best_rows).to_csv(best_path, index=False)
    pd.DataFrame(influence_rows(frame)).to_csv(influence_path, index=False)
    write_report(report_path, summaries, n_ok, n_unresolved, n_attempted, requested)

    audit = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "status": "complete" if n_ok == requested and n_unresolved == 0 else "interim",
        "successful_unique_points": n_ok,
        "unresolved_points": n_unresolved,
        "attempted_unique_points": n_attempted,
        "requested_points": requested,
        "input_points_file": str(points_path),
        "input_points_sha256": sha256(points_path),
        "source_config": str(config_path),
        "source_config_sha256": sha256(config_path),
        "programme": str(Path(__file__).resolve()),
        "programme_sha256": sha256(Path(__file__).resolve()),
    }
    audit_path = cli.outdir / "su2_bbn_lithium_gate_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")

    for source in (summary_path, best_path, influence_path, report_path, audit_path):
        destination = cli.outputs / f"{source.stem}_{RUN_DATE}{source.suffix}"
        shutil.copy2(source, destination)
        print(f"Saved: {destination}")


if __name__ == "__main__":
    main()
