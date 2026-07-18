#!/usr/bin/env python3
r"""Post hoc empirical-null audit of the locked group power calibration."""

from __future__ import annotations

import hashlib
import json
import math
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
RUN_DATE = "2026-07-18"
RUN_DIR = ROOT / "plamb_runs" / "diagnostics" / "sparc_fr_group_power_calibration_20260718"
EXPORT_DIR = (
    ROOT
    / "github_export"
    / "results"
    / "2026-07-18"
    / "fr_group_power_calibration"
)
LOCKED_PROGRAMME = (
    ROOT
    / "github_export"
    / "code"
    / "sparc"
    / "run_sparc_fr_group_power_calibration_2026-07-18.py"
)
PRIMARY_FRAME = (
    ROOT
    / "github_export"
    / "results"
    / "2026-07-18"
    / "fr_group_catalogue_replication"
    / "sparc_fr_group_catalogue_primary_frame.csv"
)
GROUP_PROGRAMME = (
    ROOT
    / "github_export"
    / "code"
    / "sparc"
    / "run_sparc_fr_group_catalogue_replication_2026-07-18.py"
)
GROUP_PREREGISTRATION = PRIMARY_FRAME.parent / "sparc_fr_group_catalogue_preregistration.json"
GROUP_GATE = PRIMARY_FRAME.parent / "sparc_fr_group_catalogue_gate.json"
GROUP_SUMMARY = PRIMARY_FRAME.parent / "sparc_fr_group_catalogue_test_summary.csv"
PRIOR_ENVIRONMENT_SUMMARY = (
    ROOT
    / "github_export"
    / "results"
    / "2026-07-18"
    / "fr_environment_asymmetry"
    / "sparc_fr_environment_test_summary.csv"
)
PRIOR_EFFECT_ANCHOR = 0.307254


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


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


def make_plot(frame: pd.DataFrame, path: Path) -> None:
    x = frame["target_abs_residual_correlation"].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    ax.plot(
        x,
        frame["locked_joint_primary_recovery"],
        marker="o",
        linewidth=2.0,
        color="#2f6f6d",
        label="locked joint recovery",
    )
    ax.plot(
        x,
        frame["empirical_null_joint_recovery"],
        marker="s",
        linewidth=2.0,
        color="#b54435",
        label="empirical-null joint recovery",
    )
    ax.plot(
        x,
        frame["empirical_null_p_recovery"],
        marker="^",
        linewidth=1.6,
        color="#586994",
        label="empirical-null p component",
    )
    ax.axhline(0.80, color="#777777", linestyle="--", linewidth=1.0)
    ax.axvline(PRIOR_EFFECT_ANCHOR, color="#777777", linestyle=":", linewidth=1.0)
    ax.set_xlim(-0.01, max(x) + 0.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Injected absolute host-residual correlation")
    ax.set_ylabel("Recovery probability")
    ax.set_title("SPARC Group Power: Empirical-Null Audit")
    ax.grid(alpha=0.18)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_manifest() -> None:
    rows: list[dict[str, Any]] = []
    for path in sorted(EXPORT_DIR.iterdir()):
        if not path.is_file() or path.name == "manifest.csv":
            continue
        rows.append(
            {
                "path": path.name,
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
                "role": "committed power-calibration product",
                "tracked_in_git": True,
            }
        )
    for path, role in [
        (LOCKED_PROGRAMME, "locked power-calibration programme"),
        (Path(__file__).resolve(), "post hoc empirical-null audit programme"),
        (PRIMARY_FRAME, "locked 65-galaxy primary frame"),
        (GROUP_PROGRAMME, "locked group-replication programme"),
        (GROUP_PREREGISTRATION, "locked group-replication preregistration"),
        (GROUP_GATE, "locked failed group-replication gate"),
        (GROUP_SUMMARY, "locked group-replication summary"),
        (PRIOR_ENVIRONMENT_SUMMARY, "earlier 2MRS environment summary"),
    ]:
        rows.append(
            {
                "path": str(path.relative_to(ROOT)).replace("\\", "/"),
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
                "role": role,
                "tracked_in_git": True,
            }
        )
    pd.DataFrame(rows).to_csv(EXPORT_DIR / "manifest.csv", index=False)


def main() -> None:
    config = json.loads(
        (RUN_DIR / "sparc_fr_group_power_preregistration.json").read_text(encoding="utf-8")
    )
    if sha256(LOCKED_PROGRAMME) != config["programme_sha256"]:
        raise RuntimeError("Locked power programme no longer matches its preregistration")
    draws = pd.read_csv(RUN_DIR / "sparc_fr_group_power_simulation_draws.csv.gz")
    locked = pd.read_csv(RUN_DIR / "sparc_fr_group_power_summary.csv")
    null = draws[np.isclose(draws["target_abs_residual_correlation"], 0.0)].copy()
    empirical_critical = float(
        np.quantile(np.abs(null["partial_spearman_rho"]), 0.99, method="higher")
    )
    draws["empirical_null_p_passed"] = (
        np.abs(draws["partial_spearman_rho"]) >= empirical_critical
    )
    draws["empirical_null_joint_passed"] = (
        draws["empirical_null_p_passed"]
        & draws["cv_gate_passed"].astype(bool)
        & draws["sign_gate_passed"].astype(bool)
    )
    rows: list[dict[str, Any]] = []
    for target, selected in draws.groupby("target_abs_residual_correlation", sort=True):
        locked_row = locked[
            np.isclose(locked["target_abs_residual_correlation"], target)
        ].iloc[0]
        p_rate = float(selected["empirical_null_p_passed"].mean())
        joint_rate = float(selected["empirical_null_joint_passed"].mean())
        count = len(selected)
        rows.append(
            {
                "target_abs_residual_correlation": float(target),
                "replicates": count,
                "locked_p_recovery": float(locked_row["p_component_recovery"]),
                "empirical_null_p_recovery": p_rate,
                "empirical_null_p_mc_se": math.sqrt(p_rate * (1.0 - p_rate) / count),
                "locked_joint_primary_recovery": float(
                    locked_row["joint_primary_recovery"]
                ),
                "empirical_null_joint_recovery": joint_rate,
                "empirical_null_joint_mc_se": math.sqrt(
                    joint_rate * (1.0 - joint_rate) / count
                ),
            }
        )
    recalibrated = pd.DataFrame(rows)
    anchor = recalibrated[
        np.isclose(recalibrated["target_abs_residual_correlation"], PRIOR_EFFECT_ANCHOR)
    ].iloc[0]
    qualifying = recalibrated[recalibrated["empirical_null_joint_recovery"] >= 0.80]
    threshold = (
        float(qualifying["target_abs_residual_correlation"].min())
        if len(qualifying)
        else float("nan")
    )
    recalibrated.to_csv(
        RUN_DIR / "sparc_fr_group_power_empirical_null_recalibration.csv", index=False
    )
    make_plot(
        recalibrated,
        RUN_DIR / "sparc_fr_group_power_empirical_null_recalibration.png",
    )
    columns = [
        "target_abs_residual_correlation",
        "locked_p_recovery",
        "empirical_null_p_recovery",
        "locked_joint_primary_recovery",
        "empirical_null_joint_recovery",
    ]
    report = [
        "# SPARC Group-Power Empirical-Null Audit",
        "",
        f"Date: {RUN_DATE}",
        f"Completed: `{utc_now()}`",
        "",
        "## Status",
        "",
        "This is a post hoc calibration diagnostic. It preserves the locked simulation and cannot alter its preregistered result.",
        "",
        "## Finding",
        "",
        f"The locked rank-permutation threshold was |rho|>=`{float(json.loads((RUN_DIR / 'sparc_fr_group_power_diagnostics.json').read_text(encoding='utf-8'))['critical_absolute_spearman_rho_for_p_le_0p01']):.6g}`. Under the complete empirical-residual null pipeline, its p-component false-positive rate was `{100.0 * float(locked.iloc[0]['p_component_recovery']):.2f}` per cent rather than 1 per cent.",
        "",
        f"The empirical 1 per cent threshold is |rho|>=`{empirical_critical:.6g}`. At the earlier-effect anchor, p-component recovery changes from `{100.0 * float(anchor['locked_p_recovery']):.2f}` to `{100.0 * float(anchor['empirical_null_p_recovery']):.2f}` per cent, and joint recovery changes from `{100.0 * float(anchor['locked_joint_primary_recovery']):.2f}` to `{100.0 * float(anchor['empirical_null_joint_recovery']):.2f}` per cent.",
        "",
        f"The smallest tested effect with at least 80 per cent empirically recalibrated joint recovery remains `{threshold:.6g}`.",
        "",
        "## Recalibrated Grid",
        "",
        markdown_table(recalibrated, columns),
        "",
        "## Interpretation",
        "",
        "The rank-permutation p component is anti-conservative for this empirical-noise and cross-fitting pipeline. The full locked joint gate remains conservative under the null because it also requires a 5 per cent predictive improvement: its locked null pass rate is 0.54 per cent, while the empirically recalibrated joint null pass rate is 0.32 per cent.",
        "",
        "The correction strengthens the conclusion that the 65-galaxy group-richness design is underpowered at the earlier-effect anchor. It does not change the observed catalogue result, whose p-value was 0.704, and it creates no evidence for FR or antimatter.",
        "",
    ]
    report_path = RUN_DIR / "sparc_fr_group_power_empirical_null_audit.md"
    report_path.write_text("\n".join(report), encoding="utf-8")
    for name in [
        "sparc_fr_group_power_empirical_null_recalibration.csv",
        "sparc_fr_group_power_empirical_null_recalibration.png",
        "sparc_fr_group_power_empirical_null_audit.md",
    ]:
        shutil.copy2(RUN_DIR / name, EXPORT_DIR / name)
    write_manifest()
    print(f"Saved empirical-null audit: {report_path}")
    print(f"Empirical p critical |rho|: {empirical_critical:.6f}")
    print(f"Anchor empirical joint recovery: {float(anchor['empirical_null_joint_recovery']):.6f}")


if __name__ == "__main__":
    main()
