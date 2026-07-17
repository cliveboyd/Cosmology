#!/usr/bin/env python3
"""Quantify the SU2-compatible BBN lithium near-miss frontier."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
DEFAULT_OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "su2_bbn_lithium_frontier_20260717"
DEFAULT_OUTPUTS = TASK_ROOT / "outputs"
RUN_DATE = "2026-07-17"
REQUESTED_POINTS = 30_503

OBSERVED = {
    "D_H_1e5": (2.508, 0.029),
    "Yp_mass": (0.245, 0.003),
    "Li7_H_1e10": (1.45, 0.25),
    "eta_fac": (1.0, 0.04 / 6.12),
    "tau_n_fac": (1.0, 0.6 / 879.4),
}

NUMERIC_COLUMNS = (
    "eta_fac",
    "clock_factor",
    "delta_neff_initial",
    "tau_n_fac",
    "q_pull_chi2",
    "D_H_1e5",
    "Yp_mass",
    "Li7_H_1e10",
    "chi2_total_li_measurement",
    "chi2_total_li_lower_bound",
)


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
    unresolved_ids = set(raw.loc[raw["status"].eq("failed"), "point_id"].dropna().astype(str)) - ok_ids
    return ok.drop(columns="_row_order"), len(ok_ids), len(unresolved_ids), len(attempted_ids)


def add_diagnostics(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    for column in NUMERIC_COLUMNS:
        result[column] = pd.to_numeric(result[column], errors="coerce")
    for column in ("D_H_1e5", "Yp_mass", "Li7_H_1e10", "eta_fac", "tau_n_fac"):
        mean, sigma = OBSERVED[column]
        result[f"z_{column}"] = (result[column] - mean) / sigma
    result["chi2_D_He"] = result["z_D_H_1e5"] ** 2 + result["z_Yp_mass"] ** 2
    result["max_abs_z_D_He"] = result[["z_D_H_1e5", "z_Yp_mass"]].abs().max(axis=1)
    result["li_depletion_to_centre"] = result["Li7_H_1e10"] / OBSERVED["Li7_H_1e10"][0]
    result["li_suppression_to_2sigma_upper"] = result["Li7_H_1e10"] / (
        OBSERVED["Li7_H_1e10"][0] + 2.0 * OBSERVED["Li7_H_1e10"][1]
    )
    return result


def scenario_frames(frame: pd.DataFrame) -> dict[str, pd.DataFrame]:
    cmb = frame["z_eta_fac"].abs().le(2.0)
    tau = frame["z_tau_n_fac"].abs().le(2.0)
    positive_expansion = frame["clock_factor"].ge(1.0)
    base = cmb & tau & positive_expansion
    return {
        "su2_expansion_only": frame.loc[base & frame["q_pull_chi2"].le(1e-12)].copy(),
        "su2_plus_modest_rate_controls": frame.loc[base & frame["q_pull_chi2"].le(9.0)].copy(),
        "su2_plus_scanned_rate_controls": frame.loc[base].copy(),
        "all_scanned_controls": frame.copy(),
    }


def finite(value: Any) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def tolerance_profiles(scenarios: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries: list[dict[str, Any]] = []
    best_rows: list[dict[str, Any]] = []
    for scenario, subset in scenarios.items():
        for tolerance in (1.0, 2.0, 3.0):
            accepted = subset.loc[subset["max_abs_z_D_He"].le(tolerance)].copy()
            best = None if accepted.empty else accepted.loc[accepted["Li7_H_1e10"].idxmin()]
            minimum_li = finite(None if best is None else best["Li7_H_1e10"])
            summaries.append(
                {
                    "scenario": scenario,
                    "D_He_each_sigma_limit": tolerance,
                    "n_scenario_rows": len(subset),
                    "n_D_He_pass": len(accepted),
                    "minimum_Li7_H_1e10": minimum_li,
                    "minimum_Li_z": (
                        (minimum_li - OBSERVED["Li7_H_1e10"][0]) / OBSERVED["Li7_H_1e10"][1]
                        if math.isfinite(minimum_li)
                        else math.nan
                    ),
                    "minimum_depletion_to_plateau_centre": finite(
                        None if best is None else best["li_depletion_to_centre"]
                    ),
                    "minimum_suppression_to_2sigma_upper": finite(
                        None if best is None else best["li_suppression_to_2sigma_upper"]
                    ),
                    "best_delta_neff_initial": finite(None if best is None else best["delta_neff_initial"]),
                    "best_clock_factor": finite(None if best is None else best["clock_factor"]),
                    "best_eta_fac": finite(None if best is None else best["eta_fac"]),
                    "best_q_pull_chi2": finite(None if best is None else best["q_pull_chi2"]),
                }
            )
            if best is not None:
                best_rows.append(
                    {
                        "scenario": scenario,
                        "D_He_each_sigma_limit": tolerance,
                        **{
                            key: best.get(key)
                            for key in (
                                "point_id",
                                "family",
                                "eta_fac",
                                "clock_factor",
                                "delta_neff_initial",
                                "tau_n_fac",
                                "q_pull_chi2",
                                "D_H_1e5",
                                "Yp_mass",
                                "Li7_H_1e10",
                                "z_D_H_1e5",
                                "z_Yp_mass",
                                "z_Li7_H_1e10",
                                "li_depletion_to_centre",
                                "li_suppression_to_2sigma_upper",
                            )
                        },
                    }
                )
    return pd.DataFrame(summaries), pd.DataFrame(best_rows)


def delta_neff_profiles(scenarios: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    edges = np.arange(0.0, 1.2001, 0.1)
    for scenario in (
        "su2_expansion_only",
        "su2_plus_modest_rate_controls",
        "su2_plus_scanned_rate_controls",
    ):
        subset = scenarios[scenario]
        for low, high in zip(edges[:-1], edges[1:], strict=True):
            in_bin = subset.loc[
                subset["delta_neff_initial"].ge(low)
                & subset["delta_neff_initial"].lt(high if high < edges[-1] else high + 1e-12)
            ]
            accepted = in_bin.loc[in_bin["max_abs_z_D_He"].le(2.0)]
            rows.append(
                {
                    "scenario": scenario,
                    "delta_neff_low": low,
                    "delta_neff_high": high,
                    "n_rows": len(in_bin),
                    "n_D_He_2sigma": len(accepted),
                    "minimum_chi2_D_He": finite(in_bin["chi2_D_He"].min()),
                    "minimum_Li7_H_1e10_given_D_He": finite(accepted["Li7_H_1e10"].min()),
                    "median_Li7_H_1e10_given_D_He": finite(accepted["Li7_H_1e10"].median()),
                    "minimum_q_pull_chi2_given_D_He": finite(accepted["q_pull_chi2"].min()),
                }
            )
    return pd.DataFrame(rows)


def pareto_frontier(subset: pd.DataFrame, scenario: str) -> pd.DataFrame:
    ordered = subset.sort_values(["chi2_D_He", "Li7_H_1e10", "q_pull_chi2"])
    selected: list[pd.Series] = []
    best_li = math.inf
    for _, row in ordered.iterrows():
        lithium = finite(row["Li7_H_1e10"])
        if lithium < best_li:
            selected.append(row)
            best_li = lithium
    if not selected:
        return pd.DataFrame()
    frontier = pd.DataFrame(selected)
    keep = [
        "point_id",
        "family",
        "eta_fac",
        "clock_factor",
        "delta_neff_initial",
        "tau_n_fac",
        "q_pull_chi2",
        "D_H_1e5",
        "Yp_mass",
        "Li7_H_1e10",
        "z_D_H_1e5",
        "z_Yp_mass",
        "z_Li7_H_1e10",
        "chi2_D_He",
        "max_abs_z_D_He",
        "li_depletion_to_centre",
        "li_suppression_to_2sigma_upper",
    ]
    frontier = frontier[keep].copy()
    frontier.insert(0, "scenario", scenario)
    frontier.insert(1, "frontier_rank", np.arange(1, len(frontier) + 1))
    return frontier


def make_plot(frame: pd.DataFrame, scenarios: dict[str, pd.DataFrame], path: Path) -> None:
    subset = scenarios["su2_plus_scanned_rate_controls"]
    accepted = subset["max_abs_z_D_He"].le(2.0)
    colour = np.log10(1.0 + subset["q_pull_chi2"].clip(lower=0.0))
    li_mean, li_sigma = OBSERVED["Li7_H_1e10"]

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), constrained_layout=True)
    axes[0].scatter(
        subset["delta_neff_initial"],
        subset["Li7_H_1e10"],
        c=colour,
        cmap="viridis",
        s=9,
        alpha=0.35,
        linewidths=0,
    )
    axes[0].scatter(
        subset.loc[accepted, "delta_neff_initial"],
        subset.loc[accepted, "Li7_H_1e10"],
        facecolors="none",
        edgecolors="#c43c39",
        s=22,
        linewidths=0.7,
        label="D and He within 2 sigma",
    )
    axes[0].axhspan(li_mean - 2 * li_sigma, li_mean + 2 * li_sigma, color="#e6b84a", alpha=0.25)
    axes[0].axhline(li_mean, color="#8b5d00", linewidth=1.2, label="Lithium plateau")
    axes[0].set_xlabel(r"Effective $\Delta N_{\rm eff}$ at BBN")
    axes[0].set_ylabel(r"$^{7}{\rm Li}/{\rm H}\;[10^{-10}]$")
    axes[0].set_title("Expansion and lithium")
    axes[0].legend(frameon=False, fontsize=8)

    scatter = axes[1].scatter(
        subset["max_abs_z_D_He"],
        subset["Li7_H_1e10"],
        c=colour,
        cmap="viridis",
        s=11,
        alpha=0.5,
        linewidths=0,
    )
    axes[1].axvline(2.0, color="#c43c39", linewidth=1.2, linestyle="--", label="D and He gate")
    axes[1].axhspan(li_mean - 2 * li_sigma, li_mean + 2 * li_sigma, color="#e6b84a", alpha=0.25)
    axes[1].axhline(li_mean, color="#8b5d00", linewidth=1.2)
    axes[1].set_xlabel(r"$\max(|z_{\rm D}|, |z_{\rm He}|)$")
    axes[1].set_ylabel(r"$^{7}{\rm Li}/{\rm H}\;[10^{-10}]$")
    axes[1].set_title("Abundance near-miss frontier")
    axes[1].set_xlim(left=0)
    axes[1].legend(frameon=False, fontsize=8)
    colour_bar = fig.colorbar(scatter, ax=axes, shrink=0.86)
    colour_bar.set_label(r"$\log_{10}(1+\sum q_i^2)$")
    fig.suptitle("SU2-compatible effective BBN lithium diagnostic", fontsize=13)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def fmt(value: Any, digits: int = 4) -> str:
    number = finite(value)
    return f"{number:.{digits}g}" if math.isfinite(number) else "NA"


def write_report(
    path: Path,
    tolerance: pd.DataFrame,
    n_ok: int,
    n_unresolved: int,
    n_attempted: int,
) -> None:
    rows_2sigma = tolerance.loc[tolerance["D_He_each_sigma_limit"].eq(2.0)]
    by_name = rows_2sigma.set_index("scenario")
    complete = n_ok == REQUESTED_POINTS and n_unresolved == 0
    lines = [
        "# SU2 BBN Lithium Near-Miss Frontier",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        f"Status: **{'complete' if complete else 'interim'}**; successful unique points `{n_ok} / {REQUESTED_POINTS}`, attempted `{n_attempted}`, unresolved `{n_unresolved}`.",
        "",
        "## Scope",
        "",
        "This is an exploratory diagnostic downstream of the separately registered SU2 completion gate. It does not alter that gate. It asks how close the sampled effective expansion model can move lithium while retaining deuterium, helium, Cosmic Microwave Background baryon-density and neutron-lifetime consistency.",
        "",
        "The standardised residuals and frontier quantities are",
        "",
        r"\[",
        r"\begin{aligned}",
        r"z_i                         &= (X_i-\mu_i)/\sigma_i, \\",
        r"\chi^2_{\rm D+He}           &= z_{\rm D}^{2}+z_{\rm He}^{2}, \\",
        r"f_{\rm dep}                 &= ({}^7{\rm Li/H})_{\rm BBN}/\mu_{\rm Li}, \\",
        r"f_{2\sigma}                 &= ({}^7{\rm Li/H})_{\rm BBN}/(\mu_{\rm Li}+2\sigma_{\rm Li}).",
        r"\end{aligned}",
        r"\]",
        "",
        "The abundance tolerances below require D/H and He-4 each, rather than only their combined chi-squared, to lie inside the stated interval. The uncertainties are the fixed observational errors used by the registered LINX scan; this table does not add theory-error floors.",
        "",
        "## Two-Sigma Frontier",
        "",
        "| Scenario | Rows passing D+He | Minimum Li/H x 1e10 | Li tension | Depletion to centre | Suppression to 2-sigma upper |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    labels = {
        "su2_expansion_only": "Expansion only",
        "su2_plus_modest_rate_controls": "Expansion + modest rate controls",
        "su2_plus_scanned_rate_controls": "Expansion + all scanned rate controls",
        "all_scanned_controls": "All scanned controls",
    }
    for scenario, row in by_name.iterrows():
        lines.append(
            f"| {labels[scenario]} | {int(row['n_D_He_pass'])} | {fmt(row['minimum_Li7_H_1e10'])} | "
            f"{fmt(row['minimum_Li_z'])} sigma | {fmt(row['minimum_depletion_to_plateau_centre'])} | "
            f"{fmt(row['minimum_suppression_to_2sigma_upper'])} |"
        )

    expansion = by_name.loc["su2_expansion_only"]
    modest = by_name.loc["su2_plus_modest_rate_controls"]
    full = by_name.loc["su2_plus_scanned_rate_controls"]
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"- Expansion only bottoms out at Li/H x 1e10 = `{fmt(expansion['minimum_Li7_H_1e10'])}` after D and He pass, still `{fmt(expansion['minimum_Li_z'])}` observational sigma above the plateau centre.",
            f"- Modest selected-rate controls lower that floor to `{fmt(modest['minimum_Li7_H_1e10'])}`, requiring a residual depletion factor `{fmt(modest['minimum_depletion_to_plateau_centre'])}` to reach the plateau centre.",
            f"- The unrestricted selected-rate scan reaches `{fmt(full['minimum_Li7_H_1e10'])}`, but still requires suppression by a factor `{fmt(full['minimum_suppression_to_2sigma_upper'])}` merely to enter the upper two-sigma lithium interval.",
            "- A point on this frontier is a near miss, not a model preference. Nuclear-rate pull improvements remain controls unless a specified SU2 microphysical model predicts their signs and magnitudes.",
            "",
            "## Claim Boundary",
            "",
            "The result tests an effective non-negative expansion shift only. A first-principles SU2 claim would need a self-consistent temperature history and predicted changes to weak rates, binding energies, mass splittings and nuclear Q-values, followed by a new network calculation. The present frontier identifies the size and direction of that burden; it does not substitute for it.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    cli = parse_args()
    cli.outdir.mkdir(parents=True, exist_ok=True)
    cli.outputs.mkdir(parents=True, exist_ok=True)
    points_path = cli.scan_dir / "bbn_lithium_linx_fr_network_points.csv"
    frame, n_ok, n_unresolved, n_attempted = latest_status(points_path)
    frame = add_diagnostics(frame)
    scenarios = scenario_frames(frame)

    tolerance, best_rows = tolerance_profiles(scenarios)
    delta_profiles = delta_neff_profiles(scenarios)
    frontiers = pd.concat(
        [
            pareto_frontier(scenarios["su2_plus_modest_rate_controls"], "su2_plus_modest_rate_controls"),
            pareto_frontier(scenarios["su2_plus_scanned_rate_controls"], "su2_plus_scanned_rate_controls"),
        ],
        ignore_index=True,
    )

    files = {
        "tolerance": cli.outdir / "su2_bbn_lithium_tolerance_profile.csv",
        "best": cli.outdir / "su2_bbn_lithium_frontier_best_rows.csv",
        "delta": cli.outdir / "su2_bbn_lithium_delta_neff_profile.csv",
        "pareto": cli.outdir / "su2_bbn_lithium_pareto_frontier.csv",
        "plot": cli.outdir / "su2_bbn_lithium_frontier.png",
        "report": cli.outdir / "su2_bbn_lithium_frontier_report.md",
        "audit": cli.outdir / "su2_bbn_lithium_frontier_audit.json",
    }
    tolerance.to_csv(files["tolerance"], index=False)
    best_rows.to_csv(files["best"], index=False)
    delta_profiles.to_csv(files["delta"], index=False)
    frontiers.to_csv(files["pareto"], index=False)
    make_plot(frame, scenarios, files["plot"])
    write_report(files["report"], tolerance, n_ok, n_unresolved, n_attempted)

    audit = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "status": "complete" if n_ok == REQUESTED_POINTS and n_unresolved == 0 else "interim",
        "successful_unique_points": n_ok,
        "unresolved_points": n_unresolved,
        "attempted_unique_points": n_attempted,
        "requested_points": REQUESTED_POINTS,
        "input_points_file": str(points_path),
        "input_points_sha256": sha256(points_path),
        "programme": str(Path(__file__).resolve()),
        "programme_sha256": sha256(Path(__file__).resolve()),
        "analysis_role": "exploratory downstream diagnostic; registered gate unchanged",
    }
    files["audit"].write_text(json.dumps(audit, indent=2), encoding="utf-8")

    for source in files.values():
        destination = cli.outputs / f"{source.stem}_{RUN_DATE}{source.suffix}"
        shutil.copy2(source, destination)
        print(f"Saved: {destination}")


if __name__ == "__main__":
    main()
