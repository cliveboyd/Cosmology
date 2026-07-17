#!/usr/bin/env python3
"""Run and analyse the preregistered FR/LINX BBN validation matrix."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import sys
import time
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
DEFAULT_LINX_PATH = REPO_ROOT / "plamb_runs" / "deps" / "LINX"
DEFAULT_OUTDIR = (
    REPO_ROOT
    / "plamb_runs"
    / "diagnostics"
    / "bbn_lithium_linx_validation_matrix_20260717"
)
NETWORKS = ("key", "small", "full")
D_RATES = ("dpHe3g", "ddHe3n", "ddtp")
OMITTED_LI_BE_RATES = ("Li7paag", "Be7naa", "Be7daap", "Be7pB8g", "Li7daan")
STRICT_RTOL = 1e-7
STRICT_ATOL = 1e-11
PRIMARY_SAMPLING = 300
CHECK_SAMPLING = 600
TAU_REFERENCE_SECONDS = 879.4
PDG_TAU_FACTOR = 878.4 / TAU_REFERENCE_SECONDS
REGISTERED = {"D": (2.508, 0.029), "He": (0.245, 0.003), "Li": (1.45, 0.25), "tau": (879.4, 0.6)}
CURRENT = {"D": (2.533, 0.024), "He": (0.2458, 0.0013), "Li": (1.45, 0.25), "tau": (878.4, 0.5)}
ABUNDANCE_COLUMNS = ("D_H_1e5", "Yp_mass", "Li7_H_1e10")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("prepare", "run", "analyse"), required=True)
    parser.add_argument("--network", choices=NETWORKS)
    parser.add_argument("--scan-dir", type=Path, default=SCAN_DIR)
    parser.add_argument("--linx-path", type=Path, default=DEFAULT_LINX_PATH)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_scan_module():
    spec = importlib.util.spec_from_file_location("fr_linx_scan_for_validation", SCAN_PROGRAM)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load {SCAN_PROGRAM}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def latest_successes(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, low_memory=False)
    raw["_row_order"] = np.arange(len(raw))
    frame = raw.loc[raw["status"].eq("ok")].sort_values("_row_order")
    frame = frame.drop_duplicates("point_id", keep="last").drop(columns="_row_order")
    numeric_columns = [
        "eta_fac",
        "clock_factor",
        "tau_n_fac",
        "D_H_1e5",
        "Yp_mass",
        "Li7_H_1e10",
        "q_pull_chi2",
        "chi2_D",
        "chi2_He4",
        "chi2_Li_measurement",
        "chi2_eta_cmb",
        "chi2_tau",
        "chi2_total_li_measurement",
        "chi2_total_li_lower_bound",
    ]
    numeric_columns.extend(column for column in frame if column.startswith("q_"))
    for column in numeric_columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def state_from_row(label: str, row: pd.Series) -> dict[str, Any]:
    q_values = {
        column[2:]: float(row[column])
        for column in row.index
        if column.startswith("q_")
        and column != "q_pull_chi2"
        and math.isfinite(float(row[column]))
        and abs(float(row[column])) > 1e-12
    }
    return {
        "state_label": label,
        "source_point_id": str(row["point_id"]),
        "eta_fac": float(row["eta_fac"]),
        "clock_factor": float(row["clock_factor"]),
        "tau_n_fac": float(row["tau_n_fac"]),
        "q_values": q_values,
        "source_D_H_1e5": float(row["D_H_1e5"]),
        "source_Yp_mass": float(row["Yp_mass"]),
        "source_Li7_H_1e10": float(row["Li7_H_1e10"]),
    }


def select_states(frame: pd.DataFrame) -> list[dict[str, Any]]:
    eta_gate = (frame["eta_fac"] - 1.0).abs().le(2.0 * 0.04 / 6.12)
    tau_gate = (frame["tau_n_fac"] - 1.0).abs().le(2.0 * 0.6 / TAU_REFERENCE_SECONDS)
    su2_base = frame["clock_factor"].ge(1.0) & eta_gate & tau_gate
    d_gate = (frame["D_H_1e5"] - REGISTERED["D"][0]).abs().le(2.0 * REGISTERED["D"][1])
    he_gate = (frame["Yp_mass"] - REGISTERED["He"][0]).abs().le(2.0 * REGISTERED["He"][1])

    base_without_updated = (
        frame["chi2_total_li_measurement"] - frame["chi2_D"] - frame["chi2_He4"] - frame["chi2_tau"]
    )
    tau_mean = CURRENT["tau"][0] / TAU_REFERENCE_SECONDS
    tau_sigma = CURRENT["tau"][1] / TAU_REFERENCE_SECONDS
    combined_score = (
        base_without_updated
        + ((frame["D_H_1e5"] - CURRENT["D"][0]) / CURRENT["D"][1]) ** 2
        + ((frame["Yp_mass"] - CURRENT["He"][0]) / CURRENT["He"][1]) ** 2
        + ((frame["tau_n_fac"] - tau_mean) / tau_sigma) ** 2
    )

    selectors = [
        ("baseline", frame["family"].eq("baseline")),
        ("best_registered_direct", frame["chi2_total_li_measurement"].eq(frame["chi2_total_li_measurement"].min())),
        ("best_registered_lower_bound", frame["chi2_total_li_lower_bound"].eq(frame["chi2_total_li_lower_bound"].min())),
        (
            "best_su2_expansion_only",
            su2_base & frame["q_pull_chi2"].le(1e-12),
        ),
        (
            "best_su2_modest_rates",
            su2_base & frame["q_pull_chi2"].le(9.0),
        ),
        ("best_su2_all_rates", su2_base),
        ("minimum_li_given_registered_D_He", d_gate & he_gate),
        ("best_combined_current_reweight", pd.Series(True, index=frame.index)),
    ]

    states = []
    for label, mask in selectors:
        subset = frame.loc[mask]
        if subset.empty:
            raise RuntimeError(f"No source row for {label}")
        if label == "minimum_li_given_registered_D_He":
            row = subset.loc[subset["Li7_H_1e10"].idxmin()]
        elif label == "best_combined_current_reweight":
            row = frame.loc[combined_score.idxmin()]
        elif label == "best_registered_lower_bound":
            row = subset.loc[subset["chi2_total_li_lower_bound"].idxmin()]
        else:
            row = subset.loc[subset["chi2_total_li_measurement"].idxmin()]
        states.append(state_from_row(label, row))

    baseline = next(state for state in states if state["state_label"] == "baseline")
    best = next(state for state in states if state["state_label"] == "best_registered_direct")
    for label, source in (("baseline_pdg_tau", baseline), ("best_registered_direct_pdg_tau", best)):
        derived = dict(source)
        derived["q_values"] = dict(source["q_values"])
        derived.update(
            {
                "state_label": label,
                "source_point_id": f"derived_from:{source['source_point_id']}",
                "tau_n_fac": PDG_TAU_FACTOR,
                "source_D_H_1e5": math.nan,
                "source_Yp_mass": math.nan,
                "source_Li7_H_1e10": math.nan,
            }
        )
        states.append(derived)
    return states


def matrix_row(
    group: str,
    case_label: str,
    state: dict[str, Any],
    network: str,
    sampling: int,
    q_values: dict[str, float] | None = None,
    reaction: str = "",
    sign: int = 0,
) -> dict[str, Any]:
    q = dict(state["q_values"] if q_values is None else q_values)
    matrix_id = f"{group}|{case_label}|{network}|sampling{sampling}"
    return {
        "matrix_id": matrix_id,
        "case_group": group,
        "case_label": case_label,
        "state_label": state["state_label"],
        "source_point_id": state["source_point_id"],
        "network": network,
        "sampling_nTOp": sampling,
        "rtol": STRICT_RTOL,
        "atol": STRICT_ATOL,
        "eta_fac": state["eta_fac"],
        "clock_factor": state["clock_factor"],
        "tau_n_fac": state["tau_n_fac"],
        "reaction": reaction,
        "perturbation_sign": sign,
        "q_values_json": json.dumps(q, sort_keys=True, separators=(",", ":")),
        "source_D_H_1e5": state["source_D_H_1e5"],
        "source_Yp_mass": state["source_Yp_mass"],
        "source_Li7_H_1e10": state["source_Li7_H_1e10"],
    }


def build_matrix(states: list[dict[str, Any]]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for state in states:
        for network in NETWORKS:
            rows.append(matrix_row("core", state["state_label"], state, network, PRIMARY_SAMPLING))

    response_states = [
        next(state for state in states if state["state_label"] == "baseline"),
        next(state for state in states if state["state_label"] == "best_registered_direct"),
    ]
    for state in response_states:
        for reaction in D_RATES:
            for sign in (-1, 1):
                q = dict(state["q_values"])
                q[reaction] = float(sign)
                label = f"{state['state_label']}__{reaction}{sign:+d}"
                for network in NETWORKS:
                    rows.append(
                        matrix_row(
                            "d_rate_response", label, state, network, PRIMARY_SAMPLING, q, reaction, sign
                        )
                    )

    for state in response_states:
        for reaction in OMITTED_LI_BE_RATES:
            for sign in (-1, 1):
                q = dict(state["q_values"])
                q[reaction] = float(sign)
                label = f"{state['state_label']}__{reaction}{sign:+d}"
                rows.append(
                    matrix_row(
                        "full_omitted_rate_response",
                        label,
                        state,
                        "full",
                        PRIMARY_SAMPLING,
                        q,
                        reaction,
                        sign,
                    )
                )

    for state in response_states:
        for network in NETWORKS:
            rows.append(matrix_row("sampling_check", state["state_label"], state, network, CHECK_SAMPLING))

    matrix = pd.DataFrame(rows)
    if len(matrix) != 92 or matrix["matrix_id"].duplicated().any():
        raise RuntimeError(f"Validation matrix invariant failed: rows={len(matrix)}")
    return matrix


def write_registration(path: Path, matrix: pd.DataFrame) -> None:
    counts = matrix.groupby(["network", "case_group"]).size().unstack(fill_value=0)
    lines = [
        "# FR/LINX BBN Targeted Validation Matrix Registration",
        "",
        "Date: 17 July 2026",
        "Status: fixed before validation abundances were evaluated",
        "",
        "## Objective",
        "",
        "Validate the interpretation-bearing rows of the completed 30,503-point key-network scan against tighter numerical settings, larger networks, released deuterium-rate uncertainties, omitted lithium/beryllium channels and the updated neutron-lifetime centre.",
        "",
        "## Fixed Numerical Settings",
        "",
        r"\[",
        r"\begin{aligned}",
        rf"\mathrm{{rtol}}             &= {STRICT_RTOL:.0e},\\",
        rf"\mathrm{{atol}}             &= {STRICT_ATOL:.0e},\\",
        rf"N_{{\rm interp,primary}}    &= {PRIMARY_SAMPLING},\\",
        rf"N_{{\rm interp,check}}      &= {CHECK_SAMPLING}.",
        r"\end{aligned}",
        r"\]",
        "",
        "The networks are `key_PRIMAT_2023`, `small_PRIMAT_2023` and `full_PRIMAT_2023`. The full network retains the LINX small-to-full switch at its upstream default temperature.",
        "",
        "## Matrix",
        "",
        f"Total evaluations: **{len(matrix)}**.",
        "",
        "| Network | Core | D-rate response | Omitted Li/Be response | Sampling check | Total |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for network in NETWORKS:
        row = counts.loc[network]
        values = [int(row.get(group, 0)) for group in ("core", "d_rate_response", "full_omitted_rate_response", "sampling_check")]
        lines.append(f"| {network} | {values[0]} | {values[1]} | {values[2]} | {values[3]} | {sum(values)} |")
    lines.extend(
        [
            "",
            "The ten core states are the baseline; registered direct-lithium and lower-bound optima; expansion-only, modest-rate and unrestricted SU2-compatible optima; minimum lithium after the registered D+He gate; combined-current-anchor optimum; and baseline/direct-optimum states at the 2025 PDG neutron-lifetime centre.",
            "",
            "The D-rate response uses independent unit pulls for `dpHe3g`, `ddHe3n` and `ddtp` around both the baseline and registered direct-lithium optimum in every network. The full-only completeness response uses unit pulls for `Li7paag`, `Be7naa`, `Be7daap`, `Be7pB8g` and `Li7daan`.",
            "",
            "## Readout Rules",
            "",
            "1. Numerical convergence compares strict key-network core rows with their stored `rtol=1e-5`, `atol=1e-9`, sampling-150 values.",
            "2. Interpolation convergence compares sampling 300 and 600 at the baseline and registered direct-lithium optimum in each network.",
            "3. Network sensitivity compares key, small and full predictions at identical strict settings.",
            "4. Nuclear theory sensitivity is estimated from centred unit-pull responses. It is descriptive and does not assume that the response is globally linear.",
            "5. Registered and combined-current abundance gates are reported, but this designed matrix is not treated as a new parameter search.",
            "6. No FR or SU2 mechanism is inferred from a nuisance-rate response unless an explicit model predicts that coupled pull pattern.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def prepare(cli: argparse.Namespace) -> None:
    cli.outdir.mkdir(parents=True, exist_ok=True)
    points_path = cli.scan_dir / "bbn_lithium_linx_fr_network_points.csv"
    frame = latest_successes(points_path)
    states = select_states(frame)
    matrix = build_matrix(states)
    matrix_path = cli.outdir / "validation_matrix.csv"
    matrix.to_csv(matrix_path, index=False)
    registration_path = cli.outdir / "validation_matrix_registration.md"
    write_registration(registration_path, matrix)
    config = {
        "generated_utc": datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "status": "registered_before_validation_evaluation",
        "programme": str(SCRIPT_PATH),
        "programme_sha256": sha256(SCRIPT_PATH),
        "source_catalogue": str(points_path),
        "source_catalogue_sha256": sha256(points_path),
        "matrix_sha256": sha256(matrix_path),
        "linx_path": str(cli.linx_path),
        "networks": list(NETWORKS),
        "strict_rtol": STRICT_RTOL,
        "strict_atol": STRICT_ATOL,
        "primary_sampling": PRIMARY_SAMPLING,
        "check_sampling": CHECK_SAMPLING,
        "n_evaluations": len(matrix),
    }
    (cli.outdir / "validation_matrix_config.json").write_text(
        json.dumps(config, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Prepared {len(matrix)} evaluations: {matrix_path}")
    print(f"Registration: {registration_path}")


RESULT_FIELDS = [
    "matrix_id",
    "case_group",
    "case_label",
    "state_label",
    "source_point_id",
    "network",
    "sampling_nTOp",
    "rtol",
    "atol",
    "eta_fac",
    "clock_factor",
    "tau_n_fac",
    "reaction",
    "perturbation_sign",
    "q_values_json",
    "q_pull_chi2",
    "source_D_H_1e5",
    "source_Yp_mass",
    "source_Li7_H_1e10",
    "neff_final",
    "D_H_1e5",
    "Yp_mass",
    "Li7_H_1e10",
    "He3_H_1e5",
    "Li6_H_1e14",
    "chi2_registered",
    "chi2_combined_current",
    "registered_D_He_Li_2sigma",
    "combined_current_D_He_Li_2sigma",
    "elapsed_s",
    "status",
    "message",
    "run_utc",
]


def append_result(path: Path, row: dict[str, Any]) -> None:
    exists = path.exists()
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS)
        if not exists:
            writer.writeheader()
        writer.writerow({field: row.get(field, "") for field in RESULT_FIELDS})
        handle.flush()


def score_prediction(pred: dict[str, float], eta: float, tau: float, q_penalty: float, anchors: dict[str, tuple[float, float]]) -> float:
    tau_mean = anchors["tau"][0] / TAU_REFERENCE_SECONDS
    tau_sigma = anchors["tau"][1] / TAU_REFERENCE_SECONDS
    return float(
        ((pred["D_H_1e5"] - anchors["D"][0]) / anchors["D"][1]) ** 2
        + ((pred["Yp_mass"] - anchors["He"][0]) / anchors["He"][1]) ** 2
        + ((pred["Li7_H_1e10"] - anchors["Li"][0]) / anchors["Li"][1]) ** 2
        + ((eta - 1.0) / (0.04 / 6.12)) ** 2
        + ((tau - tau_mean) / tau_sigma) ** 2
        + q_penalty
    )


def joint_gate(pred: dict[str, float], anchors: dict[str, tuple[float, float]]) -> bool:
    return all(
        abs(pred[column] - anchors[key][0]) <= 2.0 * anchors[key][1]
        for column, key in (("D_H_1e5", "D"), ("Yp_mass", "He"), ("Li7_H_1e10", "Li"))
    )


def run_network(cli: argparse.Namespace) -> None:
    if cli.network is None:
        raise ValueError("--network is required in run mode")
    matrix = pd.read_csv(cli.outdir / "validation_matrix.csv")
    subset = matrix.loc[matrix["network"].eq(cli.network)].copy()
    result_path = cli.outdir / f"validation_results_{cli.network}.csv"
    done: set[str] = set()
    if cli.resume and result_path.exists():
        previous = pd.read_csv(result_path, low_memory=False)
        done = set(previous.loc[previous["status"].eq("ok"), "matrix_id"].astype(str))

    scan = load_scan_module()
    runner = scan.LinxRunner(
        cli.linx_path,
        cli.network,
        STRICT_RTOL,
        STRICT_ATOL,
        PRIMARY_SAMPLING,
        PRIMARY_SAMPLING,
        PRIMARY_SAMPLING,
    )
    available = set(runner.reaction_names)
    print(
        f"[{datetime.now().isoformat(timespec='seconds')}] {cli.network}: "
        f"{len(subset)} evaluations, {len(done)} already complete",
        flush=True,
    )
    for sequence, (_, item) in enumerate(subset.iterrows(), start=1):
        matrix_id = str(item["matrix_id"])
        if matrix_id in done:
            continue
        q_values = {key: float(value) for key, value in json.loads(item["q_values_json"]).items()}
        unavailable = sorted(set(q_values) - available)
        start = time.time()
        base = item.to_dict()
        base["run_utc"] = datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z")
        try:
            if unavailable:
                raise ValueError(f"Unavailable reactions in {cli.network}: {unavailable}")
            sampling = int(item["sampling_nTOp"])
            runner.sampling_key = sampling
            runner.sampling_small = sampling
            runner.sampling_full = sampling
            point = scan.ScanPoint(
                point_id=matrix_id,
                family=str(item["case_group"]),
                eta_fac=float(item["eta_fac"]),
                clock_factor=float(item["clock_factor"]),
                tau_n_fac=float(item["tau_n_fac"]),
                q_values=q_values,
            )
            pred = runner.run_point(point)
            q_penalty = sum(value * value for value in q_values.values())
            row = {
                **base,
                **pred,
                "q_pull_chi2": q_penalty,
                "chi2_registered": score_prediction(
                    pred, point.eta_fac, point.tau_n_fac, q_penalty, REGISTERED
                ),
                "chi2_combined_current": score_prediction(
                    pred, point.eta_fac, point.tau_n_fac, q_penalty, CURRENT
                ),
                "registered_D_He_Li_2sigma": joint_gate(pred, REGISTERED),
                "combined_current_D_He_Li_2sigma": joint_gate(pred, CURRENT),
                "status": "ok",
                "message": "",
            }
        except Exception as exc:
            row = {**base, "status": "failed", "message": repr(exc)}
        row["elapsed_s"] = time.time() - start
        append_result(result_path, row)
        print(
            f"[{datetime.now().isoformat(timespec='seconds')}] {cli.network} "
            f"{sequence}/{len(subset)} {item['case_group']} {item['case_label']} "
            f"status={row['status']} elapsed={row['elapsed_s']:.2f}s",
            flush=True,
        )


def latest_results(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, low_memory=False)
    raw["_row_order"] = np.arange(len(raw))
    return raw.sort_values("_row_order").drop_duplicates("matrix_id", keep="last").drop(columns="_row_order")


def paired_differences(results: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    core = results.loc[(results["case_group"].eq("core")) & (results["sampling_nTOp"].eq(PRIMARY_SAMPLING))].copy()
    numerical_rows = []
    for _, row in core.loc[core["network"].eq("key")].iterrows():
        if not math.isfinite(float(row["source_D_H_1e5"])):
            continue
        out = {"state_label": row["state_label"]}
        for column in ABUNDANCE_COLUMNS:
            source = float(row[f"source_{column}"])
            out[f"delta_{column}"] = float(row[column]) - source
        numerical_rows.append(out)

    sampling_rows = []
    for (network, state), group in results.loc[
        results["state_label"].isin(("baseline", "best_registered_direct"))
        & results["case_group"].isin(("core", "sampling_check"))
    ].groupby(["network", "state_label"]):
        by_sampling = group.set_index("sampling_nTOp")
        if PRIMARY_SAMPLING not in by_sampling.index or CHECK_SAMPLING not in by_sampling.index:
            continue
        out = {"network": network, "state_label": state}
        for column in ABUNDANCE_COLUMNS:
            out[f"delta_{column}_600_minus_300"] = float(by_sampling.at[CHECK_SAMPLING, column]) - float(
                by_sampling.at[PRIMARY_SAMPLING, column]
            )
        sampling_rows.append(out)

    network_rows = []
    for state, group in core.groupby("state_label"):
        by_network = group.set_index("network")
        if not set(NETWORKS).issubset(by_network.index):
            continue
        for network in ("small", "full"):
            out = {"state_label": state, "comparison": f"{network}_minus_key"}
            for column in ABUNDANCE_COLUMNS:
                out[f"delta_{column}"] = float(by_network.at[network, column]) - float(by_network.at["key", column])
            network_rows.append(out)

    response_rows = []
    for group_name in ("d_rate_response", "full_omitted_rate_response"):
        subset = results.loc[results["case_group"].eq(group_name)]
        for (network, state, reaction), group in subset.groupby(["network", "state_label", "reaction"]):
            by_sign = group.set_index("perturbation_sign")
            if -1 not in by_sign.index or 1 not in by_sign.index:
                continue
            out = {
                "case_group": group_name,
                "network": network,
                "state_label": state,
                "reaction": reaction,
            }
            for column in ABUNDANCE_COLUMNS:
                out[f"centred_response_{column}"] = 0.5 * (
                    float(by_sign.at[1, column]) - float(by_sign.at[-1, column])
                )
            response_rows.append(out)
    return tuple(pd.DataFrame(rows) for rows in (numerical_rows, sampling_rows, network_rows, response_rows))


def max_abs(frame: pd.DataFrame, column: str) -> float:
    return float(frame[column].abs().max()) if not frame.empty and column in frame else math.nan


def analyse(cli: argparse.Namespace) -> None:
    matrix = pd.read_csv(cli.outdir / "validation_matrix.csv")
    frames = []
    for network in NETWORKS:
        path = cli.outdir / f"validation_results_{network}.csv"
        if path.exists():
            frames.append(latest_results(path))
    if not frames:
        raise FileNotFoundError("No validation result files found")
    results = pd.concat(frames, ignore_index=True)
    for column in [*ABUNDANCE_COLUMNS, "sampling_nTOp", "elapsed_s", "source_D_H_1e5", "source_Yp_mass", "source_Li7_H_1e10"]:
        results[column] = pd.to_numeric(results[column], errors="coerce")
    combined_path = cli.outdir / "validation_results_combined.csv"
    results.to_csv(combined_path, index=False)
    numerical, sampling, network, responses = paired_differences(results.loc[results["status"].eq("ok")])
    numerical.to_csv(cli.outdir / "validation_numerical_convergence.csv", index=False)
    sampling.to_csv(cli.outdir / "validation_sampling_convergence.csv", index=False)
    network.to_csv(cli.outdir / "validation_network_differences.csv", index=False)
    responses.to_csv(cli.outdir / "validation_rate_responses.csv", index=False)

    ok = results.loc[results["status"].eq("ok")]
    missing = set(matrix["matrix_id"]) - set(ok["matrix_id"])
    baseline = ok.loc[
        ok["case_group"].eq("core")
        & ok["state_label"].eq("baseline")
        & ok["sampling_nTOp"].eq(PRIMARY_SAMPLING)
    ].sort_values("network")
    d_response = responses.loc[responses["case_group"].eq("d_rate_response")]
    d_budget_rows = []
    for (network_name, state), group in d_response.groupby(["network", "state_label"]):
        d_budget_rows.append(
            {
                "network": network_name,
                "state_label": state,
                **{
                    f"quadrature_sigma_{column}": float(
                        np.sqrt(np.sum(np.square(group[f"centred_response_{column}"])))
                    )
                    for column in ABUNDANCE_COLUMNS
                },
            }
        )
    d_budget = pd.DataFrame(d_budget_rows)
    d_budget.to_csv(cli.outdir / "validation_d_rate_theory_budget.csv", index=False)

    lines = [
        "# FR/LINX BBN Targeted Validation Matrix Report",
        "",
        f"Generated: {datetime.now(UTC).isoformat(timespec='seconds').replace('+00:00', 'Z')}",
        "",
        "## Completion",
        "",
        f"- Registered evaluations: `{len(matrix)}`.",
        f"- Successful evaluations: `{len(ok)}`.",
        f"- Missing or unresolved evaluations: `{len(missing)}`.",
        "",
        "## Strict Baseline",
        "",
        "| Network | D/H x 1e5 | Yp | Li-7/H x 1e10 | Runtime (s) |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in baseline.itertuples(index=False):
        lines.append(
            f"| {row.network} | {row.D_H_1e5:.7g} | {row.Yp_mass:.7g} | "
            f"{row.Li7_H_1e10:.7g} | {row.elapsed_s:.3g} |"
        )
    lines.extend(
        [
            "",
            "## Maximum Absolute Differences",
            "",
            "| Check | D/H x 1e5 | Yp | Li-7/H x 1e10 |",
            "|---|---:|---:|---:|",
            f"| Strict key minus stored loose key | {max_abs(numerical, 'delta_D_H_1e5'):.4g} | {max_abs(numerical, 'delta_Yp_mass'):.4g} | {max_abs(numerical, 'delta_Li7_H_1e10'):.4g} |",
            f"| Sampling 600 minus 300 | {max_abs(sampling, 'delta_D_H_1e5_600_minus_300'):.4g} | {max_abs(sampling, 'delta_Yp_mass_600_minus_300'):.4g} | {max_abs(sampling, 'delta_Li7_H_1e10_600_minus_300'):.4g} |",
            f"| Small minus key | {max_abs(network.loc[network['comparison'].eq('small_minus_key')], 'delta_D_H_1e5'):.4g} | {max_abs(network.loc[network['comparison'].eq('small_minus_key')], 'delta_Yp_mass'):.4g} | {max_abs(network.loc[network['comparison'].eq('small_minus_key')], 'delta_Li7_H_1e10'):.4g} |",
            f"| Full minus key | {max_abs(network.loc[network['comparison'].eq('full_minus_key')], 'delta_D_H_1e5'):.4g} | {max_abs(network.loc[network['comparison'].eq('full_minus_key')], 'delta_Yp_mass'):.4g} | {max_abs(network.loc[network['comparison'].eq('full_minus_key')], 'delta_Li7_H_1e10'):.4g} |",
            "",
            "## D-Rate Theory Response",
            "",
            "Independent centred unit-pull responses were combined in quadrature. These are local response estimates, not a complete marginal likelihood.",
            "",
            "| Network | State | sigma(D/H x 1e5) | sigma(Yp) | sigma(Li x 1e10) |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for row in d_budget.itertuples(index=False):
        lines.append(
            f"| {row.network} | {row.state_label} | {row.quadrature_sigma_D_H_1e5:.4g} | "
            f"{row.quadrature_sigma_Yp_mass:.4g} | {row.quadrature_sigma_Li7_H_1e10:.4g} |"
        )
    omitted = responses.loc[responses["case_group"].eq("full_omitted_rate_response")]
    lines.extend(
        [
            "",
            "## Gate and Claim Boundary",
            "",
            f"- Registered D+He+Li two-sigma passes among successful designed evaluations: `{int(ok['registered_D_He_Li_2sigma'].astype(str).str.lower().eq('true').sum())}`.",
            f"- Combined-current D+He+Li two-sigma passes: `{int(ok['combined_current_D_He_Li_2sigma'].astype(str).str.lower().eq('true').sum())}`.",
            f"- Largest absolute full-network omitted-channel Li response: `{max_abs(omitted, 'centred_response_Li7_H_1e10'):.4g}` in Li/H x 1e10 units.",
            "",
            "This matrix validates numerical and selected-network stability and quantifies local nuclear-rate sensitivity. It remains a standard-background LINX calculation. It does not validate a non-expanding FR thermal history or identify nuclear pulls with FR or SU2 dynamics.",
        ]
    )
    report_path = cli.outdir / "validation_matrix_report.md"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    audit = {
        "generated_utc": datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "expected": len(matrix),
        "successful": len(ok),
        "missing_matrix_ids": sorted(missing),
        "matrix_sha256": sha256(cli.outdir / "validation_matrix.csv"),
        "combined_results_sha256": sha256(combined_path),
        "programme_sha256": sha256(SCRIPT_PATH),
    }
    (cli.outdir / "validation_matrix_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Saved report: {report_path}")
    print(f"Successful: {len(ok)}/{len(matrix)}; missing: {len(missing)}")


def main() -> None:
    cli = parse_args()
    cli.outdir.mkdir(parents=True, exist_ok=True)
    if cli.mode == "prepare":
        prepare(cli)
    elif cli.mode == "run":
        run_network(cli)
    else:
        analyse(cli)


if __name__ == "__main__":
    main()
