r"""Follow-up gates for FR c(z)=c0(1+z) alpha=0 versus LCDM.

This wrapper runs the already-defined alpha=0 FR no-loss versus LCDM ladder
programme across:

* leave-one-dataset-out and two-dataset combinations;
* redshift bands; and
* calibration-pull summaries.

It deliberately reuses
`analyze_fr_c1pz_noloss_lcdm_calibration_ladder_2026-07-16.py` so the
likelihood, covariance handling, frame choices, calibration priors, and LCDM
comparison are identical to the preceding promotion-gate run.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
LADDER_SCRIPT = SCRIPT_PATH.with_name(
    "analyze_fr_c1pz_noloss_lcdm_calibration_ladder_2026-07-16.py"
)
DEFAULT_OUTDIR = (
    REPO_ROOT
    / "plamb_runs"
    / "diagnostics"
    / "fr_c1pz_noloss_followup_gates_20260716"
)

PRIOR_CONFIGS = (
    "budget_010mmag_ds020:0.020:0.010,"
    "budget_016mmag_ds032:0.032:0.016,"
    "budget_025mmag_ds050:0.050:0.025,"
    "budget_050mmag_ds100_stress:0.100:0.050"
)


@dataclass(frozen=True)
class RunCase:
    gate: str
    label: str
    datasets: str
    frames: str
    z_min: float
    max_z: float | None
    min_n: int


DATASET_CASES = [
    RunCase("dataset_combo", "pantheon_only", "pantheon", "HD,HEL,CMB_PANTHEON_ONLY", 0.01, None, 10),
    RunCase("dataset_combo", "des_only", "des_raw", "HD,HEL", 0.01, None, 10),
    RunCase("dataset_combo", "union31_only", "union31", "HD", 0.01, None, 10),
    RunCase("dataset_combo", "pantheon_des", "pantheon,des_raw", "HD,HEL,CMB_PANTHEON_ONLY", 0.01, None, 10),
    RunCase("dataset_combo", "pantheon_union31", "pantheon,union31", "HD,HEL,CMB_PANTHEON_ONLY", 0.01, None, 10),
    RunCase("dataset_combo", "des_union31", "des_raw,union31", "HD,HEL", 0.01, None, 10),
    RunCase("dataset_combo", "all_three", "pantheon,des_raw,union31", "HD,HEL,CMB_PANTHEON_ONLY", 0.01, None, 10),
]

REDSHIFT_CASES = [
    RunCase("redshift_band", "all_z", "pantheon,des_raw,union31", "HD,HEL,CMB_PANTHEON_ONLY", 0.01, None, 10),
    RunCase("redshift_band", "low_z_0p01_0p10", "pantheon,des_raw,union31", "HD,HEL,CMB_PANTHEON_ONLY", 0.01, 0.10, 1),
    RunCase("redshift_band", "mid_z_0p10_0p50", "pantheon,des_raw,union31", "HD,HEL,CMB_PANTHEON_ONLY", 0.10, 0.50, 1),
    RunCase("redshift_band", "high_z_0p50_plus", "pantheon,des_raw,union31", "HD,HEL,CMB_PANTHEON_ONLY", 0.50, None, 1),
]


def safe_label(value: str) -> str:
    return "".join(c if c.isalnum() or c in {"_", "-", "."} else "_" for c in value)


def run_ladder(case: RunCase, outdir: Path, force: bool) -> tuple[Path, str]:
    case_dir = outdir / "runs" / safe_label(f"{case.gate}_{case.label}")
    comparison_path = case_dir / "fr_c1pz_noloss_lcdm_ladder_comparison.csv"
    if comparison_path.exists() and not force:
        return case_dir, "reused"

    case_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(LADDER_SCRIPT),
        "--outdir",
        str(case_dir),
        "--datasets",
        case.datasets,
        "--frames",
        case.frames,
        "--offset-modes",
        "none,dataset+idsurvey",
        "--prior-configs",
        PRIOR_CONFIGS,
        "--z-min",
        str(case.z_min),
        "--min-n",
        str(case.min_n),
    ]
    if case.max_z is not None:
        cmd.extend(["--max-z", str(case.max_z)])

    completed = subprocess.run(
        cmd,
        cwd=str(REPO_ROOT),
        text=True,
        capture_output=True,
        check=False,
    )
    (case_dir / "stdout.log").write_text(completed.stdout, encoding="utf-8")
    (case_dir / "stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"Ladder case {case.gate}/{case.label} failed with code {completed.returncode}. "
            f"See {case_dir / 'stderr.log'}"
        )
    return case_dir, "ran"


def read_case_outputs(case: RunCase, case_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    comp = pd.read_csv(case_dir / "fr_c1pz_noloss_lcdm_ladder_comparison.csv")
    summary = pd.read_csv(case_dir / "fr_c1pz_noloss_lcdm_ladder_summary.csv")
    offsets = pd.read_csv(case_dir / "fr_c1pz_noloss_lcdm_ladder_offsets.csv")
    for df in (comp, summary, offsets):
        df.insert(0, "case_label", case.label)
        df.insert(0, "gate", case.gate)
        df["case_datasets"] = case.datasets
        df["case_z_min"] = case.z_min
        df["case_max_z"] = np.nan if case.max_z is None else case.max_z
    return comp, summary, offsets


def aggregate_pull_summary(offsets: pd.DataFrame) -> pd.DataFrame:
    keep = offsets[
        (offsets["offset_mode"] == "dataset+idsurvey")
        & (offsets["model"] == "FR_C1PZ_ALPHA0_H0fixed675")
        & np.isfinite(offsets["pull_sigma"])
    ].copy()
    if keep.empty:
        return pd.DataFrame()
    keep["abs_pull_sigma"] = keep["pull_sigma"].abs()
    keep["offset_dataset"] = keep["offset_label"].astype(str).str.split(":").str[0]
    rows = []
    group_cols = ["gate", "case_label", "frame", "prior_config", "offset_dataset", "offset_type"]
    for key, grp in keep.groupby(group_cols):
        rows.append(
            {
                "gate": key[0],
                "case_label": key[1],
                "frame": key[2],
                "prior_config": key[3],
                "offset_dataset": key[4],
                "offset_type": key[5],
                "n_offsets": int(len(grp)),
                "max_abs_pull_sigma": float(grp["abs_pull_sigma"].max()),
                "n_abs_pull_ge_3": int((grp["abs_pull_sigma"] >= 3.0).sum()),
                "n_abs_pull_ge_5": int((grp["abs_pull_sigma"] >= 5.0).sum()),
                "median_abs_pull_sigma": float(grp["abs_pull_sigma"].median()),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["gate", "case_label", "frame", "prior_config", "max_abs_pull_sigma"],
        ascending=[True, True, True, True, False],
    )


def best_hd_dataset_rows(comp: pd.DataFrame) -> pd.DataFrame:
    keep = comp[
        (comp["gate"] == "dataset_combo")
        & (comp["frame"] == "HD")
        & (comp["offset_mode"] == "dataset+idsurvey")
    ].copy()
    if keep.empty:
        return keep
    keep["fr_beats_lcdm_bic"] = keep["delta_BIC_FR_minus_LCDMfree"] < 0
    keep["fr_beats_lcdm_posterior"] = keep["delta_posterior_FR_minus_LCDMfree"] < 0
    return keep.sort_values(["case_label", "prior_config"])


def redshift_hd_rows(comp: pd.DataFrame) -> pd.DataFrame:
    keep = comp[
        (comp["gate"] == "redshift_band")
        & (comp["frame"] == "HD")
        & (comp["offset_mode"] == "dataset+idsurvey")
    ].copy()
    if keep.empty:
        return keep
    keep["fr_beats_lcdm_bic"] = keep["delta_BIC_FR_minus_LCDMfree"] < 0
    keep["fr_beats_lcdm_posterior"] = keep["delta_posterior_FR_minus_LCDMfree"] < 0
    return keep.sort_values(["case_label", "prior_config"])


def markdown_table(df: pd.DataFrame, cols: list[str], max_rows: int | None = None) -> str:
    if df.empty:
        return "_No rows._"
    out = df[cols].copy()
    if max_rows is not None:
        out = out.head(max_rows)
    return out.to_markdown(index=False, floatfmt=".6g")


def write_report(
    path: Path,
    comp: pd.DataFrame,
    summary: pd.DataFrame,
    offsets: pd.DataFrame,
    pull_summary: pd.DataFrame,
    run_rows: list[dict[str, object]],
) -> None:
    lines: list[str] = []
    lines.append("# FR c(z)=c0(1+z) No-Loss Follow-Up Gates")
    lines.append("")
    lines.append(f"Date: {datetime.now():%Y-%m-%d %H:%M:%S}")
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append("This follow-up reuses the identical-offset FR alpha=0 versus LCDM ladder programme and asks whether the apparent FR advantage is broad across catalogues, redshift bands, and calibration-pull structure.")
    lines.append("")
    lines.append("Negative `delta_BIC_FR_minus_LCDMfree` favours FR. Negative `delta_posterior_FR_minus_LCDMfree` favours FR before the base-parameter BIC penalty.")
    lines.append("")

    run_df = pd.DataFrame(run_rows)
    lines.append("## Executed Cases")
    lines.append("")
    lines.append(markdown_table(run_df, ["gate", "case_label", "datasets", "frames", "z_min", "max_z", "status"]))
    lines.append("")

    lines.append("## Leave-One-Dataset-Out: HD Dataset+Survey")
    lines.append("")
    ds = best_hd_dataset_rows(comp)
    lines.append(markdown_table(
        ds,
        [
            "case_label",
            "prior_config",
            "N",
            "FR_BIC_base",
            "LCDMfree_Om",
            "LCDMfree_BIC_base",
            "delta_BIC_FR_minus_LCDMfree",
            "delta_posterior_FR_minus_LCDMfree",
            "fr_beats_lcdm_bic",
            "fr_beats_lcdm_posterior",
        ],
    ))
    lines.append("")

    if not ds.empty:
        tally = ds.groupby("case_label").agg(
            bic_wins=("fr_beats_lcdm_bic", "sum"),
            posterior_wins=("fr_beats_lcdm_posterior", "sum"),
            n_rows=("fr_beats_lcdm_bic", "size"),
        ).reset_index()
        lines.append("### Dataset-Combination Tally")
        lines.append("")
        lines.append(markdown_table(tally, ["case_label", "bic_wins", "posterior_wins", "n_rows"]))
        lines.append("")

    lines.append("## Redshift-Band Check: HD Dataset+Survey")
    lines.append("")
    rz = redshift_hd_rows(comp)
    lines.append(markdown_table(
        rz,
        [
            "case_label",
            "prior_config",
            "N",
            "FR_BIC_base",
            "LCDMfree_Om",
            "LCDMfree_BIC_base",
            "delta_BIC_FR_minus_LCDMfree",
            "delta_posterior_FR_minus_LCDMfree",
            "fr_beats_lcdm_bic",
            "fr_beats_lcdm_posterior",
        ],
    ))
    lines.append("")

    if not rz.empty:
        tally = rz.groupby("case_label").agg(
            bic_wins=("fr_beats_lcdm_bic", "sum"),
            posterior_wins=("fr_beats_lcdm_posterior", "sum"),
            n_rows=("fr_beats_lcdm_bic", "size"),
        ).reset_index()
        lines.append("### Redshift-Band Tally")
        lines.append("")
        lines.append(markdown_table(tally, ["case_label", "bic_wins", "posterior_wins", "n_rows"]))
        lines.append("")

    lines.append("## Calibration-Pull Audit")
    lines.append("")
    if pull_summary.empty:
        lines.append("_No finite pull-summary rows._")
    else:
        hd_pulls = pull_summary[
            (pull_summary["gate"].isin(["dataset_combo", "redshift_band"]))
            & (pull_summary["frame"] == "HD")
        ].sort_values("max_abs_pull_sigma", ascending=False)
        lines.append(markdown_table(
            hd_pulls,
            [
                "gate",
                "case_label",
                "prior_config",
                "offset_dataset",
                "offset_type",
                "n_offsets",
                "max_abs_pull_sigma",
                "n_abs_pull_ge_3",
                "n_abs_pull_ge_5",
                "median_abs_pull_sigma",
            ],
            max_rows=40,
        ))
    lines.append("")

    lines.append("## Readout Boundary")
    lines.append("")
    lines.append("This is still a promotion-gate diagnostic, not a discovery claim. A robust FR result should not be carried by a single catalogue, a single redshift band, or implausibly large calibration pulls. The output tables are intended to identify exactly which gate is currently strongest or weakest.")
    lines.append("")
    lines.append("## Output Files")
    lines.append("")
    lines.append("- `fr_c1pz_followup_dataset_comparison.csv`")
    lines.append("- `fr_c1pz_followup_redshift_comparison.csv`")
    lines.append("- `fr_c1pz_followup_all_summary.csv`")
    lines.append("- `fr_c1pz_followup_all_offsets.csv`")
    lines.append("- `fr_c1pz_followup_pull_summary.csv`")
    lines.append("- `fr_c1pz_followup_run_manifest.csv`")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_readout(path: Path, comp: pd.DataFrame, pull_summary: pd.DataFrame) -> None:
    ds = best_hd_dataset_rows(comp)
    rz = redshift_hd_rows(comp)
    lines: list[str] = []
    lines.append("# FR No-Loss Follow-Up Gates Readout")
    lines.append("")
    lines.append("Date: 2026-07-16")
    lines.append("")
    lines.append("## Executive readout")
    lines.append("")
    if not ds.empty:
        ds_tally = ds.groupby("case_label").agg(
            bic_wins=("fr_beats_lcdm_bic", "sum"),
            posterior_wins=("fr_beats_lcdm_posterior", "sum"),
            n_rows=("fr_beats_lcdm_bic", "size"),
        ).reset_index()
        lines.append("Leave-one-dataset-out gate, HD dataset+survey rows:")
        lines.append("")
        lines.append(markdown_table(ds_tally, ["case_label", "bic_wins", "posterior_wins", "n_rows"]))
        lines.append("")
    if not rz.empty:
        rz_tally = rz.groupby("case_label").agg(
            bic_wins=("fr_beats_lcdm_bic", "sum"),
            posterior_wins=("fr_beats_lcdm_posterior", "sum"),
            n_rows=("fr_beats_lcdm_bic", "size"),
        ).reset_index()
        lines.append("Redshift-band gate, HD dataset+survey rows:")
        lines.append("")
        lines.append(markdown_table(rz_tally, ["case_label", "bic_wins", "posterior_wins", "n_rows"]))
        lines.append("")

    lines.append("## Largest calibration-pull structures")
    lines.append("")
    if pull_summary.empty:
        lines.append("_No finite pull-summary rows._")
    else:
        rows = pull_summary[pull_summary["frame"] == "HD"].sort_values("max_abs_pull_sigma", ascending=False)
        lines.append(markdown_table(
            rows,
            [
                "gate",
                "case_label",
                "prior_config",
                "offset_dataset",
                "offset_type",
                "max_abs_pull_sigma",
                "n_abs_pull_ge_3",
                "n_abs_pull_ge_5",
            ],
            max_rows=15,
        ))
    lines.append("")
    lines.append("## Claim boundary")
    lines.append("")
    lines.append("The right reading is gate-based: a favourable FR-vs-LCDM delta is only persuasive if it persists across catalogue combinations and redshift bands without relying on implausibly large calibration pulls.")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--force", action="store_true")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cases = DATASET_CASES + REDSHIFT_CASES
    comp_frames: list[pd.DataFrame] = []
    summary_frames: list[pd.DataFrame] = []
    offset_frames: list[pd.DataFrame] = []
    run_rows: list[dict[str, object]] = []

    for case in cases:
        print(f"[case] {case.gate}/{case.label}", flush=True)
        case_dir, status = run_ladder(case, args.outdir, args.force)
        comp, summary, offsets = read_case_outputs(case, case_dir)
        comp_frames.append(comp)
        summary_frames.append(summary)
        offset_frames.append(offsets)
        run_rows.append(
            {
                "gate": case.gate,
                "case_label": case.label,
                "datasets": case.datasets,
                "frames": case.frames,
                "z_min": case.z_min,
                "max_z": "" if case.max_z is None else case.max_z,
                "min_n": case.min_n,
                "run_dir": str(case_dir),
                "status": status,
            }
        )

    comp_all = pd.concat(comp_frames, ignore_index=True)
    summary_all = pd.concat(summary_frames, ignore_index=True)
    offsets_all = pd.concat(offset_frames, ignore_index=True)
    pull_summary = aggregate_pull_summary(offsets_all)
    run_manifest = pd.DataFrame(run_rows)

    ds_comp = comp_all[comp_all["gate"] == "dataset_combo"].copy()
    rz_comp = comp_all[comp_all["gate"] == "redshift_band"].copy()

    ds_comp.to_csv(args.outdir / "fr_c1pz_followup_dataset_comparison.csv", index=False)
    rz_comp.to_csv(args.outdir / "fr_c1pz_followup_redshift_comparison.csv", index=False)
    comp_all.to_csv(args.outdir / "fr_c1pz_followup_all_comparison.csv", index=False)
    summary_all.to_csv(args.outdir / "fr_c1pz_followup_all_summary.csv", index=False)
    offsets_all.to_csv(args.outdir / "fr_c1pz_followup_all_offsets.csv", index=False)
    pull_summary.to_csv(args.outdir / "fr_c1pz_followup_pull_summary.csv", index=False)
    run_manifest.to_csv(args.outdir / "fr_c1pz_followup_run_manifest.csv", index=False)

    write_report(
        args.outdir / "fr_c1pz_followup_report.md",
        comp_all,
        summary_all,
        offsets_all,
        pull_summary,
        run_rows,
    )
    write_readout(args.outdir / "fr_c1pz_followup_readout.md", comp_all, pull_summary)
    (args.outdir / "fr_c1pz_followup_config.json").write_text(
        json.dumps(
            {
                "script": str(SCRIPT_PATH),
                "ladder_script": str(LADDER_SCRIPT),
                "outdir": str(args.outdir),
                "prior_configs": PRIOR_CONFIGS,
                "n_cases": len(cases),
                "claim_boundary": "promotion-gate diagnostic only",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Saved readout: {args.outdir / 'fr_c1pz_followup_readout.md'}")
    print(f"Saved report: {args.outdir / 'fr_c1pz_followup_report.md'}")


if __name__ == "__main__":
    main()
