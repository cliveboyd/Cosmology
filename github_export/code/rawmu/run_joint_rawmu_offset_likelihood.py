r"""
Joint DES/Pantheon+/Union3.1 raw-MU nuisance likelihood.

This runner is deliberately narrow. It tests the paper recommendation:

* DES-Dovekie raw + Pantheon+ + Union3.1 in one likelihood;
* full covariance / precision for each input, combined as independent
  block-diagonal covariance blocks;
* shared and survey-specific magnitude offsets;
* per-dataset redshift-column mapping, because Union3.1 compressed nodes use
  `z` while DES/Pantheon use `zHD`, `zHEL`, and `zCMB`.

Outputs:

    plamb_runs/diagnostics/joint_rawmu_offset_likelihood_YYYYMMDD/
        joint_rawmu_offset_summary.csv
        joint_rawmu_offset_offsets.csv
        joint_rawmu_offset_blocks.csv
        joint_rawmu_offset_report.md
        joint_rawmu_offset_config.json
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from fit_plamb_rawmu_nuisance import (
    PRESETS,
    ModelChoice,
    combine_blocks,
    design_matrix,
    fit_model,
    load_dataset,
    offset_rows,
    row_for_fit,
    select_block,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTDIR = ROOT / "plamb_runs" / "diagnostics" / f"joint_rawmu_offset_likelihood_{datetime.now():%Y%m%d}"


FRAME_MAPS = {
    "HD": {
        "PantheonPlusSH0ES": "zHD",
        "DES_Dovekie_raw": "zHD",
        "Union3p1_UNITY1p8": "z",
    },
    "HEL": {
        "PantheonPlusSH0ES": "zHEL",
        "DES_Dovekie_raw": "zHEL",
        "Union3p1_UNITY1p8": "z",
    },
    "CMB_PANTHEON_ONLY": {
        "PantheonPlusSH0ES": "zCMB",
        "DES_Dovekie_raw": "zHD",
        "Union3p1_UNITY1p8": "z",
    },
}


MODELS = [
    ModelChoice("FR_BETA05_H0fixed675", 67.5, 0.5),
    ModelChoice("FR_BETAfree_H0fixed675", 67.5, None),
    ModelChoice("FR_BETA05_H0free", None, 0.5),
    ModelChoice("FR_BETAfree_H0free", None, None),
]


def parse_list(raw: str) -> list[str]:
    return [item.strip() for item in str(raw).split(",") if item.strip()]


def analysis_args(cli: argparse.Namespace) -> SimpleNamespace:
    return SimpleNamespace(
        mu_col=None,
        calibrator_col=None,
        z_min=cli.z_min,
        max_z=cli.max_z,
        min_n=cli.min_n,
        keep_calibrators=cli.keep_calibrators,
        min_survey_n=cli.min_survey_n,
        allow_precision_submatrix=cli.allow_precision_submatrix,
        h0_min=cli.h0_min,
        h0_max=cli.h0_max,
        beta_min=cli.beta_min,
        beta_max=cli.beta_max,
        grid_steps=cli.grid_steps,
    )


def build_blocks(frame: str, cli: argparse.Namespace):
    specs = [PRESETS["pantheon"], PRESETS["des_raw"], PRESETS["union31"]]
    z_map = FRAME_MAPS[frame]
    args = analysis_args(cli)
    blocks = []
    block_rows = []
    for spec in specs:
        loaded = load_dataset(spec, args)
        z_col = z_map[loaded.label]
        block = select_block(loaded, z_col, args)
        blocks.append(block)
        block_rows.append(
            {
                "frame": frame,
                "dataset": block.label,
                "z_col": block.z_col,
                "N_source": block.n_source,
                "N_used": block.n_used,
                "cov_note": block.cov_note,
                "subset_note": block.subset_note,
                "source_path": str(block.source_path),
                "cov_path": str(block.cov_path) if block.cov_path else "",
            }
        )
    return blocks, block_rows


def run_one(frame: str, offset_mode: str, model: ModelChoice, cli: argparse.Namespace):
    blocks, block_rows = build_blocks(frame, cli)
    z, mu, precision = combine_blocks(blocks)
    x, labels, design_note = design_matrix(blocks, offset_mode)
    run_label = f"joint_des_pantheon_union31_{frame}_{offset_mode}_{model.name}"
    fit, offsets, _profiled = fit_model(z, mu, precision, x, model, analysis_args(cli))
    row = row_for_fit(run_label, blocks, frame, offset_mode, x, labels, design_note, model, fit, analysis_args(cli))
    row["frame"] = frame
    row["frame_z_map"] = json.dumps(FRAME_MAPS[frame], sort_keys=True)
    off = offset_rows(run_label, blocks, frame, offset_mode, model.name, labels, offsets)
    for item in off:
        item["frame"] = frame
    return row, off, block_rows


def write_report(path: Path, summary: pd.DataFrame, offsets: pd.DataFrame, blocks: pd.DataFrame, cli: argparse.Namespace) -> None:
    lines: list[str] = []
    lines.append("# Joint Raw-MU Offset Likelihood")
    lines.append("")
    lines.append("Datasets: DES-Dovekie raw, Pantheon+SH0ES, and Union3.1 UNITY1.8 compressed MU nodes.")
    lines.append("")
    lines.append("Covariance treatment: each dataset uses its available full covariance or inverse covariance. Dataset blocks are combined as an independent block-diagonal precision matrix; no cross-dataset covariance is available.")
    lines.append("")
    lines.append("Core model: `MU = 5 log10((c/H0) z (1 + beta z)) + 25 + offsets`.")
    lines.append("")
    lines.append("## Configuration")
    lines.append("")
    lines.append(f"- frames: `{cli.frames}`")
    lines.append(f"- offset modes: `{cli.offset_modes}`")
    lines.append(f"- z_min: `{cli.z_min}`")
    lines.append(f"- H0/beta bounds: H0 `{cli.h0_min}` to `{cli.h0_max}`, beta `{cli.beta_min}` to `{cli.beta_max}`")
    lines.append(f"- min IDSURVEY/group size before OTHER pooling: `{cli.min_survey_n}`")
    lines.append("")
    lines.append("Offset-mode meaning:")
    lines.append("")
    lines.append("- `global`: one shared magnitude offset across all datasets.")
    lines.append("- `dataset`: one magnitude offset per dataset.")
    lines.append("- `dataset+idsurvey`: IDSURVEY offsets where available, with one dataset-level offset for compressed Union3.1.")
    lines.append("")
    lines.append("## Ranked Fits")
    lines.append("")
    if summary.empty:
        lines.append("_No successful rows._")
    else:
        show = summary.sort_values(["frame", "BIC", "AIC"]).copy()
        cols = [
            "frame",
            "offset_mode",
            "model",
            "N",
            "k_eff",
            "H0",
            "beta",
            "chi2",
            "chi2_dof",
            "AIC",
            "BIC",
            "note",
        ]
        lines.append(show[[c for c in cols if c in show.columns]].to_markdown(index=False, floatfmt=".6g"))
    lines.append("")
    lines.append("## Dataset Blocks")
    lines.append("")
    if blocks.empty:
        lines.append("_No block rows._")
    else:
        lines.append(blocks.to_markdown(index=False))
    lines.append("")
    lines.append("## Offset Rows")
    lines.append("")
    if offsets.empty:
        lines.append("_No profiled offsets._")
    else:
        show_offsets = offsets.copy()
        if len(show_offsets) > 80:
            show_offsets = show_offsets.head(80)
            lines.append("Showing first 80 offset rows; see CSV for all rows.")
            lines.append("")
        lines.append(show_offsets.to_markdown(index=False, floatfmt=".6g"))
    lines.append("")
    lines.append("## Interpretation Notes")
    lines.append("")
    lines.append("- Stable beta near 0.5 after `global` and `dataset+idsurvey` offsets would support the raw-MU shape.")
    lines.append("- Large improvement only in high-offset modes should be treated as calibration/survey structure until penalized by BIC or a hierarchical model.")
    lines.append("- H0-free rows with magnitude offsets are scale-degenerate; beta and information criteria are more meaningful than H0 there.")
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--frames", default="HD,HEL")
    parser.add_argument("--offset-modes", default="none,global,dataset,dataset+idsurvey")
    parser.add_argument("--z-min", type=float, default=0.01)
    parser.add_argument("--max-z", type=float, default=None)
    parser.add_argument("--min-n", type=int, default=10)
    parser.add_argument("--min-survey-n", type=int, default=10)
    parser.add_argument("--keep-calibrators", action="store_true")
    parser.add_argument("--allow-precision-submatrix", action="store_true")
    parser.add_argument("--include-degenerate", action="store_true", help="Keep H0-free models when magnitude offsets are profiled.")
    parser.add_argument("--h0-min", type=float, default=45.0)
    parser.add_argument("--h0-max", type=float, default=90.0)
    parser.add_argument("--beta-min", type=float, default=-0.25)
    parser.add_argument("--beta-max", type=float, default=1.25)
    parser.add_argument("--grid-steps", type=int, default=301)
    return parser


def main() -> None:
    cli = build_parser().parse_args()
    cli.outdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    offsets: list[dict[str, object]] = []
    block_rows: list[dict[str, object]] = []
    errors: list[dict[str, object]] = []

    frames = parse_list(cli.frames)
    offset_modes = parse_list(cli.offset_modes)
    for frame in frames:
        if frame not in FRAME_MAPS:
            errors.append({"frame": frame, "offset_mode": "", "model": "", "error": f"Unknown frame {frame}"})
            continue
        for offset_mode in offset_modes:
            for model in MODELS:
                if offset_mode != "none" and model.h0_fixed is None and not cli.include_degenerate:
                    continue
                try:
                    row, off, blocks = run_one(frame, offset_mode, model, cli)
                    rows.append(row)
                    offsets.extend(off)
                    block_rows.extend(blocks)
                    print(
                        f"{frame} {offset_mode} {model.name}: "
                        f"chi2={float(row['chi2']):.3f} beta={float(row['beta']):.5f} BIC={float(row['BIC']):.3f}",
                        flush=True,
                    )
                except Exception as exc:
                    errors.append({"frame": frame, "offset_mode": offset_mode, "model": model.name, "error": str(exc)})
                    print(f"ERROR {frame} {offset_mode} {model.name}: {exc}", flush=True)

    summary = pd.DataFrame(rows)
    offsets_df = pd.DataFrame(offsets)
    blocks_df = pd.DataFrame(block_rows).drop_duplicates() if block_rows else pd.DataFrame()
    errors_df = pd.DataFrame(errors)

    summary_path = cli.outdir / "joint_rawmu_offset_summary.csv"
    offsets_path = cli.outdir / "joint_rawmu_offset_offsets.csv"
    blocks_path = cli.outdir / "joint_rawmu_offset_blocks.csv"
    errors_path = cli.outdir / "joint_rawmu_offset_errors.csv"
    report_path = cli.outdir / "joint_rawmu_offset_report.md"
    config_path = cli.outdir / "joint_rawmu_offset_config.json"

    summary.to_csv(summary_path, index=False)
    offsets_df.to_csv(offsets_path, index=False)
    blocks_df.to_csv(blocks_path, index=False)
    errors_df.to_csv(errors_path, index=False)
    config_path.write_text(json.dumps(vars(cli), indent=2, default=str), encoding="utf-8")
    write_report(report_path, summary, offsets_df, blocks_df, cli)

    print(f"Saved summary: {summary_path}")
    print(f"Saved offsets: {offsets_path}")
    print(f"Saved blocks: {blocks_path}")
    print(f"Saved report: {report_path}")
    if not errors_df.empty:
        print(f"Saved errors: {errors_path}")


if __name__ == "__main__":
    main()
