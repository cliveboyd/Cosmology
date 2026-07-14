r"""
Run the raw-MU FR/no-expansion diagnostic across the available SN datasets.

This is the project-level comparison wrapper for Peter Lamb's updated
distance-modulus proposal. It reuses `diagnose_pantheon_rawmu_fr.py` and
produces one combined evidence table for:

    * Pantheon+SH0ES
    * DES SN5YR Dovekie raw Hubble diagram
    * DES SN5YR Dovekie PLAMB-prepared table
    * Union3 compressed MU nodes, where available

The key question is whether the direct-MU distance relation prefers Peter's
travel-time correction

    DL_proxy = 10 ** ((MU - 25) / 5) / (1 + z/2)

and whether the implied H0 remains close to 67.5 km/s/Mpc across datasets.

Outputs:
    plamb_runs/diagnostics/rawmu_fr_multidataset/rawmu_fr_multidataset_fit_summary.csv
    plamb_runs/diagnostics/rawmu_fr_multidataset/rawmu_fr_multidataset_beta_scan.csv
    plamb_runs/diagnostics/rawmu_fr_multidataset/rawmu_fr_multidataset_key_results.csv
    plamb_runs/diagnostics/rawmu_fr_multidataset/rawmu_fr_multidataset_summary.md
"""

from __future__ import annotations

import argparse
from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd

from diagnose_pantheon_rawmu_fr import (
    DEFAULT_PANTHEON,
    DatasetSpec,
    analyse_dataset,
    make_plots,
    result_to_dict,
    sanitize_label,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "rawmu_fr_multidataset"

DEFAULT_DATASETS = [
    DatasetSpec("PantheonPlusSH0ES", DEFAULT_PANTHEON),
    DatasetSpec(
        "DES_Dovekie_raw",
        ROOT
        / "external_datasets"
        / "current_cosmology_sources"
        / "DES_SN5YR_Dovekie"
        / "DES-Dovekie_HD.csv",
    ),
    DatasetSpec(
        "DES_Dovekie_PLamb",
        ROOT
        / "external_datasets"
        / "current_cosmology_sources"
        / "prepared_for_plamb"
        / "DES_SN5YR_Dovekie_PLamb.dat",
    ),
    DatasetSpec(
        "Union3_compressed",
        ROOT
        / "external_datasets"
        / "current_cosmology_sources"
        / "Union3_release"
        / "mu_mat_union3_cosmo=2_mu.fits",
    ),
    DatasetSpec(
        "Union3p1_UNITY1p8",
        ROOT
        / "external_datasets"
        / "current_cosmology_sources"
        / "Union3_release"
        / "mu_mat_union3.1_UNITY1.8_template_cosmo=2_0_mu.fits",
    ),
]


def parse_dataset(raw: str) -> DatasetSpec:
    if "=" in raw:
        label, path = raw.split("=", 1)
        return DatasetSpec(sanitize_label(label), Path(path))
    path = Path(raw)
    return DatasetSpec(sanitize_label(path.stem), path)


def build_base_args(args: argparse.Namespace, outdir: Path) -> Namespace:
    return Namespace(
        outdir=outdir,
        z_cols=args.z_cols,
        mu_col=None,
        mu_err_col=None,
        calibrator_col=None,
        hf_col=None,
        prob_col=None,
        prob_min=args.prob_min,
        zmin=args.zmin,
        lowz_max=args.lowz_max,
        midz_max=args.midz_max,
        min_n=args.min_n,
        beta_min=args.beta_min,
        beta_max=args.beta_max,
        beta_steps=args.beta_steps,
        weight_mode=args.weight_mode,
        no_plots=args.no_plots,
        plot_max_points=args.plot_max_points,
    )


def add_offsets(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["H0_minus_67p5"] = out["H0_km_s_Mpc"] - 67.5
    out["abs_H0_minus_67p5"] = out["H0_minus_67p5"].abs()
    out["beta_minus_0p5"] = out["beta"] - 0.5
    out["abs_beta_minus_0p5"] = out["beta_minus_0p5"].abs()
    return out


def key_result_rows(fit_df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    headline_subsets = [
        "SH0ES_HF_flag",
        f"lowz_{args.zmin:g}_{args.lowz_max:g}_noncal",
        f"full_noncal_zgt{args.zmin:g}",
        "all_finite_zpos",
    ]
    headline_corrections = [
        "raw_DL_from_MU",
        "DL_div_1plus_z_over_2",
        "DL_div_1plus_beta_z_best",
    ]
    keep = fit_df[
        (fit_df["fit_type"] == "through_origin")
        & (fit_df["subset"].isin(headline_subsets))
        & (fit_df["correction"].isin(headline_corrections))
    ].copy()
    if keep.empty:
        return keep

    keep = add_offsets(keep)
    sort_cols = [
        "dataset",
        "z_col",
        "subset",
        "weighting",
        "correction",
    ]
    return keep.sort_values(sort_cols)


def best_by_dataset(fit_df: pd.DataFrame) -> pd.DataFrame:
    best = fit_df[
        (fit_df["fit_type"] == "through_origin")
        & (fit_df["correction"] == "DL_div_1plus_beta_z_best")
        & (fit_df["weighting"] == "weighted_diag_mu")
    ].copy()
    if best.empty:
        return best
    best = add_offsets(best)
    best["score"] = best["abs_H0_minus_67p5"] + 3.0 * best["abs_beta_minus_0p5"]
    idx = best.groupby("dataset")["score"].idxmin()
    return best.loc[idx].sort_values(["dataset", "score"])


def write_multidataset_summary(
    path: Path,
    fit_df: pd.DataFrame,
    beta_df: pd.DataFrame,
    key_df: pd.DataFrame,
    readiness_df: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    lines: list[str] = []
    lines.append("# Raw-MU FR Multidataset Diagnostic")
    lines.append("")
    lines.append("This report reruns Peter Lamb's updated raw distance-modulus test across the available supernova distance datasets.")
    lines.append("")
    lines.append("Core rule under test:")
    lines.append("")
    lines.append("`DL_proxy = 10 ** ((MU - 25) / 5) / (1 + z/2)`")
    lines.append("")
    lines.append("The report also tests raw `DL`, `/ (1+z)`, and fitted `/ (1+beta z)` variants.")
    lines.append("")

    if not readiness_df.empty:
        lines.append("## Dataset Readiness")
        lines.append("")
        lines.append(readiness_df.to_markdown(index=False))
        lines.append("")

    if not key_df.empty:
        preferred_cols = [
            "dataset",
            "z_col",
            "subset",
            "weighting",
            "correction",
            "N",
            "beta",
            "H0_km_s_Mpc",
            "H0_minus_67p5",
            "chi2_dof",
            "frac_rms",
        ]
        lines.append("## Key Results")
        lines.append("")
        lines.append(key_df[preferred_cols].to_markdown(index=False, floatfmt=".5g"))
        lines.append("")

    best = best_by_dataset(fit_df)
    if not best.empty:
        cols = [
            "dataset",
            "z_col",
            "subset",
            "N",
            "beta",
            "beta_minus_0p5",
            "H0_km_s_Mpc",
            "H0_minus_67p5",
            "chi2_dof",
            "frac_rms",
        ]
        lines.append("## Closest Weighted Fitted-Beta Rows")
        lines.append("")
        lines.append("Rows are selected by closeness to `H0=67.5` and `beta=0.5`; this is a diagnostic ranking, not a likelihood preference.")
        lines.append("")
        lines.append(best[cols].to_markdown(index=False, floatfmt=".5g"))
        lines.append("")

    lines.append("## Interpretation Notes")
    lines.append("")
    lines.append("- DES SN5YR is the cleanest test of Peter's quoted DES behaviour because it supplies `MU`, `MUERR`, `zHD`, and `zHEL` directly.")
    lines.append("- Pantheon+ also contains `MU_SH0ES`, but its calibration and weighting choices move the result away from the simple low-z unweighted `H0 ~ 67.5` behaviour.")
    lines.append("- Union3 compressed nodes are not individual supernovae; treat them as a compressed distance-product cross-check, not a raw SN Hubble diagram.")
    lines.append("- These searches use diagonal distance-modulus errors or unweighted fits, not the full covariance. A full likelihood branch is the next step if this diagnostic remains stable.")
    lines.append("")
    lines.append("## Output Files")
    lines.append("")
    lines.append(f"- Fit summary rows: `{len(fit_df)}`")
    lines.append(f"- Beta scan rows: `{len(beta_df)}`")
    lines.append(f"- Key result rows: `{len(key_df)}`")
    lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--dataset", action="append", default=None, help="Extra or replacement dataset as LABEL=PATH.")
    parser.add_argument("--only-defaults", action="store_true", help="Ignore --dataset extras and use built-in defaults only.")
    parser.add_argument("--replace-defaults", action="store_true", help="Use only datasets passed via --dataset.")
    parser.add_argument("--z-cols", default="auto")
    parser.add_argument("--prob-min", type=float, default=None)
    parser.add_argument("--zmin", type=float, default=0.01)
    parser.add_argument("--lowz-max", type=float, default=0.15)
    parser.add_argument("--midz-max", type=float, default=0.7)
    parser.add_argument("--min-n", type=int, default=5)
    parser.add_argument("--beta-min", type=float, default=-0.25)
    parser.add_argument("--beta-max", type=float, default=1.25)
    parser.add_argument("--beta-steps", type=int, default=301)
    parser.add_argument("--weight-mode", choices=["weighted", "unweighted", "both"], default="both")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--plot-max-points", type=int, default=5000)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    if args.replace_defaults:
        datasets = []
    else:
        datasets = list(DEFAULT_DATASETS)
    if args.dataset and not args.only_defaults:
        datasets.extend(parse_dataset(raw) for raw in args.dataset)

    base_args = build_base_args(args, outdir)
    all_results = []
    all_beta_rows: list[dict[str, object]] = []
    all_point_rows: list[dict[str, object]] = []
    readiness_rows: list[dict[str, object]] = []

    for dataset in datasets:
        if not dataset.path.exists():
            readiness_rows.append(
                {
                    "dataset": dataset.label,
                    "path": str(dataset.path),
                    "status": "missing",
                    "note": "File not found; skipped.",
                }
            )
            continue
        try:
            results, beta_rows, point_rows = analyse_dataset(dataset, base_args, outdir)
            all_results.extend(results)
            all_beta_rows.extend(beta_rows)
            all_point_rows.extend(point_rows)
            readiness_rows.append(
                {
                    "dataset": dataset.label,
                    "path": str(dataset.path),
                    "status": "ok",
                    "note": f"{len(results)} fit rows, {len(point_rows)} point rows.",
                }
            )
        except Exception as exc:
            readiness_rows.append(
                {
                    "dataset": dataset.label,
                    "path": str(dataset.path),
                    "status": "error",
                    "note": str(exc),
                }
            )

    fit_df = pd.DataFrame([result_to_dict(result) for result in all_results])
    beta_df = pd.DataFrame(all_beta_rows)
    point_df = pd.DataFrame(all_point_rows)
    readiness_df = pd.DataFrame(readiness_rows)
    key_df = key_result_rows(fit_df, args) if not fit_df.empty else pd.DataFrame()

    fit_path = outdir / "rawmu_fr_multidataset_fit_summary.csv"
    beta_path = outdir / "rawmu_fr_multidataset_beta_scan.csv"
    point_path = outdir / "rawmu_fr_multidataset_point_table.csv"
    readiness_path = outdir / "rawmu_fr_multidataset_readiness.csv"
    key_path = outdir / "rawmu_fr_multidataset_key_results.csv"
    summary_path = outdir / "rawmu_fr_multidataset_summary.md"

    fit_df.to_csv(fit_path, index=False)
    beta_df.to_csv(beta_path, index=False)
    point_df.to_csv(point_path, index=False)
    readiness_df.to_csv(readiness_path, index=False)
    key_df.to_csv(key_path, index=False)
    write_multidataset_summary(summary_path, fit_df, beta_df, key_df, readiness_df, args)
    if not point_df.empty and not fit_df.empty:
        for path in make_plots(outdir, fit_df, point_df, base_args):
            print(f"Saved plot: {path}")

    print(f"Saved fit summary: {fit_path}")
    print(f"Saved beta scan: {beta_path}")
    print(f"Saved point table: {point_path}")
    print(f"Saved readiness: {readiness_path}")
    print(f"Saved key results: {key_path}")
    print(f"Saved summary: {summary_path}")

    if not key_df.empty:
        view = key_df[
            (key_df["weighting"] == "weighted_diag_mu")
            & (key_df["correction"].isin(["DL_div_1plus_z_over_2", "DL_div_1plus_beta_z_best"]))
        ].copy()
        view = view.sort_values(["abs_H0_minus_67p5", "abs_beta_minus_0p5"]).head(12)
        print("\nClosest weighted raw-MU FR rows:")
        for _, row in view.iterrows():
            print(
                f"{row['dataset']} {row['z_col']} {row['subset']} {row['correction']}: "
                f"N={int(row['N'])} beta={row['beta']:.3f} "
                f"H0={row['H0_km_s_Mpc']:.3f} "
                f"dH0={row['H0_minus_67p5']:+.3f} "
                f"chi2/dof={row['chi2_dof']:.3f}"
            )


if __name__ == "__main__":
    main()

