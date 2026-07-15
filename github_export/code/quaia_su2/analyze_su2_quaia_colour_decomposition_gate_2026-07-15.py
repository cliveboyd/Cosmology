from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_colour_decomposition_gate_20260715"
BASE_SCRIPT_NAME = "analyze_su2_quaia_matched_quality_gate_2026-07-15.py"


def load_base_module() -> Any:
    candidates = [
        Path(__file__).with_name(BASE_SCRIPT_NAME),
        ROOT / "github_export" / "code" / "quaia_su2" / BASE_SCRIPT_NAME,
        Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs") / BASE_SCRIPT_NAME,
    ]
    for path in candidates:
        if path.exists():
            spec = importlib.util.spec_from_file_location("su2_matched_quality_base", path)
            if spec is None or spec.loader is None:
                continue
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            return module
    raise FileNotFoundError(f"Could not find {BASE_SCRIPT_NAME}")


BASE = load_base_module()


def quality_groups() -> list[dict[str, Any]]:
    return [
        {"group": "zerr_only", "variables": ["zerr"]},
        {"group": "gaia_g_only", "variables": ["g"]},
        {"group": "gaia_bp_rp_only", "variables": ["bp_rp"]},
        {"group": "gaia_g_plus_bp_rp", "variables": ["g", "bp_rp"]},
        {"group": "wise_w1_only", "variables": ["w1"]},
        {"group": "wise_w2_only", "variables": ["w2"]},
        {"group": "wise_w1_w2_only", "variables": ["w1_w2"]},
        {"group": "wise_w1_plus_w2", "variables": ["w1", "w2"]},
        {"group": "wise_w1_w2_plus_mags", "variables": ["w1", "w2", "w1_w2"]},
        {"group": "gaia_wise_colour_cross", "variables": ["bp_rp", "w1_w2"]},
        {"group": "zerr_x_gaia_colour", "variables": ["zerr", "bp_rp"]},
        {"group": "zerr_x_wise_colour", "variables": ["zerr", "w1_w2"]},
        {"group": "zerr_x_gaia_wise_colour", "variables": ["zerr", "bp_rp", "w1_w2"]},
        {"group": "colour_redshift_estimator_compact", "variables": ["zerr", "g", "bp_rp", "w1", "w2", "w1_w2"]},
        {"group": "pm_error_quality", "variables": ["pmra_error", "pmdec_error"]},
        {"group": "all_catalogue_quality", "variables": ["zerr", "g", "bp_rp", "w1", "w2", "w1_w2", "pmra_error", "pmdec_error"]},
    ]


GROUP_NOTES = {
    "zerr_only": "Quaia redshift-uncertainty control.",
    "gaia_g_only": "Gaia optical magnitude only.",
    "gaia_bp_rp_only": "Gaia optical colour only.",
    "gaia_g_plus_bp_rp": "Gaia optical magnitude plus colour.",
    "wise_w1_only": "WISE W1 infrared magnitude only.",
    "wise_w2_only": "WISE W2 infrared magnitude only.",
    "wise_w1_w2_only": "WISE infrared colour only.",
    "wise_w1_plus_w2": "WISE infrared magnitudes without explicit colour column.",
    "wise_w1_w2_plus_mags": "WISE infrared magnitude-colour block.",
    "gaia_wise_colour_cross": "Gaia and WISE colour cross-block; polynomial features include their interaction.",
    "zerr_x_gaia_colour": "Redshift-uncertainty by Gaia-colour interaction block.",
    "zerr_x_wise_colour": "Redshift-uncertainty by WISE-colour interaction block.",
    "zerr_x_gaia_wise_colour": "Compact redshift-estimator plus optical/infrared colour interaction block.",
    "colour_redshift_estimator_compact": "All optical/infrared colour-magnitude variables plus redshift uncertainty, excluding proper-motion errors.",
    "pm_error_quality": "Proper-motion error control.",
    "all_catalogue_quality": "Full catalogue-quality block used as the current promotion gate.",
}


def fmt(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.6g}" if math.isfinite(value) else ""
    return str(value)


def markdown_table(rows: list[dict[str, Any]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines)


def pvalue_lookup(pvalues: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {row["fit_group"]: row for row in pvalues}


def ranked_suppressors(summaries: list[dict[str, Any]], pvalues: list[dict[str, Any]]) -> list[dict[str, Any]]:
    pvals = pvalue_lookup(pvalues)
    rows: list[dict[str, Any]] = []
    for row in summaries:
        group = row["fit_group"]
        if group == "baseline_unweighted":
            continue
        p = pvals.get(group, {})
        rows.append(
            {
                "rank_by_amp_suppression": 0,
                "fit_group": group,
                "variables": ",".join(next(item["variables"] for item in quality_groups() if item["group"] == group)),
                "hypothesis": GROUP_NOTES.get(group, ""),
                "amp_ratio_vs_baseline": row.get("amp_ratio_vs_baseline", float("nan")),
                "max_direction_sep_vs_baseline_deg": row.get("max_direction_sep_vs_baseline_deg", float("nan")),
                "snr_max": row.get("snr_max", float("nan")),
                "max_abs_smd_after": row.get("max_abs_smd_after", float("nan")),
                "max_weighted_ks_after": row.get("max_weighted_ks_after", float("nan")),
                "p_joint_snr_and_pair_sep": p.get("p_joint_snr_and_pair_sep", float("nan")),
                "readout": row.get("readout", ""),
            }
        )
    rows.sort(key=lambda item: (float(item["amp_ratio_vs_baseline"]), -float(item["max_direction_sep_vs_baseline_deg"])))
    for idx, row in enumerate(rows, start=1):
        row["rank_by_amp_suppression"] = idx
    return rows


def compact_rows(rows: list[dict[str, Any]], keep: int = 10) -> list[dict[str, Any]]:
    return rows[:keep]


def write_report(path: Path, summaries: list[dict[str, Any]], pvalues: list[dict[str, Any]], ranked: list[dict[str, Any]], args: argparse.Namespace) -> None:
    nonbase = [row for row in summaries if row["fit_group"] != "baseline_unweighted"]
    lines = [
        "# SU2 / Quaia Colour-Decomposition Gate",
        "",
        "Date: July 15, 2026",
        "",
        "## Purpose",
        "",
        "This follow-up decomposes the catalogue-quality suppression found in the matched-quality gate. The locked target remains the Quaia `0.95 <= z < 1.45` angular mode, with the same latitude cuts and redshift-dipole model. Only the matched feature groups change.",
        "",
        "The aim is to distinguish four explanations:",
        "",
        "- Gaia optical magnitude or colour.",
        "- WISE infrared magnitude or colour.",
        "- Quaia redshift-uncertainty structure.",
        "- Interactions between colour space and the Quaia redshift estimator.",
        "",
        "## Model",
        "",
        "```text",
        "z_i       = c0 + d . n_i + epsilon_i",
        "A_dipole  = ||d||",
        "SNR_d     = sqrt(d^T Cov(d)^(-1) d)",
        "p_i       = P(H_i = 1 | x_i)",
        "w_i       = 1 - p_i,  H_i = 1",
        "w_i       = p_i,      H_i = 0",
        "```",
        "",
        "Here `H_i` is the hemisphere defined by the baseline dipole axis for the relevant latitude cut, and `x_i` is the tested catalogue-quality feature block. Polynomial degree-2 propensity features are used, so multi-variable blocks include squares and pairwise interactions.",
        "",
        "## Locked Configuration",
        "",
        f"- redshift window: `{args.z_min} <= z < {args.z_max}`",
        f"- latitude cuts: `{', '.join(str(x) for x in args.bcuts)}`",
        f"- matched-feature shuffle mocks per group: `{args.n_mocks}`",
        f"- propensity strata for redshift shuffles: `{args.n_strata}`",
        "",
        "## Ranked Suppressors",
        "",
        markdown_table(
            compact_rows(ranked, 12),
            [
                "rank_by_amp_suppression",
                "fit_group",
                "amp_ratio_vs_baseline",
                "max_direction_sep_vs_baseline_deg",
                "snr_max",
                "p_joint_snr_and_pair_sep",
                "readout",
            ],
        ),
        "",
        "## Full Summary",
        "",
        markdown_table(
            nonbase,
            [
                "fit_group",
                "amp_ratio_vs_baseline",
                "max_direction_sep_vs_baseline_deg",
                "snr_max",
                "max_abs_smd_after",
                "max_weighted_ks_after",
                "readout",
            ],
        ),
        "",
        "## Matched-Feature Shuffle P-Values",
        "",
        markdown_table(
            pvalues,
            [
                "fit_group",
                "n_mocks",
                "p_snr_max_ge_observed",
                "p_pair_sep_le_observed",
                "p_coherence_ge_observed",
                "p_joint_snr_and_pair_sep",
            ],
        ),
        "",
        "## Interpretation Rule",
        "",
        "A feature block is treated as a plausible catalogue explanation when it both balances successfully and strongly suppresses the mode. In practice, the most informative rows are those with low amplitude ratio, low SNR, and large direction rotation after good post-weight balance.",
        "",
        "## Outputs",
        "",
        f"- ranked suppressors: `{args.out / 'su2_quaia_colour_decomposition_ranked_suppressors.csv'}`",
        f"- summary: `{args.out / 'su2_quaia_colour_decomposition_summary.csv'}`",
        f"- fit rows: `{args.out / 'su2_quaia_colour_decomposition_fit_rows.csv'}`",
        f"- balance diagnostics: `{args.out / 'su2_quaia_colour_decomposition_balance.csv'}`",
        f"- p-values: `{args.out / 'su2_quaia_colour_decomposition_pvalues.csv'}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_readout(path: Path, ranked: list[dict[str, Any]], pvalues: list[dict[str, Any]]) -> None:
    top = ranked[0] if ranked else None
    compact = next((row for row in ranked if row["fit_group"] == "colour_redshift_estimator_compact"), None)
    all_quality = next((row for row in ranked if row["fit_group"] == "all_catalogue_quality"), None)
    if top is None:
        bottom = "The colour-decomposition gate did not produce ranked rows."
    else:
        bottom = (
            f"The strongest suppressor is `{top['fit_group']}` with amplitude ratio "
            f"`{fmt(top['amp_ratio_vs_baseline'])}` and direction rotation "
            f"`{fmt(top['max_direction_sep_vs_baseline_deg'])}` deg."
        )
    lines = [
        "# SU2 / Quaia Colour-Decomposition Gate Readout",
        "",
        "Date: July 15, 2026",
        "",
        "## Bottom Line",
        "",
        bottom,
        "",
        "## Ranked Suppressors",
        "",
        markdown_table(
            compact_rows(ranked, 12),
            [
                "rank_by_amp_suppression",
                "fit_group",
                "amp_ratio_vs_baseline",
                "max_direction_sep_vs_baseline_deg",
                "snr_max",
                "p_joint_snr_and_pair_sep",
                "readout",
            ],
        ),
    ]
    if compact is not None:
        lines.extend(
            [
                "",
                "## Compact Colour-Redshift Estimator Block",
                "",
                f"- amplitude ratio vs baseline: `{fmt(compact['amp_ratio_vs_baseline'])}`",
                f"- maximum direction rotation: `{fmt(compact['max_direction_sep_vs_baseline_deg'])}` deg",
                f"- SNR max: `{fmt(compact['snr_max'])}`",
                f"- matched-feature joint p-value: `{fmt(compact['p_joint_snr_and_pair_sep'])}`",
            ]
        )
    if all_quality is not None:
        lines.extend(
            [
                "",
                "## Full Catalogue-Quality Block",
                "",
                f"- amplitude ratio vs baseline: `{fmt(all_quality['amp_ratio_vs_baseline'])}`",
                f"- maximum direction rotation: `{fmt(all_quality['max_direction_sep_vs_baseline_deg'])}` deg",
                f"- SNR max: `{fmt(all_quality['snr_max'])}`",
                f"- matched-feature joint p-value: `{fmt(all_quality['p_joint_snr_and_pair_sep'])}`",
            ]
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Colour-decomposition follow-up for the SU2/Quaia matched-quality gate.")
    parser.add_argument("--quaia-fits", type=Path, default=BASE.DEFAULT_QUAIA)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--z-min", type=float, default=0.95)
    parser.add_argument("--z-max", type=float, default=1.45)
    parser.add_argument("--bcuts", type=float, nargs="+", default=[10, 15, 20, 25, 30, 35])
    parser.add_argument("--n-mocks", type=int, default=300)
    parser.add_argument("--n-strata", type=int, default=10)
    parser.add_argument("--seed", type=int, default=170716)
    parser.add_argument("--poly-degree", type=int, default=2)
    parser.add_argument("--progress-every", type=int, default=50)
    parser.add_argument("--pass-amp-ratio", type=float, default=0.8)
    parser.add_argument("--pass-direction-deg", type=float, default=30.0)
    parser.add_argument("--pass-smd", type=float, default=0.1)
    parser.add_argument("--pass-ks", type=float, default=0.1)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    BASE.quality_groups = quality_groups
    fit_rows, summaries, balance, mocks, pvalues, config = BASE.run(args)
    config.update(
        {
            "analysis": "su2_quaia_colour_decomposition_gate",
            "purpose": "Decompose the matched catalogue-quality suppression into Gaia colour, WISE colour, redshift-uncertainty, and interaction blocks.",
            "group_notes": GROUP_NOTES,
        }
    )
    ranked = ranked_suppressors(summaries, pvalues)
    BASE.write_csv(args.out / "su2_quaia_colour_decomposition_fit_rows.csv", fit_rows)
    BASE.write_csv(args.out / "su2_quaia_colour_decomposition_summary.csv", summaries)
    BASE.write_csv(args.out / "su2_quaia_colour_decomposition_balance.csv", balance)
    BASE.write_csv(args.out / "su2_quaia_colour_decomposition_mock_summary.csv", mocks)
    BASE.write_csv(args.out / "su2_quaia_colour_decomposition_pvalues.csv", pvalues)
    BASE.write_csv(args.out / "su2_quaia_colour_decomposition_ranked_suppressors.csv", ranked)
    (args.out / "su2_quaia_colour_decomposition_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    write_report(args.out / "su2_quaia_colour_decomposition_report.md", summaries, pvalues, ranked, args)
    write_readout(args.out / "su2_quaia_colour_decomposition_readout.md", ranked, pvalues)
    print(f"Saved ranked suppressors: {args.out / 'su2_quaia_colour_decomposition_ranked_suppressors.csv'}", flush=True)
    print(f"Saved summary: {args.out / 'su2_quaia_colour_decomposition_summary.csv'}", flush=True)
    print(f"Saved report: {args.out / 'su2_quaia_colour_decomposition_report.md'}", flush=True)
    print(f"Saved readout: {args.out / 'su2_quaia_colour_decomposition_readout.md'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
