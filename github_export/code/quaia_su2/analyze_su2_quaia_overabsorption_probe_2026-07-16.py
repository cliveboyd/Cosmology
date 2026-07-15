from __future__ import annotations

import csv
import hashlib
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
OUTPUTS = Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs")
RESULT_DIR = ROOT / "github_export" / "results" / "2026-07-16" / "su2_quaia"
CODE_PATH = ROOT / "github_export" / "code" / "quaia_su2" / "analyze_su2_quaia_overabsorption_probe_2026-07-16.py"
TRIALS = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_injection_recovery_overnight_20260715" / "su2_quaia_injection_recovery_trials.csv"

DATE = "2026-07-16"
PREFIX = "su2_quaia_overabsorption_probe_2026-07-16"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin(["true", "1", "yes"])


def q(series: pd.Series, prob: float) -> float:
    values = pd.to_numeric(series, errors="coerce")
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return float("nan")
    return float(values.quantile(prob))


def rate(series: pd.Series) -> float:
    if len(series) == 0:
        return float("nan")
    return float(np.mean(series.astype(float)))


def fmt(value: object) -> str:
    if isinstance(value, (float, np.floating)):
        if not math.isfinite(float(value)):
            return ""
        return f"{float(value):.6g}"
    return str(value)


def markdown_table(rows: list[dict[str, object]], columns: list[str]) -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines)


def category_for(label: str) -> str:
    if label.startswith("locked"):
        return "locked"
    if label.startswith("random"):
        return "random"
    if label in {"cmb", "anti_cmb"}:
        return "cmb_or_anti"
    return "other"


def summarise_group(group: pd.DataFrame) -> pd.Series:
    strict = group["scanlaw_strict_physical_pass"]
    loose = group["scanlaw_loose_physical_pass"]
    bic = group["scanlaw_bic_accept"]
    return pd.Series(
        {
            "n_trials": int(len(group)),
            "raw_pass_rate": rate(group["raw_pass"]),
            "matched_pass_rate": rate(group["matched_pass"]),
            "scanlaw_strict_physical_rate": rate(strict),
            "scanlaw_loose_physical_rate": rate(loose),
            "scanlaw_bic_accept_rate": rate(bic),
            "scanlaw_colour_pass_rate": rate(group["scanlaw_colour_pass"]),
            "strict_physical_bic_rejected_rate": rate(strict & ~bic),
            "matched_to_strict_physical_loss_rate": rate(group["matched_pass"] & ~strict),
            "high_snr_rate": rate(group["scanlaw_colour_dipole_snr"] >= 3.0),
            "amp_too_high_rate": rate(group["scanlaw_colour_amp_ratio"] > 1.5),
            "amp_too_low_rate": rate(group["scanlaw_colour_amp_ratio"] < 0.5),
            "sep_too_high_rate": rate(group["scanlaw_colour_sep_from_injected_deg"] > 30.0),
            "delta_bic_p05": q(group["scanlaw_colour_delta_bic"], 0.05),
            "delta_bic_p50": q(group["scanlaw_colour_delta_bic"], 0.50),
            "delta_bic_p95": q(group["scanlaw_colour_delta_bic"], 0.95),
            "scanlaw_snr_p50": q(group["scanlaw_colour_dipole_snr"], 0.50),
            "scanlaw_sep_deg_p50": q(group["scanlaw_colour_sep_from_injected_deg"], 0.50),
            "scanlaw_amp_ratio_p50": q(group["scanlaw_colour_amp_ratio"], 0.50),
        }
    )


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def build_outputs() -> list[Path]:
    df = pd.read_csv(TRIALS)
    for col in ["raw_pass", "matched_pass", "scanlaw_colour_pass"]:
        df[col] = as_bool(df[col])

    nonzero = df[df["amp_scale"] > 0.0].copy()
    nonzero["category"] = nonzero["direction_label"].astype(str).map(category_for)
    nonzero["scanlaw_strict_physical_pass"] = (
        nonzero["scanlaw_colour_amp_ratio"].between(0.5, 1.5)
        & (nonzero["scanlaw_colour_sep_from_injected_deg"] <= 30.0)
        & (nonzero["scanlaw_colour_dipole_snr"] >= 3.0)
    )
    nonzero["scanlaw_loose_physical_pass"] = (
        nonzero["scanlaw_colour_amp_ratio"].between(0.25, 2.5)
        & (nonzero["scanlaw_colour_sep_from_injected_deg"] <= 45.0)
        & (nonzero["scanlaw_colour_dipole_snr"] >= 3.0)
    )
    nonzero["scanlaw_bic_accept"] = nonzero["scanlaw_colour_delta_bic"] < -10.0

    metric_cols = [
        "raw_pass",
        "matched_pass",
        "scanlaw_colour_pass",
        "scanlaw_colour_amp_ratio",
        "scanlaw_colour_sep_from_injected_deg",
        "scanlaw_colour_dipole_snr",
        "scanlaw_colour_delta_bic",
        "scanlaw_strict_physical_pass",
        "scanlaw_loose_physical_pass",
        "scanlaw_bic_accept",
    ]
    amp_summary = (
        nonzero.groupby("amp_scale", sort=True)[metric_cols]
        .apply(summarise_group)
        .reset_index()
    )
    category_amp_summary = (
        nonzero.groupby(["category", "amp_scale"], sort=True)[metric_cols]
        .apply(summarise_group)
        .reset_index()
    )
    one_x = nonzero[nonzero["amp_scale"] == 1.0].copy()
    condition_summary = (
        one_x.groupby(["direction_label", "category", "bcut_deg"], sort=True)[metric_cols]
        .apply(summarise_group)
        .reset_index()
    )

    top_false_negative = condition_summary.sort_values(
        ["strict_physical_bic_rejected_rate", "scanlaw_strict_physical_rate"],
        ascending=[False, False],
    ).head(12)

    locked_rows = condition_summary[condition_summary["category"] == "locked"].copy()
    category_one_x = category_amp_summary[category_amp_summary["amp_scale"] == 1.0].copy()

    output_paths: list[Path] = []
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    csv_tables = {
        f"{PREFIX}_amplitude_summary.csv": amp_summary,
        f"{PREFIX}_category_amplitude_summary.csv": category_amp_summary,
        f"{PREFIX}_condition_summary.csv": condition_summary,
        f"{PREFIX}_top_bic_false_negatives.csv": top_false_negative,
    }
    for name, table in csv_tables.items():
        out_path = OUTPUTS / name
        table.to_csv(out_path, index=False)
        output_paths.append(out_path)
        shutil.copy2(out_path, RESULT_DIR / name)

    one_x_all = summarise_group(one_x).to_dict()
    matched = one_x_all["matched_pass_rate"]
    strict = one_x_all["scanlaw_strict_physical_rate"]
    loose = one_x_all["scanlaw_loose_physical_rate"]
    colour_pass = one_x_all["scanlaw_colour_pass_rate"]
    bic_reject = one_x_all["strict_physical_bic_rejected_rate"]
    loss = one_x_all["matched_to_strict_physical_loss_rate"]

    bottom = (
        "The scan-law+colour gate is over-rejecting injected angular modes at the current 1x "
        "amplitude. The strict physical recovery rate is "
        f"{strict:.4f}, but the final scan-law+colour pass rate is only {colour_pass:.4f}; "
        f"{bic_reject:.4f} of all trials are physically recovered and then rejected only by "
        "the BIC gate. The locked observed-direction injections show a separate failure mode: "
        "they are usually high-SNR after controls, but are rotated or amplitude-inflated enough "
        "to fail the strict recovery cuts."
    )

    short_cols = [
        "category",
        "amp_scale",
        "matched_pass_rate",
        "scanlaw_strict_physical_rate",
        "scanlaw_loose_physical_rate",
        "scanlaw_colour_pass_rate",
        "strict_physical_bic_rejected_rate",
        "matched_to_strict_physical_loss_rate",
        "delta_bic_p50",
        "scanlaw_snr_p50",
        "scanlaw_sep_deg_p50",
        "scanlaw_amp_ratio_p50",
    ]
    amp_cols = [
        "amp_scale",
        "matched_pass_rate",
        "scanlaw_strict_physical_rate",
        "scanlaw_loose_physical_rate",
        "scanlaw_colour_pass_rate",
        "strict_physical_bic_rejected_rate",
        "matched_to_strict_physical_loss_rate",
        "delta_bic_p50",
    ]
    locked_cols = [
        "direction_label",
        "bcut_deg",
        "matched_pass_rate",
        "scanlaw_strict_physical_rate",
        "scanlaw_loose_physical_rate",
        "scanlaw_colour_pass_rate",
        "strict_physical_bic_rejected_rate",
        "high_snr_rate",
        "amp_too_high_rate",
        "sep_too_high_rate",
        "delta_bic_p50",
        "scanlaw_snr_p50",
        "scanlaw_sep_deg_p50",
        "scanlaw_amp_ratio_p50",
    ]

    readout = [
        "# SU2 / Quaia Scan-Law + Colour Over-Absorption Probe",
        "",
        f"Date: {DATE}",
        "",
        "## Bottom Line",
        "",
        bottom,
        "",
        "## 1x Gate Decomposition",
        "",
        f"- matched-quality recovery rate: `{matched:.4f}`",
        f"- strict scan-law+colour physical recovery rate: `{strict:.4f}`",
        f"- loose scan-law+colour physical recovery rate: `{loose:.4f}`",
        f"- final scan-law+colour pass rate: `{colour_pass:.4f}`",
        f"- strict physical recovery rejected only by BIC: `{bic_reject:.4f}`",
        f"- matched recovery lost before strict physical scan-law recovery: `{loss:.4f}`",
        "",
        "## 1x By Direction Family",
        "",
        markdown_table(category_one_x[short_cols].to_dict("records"), short_cols),
    ]
    readout_path = OUTPUTS / f"{PREFIX}_readout.md"
    readout_path.write_text("\n".join(readout) + "\n", encoding="utf-8")
    output_paths.append(readout_path)
    shutil.copy2(readout_path, RESULT_DIR / readout_path.name)

    report = [
        "# SU2 / Quaia Scan-Law + Colour Over-Absorption Probe",
        "",
        f"Date: {DATE}",
        "",
        "## Purpose",
        "",
        "This post-processes the completed 120k-row injection-recovery table to separate three effects:",
        "",
        "1. catalogue-quality matched recovery,",
        "2. physical recovery after scan-law+colour controls,",
        "3. final BIC-gated scan-law+colour acceptance.",
        "",
        "Strict physical recovery uses the original non-zero injection rule: SNR >= 3, angular separation <= 30 deg, and recovered/injected amplitude ratio in [0.5, 1.5].",
        "",
        "Loose physical recovery is a diagnostic only: SNR >= 3, angular separation <= 45 deg, and recovered/injected amplitude ratio in [0.25, 2.5]. It is not proposed as a promotion gate.",
        "",
        "The BIC gate is the current likelihood requirement Delta BIC = BIC(full dipole + controls) - BIC(controls only) < -10.",
        "",
        "## Bottom Line",
        "",
        bottom,
        "",
        "## Amplitude Summary",
        "",
        markdown_table(amp_summary[amp_cols].to_dict("records"), amp_cols),
        "",
        "## Direction-Family Summary",
        "",
        markdown_table(category_amp_summary[short_cols].to_dict("records"), short_cols),
        "",
        "## Locked-Direction 1x Breakdown",
        "",
        markdown_table(locked_rows[locked_cols].to_dict("records"), locked_cols),
        "",
        "## Largest 1x BIC False-Negative Conditions",
        "",
        markdown_table(top_false_negative[locked_cols].to_dict("records"), locked_cols),
        "",
        "## Interpretation",
        "",
        "For random injected directions at 1x amplitude, scan-law+colour often still recovers a recognisable physical dipole, but the BIC term rejects almost all of those recoveries. This says the current BIC threshold is too stringent for a promotion gate unless calibrated by injection-recovery power.",
        "",
        "For locked observed-direction injections at 1x amplitude, the problem is not mainly a BIC penalty. The full controlled fit is high-SNR, but the direction and amplitude are distorted: median separations are commonly just above 30 deg and median amplitude ratios are often above 1.5. That points to collinearity or template leakage between the locked angular mode and the scan-law+colour controls.",
        "",
        "The follow-up gate should therefore use a two-part diagnostic: first calibrate detection power and false-negative rate with injections, then inspect template collinearity for the locked z ~ 1-1.5 mode before treating a non-pass as a physical rejection.",
    ]
    report_path = OUTPUTS / f"{PREFIX}_report.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    output_paths.append(report_path)
    shutil.copy2(report_path, RESULT_DIR / report_path.name)

    manifest_rows: list[dict[str, object]] = []
    for path in sorted(output_paths):
        manifest_rows.append(
            {
                "artifact": str((RESULT_DIR / path.name).relative_to(ROOT)).replace("\\", "/"),
                "bytes": (RESULT_DIR / path.name).stat().st_size,
                "sha256": sha256_file(RESULT_DIR / path.name),
                "source": str(TRIALS.relative_to(ROOT)).replace("\\", "/"),
                "note": "SU2/Quaia scan-law+colour over-absorption probe result artefact.",
            }
        )
    manifest_rows.append(
        {
            "artifact": str(TRIALS.relative_to(ROOT)).replace("\\", "/"),
            "bytes": TRIALS.stat().st_size,
            "sha256": sha256_file(TRIALS),
            "source": str(TRIALS.relative_to(ROOT)).replace("\\", "/"),
            "note": "Large injection-recovery trial table retained locally and referenced for reproducibility.",
        }
    )
    manifest_rows.append(
        {
            "artifact": str(CODE_PATH.relative_to(ROOT)).replace("\\", "/"),
            "bytes": CODE_PATH.stat().st_size,
            "sha256": sha256_file(CODE_PATH),
            "source": str(CODE_PATH.relative_to(ROOT)).replace("\\", "/"),
            "note": "Reproduction script for the over-absorption probe.",
        }
    )

    manifest_name = f"{PREFIX}_manifest.csv"
    manifest_path = OUTPUTS / manifest_name
    write_csv(manifest_path, manifest_rows)
    shutil.copy2(manifest_path, RESULT_DIR / manifest_name)
    shutil.copy2(manifest_path, ROOT / "github_export" / "manifests" / manifest_name)
    write_csv(RESULT_DIR / "manifest.csv", manifest_rows)

    return [*output_paths, manifest_path]


def main() -> int:
    paths = build_outputs()
    for path in paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
