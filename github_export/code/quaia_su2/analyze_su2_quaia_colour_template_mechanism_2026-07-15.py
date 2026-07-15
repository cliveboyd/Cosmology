from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
OUTPUTS = Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs")
CODE_DIR = ROOT / "github_export" / "code" / "quaia_su2"
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_colour_template_mechanism_20260715"
EXT_SCRIPT = "analyze_su2_quaia_external_dust_gate_2026-07-15.py"
WISE_SCRIPT = "analyze_su2_quaia_wise_stellar_gate_2026-07-15.py"


def load_module(name: str, candidates: list[Path]) -> Any:
    for path in candidates:
        if path.exists():
            spec = importlib.util.spec_from_file_location(name, path)
            if spec is None or spec.loader is None:
                continue
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            return module
    raise FileNotFoundError(f"Could not find module source for {name}")


EXT = load_module("su2_external_dust_gate", [OUTPUTS / EXT_SCRIPT, CODE_DIR / EXT_SCRIPT])
sys.modules["analyze_su2_quaia_external_dust_gate"] = EXT
WISE = load_module("su2_wise_stellar_gate", [OUTPUTS / WISE_SCRIPT, CODE_DIR / WISE_SCRIPT])


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.6g}" if math.isfinite(value) else ""
    return str(value)


def markdown_table(rows: list[dict[str, Any]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines)


def standardize(values: np.ndarray, mask: np.ndarray | None = None) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    subset = values if mask is None else values[mask]
    med = float(np.nanmedian(subset))
    scale = float(np.nanstd(subset))
    if not math.isfinite(scale) or scale <= 0.0:
        scale = 1.0
    return (values - med) / scale


def finite_mask(*arrays: np.ndarray) -> np.ndarray:
    mask = np.ones(len(arrays[0]), dtype=bool)
    for arr in arrays:
        mask &= np.isfinite(arr)
    return mask


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    mask = finite_mask(x, y)
    if int(np.sum(mask)) < 3:
        return float("nan")
    return float(np.corrcoef(x[mask], y[mask])[0, 1])


def unit_from_lb(l_deg: float, b_deg: float) -> np.ndarray:
    l = math.radians(l_deg)
    b = math.radians(b_deg)
    cb = math.cos(b)
    return np.array([cb * math.cos(l), cb * math.sin(l), math.sin(b)], dtype=float)


def angular_sep_deg(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 0.0:
        return float("nan")
    dot = float(np.clip(np.dot(a, b) / denom, -1.0, 1.0))
    return math.degrees(math.acos(dot))


def fit_response(response: np.ndarray, vec: np.ndarray, controls: list[np.ndarray] | None = None) -> dict[str, Any]:
    columns = [np.ones(len(response)), vec[:, 0], vec[:, 1], vec[:, 2]]
    if controls:
        columns.extend(standardize(control) for control in controls)
    x = np.column_stack(columns)
    beta, *_ = np.linalg.lstsq(x, response, rcond=None)
    fitted = x @ beta
    resid = response - fitted
    sst = float(np.dot(response - np.mean(response), response - np.mean(response)))
    sse = float(np.dot(resid, resid))
    dof = max(len(response) - x.shape[1], 1)
    sigma2 = sse / dof
    dvec = beta[1:4]
    amp = float(np.linalg.norm(dvec))
    l_deg, b_deg = EXT.unit_to_lb(dvec)
    try:
        cov = sigma2 * np.linalg.pinv(x.T @ x)
        cov_d = cov[1:4, 1:4]
        snr = float(math.sqrt(max(float(dvec @ np.linalg.pinv(cov_d) @ dvec), 0.0)))
    except Exception:
        snr = float("nan")
    return {
        "amp": amp,
        "l_deg": l_deg,
        "b_deg": b_deg,
        "formal_snr": snr,
        "r2": 1.0 - sse / sst if sst > 0.0 else float("nan"),
        "resid": resid,
        "n_params": int(x.shape[1]),
        "rms_resid": float(math.sqrt(max(sigma2, 0.0))),
    }


def residualize_response(response: np.ndarray, controls: list[np.ndarray]) -> dict[str, Any]:
    if not controls:
        resid = response - np.mean(response)
        return {"resid": resid, "r2_templates": 0.0, "n_template_params": 1}
    x = np.column_stack([np.ones(len(response))] + [standardize(control) for control in controls])
    beta, *_ = np.linalg.lstsq(x, response, rcond=None)
    fitted = x @ beta
    resid = response - fitted
    sst = float(np.dot(response - np.mean(response), response - np.mean(response)))
    sse = float(np.dot(resid, resid))
    return {
        "resid": resid,
        "r2_templates": 1.0 - sse / sst if sst > 0.0 else float("nan"),
        "n_template_params": int(x.shape[1]),
    }


def smd(values: np.ndarray, hemisphere: np.ndarray) -> float:
    h1 = values[hemisphere]
    h0 = values[~hemisphere]
    if len(h1) < 2 or len(h0) < 2:
        return float("nan")
    pooled = math.sqrt(max((float(np.var(h1)) + float(np.var(h0))) / 2.0, 1e-30))
    return (float(np.mean(h1)) - float(np.mean(h0))) / pooled


def mean_direction(rows: list[dict[str, Any]]) -> np.ndarray:
    dirs = np.array([unit_from_lb(row["l_deg"], row["b_deg"]) for row in rows], dtype=float)
    vec = np.sum(dirs, axis=0)
    norm = float(np.linalg.norm(vec))
    return vec / norm if norm > 0.0 else vec


def summarise_fit_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    raw_amp = np.array([row["raw_colour_amp"] for row in rows], dtype=float)
    resid_amp = np.array([row["resid_colour_amp"] for row in rows], dtype=float)
    r2 = np.array([row["template_r2"] for row in rows], dtype=float)
    raw_smd = np.array([abs(row["raw_colour_smd_on_z_axis"]) for row in rows], dtype=float)
    resid_smd = np.array([abs(row["resid_colour_smd_on_z_axis"]) for row in rows], dtype=float)
    direction_sep = np.array([row["resid_direction_sep_from_raw_deg"] for row in rows], dtype=float)
    alignment_raw = np.array([row["raw_colour_sep_from_z_dipole_deg"] for row in rows], dtype=float)
    alignment_resid = np.array([row["resid_colour_sep_from_z_dipole_deg"] for row in rows], dtype=float)
    raw_amp_mean = float(np.nanmean(raw_amp))
    resid_amp_mean = float(np.nanmean(resid_amp))
    raw_smd_mean = float(np.nanmean(raw_smd))
    resid_smd_mean = float(np.nanmean(resid_smd))
    amp_ratio = resid_amp_mean / raw_amp_mean if raw_amp_mean > 0.0 else float("nan")
    smd_ratio = resid_smd_mean / raw_smd_mean if raw_smd_mean > 0.0 else float("nan")
    dipole_absorption = 1.0 - amp_ratio if math.isfinite(amp_ratio) else float("nan")
    smd_absorption = 1.0 - smd_ratio if math.isfinite(smd_ratio) else float("nan")
    score_parts = [
        0.45 * float(np.clip(dipole_absorption, 0.0, 1.0)) if math.isfinite(dipole_absorption) else 0.0,
        0.35 * float(np.clip(smd_absorption, 0.0, 1.0)) if math.isfinite(smd_absorption) else 0.0,
        0.20 * float(np.clip(float(np.nanmean(r2)), 0.0, 1.0)),
    ]
    return {
        "response_field": rows[0]["response_field"],
        "template_group": rows[0]["template_group"],
        "n_bcuts": len(rows),
        "N_min": int(min(row["N"] for row in rows)),
        "N_max": int(max(row["N"] for row in rows)),
        "template_r2_mean": float(np.nanmean(r2)),
        "template_r2_max": float(np.nanmax(r2)),
        "raw_colour_amp_mean": raw_amp_mean,
        "resid_colour_amp_mean": resid_amp_mean,
        "resid_amp_ratio_vs_raw": amp_ratio,
        "dipole_absorption_fraction": dipole_absorption,
        "raw_abs_smd_mean": raw_smd_mean,
        "resid_abs_smd_mean": resid_smd_mean,
        "smd_absorption_fraction": smd_absorption,
        "max_resid_direction_sep_from_raw_deg": float(np.nanmax(direction_sep)),
        "raw_sep_from_z_dipole_mean_deg": float(np.nanmean(alignment_raw)),
        "resid_sep_from_z_dipole_mean_deg": float(np.nanmean(alignment_resid)),
        "mechanism_score": float(sum(score_parts)),
    }


def load_colour_data(args: argparse.Namespace) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    gaiascanlaw_data_dir = args.gaiascanlaw_data_dir or EXT.resolve_gaiascanlaw_data_dir()
    q = EXT.load_data(EXT.DEFAULT_QUAIA, EXT.DEFAULT_SELECTION, EXT.DEFAULT_RANDOMS, EXT.DEFAULT_SFD_DIR, gaiascanlaw_data_dir)
    with fits.open(EXT.DEFAULT_QUAIA, memmap=True) as hdul:
        data = hdul[1].data
        q["bp_rp"] = np.asarray(data["phot_bp_mean_mag"], dtype=float) - np.asarray(data["phot_rp_mean_mag"], dtype=float)
    q["w1_w2"] = q["w1"] - q["w2"]
    q.update(WISE.load_wise_depth(args.wise_inventory, q))
    gaia_density, gaia_density_source = WISE.load_gaia_density(
        args.gaia_density,
        q,
        args.gaia_density_hips_base_url,
        args.gaia_density_hips_cache_dir,
        args.gaia_density_hips_order,
    )
    q.update(gaia_density)
    base_mask = (
        (q["z"] >= args.z_min)
        & (q["z"] < args.z_max)
        & (np.abs(q["b"]) >= min(args.b_cuts))
        & finite_mask(q["bp_rp"], q["w1_w2"])
    )
    q["bp_rp_z"] = standardize(q["bp_rp"], base_mask)
    q["w1_w2_z"] = standardize(q["w1_w2"], base_mask)
    cross = q["bp_rp_z"] * q["w1_w2_z"]
    q["bp_rp_x_w1_w2_z"] = standardize(cross, base_mask)
    meta = {
        "gaiascanlaw_data_dir": str(gaiascanlaw_data_dir) if gaiascanlaw_data_dir else None,
        "gaia_density_source": gaia_density_source,
        "gaia_density_present": bool(np.all(np.isfinite(q["gaia_dr3_density_hpx4_log1p"]))),
    }
    return q, meta


def template_groups(q: dict[str, np.ndarray]) -> dict[str, list[tuple[str, np.ndarray]]]:
    sfd = [
        ("ebv_sfd", q["ebv_sfd"]),
        ("sfd_log1p", q["sfd_log1p"]),
        ("sfd_sq", q["sfd_sq"]),
    ]
    selection = [
        ("selection_T", q["selection_T"]),
        ("random_density_log1p", q["random_density_log1p"]),
    ]
    scan = [
        ("gaia_scan_count_log1p_dr3", q["gaia_scan_count_log1p_dr3"]),
        ("gaia_scan_angle_cos2_mean_dr3", q["gaia_scan_angle_cos2_mean_dr3"]),
        ("gaia_scan_angle_sin2_mean_dr3", q["gaia_scan_angle_sin2_mean_dr3"]),
        ("gaia_scan_angle_resultant_dr3", q["gaia_scan_angle_resultant_dr3"]),
    ]
    wise_depth = [
        ("wise_w1_medcov_log1p", q["wise_w1_medcov_log1p"]),
        ("wise_w2_medcov_log1p", q["wise_w2_medcov_log1p"]),
        ("wise_w1_mincov", q["wise_w1_mincov"]),
        ("wise_w2_mincov", q["wise_w2_mincov"]),
        ("wise_w1_lowcov_log1p", q["wise_w1_lowcov_log1p"]),
        ("wise_w2_lowcov_log1p", q["wise_w2_lowcov_log1p"]),
        ("wise_w1_w2_medcov_diff", q["wise_w1_w2_medcov_diff"]),
        ("wise_w1_w2_medcov_ratio", q["wise_w1_w2_medcov_ratio"]),
    ]
    stellar = [("gaia_dr3_density_hpx4_log1p", q["gaia_dr3_density_hpx4_log1p"])]
    groups = {
        "sfd_dust": sfd,
        "quaia_selection_depth": selection,
        "gaia_scanlaw": scan,
        "wise_depth": wise_depth,
        "stellar_density": stellar,
        "dust_selection_scanlaw": sfd + selection + scan,
        "wise_depth_stellar_density": wise_depth + stellar,
        "all_external_templates": sfd + selection + scan + wise_depth + stellar,
    }
    if not np.all(np.isfinite(q["gaia_dr3_density_hpx4_log1p"])):
        groups.pop("stellar_density", None)
        groups.pop("wise_depth_stellar_density", None)
        groups["all_external_templates"] = sfd + selection + scan + wise_depth
    return groups


def run(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    q, meta = load_colour_data(args)
    groups = template_groups(q)
    responses = {
        "gaia_bp_rp_colour": q["bp_rp_z"],
        "wise_w1_w2_colour": q["w1_w2_z"],
        "gaia_wise_colour_cross": q["bp_rp_x_w1_w2_z"],
    }
    primary_mask = (
        (q["z"] >= args.z_min)
        & (q["z"] < args.z_max)
        & (np.abs(q["b"]) >= args.primary_b_cut)
        & finite_mask(q["bp_rp_z"], q["w1_w2_z"], q["bp_rp_x_w1_w2_z"])
    )

    fit_rows: list[dict[str, Any]] = []
    for response_name, response_values in responses.items():
        for b_cut in args.b_cuts:
            mask = (
                (q["z"] >= args.z_min)
                & (q["z"] < args.z_max)
                & (np.abs(q["b"]) >= b_cut)
                & finite_mask(q["z"], response_values, q["vec"][:, 0], q["vec"][:, 1], q["vec"][:, 2])
            )
            z_fit = EXT.fit_dipole(q["z"][mask], q["vec"][mask])
            z_axis = unit_from_lb(z_fit["l_deg"], z_fit["b_deg"])
            hemisphere = (q["vec"][mask] @ z_axis) >= 0.0
            response = response_values[mask]
            raw = fit_response(response, q["vec"][mask])
            raw_axis = unit_from_lb(raw["l_deg"], raw["b_deg"])
            raw_smd = smd(response, hemisphere)
            for group_name, controls_named in groups.items():
                control_mask = np.ones(int(np.sum(mask)), dtype=bool)
                controls_full = []
                for _, values in controls_named:
                    sub = values[mask]
                    control_mask &= np.isfinite(sub)
                    controls_full.append(sub)
                if not np.all(control_mask):
                    response_use = response[control_mask]
                    vec_use = q["vec"][mask][control_mask]
                    hemisphere_use = hemisphere[control_mask]
                    controls = [control[control_mask] for control in controls_full]
                else:
                    response_use = response
                    vec_use = q["vec"][mask]
                    hemisphere_use = hemisphere
                    controls = controls_full
                raw_use = fit_response(response_use, vec_use)
                raw_axis_use = unit_from_lb(raw_use["l_deg"], raw_use["b_deg"])
                raw_smd_use = smd(response_use, hemisphere_use)
                resid_info = residualize_response(response_use, controls)
                resid = resid_info["resid"]
                resid_fit = fit_response(resid, vec_use)
                resid_axis = unit_from_lb(resid_fit["l_deg"], resid_fit["b_deg"])
                fit_rows.append(
                    {
                        "response_field": response_name,
                        "template_group": group_name,
                        "b_cut_deg": b_cut,
                        "N": int(len(response_use)),
                        "n_templates": len(controls_named),
                        "template_names": ";".join(name for name, _ in controls_named),
                        "z_dipole_amp": z_fit["amp"],
                        "z_dipole_l_deg": z_fit["l_deg"],
                        "z_dipole_b_deg": z_fit["b_deg"],
                        "raw_colour_amp": raw_use["amp"],
                        "raw_colour_l_deg": raw_use["l_deg"],
                        "raw_colour_b_deg": raw_use["b_deg"],
                        "raw_colour_snr": raw_use["formal_snr"],
                        "raw_colour_smd_on_z_axis": raw_smd_use,
                        "raw_colour_sep_from_z_dipole_deg": angular_sep_deg(raw_axis_use, z_axis),
                        "template_r2": resid_info["r2_templates"],
                        "resid_colour_amp": resid_fit["amp"],
                        "resid_colour_l_deg": resid_fit["l_deg"],
                        "resid_colour_b_deg": resid_fit["b_deg"],
                        "resid_colour_snr": resid_fit["formal_snr"],
                        "resid_colour_smd_on_z_axis": smd(resid, hemisphere_use),
                        "resid_colour_sep_from_z_dipole_deg": angular_sep_deg(resid_axis, z_axis),
                        "resid_direction_sep_from_raw_deg": angular_sep_deg(resid_axis, raw_axis_use),
                        "resid_amp_ratio_vs_raw": resid_fit["amp"] / raw_use["amp"] if raw_use["amp"] > 0.0 else float("nan"),
                    }
                )

    summary_rows: list[dict[str, Any]] = []
    for response_name in responses:
        for group_name in groups:
            rows = [row for row in fit_rows if row["response_field"] == response_name and row["template_group"] == group_name]
            if rows:
                summary_rows.append(summarise_fit_rows(rows))
    summary_rows.sort(key=lambda row: (row["response_field"], -row["mechanism_score"]))
    for response_name in responses:
        rank = 1
        for row in summary_rows:
            if row["response_field"] == response_name:
                row["rank_within_response"] = rank
                rank += 1

    all_template_names: list[tuple[str, str, np.ndarray]] = []
    for group_name, controls_named in groups.items():
        for name, values in controls_named:
            if not any(existing[1] == name for existing in all_template_names):
                all_template_names.append((group_name, name, values))
    corr_rows: list[dict[str, Any]] = []
    for response_name, response_values in responses.items():
        for first_group, template_name, template_values in all_template_names:
            mask = primary_mask & finite_mask(response_values, template_values)
            corr = pearson(response_values[mask], template_values[mask])
            zcorr = pearson(q["z"][mask], template_values[mask])
            corr_rows.append(
                {
                    "response_field": response_name,
                    "template": template_name,
                    "first_template_group": first_group,
                    "N": int(np.sum(mask)),
                    "corr_response_template": corr,
                    "abs_corr_response_template": abs(corr) if math.isfinite(corr) else float("nan"),
                    "corr_z_template": zcorr,
                }
            )
    corr_rows.sort(key=lambda row: (row["response_field"], -row["abs_corr_response_template"]))

    readout_rows: list[dict[str, Any]] = []
    for response_name in responses:
        ranked = [row for row in summary_rows if row["response_field"] == response_name]
        if ranked:
            readout_rows.append(ranked[0])

    config = {
        "date": "2026-07-15",
        "analysis": "su2_quaia_colour_template_mechanism",
        "z_min": args.z_min,
        "z_max": args.z_max,
        "b_cuts": args.b_cuts,
        "primary_b_cut": args.primary_b_cut,
        "responses": list(responses.keys()),
        "template_groups": {name: [item[0] for item in values] for name, values in groups.items()},
        "mechanism_score": "0.45*dipole_absorption + 0.35*hemisphere-SMD absorption + 0.20*template R2, each clipped to [0,1] where applicable.",
        "external_inputs": {
            "wise_inventory": str(args.wise_inventory),
            "gaia_density": str(args.gaia_density),
            "gaia_density_source": meta["gaia_density_source"],
            "gaia_density_present": meta["gaia_density_present"],
            "gaiascanlaw_data_dir": meta["gaiascanlaw_data_dir"],
        },
    }
    return fit_rows, summary_rows, corr_rows, readout_rows, config


def write_reports(out_dir: Path, summary_rows: list[dict[str, Any]], corr_rows: list[dict[str, Any]], readout_rows: list[dict[str, Any]], config: dict[str, Any]) -> None:
    summary_cols = [
        "rank_within_response",
        "response_field",
        "template_group",
        "mechanism_score",
        "template_r2_mean",
        "resid_amp_ratio_vs_raw",
        "dipole_absorption_fraction",
        "smd_absorption_fraction",
        "max_resid_direction_sep_from_raw_deg",
    ]
    corr_cols = [
        "response_field",
        "template",
        "corr_response_template",
        "corr_z_template",
        "N",
    ]
    top_corr: list[dict[str, Any]] = []
    for response in sorted(set(row["response_field"] for row in corr_rows)):
        top_corr.extend([row for row in corr_rows if row["response_field"] == response][:6])

    bottom_lines = []
    for row in readout_rows:
        bottom_lines.append(
            f"- `{row['response_field']}` is best carried by `{row['template_group']}`: "
            f"score `{fmt(row['mechanism_score'])}`, R2 `{fmt(row['template_r2_mean'])}`, "
            f"residual amp/raw `{fmt(row['resid_amp_ratio_vs_raw'])}`, "
            f"SMD absorption `{fmt(row['smd_absorption_fraction'])}`."
        )
    report = [
        "# SU2 / Quaia Colour-Template Mechanism Gate",
        "",
        "Date: July 15, 2026",
        "",
        "## Purpose",
        "",
        "This mechanism gate tests whether the colour fields that suppressed the locked Quaia angular mode are spatially carried by external survey templates. The response fields are standardised Gaia `BP-RP`, standardised WISE `W1-W2`, and their standardised cross-product.",
        "",
        "For each latitude cut, the response field is regressed against a template group. The residual colour-field dipole and the residual colour imbalance across the locked redshift-dipole hemispheres are then measured.",
        "",
        "## Model",
        "",
        "```text",
        "c_i       = alpha + T_i beta + eta_i",
        "eta_i     = c_i - alpha - T_i beta",
        "c_i       = c0 + d_c . n_i + epsilon_i",
        "eta_i     = e0 + d_eta . n_i + epsilon_i",
        "SMD_c     = (mean(c | H=1) - mean(c | H=0)) / sigma_pooled",
        "```",
        "",
        "`H` is the hemisphere defined by the redshift dipole axis for that latitude cut.",
        "",
        "## Bottom Line",
        "",
        *bottom_lines,
        "",
        "## Top Mechanism Families",
        "",
        markdown_table(readout_rows, summary_cols),
        "",
        "## Full Family Ranking",
        "",
        markdown_table(summary_rows, summary_cols),
        "",
        "## Strongest Individual Template Correlations",
        "",
        markdown_table(top_corr, corr_cols),
        "",
        "## Configuration",
        "",
        "```json",
        json.dumps(config, indent=2),
        "```",
    ]
    (out_dir / "su2_quaia_colour_template_mechanism_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    readout = [
        "# SU2 / Quaia Colour-Template Mechanism Gate Readout",
        "",
        "Date: July 15, 2026",
        "",
        "## Bottom Line",
        "",
        *bottom_lines,
        "",
        "## Top Mechanism Families",
        "",
        markdown_table(readout_rows, summary_cols),
        "",
        "## Strongest Individual Template Correlations",
        "",
        markdown_table(top_corr, corr_cols),
    ]
    (out_dir / "su2_quaia_colour_template_mechanism_readout.md").write_text("\n".join(readout) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Identify external survey-template mechanisms behind SU2/Quaia colour suppression.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--z-min", type=float, default=0.95)
    parser.add_argument("--z-max", type=float, default=1.45)
    parser.add_argument("--b-cuts", type=float, nargs="+", default=[15.0, 20.0, 25.0, 30.0, 35.0])
    parser.add_argument("--primary-b-cut", type=float, default=15.0)
    parser.add_argument("--wise-inventory", type=Path, default=WISE.DEFAULT_WISE)
    parser.add_argument("--gaia-density", type=Path, default=WISE.DEFAULT_GAIA_DENSITY)
    parser.add_argument("--gaia-density-hips-base-url", default=WISE.DEFAULT_GAIA_HIPS_BASE)
    parser.add_argument("--gaia-density-hips-cache-dir", type=Path, default=WISE.DEFAULT_GAIA_HIPS_CACHE)
    parser.add_argument("--gaia-density-hips-order", type=int, default=2)
    parser.add_argument("--gaiascanlaw-data-dir", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    fit_rows, summary_rows, corr_rows, readout_rows, config = run(args)
    write_csv(args.out_dir / "su2_quaia_colour_template_mechanism_fit_rows.csv", fit_rows)
    write_csv(args.out_dir / "su2_quaia_colour_template_mechanism_summary.csv", summary_rows)
    write_csv(args.out_dir / "su2_quaia_colour_template_mechanism_template_correlations.csv", corr_rows)
    write_csv(args.out_dir / "su2_quaia_colour_template_mechanism_top_mechanisms.csv", readout_rows)
    (args.out_dir / "su2_quaia_colour_template_mechanism_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    write_reports(args.out_dir, summary_rows, corr_rows, readout_rows, config)
    print(f"Saved readout: {args.out_dir / 'su2_quaia_colour_template_mechanism_readout.md'}", flush=True)
    print(f"Saved report: {args.out_dir / 'su2_quaia_colour_template_mechanism_report.md'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
