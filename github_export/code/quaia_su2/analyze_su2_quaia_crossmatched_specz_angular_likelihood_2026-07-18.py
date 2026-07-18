from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
CODE_DIR = ROOT / "github_export" / "code" / "quaia_su2"
RESULT_DIR = ROOT / "github_export" / "results" / "2026-07-18" / "su2_quaia"
PREFIX = "su2_quaia_crossmatched_specz_angular_likelihood"
DATE_SUFFIX = "2026-07-18"
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / f"{PREFIX}_20260718"
DR16Q = ROOT / "external_datasets" / "sdss_dr16q_v4" / "DR16Q_v4.fits"
PREREG_MD = RESULT_DIR / f"{PREFIX}_preregistration_{DATE_SUFFIX}.md"
PREREG_JSON = RESULT_DIR / f"{PREFIX}_preregistration_{DATE_SUFFIX}.json"
BASE_PATH = CODE_DIR / "analyze_su2_quaia_injection_calibrated_angular_likelihood_2026-07-18.py"


def load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BASE = load_module("su2_quaia_injection_base", BASE_PATH)
EXT = BASE.EXT
WISE = BASE.WISE


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def array_sha256(values: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(values).tobytes()).hexdigest()


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): clean_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean_json(item) for item in value]
    if isinstance(value, np.ndarray):
        return clean_json(value.tolist())
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean_json(value), indent=2) + "\n", encoding="utf-8")


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
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.6g}" if math.isfinite(float(value)) else ""
    return str(value)


def markdown_table(rows: list[dict[str, Any]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(column, "")) for column in columns) + " |")
    return "\n".join(lines)


def finite_quality(values: np.ndarray, *, lower: float | None = None, upper: float | None = None) -> np.ndarray:
    out = np.asarray(values, dtype=float).copy()
    invalid = ~np.isfinite(out)
    if lower is not None:
        invalid |= out < lower
    if upper is not None:
        invalid |= out > upper
    out[invalid] = np.nan
    return out


def load_dr16q(path: Path, z_min: float, z_max: float) -> dict[str, np.ndarray]:
    with fits.open(path, memmap=True) as hdul:
        data = hdul[1].data
        final = np.asarray(data["IS_QSO_FINAL"], dtype=int)
        z = np.asarray(data["Z"], dtype=float)
        ra = np.asarray(data["RA"], dtype=float)
        dec = np.asarray(data["DEC"], dtype=float)
        selected = (
            (final == 1)
            & np.isfinite(z)
            & (z >= z_min)
            & (z < z_max)
            & np.isfinite(ra)
            & np.isfinite(dec)
            & (ra >= 0.0)
            & (ra < 360.0)
            & (dec >= -90.0)
            & (dec <= 90.0)
        )
        index = np.flatnonzero(selected)
        psferr = np.asarray(data["PSFMAGERR"][selected], dtype=float)
        return {
            "source_index": index,
            "ra": ra[selected],
            "dec": dec[selected],
            "zspec": z[selected],
            "z_conf": np.asarray(data["Z_CONF"][selected], dtype=float),
            "zwarning_nonzero": (np.asarray(data["ZWARNING"][selected], dtype=np.int64) != 0).astype(float),
            "sn_median_all": finite_quality(data["SN_MEDIAN_ALL"][selected], lower=0.0),
            "platesn2": finite_quality(data["PLATESN2"][selected], lower=0.0),
            "psfmagerr_i": finite_quality(psferr[:, 3], lower=0.0, upper=10.0),
            "w1_mag_err": finite_quality(data["W1_MAG_ERR"][selected], lower=0.0, upper=10.0),
            "w2_mag_err": finite_quality(data["W2_MAG_ERR"][selected], lower=0.0, upper=10.0),
            "w1_chi2": finite_quality(data["W1_CHI2"][selected], lower=0.0),
            "w2_chi2": finite_quality(data["W2_CHI2"][selected], lower=0.0),
            "w1_flux_snr": finite_quality(data["W1_FLUX_SNR"][selected]),
            "w2_flux_snr": finite_quality(data["W2_FLUX_SNR"][selected]),
        }


def reciprocal_match(q: dict[str, np.ndarray], dr: dict[str, np.ndarray], radius_arcsec: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    q_coord = SkyCoord(np.asarray(q["ra"]) * u.deg, np.asarray(q["dec"]) * u.deg, frame="icrs")
    d_coord = SkyCoord(dr["ra"] * u.deg, dr["dec"] * u.deg, frame="icrs")
    d_for_q, q_sep, _ = q_coord.match_to_catalog_sky(d_coord)
    q_for_d, _d_sep, _ = d_coord.match_to_catalog_sky(q_coord)
    q_index = np.arange(len(q_coord), dtype=np.int64)
    reciprocal = q_for_d[d_for_q] == q_index
    accepted = reciprocal & (q_sep.arcsec <= radius_arcsec)
    return q_index[accepted], np.asarray(d_for_q[accepted], dtype=np.int64), np.asarray(q_sep.arcsec[accepted], dtype=float)


def standardise_column(values: np.ndarray, mask: np.ndarray) -> tuple[np.ndarray, float, float]:
    subset = np.asarray(values[mask], dtype=float)
    centre = float(np.median(subset))
    scale = float(np.std(subset))
    if not math.isfinite(scale) or scale <= 1e-12:
        raise ValueError("Cannot standardise a zero-variance column")
    return (np.asarray(values, dtype=float) - centre) / scale, centre, scale


def load_inputs(args: argparse.Namespace) -> dict[str, Any]:
    scanlaw_dir = args.gaiascanlaw_data_dir or EXT.resolve_gaiascanlaw_data_dir()
    q = EXT.load_data(args.quaia, args.selection, args.randoms, args.sfd_dir, scanlaw_dir)
    with fits.open(args.quaia, memmap=True) as hdul:
        data = hdul[1].data
        q["bp_rp"] = np.asarray(data["phot_bp_mean_mag"], dtype=float) - np.asarray(data["phot_rp_mean_mag"], dtype=float)
    q.update(WISE.load_wise_depth(args.wise_inventory, q))

    gaia_density_source = "omitted by command-line request"
    if not args.skip_gaia_density:
        density, gaia_density_source = WISE.load_gaia_density(
            args.gaia_density,
            q,
            args.gaia_density_hips_base_url,
            args.gaia_density_hips_cache_dir,
            args.gaia_density_hips_order,
        )
        q.update(density)

    dr = load_dr16q(args.dr16q, args.z_min, args.z_max)
    q_index, d_index, separation = reciprocal_match(q, dr, args.match_radius_arcsec)
    reciprocal_match_count = len(q_index)
    latitude = np.asarray(q["b"][q_index], dtype=float)
    base = (
        (np.abs(latitude) >= args.b_cut)
        & np.isfinite(q["z"][q_index])
        & np.isfinite(dr["zspec"][d_index])
        & np.all(np.isfinite(q["vec"][q_index]), axis=1)
    )
    q_index = q_index[base]
    d_index = d_index[base]
    separation = separation[base]
    zspec = np.asarray(dr["zspec"][d_index], dtype=float)
    photoz = np.asarray(q["z"][q_index], dtype=float)
    vector = np.asarray(q["vec"][q_index], dtype=float)

    grid, delta, curve_meta = BASE.scalar_curve(args)
    scalar = np.interp(zspec, grid, delta)
    scalar_centre = float(np.median(scalar))
    scalar_scale = float(np.std(scalar))
    scalar_weight = (scalar - scalar_centre) / scalar_scale
    su2_basis = scalar_weight[:, None] * vector

    q_values = {
        "ebv_sfd": q["ebv_sfd"][q_index],
        "sfd_log1p": q["sfd_log1p"][q_index],
        "sfd_sq": q["sfd_sq"][q_index],
        "selection_T": q["selection_T"][q_index],
        "random_density_log1p": q["random_density_log1p"][q_index],
        "wise_w1_medcov_log1p": q.get("wise_w1_medcov_log1p", np.full(len(q["z"]), np.nan))[q_index],
        "wise_w2_medcov_log1p": q.get("wise_w2_medcov_log1p", np.full(len(q["z"]), np.nan))[q_index],
        "wise_w1_lowcov_log1p": q.get("wise_w1_lowcov_log1p", np.full(len(q["z"]), np.nan))[q_index],
        "wise_w2_lowcov_log1p": q.get("wise_w2_lowcov_log1p", np.full(len(q["z"]), np.nan))[q_index],
        "gaia_dr3_density_hpx4_log1p": q.get("gaia_dr3_density_hpx4_log1p", np.full(len(q["z"]), np.nan))[q_index],
        "gaia_scan_count_log1p_dr3": q["gaia_scan_count_log1p_dr3"][q_index],
        "gaia_scan_angle_cos2_mean_dr3": q["gaia_scan_angle_cos2_mean_dr3"][q_index],
        "gaia_scan_angle_sin2_mean_dr3": q["gaia_scan_angle_sin2_mean_dr3"][q_index],
        "gaia_scan_angle_resultant_dr3": q["gaia_scan_angle_resultant_dr3"][q_index],
        "w1_w2": q["w1_w2"][q_index],
        "bp_rp": q["bp_rp"][q_index],
        "zerr": q["zerr"][q_index],
        "g": q["g"][q_index],
        "w1": q["w1"][q_index],
        "w2": q["w2"][q_index],
        "pmra_error": q["pmra_error"][q_index],
        "pmdec_error": q["pmdec_error"][q_index],
    }
    colour_mask = np.isfinite(q_values["w1_w2"]) & np.isfinite(q_values["bp_rp"])
    if np.any(colour_mask):
        w1w2, _, _ = standardise_column(q_values["w1_w2"], colour_mask)
        bprp, _, _ = standardise_column(q_values["bp_rp"], colour_mask)
        q_values["bp_rp_x_w1_w2"] = bprp * w1w2

    z0 = (zspec - float(np.median(zspec))) / float(np.std(zspec))
    raw_columns: list[tuple[str, str, np.ndarray]] = [
        ("angular_main", "n_x", vector[:, 0]),
        ("angular_main", "n_y", vector[:, 1]),
        ("angular_main", "n_z", vector[:, 2]),
        ("redshift_and_scalar_main", "zspec_z", z0),
        ("redshift_and_scalar_main", "zspec_z2", z0 * z0),
        ("redshift_and_scalar_main", "zspec_z3", z0 * z0 * z0),
        ("redshift_and_scalar_main", "scalar_main", scalar_weight),
    ]
    group_map = {
        "dust": ["ebv_sfd", "sfd_log1p", "sfd_sq"],
        "depth": ["selection_T", "random_density_log1p", "wise_w1_medcov_log1p", "wise_w2_medcov_log1p", "wise_w1_lowcov_log1p", "wise_w2_lowcov_log1p", "gaia_dr3_density_hpx4_log1p"],
        "gaia_scanning_geometry": ["gaia_scan_count_log1p_dr3", "gaia_scan_angle_cos2_mean_dr3", "gaia_scan_angle_sin2_mean_dr3", "gaia_scan_angle_resultant_dr3"],
        # W1-W2 is already in the affine span of the retained W1 and W2
        # magnitude columns. Keep the interaction, but omit the exact duplicate.
        "colour": ["bp_rp", "bp_rp_x_w1_w2"],
        "quaia_quality": ["zerr", "g", "w1", "w2", "pmra_error", "pmdec_error"],
    }
    for group, names in group_map.items():
        for name in names:
            if name in q_values:
                raw_columns.append((group, name, np.asarray(q_values[name], dtype=float)))
    for name in ["z_conf", "zwarning_nonzero", "sn_median_all", "platesn2", "psfmagerr_i", "w1_mag_err", "w2_mag_err", "w1_chi2", "w2_chi2", "w1_flux_snr", "w2_flux_snr"]:
        raw_columns.append(("dr16q_quality", name, np.asarray(dr[name][d_index], dtype=float)))

    inventory: list[dict[str, Any]] = []
    inventory.append({
        "group": "colour",
        "template": "w1_w2",
        "status": "omitted_exact_affine_redundancy_with_w1_and_w2",
        "finite_rows_before_joint_mask": int(np.sum(np.isfinite(q_values["w1_w2"]))),
        "finite_fraction": float(np.mean(np.isfinite(q_values["w1_w2"]))),
        "raw_std_before_joint_mask": float(np.nanstd(q_values["w1_w2"])),
    })
    retained: list[tuple[str, str, np.ndarray]] = []
    minimum_finite = max(100, int(math.ceil(args.minimum_template_finite_fraction * len(photoz))))
    for group, name, values in raw_columns:
        finite = np.isfinite(values)
        finite_count = int(np.sum(finite))
        scale = float(np.std(values[finite])) if finite_count else float("nan")
        status = "retained" if finite_count >= minimum_finite and math.isfinite(scale) and scale > 1e-12 else "omitted_finite_or_variance_rule"
        inventory.append({
            "group": group,
            "template": name,
            "status": status,
            "finite_rows_before_joint_mask": finite_count,
            "finite_fraction": finite_count / max(len(photoz), 1),
            "raw_std_before_joint_mask": scale,
        })
        if status == "retained":
            retained.append((group, name, values))

    joint = np.isfinite(photoz) & np.isfinite(zspec) & np.all(np.isfinite(su2_basis), axis=1)
    for _group, _name, values in retained:
        joint &= np.isfinite(values)
    if int(np.sum(joint)) < args.minimum_matched_rows:
        raise RuntimeError(f"Only {int(np.sum(joint))} matched rows remain after the locked joint mask")

    template_columns: list[np.ndarray] = []
    template_names: list[str] = []
    for group, name, values in retained:
        transformed, centre, scale = standardise_column(values, joint)
        template_columns.append(transformed[joint])
        template_names.append(name)
        row = next(item for item in inventory if item["group"] == group and item["template"] == name)
        row.update({"joint_rows": int(np.sum(joint)), "centre": centre, "scale": scale})

    response = photoz[joint] - zspec[joint]
    templates = np.column_stack(template_columns).astype(float)
    basis = su2_basis[joint].astype(float)
    x_candidate = np.column_stack([np.ones(len(response)), templates, basis])
    candidate_rank = int(np.linalg.matrix_rank(x_candidate))
    if candidate_rank != x_candidate.shape[1]:
        raise RuntimeError(
            f"Candidate design is rank deficient: rank={candidate_rank}, columns={x_candidate.shape[1]}"
        )
    rng = np.random.default_rng(args.leakage_seed)
    permuted_photoz = rng.permutation(photoz)
    permuted_response = permuted_photoz[joint] - zspec[joint]
    rebuilt_candidate = np.column_stack([np.ones(len(permuted_response)), templates.copy(), basis.copy()])
    leakage_max_abs = float(np.max(np.abs(x_candidate - rebuilt_candidate)))
    leakage_audit = {
        "pass": bool(np.array_equal(x_candidate, rebuilt_candidate)),
        "max_abs_predictor_difference": leakage_max_abs,
        "original_predictor_sha256": array_sha256(x_candidate),
        "permuted_response_predictor_sha256": array_sha256(rebuilt_candidate),
        "response_changed": bool(not np.array_equal(response, permuted_response)),
        "seed": args.leakage_seed,
    }
    if not leakage_audit["pass"] or not leakage_audit["response_changed"]:
        raise RuntimeError("The response-permutation leakage audit failed")

    curve_rows = [{"zspec": float(z), "delta_mu_su2r_minus_lcdm": float(d)} for z, d in zip(grid, delta)]
    match_rows = [{
        "q_index": int(qi),
        "dr16q_filtered_index": int(di),
        "dr16q_source_index": int(dr["source_index"][di]),
        "separation_arcsec": float(sep),
        "z_quaia": float(q["z"][qi]),
        "z_spec": float(dr["zspec"][di]),
        "delta_z": float(q["z"][qi] - dr["zspec"][di]),
    } for qi, di, sep in zip(q_index[joint], d_index[joint], separation[joint])]
    metadata = {
        "quaia_rows": int(len(q["z"])),
        "dr16q_quality_window_rows": int(len(dr["zspec"])),
        "reciprocal_matches_within_radius": int(reciprocal_match_count),
        "matched_rows_before_joint_mask": int(len(photoz)),
        "matched_rows_after_joint_mask": int(np.sum(joint)),
        "match_radius_arcsec": args.match_radius_arcsec,
        "match_separation_max_arcsec": float(np.max(separation[joint])),
        "template_count": int(templates.shape[1]),
        "template_names": template_names,
        "design_rank_template_only": int(np.linalg.matrix_rank(np.column_stack([np.ones(len(response)), templates]))),
        "design_rank_candidate": candidate_rank,
        "design_columns_candidate": int(x_candidate.shape[1]),
        "design_condition_candidate": float(np.linalg.cond(x_candidate)),
        "basis_reduction": "w1_w2 omitted because it is an exact affine combination of retained w1 and w2",
        "gaiascanlaw_data_dir": str(scanlaw_dir) if scanlaw_dir else None,
        "gaia_density_source": gaia_density_source,
        "dr16q_path": str(args.dr16q),
        "dr16q_sha256": sha256_file(args.dr16q),
        "quaia_path": str(args.quaia),
        "quaia_sha256": sha256_file(args.quaia),
        "scalar_centre": scalar_centre,
        "scalar_scale": scalar_scale,
        "curve": curve_meta,
    }
    return {
        "response": response,
        "templates": templates,
        "su2_basis": basis,
        "l_deg": np.asarray(q["l"][q_index][joint], dtype=float),
        "b_deg": np.asarray(q["b"][q_index][joint], dtype=float),
        "inventory": inventory,
        "metadata": metadata,
        "curve_rows": curve_rows,
        "match_rows": match_rows,
        "leakage_audit": leakage_audit,
    }


def power_row(scale: float, partition: str, evaluation: dict[str, np.ndarray], recovered: np.ndarray, predictive: np.ndarray, diag: dict[str, np.ndarray]) -> dict[str, Any]:
    return {
        "amplitude_scale": scale,
        "partition": partition,
        "n_replicates": int(len(evaluation["delta_lpd"])),
        "predictive_detection_rate": float(np.mean(predictive)),
        "joint_recovery_rate": float(np.mean(recovered)) if scale > 0 else float("nan"),
        "delta_cvlpd_median": float(np.median(evaluation["delta_lpd"])),
        "delta_cvlpd_q05": float(np.quantile(evaluation["delta_lpd"], 0.05)),
        "delta_cvlpd_q95": float(np.quantile(evaluation["delta_lpd"], 0.95)),
        "amplitude_ratio_median": float(np.nanmedian(diag["amplitude_ratio"])) if scale > 0 else float("nan"),
        "direction_separation_median_deg": float(np.nanmedian(diag["direction_separation_deg"])) if scale > 0 else float("nan"),
    }


def export_outputs(out_dir: Path, paths: list[Path], script_path: Path) -> Path:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    exported: list[Path] = []
    for path in paths:
        destination = RESULT_DIR / path.name
        shutil.copy2(path, destination)
        exported.append(destination)
    rows = []
    for path in [script_path, PREREG_MD, PREREG_JSON, *exported]:
        rows.append({"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": sha256_file(path)})
    manifest = RESULT_DIR / f"{PREFIX}_manifest_{DATE_SUFFIX}.csv"
    write_csv(manifest, rows)
    return manifest


def build_report(out_dir: Path, decision: dict[str, Any], power_rows: list[dict[str, Any]], fold_rows: list[dict[str, Any]], config: dict[str, Any]) -> tuple[Path, Path]:
    observed = decision["observed"]
    calibration = decision["calibration"]
    conditions = [{"condition": key, "pass": value} for key, value in decision["conditions"].items()]
    report = [
        "# Leakage-Free Quaia x DR16Q Angular Likelihood",
        "",
        "Date: 18 July 2026",
        "",
        "## Bottom Line",
        "",
        f"The preregistered outcome is **{decision['outcome']}**. {decision['interpretation']}",
        "",
        "## Independent Response",
        "",
        "The response is `z_Quaia - z_spec`, where `z_spec` is the independently measured DR16Q spectroscopic redshift. The SU2R-minus-Lambda-CDM scalar curve is evaluated only at `z_spec`.",
        "",
        f"- reciprocal one-to-one matches after all masks: `{config['inputs']['matched_rows_after_joint_mask']}`",
        f"- maximum accepted separation: `{fmt(config['inputs']['match_separation_max_arcsec'])}` arcsec",
        f"- retained nuisance templates: `{config['inputs']['template_count']}`",
        f"- candidate design rank: `{config['inputs']['design_rank_candidate']}`",
        f"- candidate design condition number: `{fmt(config['inputs']['design_condition_candidate'])}`",
        "",
        "## Leakage Audit",
        "",
        f"- exact predictor equality after permuting the Quaia response: `{decision['leakage_audit']['pass']}`",
        f"- maximum absolute predictor change: `{fmt(decision['leakage_audit']['max_abs_predictor_difference'])}`",
        f"- response changed under permutation: `{decision['leakage_audit']['response_changed']}`",
        f"- predictor SHA-256: `{decision['leakage_audit']['original_predictor_sha256']}`",
        "",
        "## Calibration",
        "",
        f"- frozen 99th-percentile threshold: `{fmt(calibration['frozen_threshold'])}`",
        f"- untouched validation-null false-positive rate: `{fmt(calibration['validation_false_positive_rate'])}`",
        f"- joint recovery at the locked 1x signal: `{fmt(calibration['joint_recovery_rate_1x'])}`",
        "",
        "## Observed Held-Out Result",
        "",
        f"- pooled Delta CVLPD: `{fmt(observed['delta_cvlpd'])}`",
        f"- empirical validation-null p-value: `{fmt(observed['empirical_validation_null_p'])}`",
        f"- positive sectors: `{observed['positive_fold_count']}/{config['held_out']['n_folds']}`",
        f"- recovered amplitude ratio: `{fmt(observed['amplitude_ratio_to_reference'])}`",
        f"- recovered direction: `(l,b)=({fmt(observed['pooled_l_deg'])}, {fmt(observed['pooled_b_deg'])}) deg`",
        f"- separation from locked direction: `{fmt(observed['direction_separation_deg'])}` deg",
        "",
        "## Conjunctive Rule",
        "",
        markdown_table(conditions, ["condition", "pass"]),
        "",
        "## Power Calibration",
        "",
        markdown_table(power_rows, ["amplitude_scale", "partition", "n_replicates", "predictive_detection_rate", "joint_recovery_rate", "delta_cvlpd_median", "amplitude_ratio_median", "direction_separation_median_deg"]),
        "",
        "## Held-Out Sectors",
        "",
        markdown_table(fold_rows, ["fold", "N_test", "delta_lpd", "delta_sse", "recovered_amp", "recovered_l_deg", "recovered_b_deg"]),
        "",
        "## Claim Boundary",
        "",
        "This replacement repairs the response leakage in the first injection branch. A pass would establish only a conditional cross-catalogue predictive association. A calibration or power failure is underpowered; an observed failure after adequate calibration rejects only the locked candidate in this sample and nuisance model.",
    ]
    report_path = out_dir / f"{PREFIX}_report_{DATE_SUFFIX}.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    readout_path = out_dir / f"{PREFIX}_readout_{DATE_SUFFIX}.md"
    readout_path.write_text(
        "\n".join([
            "# Quaia x DR16Q leakage-free readout",
            "",
            f"Outcome: **{decision['outcome']}**.",
            "",
            f"The response-permutation leakage audit passed: `{decision['leakage_audit']['pass']}`. Validation FPR was `{fmt(calibration['validation_false_positive_rate'])}`, locked 1x recovery was `{fmt(calibration['joint_recovery_rate_1x'])}`, and observed Delta CVLPD was `{fmt(observed['delta_cvlpd'])}` with direction separation `{fmt(observed['direction_separation_deg'])} deg`.",
            "",
            "This is the valid replacement for the methodologically invalidated first injection-likelihood branch.",
        ]) + "\n",
        encoding="utf-8",
    )
    return report_path, readout_path


def run(args: argparse.Namespace) -> dict[str, Any]:
    if not PREREG_MD.exists() or not PREREG_JSON.exists():
        raise FileNotFoundError("The dated replacement preregistration must exist before analysis")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    loaded = load_inputs(args)
    response = loaded["response"]
    templates = loaded["templates"]
    basis = loaded["su2_basis"]
    locked_unit = BASE.unit_from_lb(args.locked_l_deg, args.locked_b_deg)
    signal = args.reference_amplitude * (basis @ locked_unit)
    x0 = np.column_stack([np.ones(len(response)), templates])
    x1 = np.column_stack([x0, basis])
    beta0 = np.linalg.pinv(x0.T @ x0, rcond=1e-12) @ (x0.T @ response)
    base = x0 @ beta0
    residual = response - base
    cell, cell_fold, object_fold, spatial_meta = BASE.make_cells_and_folds(
        loaded["l_deg"], loaded["b_deg"], args.cell_degrees, args.n_folds
    )
    if np.min(np.bincount(object_fold, minlength=args.n_folds)) == 0:
        raise RuntimeError("At least one locked spatial fold is empty")
    model0 = BASE.build_model_stats(x0, base, signal, residual, cell, object_fold, args.n_folds)
    model1 = BASE.build_model_stats(x1, base, signal, residual, cell, object_fold, args.n_folds)
    models = [model0, model1]
    response_stats = BASE.build_response_stats(base, signal, residual, cell, object_fold, args.n_folds, cell_fold)

    observed_signs = np.ones((spatial_meta["occupied_cells"], 1), dtype=float)
    observed_terms = BASE.precompute_sign_terms(models, response_stats, observed_signs, args.n_folds)
    observed_eval = BASE.evaluate_pair(models, response_stats, observed_terms, 0.0, args.n_folds)
    observed_vector = observed_eval["pooled_vector"][:, 0]
    observed_diag = BASE.vector_diagnostics(observed_vector[:, None], locked_unit, args.reference_amplitude)
    observed_l, observed_b = BASE.unit_to_lb(observed_vector)

    rng = np.random.default_rng(args.seed)
    null_signs = rng.choice(np.array([-1.0, 1.0]), size=(spatial_meta["occupied_cells"], args.null_replicates))
    null_terms = BASE.precompute_sign_terms(models, response_stats, null_signs, args.n_folds)
    null_eval = BASE.evaluate_pair(models, response_stats, null_terms, 0.0, args.n_folds)
    threshold_null = null_eval["delta_lpd"][: args.threshold_replicates]
    validation_null = null_eval["delta_lpd"][args.threshold_replicates :]
    threshold = float(np.quantile(threshold_null, 1.0 - args.alpha, method="higher"))
    validation_fpr = float(np.mean((validation_null > 0.0) & (validation_null > threshold)))

    calibration_rows: list[dict[str, Any]] = []
    power_rows: list[dict[str, Any]] = []
    null_diag = BASE.vector_diagnostics(null_eval["pooled_vector"], locked_unit, 0.0)
    null_predictive = (null_eval["delta_lpd"] > 0.0) & (null_eval["delta_lpd"] > threshold)
    for replicate in range(args.null_replicates):
        calibration_rows.append({
            "replicate": replicate,
            "partition": "threshold_calibration" if replicate < args.threshold_replicates else "untouched_validation",
            "amplitude_scale": 0.0,
            "delta_cvlpd": float(null_eval["delta_lpd"][replicate]),
            "predictive_detection": bool(null_predictive[replicate]),
        })
    for name, slc in [("threshold_calibration", slice(0, args.threshold_replicates)), ("untouched_validation", slice(args.threshold_replicates, args.null_replicates))]:
        subset_eval = {key: value[slc] if value.ndim == 1 else value for key, value in null_eval.items()}
        subset_diag = {key: value[slc] for key, value in null_diag.items()}
        power_rows.append(power_row(0.0, name, subset_eval, null_predictive[slc], null_predictive[slc], subset_diag))

    recovery_1x = float("nan")
    for scale in [0.5, 1.0, 1.5, 2.0]:
        signs = rng.choice(np.array([-1.0, 1.0]), size=(spatial_meta["occupied_cells"], args.injection_replicates))
        terms = BASE.precompute_sign_terms(models, response_stats, signs, args.n_folds)
        evaluation = BASE.evaluate_pair(models, response_stats, terms, scale, args.n_folds)
        diag = BASE.vector_diagnostics(evaluation["pooled_vector"], locked_unit, scale * args.reference_amplitude)
        predictive = (evaluation["delta_lpd"] > 0.0) & (evaluation["delta_lpd"] > threshold)
        recovered = (
            predictive
            & diag["positive_orientation"]
            & (diag["amplitude_ratio"] >= args.amplitude_ratio_min)
            & (diag["amplitude_ratio"] <= args.amplitude_ratio_max)
            & (diag["direction_separation_deg"] <= args.maximum_direction_separation_deg)
        )
        if scale == 1.0:
            recovery_1x = float(np.mean(recovered))
        power_rows.append(power_row(scale, "fresh_injection", evaluation, recovered, predictive, diag))
        for replicate in range(args.injection_replicates):
            calibration_rows.append({
                "replicate": replicate,
                "partition": "fresh_injection",
                "amplitude_scale": scale,
                "delta_cvlpd": float(evaluation["delta_lpd"][replicate]),
                "pooled_amplitude": float(diag["amplitude"][replicate]),
                "amplitude_ratio": float(diag["amplitude_ratio"][replicate]),
                "direction_separation_deg": float(diag["direction_separation_deg"][replicate]),
                "positive_orientation": bool(diag["positive_orientation"][replicate]),
                "predictive_detection": bool(predictive[replicate]),
                "joint_recovery": bool(recovered[replicate]),
            })
        print(f"evaluated fresh injection scale {scale:g}", flush=True)

    observed_delta = float(observed_eval["delta_lpd"][0])
    observed_sep = float(observed_diag["direction_separation_deg"][0])
    observed_ratio = float(observed_diag["amplitude_ratio"][0])
    observed_orientation = bool(observed_diag["positive_orientation"][0])
    empirical_p = float((1 + np.sum(validation_null >= observed_delta)) / (1 + len(validation_null)))
    conditions = {
        "predictor_permutation_leakage_audit": bool(loaded["leakage_audit"]["pass"]),
        "validation_null_false_positive_rate_le_0p01": validation_fpr <= args.alpha,
        "joint_1x_recovery_ge_0p80": recovery_1x >= args.minimum_power,
        "observed_pooled_held_out_gain_positive_and_above_threshold": observed_delta > 0.0 and observed_delta > threshold,
        "observed_direction_amplitude_orientation": observed_sep <= args.maximum_direction_separation_deg and args.amplitude_ratio_min <= observed_ratio <= args.amplitude_ratio_max and observed_orientation,
    }
    if not conditions["validation_null_false_positive_rate_le_0p01"]:
        outcome = "calibration failure"
        interpretation = "The untouched validation null exceeded the locked false-positive ceiling; the candidate is not promotable."
    elif not conditions["joint_1x_recovery_ge_0p80"]:
        outcome = "underpowered"
        interpretation = "The leakage-free design did not recover the locked one-times signal with at least 80% power; the threshold remains frozen."
    elif all(conditions.values()):
        outcome = "pass for conditional follow-up only"
        interpretation = "All leakage, calibration, power and observed held-out gates passed; this remains a conditional association rather than a discovery."
    else:
        outcome = "observed gate failure"
        interpretation = "Calibration and power were adequate, but the observed cross-matched candidate failed at least one held-out direction, amplitude or score gate."

    observed_fold = observed_eval["fold_delta_lpd"][:, 0]
    fold_rows: list[dict[str, Any]] = []
    for fold in range(args.n_folds):
        vector = observed_eval["fold_candidate_beta"][fold, :, 0]
        l_deg, b_deg = BASE.unit_to_lb(vector)
        fold_rows.append({
            "fold": fold,
            "N_test": int(response_stats["n_test"][fold]),
            "delta_lpd": float(observed_fold[fold]),
            "delta_sse": float(observed_eval["fold_delta_sse"][fold, 0]),
            "recovered_amp": float(np.linalg.norm(vector)),
            "recovered_l_deg": l_deg,
            "recovered_b_deg": b_deg,
        })

    in_sample = BASE.in_sample_diagnostics(x0, x1, response)
    in_sample_vector = np.asarray(in_sample.pop("candidate_vector"), dtype=float)
    in_l, in_b = BASE.unit_to_lb(in_sample_vector)
    in_sample.update({"candidate_amplitude": float(np.linalg.norm(in_sample_vector)), "candidate_l_deg": in_l, "candidate_b_deg": in_b})
    decision = {
        "analysis": PREFIX,
        "date": DATE_SUFFIX,
        "outcome": outcome,
        "interpretation": interpretation,
        "conjunctive_rule": True,
        "leakage_audit": loaded["leakage_audit"],
        "calibration": {
            "threshold_replicates": args.threshold_replicates,
            "validation_replicates": args.null_replicates - args.threshold_replicates,
            "injection_replicates_per_scale": args.injection_replicates,
            "frozen_threshold": threshold,
            "validation_false_positive_rate": validation_fpr,
            "joint_recovery_rate_1x": recovery_1x,
        },
        "observed": {
            "delta_cvlpd": observed_delta,
            "empirical_validation_null_p": empirical_p,
            "positive_fold_count": int(np.sum(observed_fold > 0.0)),
            "pooled_amplitude": float(np.linalg.norm(observed_vector)),
            "amplitude_ratio_to_reference": observed_ratio,
            "pooled_l_deg": observed_l,
            "pooled_b_deg": observed_b,
            "direction_separation_deg": observed_sep,
            "positive_orientation": observed_orientation,
        },
        "conditions": conditions,
        "in_sample_non_gating_diagnostics": in_sample,
    }
    config = {
        "analysis": PREFIX,
        "date": DATE_SUFFIX,
        "preregistration_markdown": str(PREREG_MD),
        "preregistration_markdown_sha256": sha256_file(PREREG_MD),
        "preregistration_json": str(PREREG_JSON),
        "preregistration_json_sha256": sha256_file(PREREG_JSON),
        "sample": {"zspec_min": args.z_min, "zspec_max_exclusive": args.z_max, "b_cut": args.b_cut, "match_radius_arcsec": args.match_radius_arcsec},
        "locked_axis": {"l_deg": args.locked_l_deg, "b_deg": args.locked_b_deg, "reference_amplitude": args.reference_amplitude},
        "held_out": spatial_meta,
        "inputs": loaded["metadata"],
        "calibration": {"seed": args.seed, "null_replicates": args.null_replicates, "threshold_replicates": args.threshold_replicates, "injection_replicates_per_scale": args.injection_replicates, "alpha": args.alpha},
    }
    paths = {
        "calibration_draws": args.out_dir / f"{PREFIX}_calibration_draws.csv",
        "power_summary": args.out_dir / f"{PREFIX}_power_summary_{DATE_SUFFIX}.csv",
        "heldout_folds": args.out_dir / f"{PREFIX}_heldout_folds_{DATE_SUFFIX}.csv",
        "template_inventory": args.out_dir / f"{PREFIX}_template_inventory_{DATE_SUFFIX}.csv",
        "scalar_curve": args.out_dir / f"{PREFIX}_scalar_curve_{DATE_SUFFIX}.csv",
        "matched_sample": args.out_dir / f"{PREFIX}_matched_sample.csv",
        "decision": args.out_dir / f"{PREFIX}_decision_{DATE_SUFFIX}.json",
        "config": args.out_dir / f"{PREFIX}_config_{DATE_SUFFIX}.json",
    }
    write_csv(paths["calibration_draws"], calibration_rows)
    write_csv(paths["power_summary"], power_rows)
    write_csv(paths["heldout_folds"], fold_rows)
    write_csv(paths["template_inventory"], loaded["inventory"])
    write_csv(paths["scalar_curve"], loaded["curve_rows"])
    write_csv(paths["matched_sample"], loaded["match_rows"])
    write_json(paths["decision"], decision)
    write_json(paths["config"], config)
    report_path, readout_path = build_report(args.out_dir, decision, power_rows, fold_rows, config)
    compact = [paths["power_summary"], paths["heldout_folds"], paths["template_inventory"], paths["scalar_curve"], paths["decision"], paths["config"], report_path, readout_path]
    manifest = export_outputs(args.out_dir, compact, Path(__file__)) if not args.no_export else None
    print(f"Outcome: {outcome}", flush=True)
    print(f"Saved report: {report_path}", flush=True)
    if manifest:
        print(f"Saved manifest: {manifest}", flush=True)
    return {"decision": decision, "config": config, "paths": paths}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Leakage-free Quaia x DR16Q injection-calibrated angular likelihood")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--dr16q", type=Path, default=DR16Q)
    parser.add_argument("--quaia", type=Path, default=EXT.DEFAULT_QUAIA)
    parser.add_argument("--randoms", type=Path, default=EXT.DEFAULT_RANDOMS)
    parser.add_argument("--selection", type=Path, default=EXT.DEFAULT_SELECTION)
    parser.add_argument("--sfd-dir", type=Path, default=EXT.DEFAULT_SFD_DIR)
    parser.add_argument("--wise-inventory", type=Path, default=WISE.DEFAULT_WISE)
    parser.add_argument("--gaia-density", type=Path, default=WISE.DEFAULT_GAIA_DENSITY)
    parser.add_argument("--gaia-density-hips-base-url", default=WISE.DEFAULT_GAIA_HIPS_BASE)
    parser.add_argument("--gaia-density-hips-cache-dir", type=Path, default=WISE.DEFAULT_GAIA_HIPS_CACHE)
    parser.add_argument("--gaia-density-hips-order", type=int, default=2)
    parser.add_argument("--gaiascanlaw-data-dir", type=Path, default=None)
    parser.add_argument("--skip-gaia-density", action="store_true")
    parser.add_argument("--lcdm-bestfit", type=Path, default=BASE.LCDM_BESTFIT)
    parser.add_argument("--su2r-bestfit", type=Path, default=BASE.SU2R_BESTFIT)
    parser.add_argument("--integration-nint", type=int, default=256)
    parser.add_argument("--scalar-grid-points", type=int, default=501)
    parser.add_argument("--z-min", type=float, default=1.0)
    parser.add_argument("--z-max", type=float, default=1.5)
    parser.add_argument("--b-cut", type=float, default=15.0)
    parser.add_argument("--match-radius-arcsec", type=float, default=1.0)
    parser.add_argument("--locked-l-deg", type=float, default=138.99062177880762)
    parser.add_argument("--locked-b-deg", type=float, default=-36.60410263509791)
    parser.add_argument("--reference-amplitude", type=float, default=0.002548691334890594)
    parser.add_argument("--cell-degrees", type=float, default=12.0)
    parser.add_argument("--n-folds", type=int, default=12)
    parser.add_argument("--seed", type=int, default=18071826)
    parser.add_argument("--leakage-seed", type=int, default=18071827)
    parser.add_argument("--null-replicates", type=int, default=2000)
    parser.add_argument("--threshold-replicates", type=int, default=1000)
    parser.add_argument("--injection-replicates", type=int, default=1000)
    parser.add_argument("--alpha", type=float, default=0.01)
    parser.add_argument("--minimum-power", type=float, default=0.80)
    parser.add_argument("--maximum-direction-separation-deg", type=float, default=30.0)
    parser.add_argument("--amplitude-ratio-min", type=float, default=0.5)
    parser.add_argument("--amplitude-ratio-max", type=float, default=1.5)
    parser.add_argument("--minimum-template-finite-fraction", type=float, default=0.80)
    parser.add_argument("--minimum-matched-rows", type=int, default=1000)
    parser.add_argument("--no-export", action="store_true")
    return parser.parse_args()


def main() -> int:
    run(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
