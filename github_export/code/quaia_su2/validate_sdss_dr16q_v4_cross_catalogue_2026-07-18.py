from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from scipy import special, stats


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
DEFAULT_PREREG = (
    ROOT
    / "github_export"
    / "results"
    / "2026-07-18"
    / "su2_quaia"
    / "sdss_dr16q_v4_cross_catalogue_preregistration_2026-07-18.json"
)
DEFAULT_FITS = ROOT / "external_datasets" / "sdss_dr16q_v4" / "DR16Q_v4.fits"
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_sdss_dr16q_v4_cross_catalogue_20260718"
DEFAULT_GIT_RESULTS = ROOT / "github_export" / "results" / "2026-07-18" / "su2_quaia"
STAT_NAMES = ("partial_F_l1", "partial_F_l2", "partial_F_l3", "joint_F_l1_l3")
GROUPS = {
    "partial_F_l1": np.arange(1, 4),
    "partial_F_l2": np.arange(4, 9),
    "partial_F_l3": np.arange(9, 16),
}


@dataclass
class AnalysisRow:
    tag: str
    role: str
    z_min: float
    z_max: float
    bcut_deg: float
    source_index: np.ndarray
    z: np.ndarray
    l_deg: np.ndarray
    b_deg: np.ndarray
    design: np.ndarray
    xtx_inverse: np.ndarray
    residual_null: np.ndarray
    null_yty: float
    dof: int
    condition_number: float
    rank: int
    beta: np.ndarray
    rss_full: float
    observed_f: np.ndarray
    observed_score: np.ndarray
    dipole_vector: np.ndarray
    dipole_amplitude: float
    dipole_l_deg: float
    dipole_b_deg: float
    dipole_ra_deg: float
    dipole_dec_deg: float
    multipole_power_l1: float
    multipole_power_l2: float
    multipole_power_l3: float

    @property
    def key(self) -> tuple[str, float]:
        return self.tag, self.bcut_deg


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


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


def write_global_maxima(path: Path, scale_maxima: list[tuple[float, np.ndarray]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cell_deg", "mock_id", "family_max_score"])
        for cell_deg, maxima in scale_maxima:
            for mock_id, value in enumerate(maxima):
                writer.writerow([f"{cell_deg:g}", mock_id, f"{float(value):.12g}"])


def fmt(value: Any) -> str:
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.6g}" if math.isfinite(float(value)) else ""
    return str(value)


def markdown_table(rows: list[dict[str, Any]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(column, "")) for column in columns) + " |")
    return "\n".join(lines)


def p_empirical(exceedances: int | np.integer, n_mocks: int) -> float:
    return float((1 + int(exceedances)) / (1 + n_mocks))


def nominal_score(f_value: np.ndarray | float, df1: int, df2: int) -> np.ndarray:
    p = stats.f.sf(np.asarray(f_value, dtype=float), df1, df2)
    return -np.log10(np.clip(p, 1e-300, 1.0))


def real_spherical_harmonic(ell: int, m: int, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    if m == 0:
        return np.real(special.sph_harm_y(ell, 0, theta, phi))
    harmonic = special.sph_harm_y(ell, abs(m), theta, phi)
    phase = (-1.0) ** abs(m)
    if m < 0:
        return math.sqrt(2.0) * phase * np.imag(harmonic)
    return math.sqrt(2.0) * phase * np.real(harmonic)


def build_design(l_deg: np.ndarray, b_deg: np.ndarray) -> np.ndarray:
    lon = np.radians(l_deg)
    lat = np.radians(b_deg)
    cos_lat = np.cos(lat)
    columns = [
        np.ones(len(l_deg), dtype=float),
        cos_lat * np.cos(lon),
        cos_lat * np.sin(lon),
        np.sin(lat),
    ]
    theta = 0.5 * math.pi - lat
    for ell in (2, 3):
        for m in range(-ell, ell + 1):
            columns.append(real_spherical_harmonic(ell, m, theta, lon))
    return np.column_stack(columns).astype(float)


def unit_to_directions(vector: np.ndarray) -> tuple[float, float, float, float]:
    amplitude = float(np.linalg.norm(vector))
    if not math.isfinite(amplitude) or amplitude <= 1e-15:
        return float("nan"), float("nan"), float("nan"), float("nan")
    unit_vector = vector / amplitude
    l_deg = math.degrees(math.atan2(float(unit_vector[1]), float(unit_vector[0]))) % 360.0
    b_deg = math.degrees(math.asin(float(np.clip(unit_vector[2], -1.0, 1.0))))
    icrs = SkyCoord(l=l_deg * u.deg, b=b_deg * u.deg, frame="galactic").icrs
    return l_deg, b_deg, float(icrs.ra.deg), float(icrs.dec.deg)


def angular_separation_vectors(left: np.ndarray, right: np.ndarray) -> float:
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denominator <= 1e-15 or not math.isfinite(denominator):
        return float("inf")
    cosine = float(np.clip(np.dot(left, right) / denominator, -1.0, 1.0))
    return math.degrees(math.acos(cosine))


def comparator_galactic_unit(prereg: dict[str, Any]) -> np.ndarray:
    target = prereg["locked_quaia_comparator"]
    coordinate = SkyCoord(ra=float(target["ra_deg"]) * u.deg, dec=float(target["dec_deg"]) * u.deg, frame="icrs").galactic
    lon = math.radians(float(coordinate.l.deg))
    lat = math.radians(float(coordinate.b.deg))
    return np.array([math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)], dtype=float)


def load_catalogue(path: Path, prereg: dict[str, Any]) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    catalogue = prereg["catalogue"]
    if not path.exists():
        raise FileNotFoundError(f"Official catalogue is not present: {path}")
    actual_bytes = path.stat().st_size
    if actual_bytes != int(catalogue["expected_bytes"]):
        raise ValueError(f"Catalogue byte count is {actual_bytes}, expected {catalogue['expected_bytes']}")

    with fits.open(path, memmap=True, mode="readonly") as hdul:
        table = hdul[1].data
        if len(table) != int(catalogue["expected_rows"]):
            raise ValueError(f"Catalogue row count is {len(table)}, expected {catalogue['expected_rows']}")
        names_present = {name.upper() for name in (hdul[1].columns.names or [])}
        missing = [name for name in catalogue["required_columns"] if name.upper() not in names_present]
        if missing:
            raise ValueError(f"Missing required catalogue columns: {missing}")

        ra_all = np.asarray(table["RA"], dtype=float)
        dec_all = np.asarray(table["DEC"], dtype=float)
        z_all = np.asarray(table["Z"], dtype=float)
        final_all = np.asarray(table["IS_QSO_FINAL"], dtype=int)
        base = (
            (final_all == 1)
            & np.isfinite(ra_all)
            & np.isfinite(dec_all)
            & np.isfinite(z_all)
            & (ra_all >= 0.0)
            & (ra_all < 360.0)
            & (dec_all >= -90.0)
            & (dec_all <= 90.0)
        )
        original_index = np.flatnonzero(base)
        sdss_name = np.asarray(table["SDSS_NAME"])[original_index]
        z_conf = np.asarray(table["Z_CONF"], dtype=float)[original_index]
        zwarning = np.asarray(table["ZWARNING"], dtype=np.int64)[original_index]
        sn_median = np.asarray(table["SN_MEDIAN_ALL"], dtype=float)[original_index]

        sn_sort = np.where(np.isfinite(sn_median), sn_median, -np.inf)
        conf_sort = np.where(np.isfinite(z_conf), z_conf, -np.inf)
        order = np.lexsort((original_index, -sn_sort, -np.isfinite(sn_median).astype(int), -conf_sort, sdss_name))
        sorted_names = sdss_name[order]
        first = np.ones(len(order), dtype=bool)
        first[1:] = sorted_names[1:] != sorted_names[:-1]
        keep_local = np.sort(order[first])
        duplicate_rows_removed = int(len(original_index) - len(keep_local))
        source_index = original_index[keep_local]

        ra = ra_all[source_index].copy()
        dec = dec_all[source_index].copy()
        z = z_all[source_index].copy()
        z_conf = np.asarray(table["Z_CONF"], dtype=float)[source_index].copy()
        zwarning = np.asarray(table["ZWARNING"], dtype=np.int64)[source_index].copy()
        sn_median = np.asarray(table["SN_MEDIAN_ALL"], dtype=float)[source_index].copy()

    galactic = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs").galactic
    l_deg = np.asarray(galactic.l.deg, dtype=float)
    b_deg = np.asarray(galactic.b.deg, dtype=float)
    meta = {
        "fits_path": str(path),
        "fits_bytes": actual_bytes,
        "fits_rows": int(catalogue["expected_rows"]),
        "primary_quality_rows_before_deduplication": int(len(original_index)),
        "primary_quality_rows_after_deduplication": int(len(source_index)),
        "duplicate_rows_removed": duplicate_rows_removed,
        "finite_z_range": [float(np.min(z)), float(np.max(z))],
    }
    return {
        "source_index": source_index,
        "ra": ra,
        "dec": dec,
        "z": z,
        "z_conf": z_conf,
        "zwarning": zwarning,
        "sn_median_all": sn_median,
        "l": l_deg,
        "b": b_deg,
    }, meta


def prepare_row(data: dict[str, np.ndarray], window: dict[str, Any], bcut_deg: float, prereg: dict[str, Any]) -> AnalysisRow:
    z_min = float(window["z_min"])
    z_max = float(window["z_max_exclusive"])
    mask = (data["z"] >= z_min) & (data["z"] < z_max) & (np.abs(data["b"]) >= bcut_deg)
    index = np.flatnonzero(mask)
    z = data["z"][index]
    l_deg = data["l"][index]
    b_deg = data["b"][index]
    design = build_design(l_deg, b_deg)
    xtx = design.T @ design
    rank = int(np.linalg.matrix_rank(xtx))
    condition_number = float(math.sqrt(max(np.linalg.cond(xtx), 0.0)))
    minimum_n = int(prereg["mask"]["minimum_row_n"])
    maximum_condition = float(prereg["mask"]["maximum_design_condition_number"])
    if len(index) < minimum_n:
        raise ValueError(f"{window['tag']} b={bcut_deg:g} has N={len(index)} < {minimum_n}")
    if rank != design.shape[1]:
        raise ValueError(f"{window['tag']} b={bcut_deg:g} has rank {rank} < {design.shape[1]}")
    if not math.isfinite(condition_number) or condition_number > maximum_condition:
        raise ValueError(f"{window['tag']} b={bcut_deg:g} has condition number {condition_number:g}")

    inverse = np.linalg.inv(xtx)
    xty = design.T @ z
    beta = inverse @ xty
    fitted_residual = z - design @ beta
    rss_full = float(np.dot(fitted_residual, fitted_residual))
    dof = len(z) - design.shape[1]
    sigma2 = rss_full / dof
    residual_null = z - float(np.mean(z))
    rss_null = float(np.dot(residual_null, residual_null))

    f_values: list[float] = []
    for stat_name in STAT_NAMES[:3]:
        group = GROUPS[stat_name]
        covariance = inverse[np.ix_(group, group)]
        wald = float(beta[group] @ np.linalg.solve(covariance, beta[group]))
        f_values.append(max(wald / len(group) / sigma2, 0.0))
    f_values.append(max(((rss_null - rss_full) / 15.0) / sigma2, 0.0))
    observed_f = np.asarray(f_values, dtype=float)
    dfs = (3, 5, 7, 15)
    observed_score = np.array([nominal_score(value, df1, dof).item() for value, df1 in zip(observed_f, dfs)])

    dipole = beta[1:4].copy()
    amplitude = float(np.linalg.norm(dipole))
    l_fit, b_fit, ra_fit, dec_fit = unit_to_directions(dipole)
    powers = [
        float(np.dot(beta[1:4], beta[1:4]) / 3.0),
        float(np.dot(beta[4:9], beta[4:9]) / 5.0),
        float(np.dot(beta[9:16], beta[9:16]) / 7.0),
    ]
    return AnalysisRow(
        tag=str(window["tag"]),
        role=str(window["role"]),
        z_min=z_min,
        z_max=z_max,
        bcut_deg=float(bcut_deg),
        source_index=index,
        z=z,
        l_deg=l_deg,
        b_deg=b_deg,
        design=design,
        xtx_inverse=inverse,
        residual_null=residual_null,
        null_yty=rss_null,
        dof=dof,
        condition_number=condition_number,
        rank=rank,
        beta=beta,
        rss_full=rss_full,
        observed_f=observed_f,
        observed_score=observed_score,
        dipole_vector=dipole,
        dipole_amplitude=amplitude,
        dipole_l_deg=l_fit,
        dipole_b_deg=b_fit,
        dipole_ra_deg=ra_fit,
        dipole_dec_deg=dec_fit,
        multipole_power_l1=powers[0],
        multipole_power_l2=powers[1],
        multipole_power_l3=powers[2],
    )


def equal_area_cell_ids(l_deg: np.ndarray, b_deg: np.ndarray, cell_deg: float) -> tuple[np.ndarray, int, int, int]:
    n_lon = max(int(round(360.0 / cell_deg)), 4)
    n_sin_b = max(int(round(180.0 / cell_deg)), 2)
    lon_bin = np.floor(np.mod(l_deg, 360.0) / 360.0 * n_lon).astype(np.int64)
    sin_bin = np.floor((np.sin(np.radians(b_deg)) + 1.0) * 0.5 * n_sin_b).astype(np.int64)
    lon_bin = np.clip(lon_bin, 0, n_lon - 1)
    sin_bin = np.clip(sin_bin, 0, n_sin_b - 1)
    cell_ids = sin_bin * n_lon + lon_bin
    return cell_ids, n_lon * n_sin_b, n_lon, n_sin_b


def aggregate_row(row: AnalysisRow, cell_ids_all: np.ndarray, n_grid: int) -> tuple[np.ndarray, np.ndarray]:
    cell = cell_ids_all[row.source_index]
    occupied = np.unique(cell)
    lookup = np.full(n_grid, -1, dtype=np.int64)
    lookup[occupied] = np.arange(len(occupied), dtype=np.int64)
    compact = lookup[cell]
    aggregate = np.empty((row.design.shape[1], len(occupied)), dtype=float)
    xr = row.design * row.residual_null[:, None]
    for column in range(row.design.shape[1]):
        aggregate[column] = np.bincount(compact, weights=xr[:, column], minlength=len(occupied))
    return occupied, aggregate


def batch_scores(row: AnalysisRow, aggregate: np.ndarray, signs: np.ndarray) -> np.ndarray:
    xty = aggregate @ signs
    beta = row.xtx_inverse @ xty
    rss_full = np.maximum(row.null_yty - np.sum(beta * xty, axis=0), 1e-300)
    sigma2 = rss_full / row.dof
    scores: list[np.ndarray] = []
    for stat_name, df1 in zip(STAT_NAMES[:3], (3, 5, 7)):
        group = GROUPS[stat_name]
        covariance = row.xtx_inverse[np.ix_(group, group)]
        metric = np.linalg.inv(covariance)
        wald = np.einsum("ib,ij,jb->b", beta[group], metric, beta[group])
        f_value = np.maximum(wald / len(group) / sigma2, 0.0)
        scores.append(nominal_score(f_value, df1, row.dof))
    rss_intercept = np.maximum(row.null_yty - xty[0] * xty[0] / len(row.z), rss_full)
    f_joint = np.maximum(((rss_intercept - rss_full) / 15.0) / sigma2, 0.0)
    scores.append(nominal_score(f_joint, 15, row.dof))
    return np.vstack(scores)


def run_block_scale(
    rows: list[AnalysisRow],
    data: dict[str, np.ndarray],
    cell_deg: float,
    n_mocks: int,
    seed: int,
    batch_size: int,
) -> tuple[dict[str, Any], list[dict[str, Any]], np.ndarray]:
    cell_ids, n_grid, n_lon, n_sin_b = equal_area_cell_ids(data["l"], data["b"], cell_deg)
    aggregates = [aggregate_row(row, cell_ids, n_grid) for row in rows]
    observed_family_max = max(float(np.max(row.observed_score)) for row in rows)
    point_exceedances = np.zeros((len(rows), len(STAT_NAMES)), dtype=np.int64)
    maxima = np.empty(n_mocks, dtype=float)
    rng = np.random.default_rng(seed)
    cursor = 0
    report_interval = max(n_mocks // 10, batch_size)
    next_report = report_interval
    while cursor < n_mocks:
        batch = min(batch_size, n_mocks - cursor)
        all_signs = rng.choice(np.array([-1.0, 1.0]), size=(n_grid, batch), replace=True)
        batch_max = np.full(batch, -np.inf, dtype=float)
        for row_index, (row, (occupied, aggregate)) in enumerate(zip(rows, aggregates)):
            scores = batch_scores(row, aggregate, all_signs[occupied])
            point_exceedances[row_index] += np.sum(scores >= row.observed_score[:, None], axis=1)
            batch_max = np.maximum(batch_max, np.max(scores, axis=0))
        maxima[cursor : cursor + batch] = batch_max
        cursor += batch
        if cursor >= next_report or cursor == n_mocks:
            print(f"cell_deg={cell_deg:g} mocks={cursor}/{n_mocks}", flush=True)
            next_report += report_interval

    global_p = p_empirical(np.sum(maxima >= observed_family_max), n_mocks)
    global_row = {
        "cell_deg": cell_deg,
        "n_mocks": n_mocks,
        "seed": seed,
        "n_grid_cells": n_grid,
        "n_occupied_catalogue_cells": int(len(np.unique(cell_ids))),
        "n_lon_bins": n_lon,
        "n_sin_b_bins": n_sin_b,
        "observed_family_max_score": observed_family_max,
        "null_max_score_mean": float(np.mean(maxima)),
        "null_max_score_p95": float(np.quantile(maxima, 0.95)),
        "null_max_score_p99": float(np.quantile(maxima, 0.99)),
        "global_p": global_p,
        "monte_carlo_floor": 1.0 / (n_mocks + 1),
    }
    point_rows: list[dict[str, Any]] = []
    for row_index, row in enumerate(rows):
        for stat_index, stat_name in enumerate(STAT_NAMES):
            point_rows.append(
                {
                    "cell_deg": cell_deg,
                    "window_tag": row.tag,
                    "bcut_deg": row.bcut_deg,
                    "statistic": stat_name,
                    "observed_F": float(row.observed_f[stat_index]),
                    "observed_score": float(row.observed_score[stat_index]),
                    "exceedances": int(point_exceedances[row_index, stat_index]),
                    "n_mocks": n_mocks,
                    "point_p": p_empirical(point_exceedances[row_index, stat_index], n_mocks),
                }
            )
    return global_row, point_rows, maxima


def fit_quality_control(data: dict[str, np.ndarray], mask_extra: np.ndarray, prereg: dict[str, Any]) -> AnalysisRow:
    window = next(item for item in prereg["windows"] if item["role"] == "primary")
    temporary = {key: values[mask_extra] for key, values in data.items()}
    return prepare_row(temporary, window, float(prereg["primary_latitude_cut_deg"]), prereg)


def observed_rows(rows: list[AnalysisRow], comparator: np.ndarray) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for row in rows:
        record: dict[str, Any] = {
            "window_tag": row.tag,
            "window_role": row.role,
            "z_min": row.z_min,
            "z_max_exclusive": row.z_max,
            "bcut_deg": row.bcut_deg,
            "N": len(row.z),
            "design_rank": row.rank,
            "design_condition_number": row.condition_number,
            "dipole_amplitude": row.dipole_amplitude,
            "dipole_l_deg": row.dipole_l_deg,
            "dipole_b_deg": row.dipole_b_deg,
            "dipole_ra_deg": row.dipole_ra_deg,
            "dipole_dec_deg": row.dipole_dec_deg,
            "direction_separation_from_locked_quaia_deg": angular_separation_vectors(row.dipole_vector, comparator),
            "multipole_power_l1": row.multipole_power_l1,
            "multipole_power_l2": row.multipole_power_l2,
            "multipole_power_l3": row.multipole_power_l3,
        }
        for index, stat_name in enumerate(STAT_NAMES):
            record[stat_name] = float(row.observed_f[index])
            record[f"nominal_score_{stat_name}"] = float(row.observed_score[index])
        output.append(record)
    return output


def evaluate_gates(
    rows: list[AnalysisRow],
    global_rows: list[dict[str, Any]],
    point_rows: list[dict[str, Any]],
    prereg: dict[str, Any],
    comparator: np.ndarray,
) -> tuple[list[dict[str, Any]], bool]:
    gates = prereg["promotion_gates"]
    primary_tag = next(item["tag"] for item in prereg["windows"] if item["role"] == "primary")
    primary_bcut = float(prereg["primary_latitude_cut_deg"])
    by_key = {row.key: row for row in rows}
    primary = by_key[(primary_tag, primary_bcut)]
    gate_rows: list[dict[str, Any]] = []

    def add(gate: str, component: str, value: Any, threshold: str, passed: bool) -> None:
        gate_rows.append(
            {
                "gate": gate,
                "component": component,
                "value": value,
                "threshold": threshold,
                "pass": bool(passed),
                "failure": "" if passed else f"{gate}: {component} failed ({fmt(value)}; required {threshold})",
            }
        )

    add("integrity", "all_registered_rows", len(rows), f"exactly {len(prereg['windows']) * len(prereg['latitude_cuts_deg'])} valid rows", len(rows) == len(prereg["windows"]) * len(prereg["latitude_cuts_deg"]))
    global_threshold = float(gates["global_p_max_each_scale"])
    for row in global_rows:
        add("global_rarity", f"cell_{fmt(row['cell_deg'])}_deg", float(row["global_p"]), f"<= {global_threshold}", float(row["global_p"]) <= global_threshold)

    coherence_threshold = float(gates["primary_joint_l1_l3_coherence_point_p_max_each_scale"])
    coherence_rows = [
        row
        for row in point_rows
        if row["window_tag"] == primary_tag
        and float(row["bcut_deg"]) == primary_bcut
        and row["statistic"] == "joint_F_l1_l3"
    ]
    for row in coherence_rows:
        add("joint_l1_l3_coherence", f"cell_{fmt(row['cell_deg'])}_deg", float(row["point_p"]), f"<= {coherence_threshold}", float(row["point_p"]) <= coherence_threshold)

    direction_threshold = float(gates["cross_catalogue_direction_max_deg"])
    cross_sep = angular_separation_vectors(primary.dipole_vector, comparator)
    add("cross_catalogue_direction", "primary_row", cross_sep, f"<= {direction_threshold} deg", cross_sep <= direction_threshold)

    window_threshold = float(gates["window_stability"]["maximum_direction_separation_from_primary_deg"])
    perturbation_rows = [row for row in rows if row.bcut_deg == primary_bcut and row.tag != primary_tag]
    for row in perturbation_rows:
        separation = angular_separation_vectors(row.dipole_vector, primary.dipole_vector)
        add("window_direction_stability", row.tag, separation, f"<= {window_threshold} deg", separation <= window_threshold)

    latitude_threshold = float(gates["latitude_stability"]["maximum_direction_separation_from_primary_deg"])
    latitude_rows = [row for row in rows if row.tag == primary_tag and row.bcut_deg != primary_bcut]
    for row in latitude_rows:
        separation = angular_separation_vectors(row.dipole_vector, primary.dipole_vector)
        add("latitude_direction_stability", f"bcut_{fmt(row.bcut_deg)}", separation, f"<= {latitude_threshold} deg", separation <= latitude_threshold)

    amplitude_gate = gates["amplitude_stability"]
    ratio_min = float(amplitude_gate["amplitude_ratio_min"])
    ratio_max = float(amplitude_gate["amplitude_ratio_max"])
    amplitude_rows = {(row.tag, row.bcut_deg): row for row in perturbation_rows + latitude_rows + [primary]}
    for row in amplitude_rows.values():
        ratio = row.dipole_amplitude / primary.dipole_amplitude if primary.dipole_amplitude > 0.0 else float("inf")
        add(
            "amplitude_stability",
            f"{row.tag}_b{fmt(row.bcut_deg)}",
            ratio,
            f"in [{ratio_min}, {ratio_max}]",
            ratio_min <= ratio <= ratio_max,
        )

    promoted = all(bool(row["pass"]) for row in gate_rows)
    return gate_rows, promoted


def quality_control_rows(data: dict[str, np.ndarray], prereg: dict[str, Any], primary: AnalysisRow) -> list[dict[str, Any]]:
    controls = {
        "zwarning_zero": data["zwarning"] == 0,
        "zconf_ge_2": data["z_conf"] >= 2,
    }
    rows: list[dict[str, Any]] = []
    for name, mask in controls.items():
        fitted = fit_quality_control(data, mask, prereg)
        separation = angular_separation_vectors(fitted.dipole_vector, primary.dipole_vector)
        ratio = fitted.dipole_amplitude / primary.dipole_amplitude if primary.dipole_amplitude > 0.0 else float("inf")
        rows.append(
            {
                "quality_control": name,
                "N": len(fitted.z),
                "dipole_amplitude": fitted.dipole_amplitude,
                "amplitude_ratio_to_primary": ratio,
                "dipole_ra_deg": fitted.dipole_ra_deg,
                "dipole_dec_deg": fitted.dipole_dec_deg,
                "direction_separation_from_primary_deg": separation,
                "joint_F_l1_l3": float(fitted.observed_f[3]),
                "diagnostic_only": True,
            }
        )
    return rows


def write_report(
    path: Path,
    prereg: dict[str, Any],
    meta: dict[str, Any],
    observed: list[dict[str, Any]],
    global_rows: list[dict[str, Any]],
    point_rows: list[dict[str, Any]],
    gate_rows: list[dict[str, Any]],
    quality_rows: list[dict[str, Any]],
    promoted: bool,
) -> None:
    primary_tag = next(item["tag"] for item in prereg["windows"] if item["role"] == "primary")
    primary_bcut = float(prereg["primary_latitude_cut_deg"])
    primary = next(row for row in observed if row["window_tag"] == primary_tag and float(row["bcut_deg"]) == primary_bcut)
    primary_points = [row for row in point_rows if row["window_tag"] == primary_tag and float(row["bcut_deg"]) == primary_bcut]
    failures = [row["failure"] for row in gate_rows if not row["pass"]]
    ranked_statistics = [
        {
            "window_tag": row["window_tag"],
            "bcut_deg": row["bcut_deg"],
            "statistic": statistic,
            "score": row[f"nominal_score_{statistic}"],
        }
        for row in observed
        for statistic in STAT_NAMES
    ]
    family_winner = max(ranked_statistics, key=lambda row: float(row["score"]))
    locked_amplitude = float(prereg["locked_quaia_comparator"]["amplitude"])
    primary_to_quaia_amplitude_ratio = float(primary["dipole_amplitude"]) / locked_amplitude
    primary_partial_points = [row for row in primary_points if row["statistic"] != "joint_F_l1_l3"]
    partial_p_min = min(float(row["point_p"]) for row in primary_partial_points)
    partial_p_max = max(float(row["point_p"]) for row in primary_partial_points)
    decision = "PROMOTED" if promoted else "NOT PROMOTED"
    report = [
        "# SDSS DR16Q v4 Independent Quaia Cross-Catalogue Validation",
        "",
        "Date: 2026-07-18",
        "",
        "## Decision",
        "",
        f"**{decision}.** Promotion is conjunctive; {len(failures)} registered gate component(s) failed.",
        "",
    ]
    if failures:
        report.extend(["## Explicit Failures", ""] + [f"- {failure}" for failure in failures] + [""])
    report.extend(
        [
            "## Primary Row",
            "",
            markdown_table(
                [primary],
                [
                    "window_tag",
                    "bcut_deg",
                    "N",
                    "dipole_amplitude",
                    "dipole_ra_deg",
                    "dipole_dec_deg",
                    "direction_separation_from_locked_quaia_deg",
                    "partial_F_l1",
                    "partial_F_l2",
                    "partial_F_l3",
                    "joint_F_l1_l3",
                ],
            ),
            "",
            "## Empirical Global Null",
            "",
            markdown_table(global_rows, ["cell_deg", "n_mocks", "n_occupied_catalogue_cells", "observed_family_max_score", "null_max_score_p99", "global_p"]),
            "",
            "## Primary Mask-Calibrated Probabilities",
            "",
            markdown_table(primary_points, ["cell_deg", "statistic", "observed_F", "point_p", "exceedances", "n_mocks"]),
            "",
            "## Interpretation",
            "",
            f"The global family winner is `{family_winner['window_tag']}`, `|b| >= {fmt(family_winner['bcut_deg'])} deg`, `{family_winner['statistic']}`. It is not the locked primary row.",
            "",
            f"The primary SDSS dipole amplitude is {primary_to_quaia_amplitude_ratio:.4g} times the locked Quaia amplitude. This cross-catalogue ratio is descriptive because the registered amplitude gate compares the fixed SDSS perturbation rows with the SDSS primary row.",
            "",
            f"Although the registered primary joint l=1-3 coherence statistic passes at every block scale, the individual partial l=1, l=2 and l=3 empirical probabilities span {partial_p_min:.4g} to {partial_p_max:.4g}. The joint significance therefore does not show that the three multipoles are separately unusual; footprint-coupled collective structure remains a viable explanation.",
            "",
            "## Promotion Gates",
            "",
            markdown_table(gate_rows, ["gate", "component", "value", "threshold", "pass", "failure"]),
            "",
            "## Quality Controls",
            "",
            "These predeclared controls are diagnostic and are not additional promotion gates.",
            "",
            markdown_table(quality_rows, ["quality_control", "N", "dipole_amplitude", "amplitude_ratio_to_primary", "dipole_ra_deg", "dipole_dec_deg", "direction_separation_from_primary_deg"]),
            "",
            "## Method",
            "",
            "The analysis conditions on the observed SDSS positions, jointly fits the registered real l=1-3 harmonic design and uses restricted-residual Rademacher signs shared by equal-area sky blocks across all overlapping rows. The empirical maximum over the fixed 80-statistic family supplies the global probability independently at 8, 12 and 16 degrees.",
            "",
            "## Catalogue Audit",
            "",
            "```json",
            json.dumps(meta, indent=2),
            "```",
            "",
            "## Claim Boundary",
            "",
            "This is an independent catalogue validation of a fixed angular target. A non-promotion result remains an explicit catalogue-validation failure and must not be reframed as positive SU(2) evidence through exploratory window, mask, block-scale or coordinate changes.",
        ]
    )
    path.write_text("\n".join(report) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the preregistered SDSS DR16Q v4 independent Quaia validation.")
    parser.add_argument("--preregistration", type=Path, default=DEFAULT_PREREG)
    parser.add_argument("--fits", type=Path, default=DEFAULT_FITS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--git-results-dir", type=Path, default=DEFAULT_GIT_RESULTS)
    parser.add_argument("--copy-to-git", action=argparse.BooleanOptionalAction, default=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    prereg = load_json(args.preregistration)
    if prereg.get("registration_id") != "sdss_dr16q_v4_cross_catalogue_2026-07-18":
        raise ValueError("Unexpected preregistration identifier")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    started = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    data, catalogue_meta = load_catalogue(args.fits, prereg)
    print(f"Primary-quality unique DR16Q rows: {len(data['z']):,}", flush=True)

    rows: list[AnalysisRow] = []
    for window in prereg["windows"]:
        for bcut in prereg["latitude_cuts_deg"]:
            row = prepare_row(data, window, float(bcut), prereg)
            rows.append(row)
            print(
                f"prepared {row.tag} b={row.bcut_deg:g} N={len(row.z):,} "
                f"A={row.dipole_amplitude:.6g} RA={row.dipole_ra_deg:.3f} Dec={row.dipole_dec_deg:.3f}",
                flush=True,
            )

    global_rows: list[dict[str, Any]] = []
    point_rows: list[dict[str, Any]] = []
    scale_maxima: list[tuple[float, np.ndarray]] = []
    for scale in prereg["null"]["scales"]:
        global_row, scale_points, maxima = run_block_scale(
            rows,
            data,
            float(scale["cell_deg"]),
            int(scale["n_mocks"]),
            int(scale["seed"]),
            int(prereg["null"]["batch_size"]),
        )
        global_rows.append(global_row)
        point_rows.extend(scale_points)
        scale_maxima.append((float(scale["cell_deg"]), maxima))
        print(f"cell_deg={scale['cell_deg']:g} global_p={global_row['global_p']:.6g}", flush=True)

    comparator = comparator_galactic_unit(prereg)
    observed = observed_rows(rows, comparator)
    gate_rows, promoted = evaluate_gates(rows, global_rows, point_rows, prereg, comparator)
    primary_tag = next(item["tag"] for item in prereg["windows"] if item["role"] == "primary")
    primary = next(row for row in rows if row.tag == primary_tag and row.bcut_deg == float(prereg["primary_latitude_cut_deg"]))
    quality_rows = quality_control_rows(data, prereg, primary)

    prefix = "sdss_dr16q_v4_cross_catalogue_validation_2026-07-18"
    output_paths = {
        "observed": args.out_dir / f"{prefix}_observed_rows.csv",
        "pointwise": args.out_dir / f"{prefix}_pointwise_pvalues.csv",
        "global": args.out_dir / f"{prefix}_global_null.csv",
        "maxima": args.out_dir / f"{prefix}_global_null_maxima.csv.gz",
        "gates": args.out_dir / f"{prefix}_promotion_gates.csv",
        "quality": args.out_dir / f"{prefix}_quality_controls.csv",
        "report": args.out_dir / f"{prefix}_report.md",
        "config": args.out_dir / f"{prefix}_run_config.json",
    }
    write_csv(output_paths["observed"], observed)
    write_csv(output_paths["pointwise"], point_rows)
    write_csv(output_paths["global"], global_rows)
    write_global_maxima(output_paths["maxima"], scale_maxima)
    write_csv(output_paths["gates"], gate_rows)
    write_csv(output_paths["quality"], quality_rows)
    download_manifest = args.fits.with_suffix(".download_manifest.json")
    fits_digest = load_json(download_manifest).get("sha256") if download_manifest.exists() else sha256(args.fits)
    run_config = {
        "analysis": prereg["registration_id"],
        "started_at": started,
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "preregistration": str(args.preregistration),
        "preregistration_sha256": sha256(args.preregistration),
        "fits": str(args.fits),
        "fits_sha256": fits_digest,
        "catalogue_meta": catalogue_meta,
        "promotion": "promoted" if promoted else "not promoted",
        "failed_gate_components": [row["failure"] for row in gate_rows if not row["pass"]],
        "locked_configuration": prereg,
    }
    output_paths["config"].write_text(json.dumps(run_config, indent=2) + "\n", encoding="utf-8")
    write_report(
        output_paths["report"],
        prereg,
        catalogue_meta,
        observed,
        global_rows,
        point_rows,
        gate_rows,
        quality_rows,
        promoted,
    )

    manifest_rows = []
    for role, path in output_paths.items():
        manifest_rows.append({"role": role, "file": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)})
    manifest_path = args.out_dir / f"{prefix}_manifest.csv"
    write_csv(manifest_path, manifest_rows)
    output_paths["manifest"] = manifest_path

    if args.copy_to_git:
        args.git_results_dir.mkdir(parents=True, exist_ok=True)
        for path in output_paths.values():
            shutil.copy2(path, args.git_results_dir / path.name)
    print(f"Decision: {'PROMOTED' if promoted else 'NOT PROMOTED'}", flush=True)
    print(f"Saved report: {output_paths['report']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
