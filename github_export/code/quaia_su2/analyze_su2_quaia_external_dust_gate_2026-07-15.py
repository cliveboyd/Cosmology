from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy_healpix import HEALPix
from scipy.ndimage import map_coordinates


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
DEFAULT_QUAIA = ROOT / ".spyder-py3" / "quaia_G20.0.fits"
DEFAULT_RANDOMS = ROOT / ".spyder-py3" / "random_G20.0_10x.fits"
DEFAULT_SELECTION = ROOT / ".spyder-py3" / "selection_function_NSIDE64_G20.0.fits"
DEFAULT_SFD_DIR = ROOT / "external_datasets" / "sfd"
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_external_dust_gate_20260715"

SOURCE_NOTES = {
    "sfd": "Schlegel, Finkbeiner & Davis 1998 full-sky E(B-V), Harvard Dataverse DOI 10.7910/DVN/EWCNL5.",
    "selection": "Local Quaia G<20.0 selection_function_NSIDE64 map.",
    "random_density": "Local Quaia G<20.0 10x random catalog binned to the selection HEALPix grid.",
    "gaia_scanlaw": "Nominal Gaia scanning-law FITS bundle from the gaiascanlaw package, counted to the DR3 end date.",
}

GAIA_TSTART = Time("2014-07-25T10:31:25.554960001", format="isot", scale="tcb").decimalyear
GAIA_TDR3 = Time("2017-05-28T08:46:28.954612431", format="isot", scale="tcb").decimalyear

EQ_TO_GAL = np.array(
    [
        [-0.0548755604, -0.8734370902, -0.4838350155],
        [0.4941094279, -0.4448296300, 0.7469822445],
        [-0.8676661490, -0.1980763734, 0.4559837762],
    ],
    dtype=float,
)


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


def gal_unit_array(l_deg: np.ndarray, b_deg: np.ndarray) -> np.ndarray:
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    cb = np.cos(b)
    return np.column_stack([cb * np.cos(l), cb * np.sin(l), np.sin(b)])


def eq_unit_array(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    cd = np.cos(dec)
    return np.column_stack([cd * np.cos(ra), cd * np.sin(ra), np.sin(dec)])


def equatorial_to_galactic(ra_deg: np.ndarray, dec_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    eq = eq_unit_array(ra_deg, dec_deg)
    gal = eq @ EQ_TO_GAL.T
    gal /= np.linalg.norm(gal, axis=1)[:, None]
    l_gal = (np.rad2deg(np.arctan2(gal[:, 1], gal[:, 0])) + 360.0) % 360.0
    b_gal = np.rad2deg(np.arcsin(np.clip(gal[:, 2], -1.0, 1.0)))
    return l_gal, b_gal


def unit_from_lb(l_deg: float, b_deg: float) -> np.ndarray:
    l = math.radians(l_deg)
    b = math.radians(b_deg)
    cb = math.cos(b)
    return np.array([cb * math.cos(l), cb * math.sin(l), math.sin(b)], dtype=float)


def unit_to_lb(vec: np.ndarray) -> tuple[float, float]:
    norm = float(np.linalg.norm(vec))
    if norm <= 0.0 or not math.isfinite(norm):
        return float("nan"), float("nan")
    unit = vec / norm
    l_deg = math.degrees(math.atan2(float(unit[1]), float(unit[0]))) % 360.0
    b_deg = math.degrees(math.asin(float(np.clip(unit[2], -1.0, 1.0))))
    return l_deg, b_deg


def angular_sep_deg(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 0.0:
        return float("nan")
    dot = float(np.clip(np.dot(a, b) / denom, -1.0, 1.0))
    return math.degrees(math.acos(dot))


def standardize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    med = np.nanmedian(values)
    scale = np.nanstd(values)
    if not math.isfinite(scale) or scale <= 0.0:
        scale = 1.0
    return (values - med) / scale


def finite_mask(*arrays: np.ndarray) -> np.ndarray:
    mask = np.ones(len(arrays[0]), dtype=bool)
    for arr in arrays:
        mask &= np.isfinite(arr)
    return mask


def healpix_ang2pix_ring(nside: int, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    z = np.cos(theta)
    za = np.abs(z)
    tt = (phi % (2.0 * np.pi)) / (0.5 * np.pi)
    npix = 12 * nside * nside
    ncap = 2 * nside * (nside - 1)
    ipix = np.empty(len(theta), dtype=np.int64)

    equ = za <= (2.0 / 3.0)
    if np.any(equ):
        temp1 = nside * (0.5 + tt[equ])
        temp2 = nside * z[equ] * 0.75
        jp = np.floor(temp1 - temp2).astype(np.int64)
        jm = np.floor(temp1 + temp2).astype(np.int64)
        ir = nside + 1 + jp - jm
        kshift = 1 - (ir & 1)
        ip = ((jp + jm - nside + kshift + 1) // 2) % (4 * nside)
        ipix[equ] = ncap + (ir - 1) * 4 * nside + ip

    pol = ~equ
    if np.any(pol):
        tp = tt[pol] - np.floor(tt[pol])
        tmp = nside * np.sqrt(3.0 * (1.0 - za[pol]))
        jp = np.floor(tp * tmp).astype(np.int64)
        jm = np.floor((1.0 - tp) * tmp).astype(np.int64)
        ir = jp + jm + 1
        ip = np.floor(tt[pol] * ir).astype(np.int64) % (4 * ir)
        north = z[pol] > 0.0
        out = np.empty(np.sum(pol), dtype=np.int64)
        out[north] = 2 * ir[north] * (ir[north] - 1) + ip[north]
        out[~north] = npix - 2 * ir[~north] * (ir[~north] + 1) + ip[~north]
        ipix[pol] = out
    return np.clip(ipix, 0, npix - 1)


def query_sfd(sfd_dir: Path, l_deg: np.ndarray, b_deg: np.ndarray, order: int = 1) -> np.ndarray:
    out = np.full(len(l_deg), np.nan, dtype=np.float64)
    for pole, pole_mask in {
        "ngp": b_deg >= 0.0,
        "sgp": b_deg < 0.0,
    }.items():
        if not np.any(pole_mask):
            continue
        path = sfd_dir / f"SFD_dust_4096_{pole}.fits"
        if not path.exists():
            raise FileNotFoundError(f"Missing SFD map: {path}")
        with fits.open(path, memmap=True) as hdul:
            data = np.asarray(hdul[0].data, dtype=np.float64)
            wcs = WCS(hdul[0].header)
        x, y = wcs.wcs_world2pix(l_deg[pole_mask], b_deg[pole_mask], 0)
        out[pole_mask] = map_coordinates(data, [y, x], order=order, mode="nearest")
    return out


def resolve_gaiascanlaw_data_dir() -> Path | None:
    spec = importlib.util.find_spec("gaiascanlaw")
    if spec is None or spec.origin is None:
        return None
    data_dir = Path(spec.origin).parent / "data"
    if (data_dir / "gost_simple_1.fits").exists() and (data_dir / "gost_simple_2.fits").exists():
        return data_dir
    return None


def load_gaia_scanlaw_templates(ra_deg: np.ndarray, dec_deg: np.ndarray, data_dir: Path | None) -> dict[str, np.ndarray]:
    if data_dir is None:
        nan = np.full(len(ra_deg), np.nan, dtype=float)
        return {
            "gaia_scan_count_dr3": nan.copy(),
            "gaia_scan_count_log1p_dr3": nan.copy(),
            "gaia_scan_angle_cos2_mean_dr3": nan.copy(),
            "gaia_scan_angle_sin2_mean_dr3": nan.copy(),
            "gaia_scan_angle_resultant_dr3": nan.copy(),
        }

    nside = 64
    npix = 12 * nside * nside
    healpix = HEALPix(nside=nside, order="nested")
    object_pix = healpix.lonlat_to_healpix(ra_deg * u.deg, dec_deg * u.deg)

    counts = np.zeros(npix, dtype=float)
    cos2_sum = np.zeros(npix, dtype=float)
    sin2_sum = np.zeros(npix, dtype=float)
    for filename in ("gost_simple_1.fits", "gost_simple_2.fits"):
        with fits.open(data_dir / filename, memmap=True) as hdul:
            data = hdul[1].data
            pix = np.asarray(data["healpixel"], dtype=np.int64)
            times = np.asarray(data["scan time [decimal year]"], dtype=float)
            angles = np.asarray(data["scan angle [radian]"], dtype=float)
        mask = (times > GAIA_TSTART) & (times < GAIA_TDR3) & (pix >= 0) & (pix < npix)
        pix = pix[mask]
        angles = angles[mask]
        counts += np.bincount(pix, minlength=npix).astype(float)
        cos2_sum += np.bincount(pix, weights=np.cos(2.0 * angles), minlength=npix).astype(float)
        sin2_sum += np.bincount(pix, weights=np.sin(2.0 * angles), minlength=npix).astype(float)

    denom = np.where(counts > 0.0, counts, 1.0)
    cos2_mean = cos2_sum / denom
    sin2_mean = sin2_sum / denom
    resultant = np.sqrt(cos2_sum * cos2_sum + sin2_sum * sin2_sum) / denom
    return {
        "gaia_scan_count_dr3": counts[object_pix],
        "gaia_scan_count_log1p_dr3": np.log1p(counts[object_pix]),
        "gaia_scan_angle_cos2_mean_dr3": cos2_mean[object_pix],
        "gaia_scan_angle_sin2_mean_dr3": sin2_mean[object_pix],
        "gaia_scan_angle_resultant_dr3": resultant[object_pix],
    }


def load_selection_map(path: Path) -> tuple[np.ndarray, int]:
    with fits.open(path, memmap=True) as hdul:
        nside = int(hdul[1].header.get("NSIDE", 64))
        values = np.asarray(hdul[1].data["T"], dtype=float).reshape(-1)
    return values, nside


def load_data(
    quaia_path: Path,
    selection_path: Path,
    randoms_path: Path,
    sfd_dir: Path,
    gaiascanlaw_data_dir: Path | None,
) -> dict[str, np.ndarray]:
    with fits.open(quaia_path, memmap=True) as hdul:
        data = hdul[1].data
        q = {
            "z": np.asarray(data["redshift_quaia"], dtype=float),
            "zerr": np.asarray(data["redshift_quaia_err"], dtype=float),
            "ra": np.asarray(data["ra"], dtype=float),
            "dec": np.asarray(data["dec"], dtype=float),
            "l": np.asarray(data["l"], dtype=float),
            "b": np.asarray(data["b"], dtype=float),
            "g": np.asarray(data["phot_g_mean_mag"], dtype=float),
            "w1": np.asarray(data["mag_w1_vg"], dtype=float),
            "w2": np.asarray(data["mag_w2_vg"], dtype=float),
            "pmra_error": np.asarray(data["pmra_error"], dtype=float),
            "pmdec_error": np.asarray(data["pmdec_error"], dtype=float),
        }
    q["vec"] = gal_unit_array(q["l"], q["b"])
    q["ebv_sfd"] = query_sfd(sfd_dir, q["l"], q["b"])

    selection_map, nside = load_selection_map(selection_path)
    theta = np.deg2rad(90.0 - q["b"])
    phi = np.deg2rad(q["l"])
    pix = healpix_ang2pix_ring(nside, theta, phi)
    q["selection_T"] = selection_map[pix]

    with fits.open(randoms_path, memmap=True) as hdul:
        r = hdul[1].data
        r_l, r_b = equatorial_to_galactic(np.asarray(r["RA"], dtype=float), np.asarray(r["DEC"], dtype=float))
    r_pix = healpix_ang2pix_ring(nside, np.deg2rad(90.0 - r_b), np.deg2rad(r_l))
    counts = np.bincount(r_pix, minlength=len(selection_map)).astype(float)
    q["random_density_log1p"] = np.log1p(counts[pix])
    q["w1_w2"] = q["w1"] - q["w2"]
    q["sfd_log1p"] = np.log1p(np.clip(q["ebv_sfd"], 0.0, None))
    q["sfd_sq"] = q["ebv_sfd"] ** 2
    q.update(load_gaia_scanlaw_templates(q["ra"], q["dec"], gaiascanlaw_data_dir))
    return q


def fit_dipole(z: np.ndarray, vec: np.ndarray, controls: list[np.ndarray] | None = None) -> dict[str, float]:
    columns = [np.ones(len(z)), vec[:, 0], vec[:, 1], vec[:, 2]]
    if controls:
        columns.extend(standardize(col) for col in controls)
    x = np.column_stack(columns)
    beta, *_ = np.linalg.lstsq(x, z, rcond=None)
    resid = z - x @ beta
    dof = max(len(z) - x.shape[1], 1)
    sigma2 = float(np.dot(resid, resid) / dof)
    try:
        cov = sigma2 * np.linalg.pinv(x.T @ x)
        dvec = beta[1:4]
        cov_d = cov[1:4, 1:4]
        snr = float(math.sqrt(max(float(dvec @ np.linalg.pinv(cov_d) @ dvec), 0.0)))
    except Exception:
        dvec = beta[1:4]
        snr = float("nan")
    amp = float(np.linalg.norm(dvec))
    l_deg, b_deg = unit_to_lb(dvec)
    return {
        "amp": amp,
        "l_deg": l_deg,
        "b_deg": b_deg,
        "formal_snr": snr,
        "rms_resid": float(math.sqrt(max(sigma2, 0.0))),
    }


def summarize(control_group: str, fits: list[dict[str, Any]], baseline: dict[str, Any] | None) -> dict[str, Any]:
    amps = np.array([row["amp"] for row in fits], dtype=float)
    snrs = np.array([row["formal_snr"] for row in fits], dtype=float)
    dirs = np.array([unit_from_lb(row["l_deg"], row["b_deg"]) for row in fits], dtype=float)
    max_pair_sep = 0.0
    for i in range(len(dirs)):
        for j in range(i + 1, len(dirs)):
            max_pair_sep = max(max_pair_sep, angular_sep_deg(dirs[i], dirs[j]))
    mean_dir = np.sum(dirs, axis=0)
    mean_resultant = float(np.linalg.norm(mean_dir) / max(len(dirs), 1))
    row = {
        "control_group": control_group,
        "N_min": int(min(fit["N"] for fit in fits)),
        "N_max": int(max(fit["N"] for fit in fits)),
        "amp_mean": float(np.nanmean(amps)),
        "amp_min": float(np.nanmin(amps)),
        "amp_max": float(np.nanmax(amps)),
        "snr_max": float(np.nanmax(snrs)),
        "snr_mean": float(np.nanmean(snrs)),
        "max_pair_direction_sep_deg": float(max_pair_sep),
        "mean_resultant_length": mean_resultant,
    }
    if baseline is not None:
        row["amp_ratio_vs_baseline"] = row["amp_mean"] / baseline["amp_mean"] if baseline["amp_mean"] else float("nan")
        row["snr_delta_vs_baseline"] = row["snr_max"] - baseline["snr_max"]
        row["pair_sep_delta_vs_baseline"] = row["max_pair_direction_sep_deg"] - baseline["max_pair_direction_sep_deg"]
    else:
        row["amp_ratio_vs_baseline"] = 1.0
        row["snr_delta_vs_baseline"] = 0.0
        row["pair_sep_delta_vs_baseline"] = 0.0
    return row


def corr_row(q: dict[str, np.ndarray], mask: np.ndarray, name: str, values: np.ndarray) -> dict[str, Any]:
    z = q["z"][mask]
    v = values[mask]
    m = np.isfinite(z) & np.isfinite(v)
    if np.sum(m) < 3:
        return {"template": name, "N": int(np.sum(m)), "corr_z": float("nan")}
    return {
        "template": name,
        "N": int(np.sum(m)),
        "mean": float(np.nanmean(v[m])),
        "std": float(np.nanstd(v[m])),
        "corr_z": float(np.corrcoef(z[m], v[m])[0, 1]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run locked SU2/Quaia external dust/selection template gate.")
    parser.add_argument("--quaia", type=Path, default=DEFAULT_QUAIA)
    parser.add_argument("--randoms", type=Path, default=DEFAULT_RANDOMS)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--sfd-dir", type=Path, default=DEFAULT_SFD_DIR)
    parser.add_argument("--gaiascanlaw-data-dir", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--z-min", type=float, default=0.95)
    parser.add_argument("--z-max", type=float, default=1.45)
    parser.add_argument("--b-cuts", type=float, nargs="+", default=[15.0, 20.0, 25.0, 30.0, 35.0])
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    gaiascanlaw_data_dir = args.gaiascanlaw_data_dir or resolve_gaiascanlaw_data_dir()
    q = load_data(args.quaia, args.selection, args.randoms, args.sfd_dir, gaiascanlaw_data_dir)

    groups: dict[str, list[tuple[str, np.ndarray]]] = {
        "baseline_no_template": [],
        "sfd_dust_external": [
            ("ebv_sfd", q["ebv_sfd"]),
            ("sfd_log1p", q["sfd_log1p"]),
            ("sfd_sq", q["sfd_sq"]),
        ],
        "quaia_selection_depth": [
            ("selection_T", q["selection_T"]),
            ("random_density_log1p", q["random_density_log1p"]),
        ],
        "sfd_plus_selection_depth": [
            ("ebv_sfd", q["ebv_sfd"]),
            ("sfd_log1p", q["sfd_log1p"]),
            ("sfd_sq", q["sfd_sq"]),
            ("selection_T", q["selection_T"]),
            ("random_density_log1p", q["random_density_log1p"]),
        ],
        "gaia_scanlaw_external": [
            ("gaia_scan_count_log1p_dr3", q["gaia_scan_count_log1p_dr3"]),
            ("gaia_scan_angle_cos2_mean_dr3", q["gaia_scan_angle_cos2_mean_dr3"]),
            ("gaia_scan_angle_sin2_mean_dr3", q["gaia_scan_angle_sin2_mean_dr3"]),
            ("gaia_scan_angle_resultant_dr3", q["gaia_scan_angle_resultant_dr3"]),
        ],
        "sfd_selection_gaia_external": [
            ("ebv_sfd", q["ebv_sfd"]),
            ("sfd_log1p", q["sfd_log1p"]),
            ("sfd_sq", q["sfd_sq"]),
            ("selection_T", q["selection_T"]),
            ("random_density_log1p", q["random_density_log1p"]),
            ("gaia_scan_count_log1p_dr3", q["gaia_scan_count_log1p_dr3"]),
            ("gaia_scan_angle_cos2_mean_dr3", q["gaia_scan_angle_cos2_mean_dr3"]),
            ("gaia_scan_angle_sin2_mean_dr3", q["gaia_scan_angle_sin2_mean_dr3"]),
            ("gaia_scan_angle_resultant_dr3", q["gaia_scan_angle_resultant_dr3"]),
        ],
        "local_catalog_quality_proxy": [
            ("zerr", q["zerr"]),
            ("g", q["g"]),
            ("w1", q["w1"]),
            ("w2", q["w2"]),
            ("w1_w2", q["w1_w2"]),
            ("pmra_error", q["pmra_error"]),
            ("pmdec_error", q["pmdec_error"]),
        ],
        "all_available_templates": [
            ("ebv_sfd", q["ebv_sfd"]),
            ("sfd_log1p", q["sfd_log1p"]),
            ("sfd_sq", q["sfd_sq"]),
            ("selection_T", q["selection_T"]),
            ("random_density_log1p", q["random_density_log1p"]),
            ("gaia_scan_count_log1p_dr3", q["gaia_scan_count_log1p_dr3"]),
            ("gaia_scan_angle_cos2_mean_dr3", q["gaia_scan_angle_cos2_mean_dr3"]),
            ("gaia_scan_angle_sin2_mean_dr3", q["gaia_scan_angle_sin2_mean_dr3"]),
            ("gaia_scan_angle_resultant_dr3", q["gaia_scan_angle_resultant_dr3"]),
            ("zerr", q["zerr"]),
            ("g", q["g"]),
            ("w1", q["w1"]),
            ("w2", q["w2"]),
            ("w1_w2", q["w1_w2"]),
            ("pmra_error", q["pmra_error"]),
            ("pmdec_error", q["pmdec_error"]),
        ],
    }

    fit_rows: list[dict[str, Any]] = []
    for group_name, group_controls in groups.items():
        for b_cut in args.b_cuts:
            mask = (
                (q["z"] >= args.z_min)
                & (q["z"] < args.z_max)
                & (np.abs(q["b"]) >= b_cut)
                & finite_mask(q["z"], q["vec"][:, 0], q["vec"][:, 1], q["vec"][:, 2])
            )
            controls = [values[mask] for _, values in group_controls]
            if controls:
                for control in controls:
                    mask2 = np.isfinite(control)
                    if not np.all(mask2):
                        raise RuntimeError(f"Unexpected non-finite controls in {group_name}")
            result = fit_dipole(q["z"][mask], q["vec"][mask], controls=controls)
            fit_rows.append(
                {
                    "control_group": group_name,
                    "b_cut_deg": b_cut,
                    "N": int(np.sum(mask)),
                    "n_templates": len(group_controls),
                    "template_names": ";".join(name for name, _ in group_controls),
                    **result,
                }
            )

    summary_rows: list[dict[str, Any]] = []
    baseline_summary: dict[str, Any] | None = None
    for group_name in groups:
        rows = [row for row in fit_rows if row["control_group"] == group_name]
        row = summarize(group_name, rows, baseline_summary)
        summary_rows.append(row)
        if group_name == "baseline_no_template":
            baseline_summary = row
            summary_rows[-1] = summarize(group_name, rows, None)

    primary_mask = (q["z"] >= args.z_min) & (q["z"] < args.z_max) & (np.abs(q["b"]) >= 15.0)
    corr_rows = [
        corr_row(q, primary_mask, "ebv_sfd", q["ebv_sfd"]),
        corr_row(q, primary_mask, "selection_T", q["selection_T"]),
        corr_row(q, primary_mask, "random_density_log1p", q["random_density_log1p"]),
        corr_row(q, primary_mask, "gaia_scan_count_log1p_dr3", q["gaia_scan_count_log1p_dr3"]),
        corr_row(q, primary_mask, "gaia_scan_angle_resultant_dr3", q["gaia_scan_angle_resultant_dr3"]),
        corr_row(q, primary_mask, "zerr", q["zerr"]),
        corr_row(q, primary_mask, "g", q["g"]),
        corr_row(q, primary_mask, "w1", q["w1"]),
        corr_row(q, primary_mask, "w2", q["w2"]),
        corr_row(q, primary_mask, "w1_w2", q["w1_w2"]),
    ]

    write_csv(args.out_dir / "su2_quaia_external_dust_gate_fit_rows.csv", fit_rows)
    write_csv(args.out_dir / "su2_quaia_external_dust_gate_summary.csv", summary_rows)
    write_csv(args.out_dir / "su2_quaia_external_dust_gate_template_correlations.csv", corr_rows)

    report_cols = [
        "control_group",
        "amp_ratio_vs_baseline",
        "snr_max",
        "snr_delta_vs_baseline",
        "max_pair_direction_sep_deg",
        "mean_resultant_length",
        "N_min",
        "N_max",
    ]
    corr_cols = ["template", "N", "mean", "std", "corr_z"]
    report = [
        "# SU2 / Quaia External Dust Template Gate",
        "",
        "Date: July 15, 2026",
        "",
        "## Locked Test",
        "",
        f"- Redshift window: {args.z_min:g} <= z < {args.z_max:g}",
        f"- Latitude cuts: {', '.join(f'|b|>{b:g}' for b in args.b_cuts)}",
        "- Response: Quaia redshift fitted with a Galactic dipole plus optional standardized template controls.",
        "",
        "## External / Template Inputs",
        "",
        f"- SFD dust: {SOURCE_NOTES['sfd']}",
        f"- Selection depth proxy: {SOURCE_NOTES['selection']}",
        f"- Random density proxy: {SOURCE_NOTES['random_density']}",
        f"- Gaia scan law: {SOURCE_NOTES['gaia_scanlaw']}",
        "- WISE-depth and stellar-density maps were not present as true external all-sky maps in the local tree; WISE magnitudes/colors and Gaia catalog-quality columns are reported as local proxies only.",
        "",
        "## Control Summary",
        "",
        markdown_table(summary_rows, report_cols),
        "",
        "## Template Correlations In Primary Mask",
        "",
        markdown_table(corr_rows, corr_cols),
        "",
        "## Readout",
        "",
    ]

    by_group = {row["control_group"]: row for row in summary_rows}
    sfd = by_group["sfd_dust_external"]
    selection = by_group["quaia_selection_depth"]
    combo = by_group["sfd_plus_selection_depth"]
    gaia = by_group["gaia_scanlaw_external"]
    external_combo = by_group["sfd_selection_gaia_external"]
    all_available = by_group["all_available_templates"]

    if sfd["amp_ratio_vs_baseline"] > 0.9 and sfd["max_pair_direction_sep_deg"] <= 20.0:
        dust_line = "SFD dust does not absorb the locked angular mode in this regression gate."
    else:
        dust_line = "SFD dust materially changes the locked angular mode and must remain a claim gate."
    if combo["amp_ratio_vs_baseline"] > 0.9 and combo["max_pair_direction_sep_deg"] <= 20.0:
        combo_line = "SFD plus local selection/depth proxies still leaves the locked angular mode stable."
    else:
        combo_line = "SFD plus local selection/depth proxies materially suppresses or rotates the mode."
    if external_combo["amp_ratio_vs_baseline"] > 0.9 and external_combo["max_pair_direction_sep_deg"] <= 20.0:
        external_line = "SFD plus selection/depth plus Gaia scan-law controls still leaves the locked angular mode stable."
    else:
        external_line = "SFD plus selection/depth plus Gaia scan-law controls materially suppresses or rotates the mode."
    proxy_line = (
        "The all-available proxy stack remains a stress test rather than a headline correction, "
        "because it mixes true external dust, Quaia selection products, and local catalog-quality columns."
    )
    report.extend(
        [
            f"- {dust_line}",
            f"- Selection/depth-only amp ratio: {selection['amp_ratio_vs_baseline']:.4g}; max direction spread: {selection['max_pair_direction_sep_deg']:.4g} deg.",
            f"- Gaia scan-law-only amp ratio: {gaia['amp_ratio_vs_baseline']:.4g}; max direction spread: {gaia['max_pair_direction_sep_deg']:.4g} deg.",
            f"- {combo_line}",
            f"- {external_line}",
            f"- All-available proxy amp ratio: {all_available['amp_ratio_vs_baseline']:.4g}; max direction spread: {all_available['max_pair_direction_sep_deg']:.4g} deg.",
            f"- {proxy_line}",
            "",
            "## Next Gate",
            "",
            "Acquire true external WISE-depth and stellar-density maps on a documented HEALPix grid, then rerun this exact template gate with those external maps replacing local proxy columns.",
            "",
        ]
    )
    report_path = args.out_dir / "su2_quaia_external_dust_gate_report.md"
    report_path.write_text("\n".join(report), encoding="utf-8")

    readout = [
        "# SU2 / Quaia External Dust Gate Readout",
        "",
        "Date: July 15, 2026",
        "",
        "## Bottom Line",
        "",
        dust_line,
        "",
        "SFD plus the local Quaia selection/depth proxies and bundled Gaia scan-law counts are the best external-template gate available from the current local files. True WISE-depth and stellar-density maps still need acquisition before this becomes the full external-systematics gate.",
        "",
        "## Key Numbers",
        "",
        markdown_table(summary_rows, report_cols),
        "",
    ]
    readout_path = args.out_dir / "su2_quaia_external_dust_gate_readout.md"
    readout_path.write_text("\n".join(readout), encoding="utf-8")

    config = {
        "quaia": str(args.quaia),
        "randoms": str(args.randoms),
        "selection": str(args.selection),
        "sfd_dir": str(args.sfd_dir),
        "gaiascanlaw_data_dir": str(gaiascanlaw_data_dir) if gaiascanlaw_data_dir else None,
        "out_dir": str(args.out_dir),
        "z_min": args.z_min,
        "z_max": args.z_max,
        "b_cuts": args.b_cuts,
        "source_notes": SOURCE_NOTES,
    }
    (args.out_dir / "su2_quaia_external_dust_gate_config.json").write_text(
        json.dumps(config, indent=2),
        encoding="utf-8",
    )
    print(f"Saved report: {report_path}")
    print(f"Saved readout: {readout_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
