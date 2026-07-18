#!/usr/bin/env python3
"""Derive the registered (R, l_A, omega_b h^2) Planck compression."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import zipfile
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
DATE = "2026-07-18"
PREFIX = "planck_pr3_distance_prior"
DEFAULT_ZIP = (
    ROOT
    / "external_datasets"
    / "current_cosmology_sources"
    / "Planck_PR3_chains"
    / "COM_CosmoParams_base-plikHM-TTTEEE-lowl-lowE_R3.00.zip"
)
RESULT_DIR = ROOT / "github_export" / "results" / DATE / "cmb"
PREREG_MD = RESULT_DIR / f"cmb_programme_preregistration_{DATE}.md"
PREREG_JSON = RESULT_DIR / f"cmb_programme_preregistration_{DATE}.json"
C_KM_S = 299792.458


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_paramnames(payload: bytes) -> tuple[list[str], list[str]]:
    raw_names: list[str] = []
    names: list[str] = []
    for line in payload.decode("utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        raw = stripped.split(maxsplit=1)[0]
        raw_names.append(raw)
        names.append(raw.rstrip("*"))
    if len(names) < 6:
        raise RuntimeError("The selected parameter-name file is unexpectedly short")
    return raw_names, names


def choose_chain_family(archive: zipfile.ZipFile) -> tuple[str, list[str]]:
    param_members = [name for name in archive.namelist() if name.endswith(".paramnames")]
    preferred = [
        name
        for name in param_members
        if Path(name).name == "base_plikHM_TTTEEE_lowl_lowE.paramnames"
    ]
    if len(preferred) != 1:
        raise RuntimeError(f"Expected one registered parameter-name member, found {preferred}")
    prefix = preferred[0][: -len(".paramnames")]
    chain_pattern = re.compile(re.escape(prefix) + r"_\d+\.txt$")
    chain_members = sorted(name for name in archive.namelist() if chain_pattern.fullmatch(name))
    if not chain_members:
        raise RuntimeError(f"No chain members found for {prefix}")
    return preferred[0], chain_members


def weighted_moments(values: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    total = float(np.sum(weights))
    if not math.isfinite(total) or total <= 0:
        raise ValueError("Posterior weights do not have a positive finite sum")
    mean = np.sum(values * weights[:, None], axis=0) / total
    centred = values - mean
    covariance = (centred * weights[:, None]).T @ centred / total
    covariance = 0.5 * (covariance + covariance.T)
    effective_n = total * total / float(np.sum(weights * weights))
    return mean, covariance, effective_n


def load_vectors(path: Path) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    with zipfile.ZipFile(path) as archive:
        param_member, chain_members = choose_chain_family(archive)
        param_payload = archive.read(param_member)
        raw_names, names = parse_paramnames(param_payload)
        lookup = {name: index for index, name in enumerate(names)}
        required = ["omegabh2", "omegam", "H0", "rstar", "thetastar"]
        missing = [name for name in required if name not in lookup]
        if missing:
            raise RuntimeError(f"Registered derived columns are absent: {missing}; available={names}")

        vector_blocks: list[np.ndarray] = []
        weight_blocks: list[np.ndarray] = []
        chain_audit: list[dict[str, Any]] = []
        reference_sample: dict[str, Any] | None = None
        for member in chain_members:
            payload = archive.read(member)
            array = np.loadtxt(io.BytesIO(payload), dtype=float, ndmin=2)
            if array.shape[1] != len(names) + 2:
                raise RuntimeError(
                    f"{member} has {array.shape[1]} columns; expected {len(names) + 2}"
                )
            weights = array[:, 0]
            parameters = array[:, 2:]
            theta_star = parameters[:, lookup["thetastar"]] / 100.0
            r_star = parameters[:, lookup["rstar"]]
            omega_m = parameters[:, lookup["omegam"]]
            h0 = parameters[:, lookup["H0"]]
            omega_b_h2 = parameters[:, lookup["omegabh2"]]
            l_a = np.pi / theta_star
            distance_mpc = r_star / theta_star
            shift = np.sqrt(omega_m) * h0 * distance_mpc / C_KM_S
            vectors = np.column_stack([shift, l_a, omega_b_h2])
            keep = (
                np.isfinite(weights)
                & (weights > 0)
                & np.all(np.isfinite(vectors), axis=1)
                & np.all(vectors > 0, axis=1)
            )
            vector_blocks.append(vectors[keep])
            weight_blocks.append(weights[keep])
            finite_loglike = np.isfinite(array[:, 1]) & np.all(np.isfinite(vectors), axis=1)
            if np.any(finite_loglike):
                candidates = np.flatnonzero(finite_loglike)
                row_index = int(candidates[np.argmin(array[candidates, 1])])
                candidate = {
                    "member": member,
                    "row_index_zero_based": row_index,
                    "minus_loglike": float(array[row_index, 1]),
                    "ombh2": float(parameters[row_index, lookup["omegabh2"]]),
                    "omch2": float(parameters[row_index, lookup["omegach2"]]),
                    "tau": float(parameters[row_index, lookup["tau"]]),
                    "logA": float(parameters[row_index, lookup["logA"]]),
                    "ns": float(parameters[row_index, lookup["ns"]]),
                    "H0": float(parameters[row_index, lookup["H0"]]),
                    "chain_vector": vectors[row_index].tolist(),
                }
                if reference_sample is None or candidate["minus_loglike"] < reference_sample["minus_loglike"]:
                    reference_sample = candidate
            chain_audit.append(
                {
                    "member": member,
                    "rows": int(len(array)),
                    "retained_rows": int(np.sum(keep)),
                    "sha256": sha256_bytes(payload),
                }
            )

    vectors = np.vstack(vector_blocks)
    weights = np.concatenate(weight_blocks)
    if reference_sample is None:
        raise RuntimeError("No finite reference sample was found")
    metadata = {
        "paramnames_member": param_member,
        "paramnames_sha256": sha256_bytes(param_payload),
        "raw_parameter_names": raw_names,
        "normalised_parameter_names": names,
        "chain_members": chain_audit,
        "maximum_posterior_reference_sample": reference_sample,
    }
    return vectors, weights, metadata


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chains", type=Path, default=DEFAULT_ZIP)
    parser.add_argument("--out-dir", type=Path, default=RESULT_DIR)
    args = parser.parse_args()
    if not PREREG_MD.exists() or not PREREG_JSON.exists():
        raise FileNotFoundError("The dated CMB preregistration is required")
    if not args.chains.exists():
        raise FileNotFoundError(args.chains)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    vectors, weights, chain_metadata = load_vectors(args.chains)
    mean, covariance, effective_n = weighted_moments(vectors, weights)
    eigenvalues = np.linalg.eigvalsh(covariance)
    if not np.all(eigenvalues > 0):
        raise RuntimeError(f"Distance-prior covariance is not positive definite: {eigenvalues}")
    if not (1.5 < mean[0] < 2.0 and 250 < mean[1] < 350 and 0.015 < mean[2] < 0.03):
        raise RuntimeError(f"Distance-prior mean fails broad physical checks: {mean}")

    labels = ["R", "l_A", "omega_b_h2"]
    correlation = covariance / np.sqrt(np.outer(np.diag(covariance), np.diag(covariance)))
    prior = {
        "analysis": PREFIX,
        "date": DATE,
        "vector_order": labels,
        "means": {label: float(value) for label, value in zip(labels, mean)},
        "cov": covariance.tolist(),
        "correlation": correlation.tolist(),
        "posterior_rows": int(len(vectors)),
        "posterior_weight_sum": float(np.sum(weights)),
        "posterior_effective_sample_size": effective_n,
        "source": {
            "path": str(args.chains),
            "bytes": args.chains.stat().st_size,
            "sha256": sha256_file(args.chains),
            **chain_metadata,
        },
        "preregistration": {
            "markdown_sha256": sha256_file(PREREG_MD),
            "json_sha256": sha256_file(PREREG_JSON),
        },
        "claim_boundary": "Valid only for models retaining standard early physics and recombination.",
    }
    prior_path = args.out_dir / f"{PREFIX}_{DATE}.json"
    prior_path.write_text(json.dumps(prior, indent=2) + "\n", encoding="utf-8")

    summary_rows = [
        {
            "parameter": label,
            "mean": float(mean[index]),
            "sigma": float(np.sqrt(covariance[index, index])),
        }
        for index, label in enumerate(labels)
    ]
    summary_path = args.out_dir / f"{PREFIX}_summary_{DATE}.csv"
    write_csv(summary_path, summary_rows)

    report = [
        "# Planck PR3 Three-variable Distance Prior",
        "",
        f"Date: 18 July 2026",
        "",
        "## Result",
        "",
        f"- retained posterior rows: `{len(vectors)}`",
        f"- weighted effective sample size: `{effective_n:.6g}`",
        f"- source SHA-256: `{prior['source']['sha256']}`",
        "",
        "| parameter | mean | sigma |",
        "| --- | ---: | ---: |",
        *[
            f"| {row['parameter']} | {row['mean']:.10g} | {row['sigma']:.10g} |"
            for row in summary_rows
        ],
        "",
        "## Correlation Matrix",
        "",
        "```text",
        np.array2string(correlation, precision=8),
        "```",
        "",
        "## Claim Boundary",
        "",
        "This posterior compression is a convenience likelihood for models with standard early physics and recombination. It is not valid evidence for an SU(2), SU2R or FR branch that changes photon thermodynamics, radiation density, sound speed or recombination; those branches require the full CMB likelihood.",
    ]
    report_path = args.out_dir / f"{PREFIX}_report_{DATE}.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Saved prior: {prior_path}")
    print(f"Saved report: {report_path}")


if __name__ == "__main__":
    main()
