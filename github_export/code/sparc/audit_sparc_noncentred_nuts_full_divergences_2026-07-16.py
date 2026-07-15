#!/usr/bin/env python3
"""Audit every divergent draw in the preregistered full SPARC NUTS run."""

from __future__ import annotations

import argparse
import hashlib
import shutil
from collections import Counter
from pathlib import Path
from typing import Any

import arviz as az
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
TASK_ROOT = Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu")
DEFAULT_RUN_DIR = ROOT / "plamb_runs" / "diagnostics" / "sparc_noncentred_nuts_full_20260716"
DEFAULT_OUTPUTS = TASK_ROOT / "outputs"
PREFIX = "sparc_noncentred_nuts_full"
MODELS = {
    "RAR": {
        "trace": f"{PREFIX}_RAR_trace.nc",
        "globals": ("ydisk", "ybul", "log10_gdagger"),
    },
    "PLAMB_OPTICAL_DEPTH_KAPPA_P": {
        "trace": f"{PREFIX}_PLAMB_OPTICAL_DEPTH_KAPPA_P_trace.nc",
        "globals": ("ydisk", "ybul", "log10_kappa", "bridge_exponent"),
    },
}


def scalar_at(dataset: Any, name: str, chain: int, draw: int) -> Any:
    value = np.asarray(dataset[name])[chain, draw]
    return value.item() if value.ndim == 0 else value


def finite_or_text(value: Any) -> Any:
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return str(value)
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def refresh_manifest(run_dir: Path) -> Path:
    files = sorted(
        path for path in run_dir.iterdir()
        if path.is_file() and path.name != "manifest.csv"
    )
    rows = [
        {"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in files
    ]
    manifest_path = run_dir / "manifest.csv"
    pd.DataFrame(rows).to_csv(manifest_path, index=False)
    return manifest_path


def audit_model(
    model: str,
    config: dict[str, Any],
    run_dir: Path,
    galaxies: list[str],
) -> tuple[list[dict[str, Any]], dict[str, Any], list[dict[str, Any]]]:
    trace_path = run_dir / config["trace"]
    idata = az.from_netcdf(trace_path)
    posterior = idata.posterior
    stats = idata.sample_stats
    divergent = np.asarray(stats["diverging"], dtype=bool)
    positions = np.argwhere(divergent)
    z_logd = np.asarray(posterior["z_logd"])
    z_logeta = np.asarray(posterior["z_logeta"])
    max_abs_z_logd = np.max(np.abs(z_logd), axis=2)
    max_abs_z_logeta = np.max(np.abs(z_logeta), axis=2)
    latent_boundary_counts = (
        np.sum(np.abs(z_logd) >= 3.8, axis=2)
        + np.sum(np.abs(z_logeta) >= 3.8, axis=2)
    )

    rows: list[dict[str, Any]] = []
    for chain, draw in positions:
        z_d = z_logd[chain, draw]
        z_eta = z_logeta[chain, draw]
        i_d = int(np.argmax(np.abs(z_d)))
        i_eta = int(np.argmax(np.abs(z_eta)))
        row: dict[str, Any] = {
            "model": model,
            "chain": int(chain) + 1,
            "draw": int(draw) + 1,
            "depth": int(scalar_at(stats, "depth", chain, draw)),
            "n_steps": int(scalar_at(stats, "n_steps", chain, draw)),
            "logp": float(scalar_at(stats, "logp", chain, draw)),
            "energy_error": float(scalar_at(stats, "energy_error", chain, draw)),
            "divergence_energy_error": float(scalar_at(stats, "divergence_energy_error", chain, draw)),
            "divergence_message": str(scalar_at(stats, "divergence_message", chain, draw)),
            "max_abs_z_logd": float(abs(z_d[i_d])),
            "max_abs_z_logd_galaxy": galaxies[i_d],
            "z_logd_at_max": float(z_d[i_d]),
            "max_abs_z_logeta": float(abs(z_eta[i_eta])),
            "max_abs_z_logeta_galaxy": galaxies[i_eta],
            "z_logeta_at_max": float(z_eta[i_eta]),
            "n_latents_abs_ge_3p5": int(np.sum(np.abs(z_d) >= 3.5) + np.sum(np.abs(z_eta) >= 3.5)),
            "n_latents_abs_ge_3p8": int(np.sum(np.abs(z_d) >= 3.8) + np.sum(np.abs(z_eta) >= 3.8)),
        }
        for name in config["globals"]:
            row[name] = float(scalar_at(posterior, name, chain, draw))
        rows.append(row)

    total_draws = int(divergent.size)
    summary = {
        "model": model,
        "retained_draws": total_draws,
        "divergences": int(divergent.sum()),
        "divergence_rate": float(divergent.mean()),
        "divergences_by_chain": ";".join(str(int(value)) for value in divergent.sum(axis=1)),
        "divergent_max_abs_z_logd_min": float(np.min(max_abs_z_logd[divergent])),
        "divergent_max_abs_z_logd_max": float(np.max(max_abs_z_logd[divergent])),
        "all_max_abs_z_logd_p99": float(np.quantile(max_abs_z_logd, 0.99)),
        "divergent_max_abs_z_logeta_min": float(np.min(max_abs_z_logeta[divergent])),
        "divergent_max_abs_z_logeta_max": float(np.max(max_abs_z_logeta[divergent])),
        "all_max_abs_z_logeta_p99": float(np.quantile(max_abs_z_logeta, 0.99)),
        "divergent_draws_with_latent_abs_ge_3p8": int(
            np.sum((max_abs_z_logd[divergent] >= 3.8) | (max_abs_z_logeta[divergent] >= 3.8))
        ),
        "all_draws_with_latent_abs_ge_3p8_fraction": float(np.mean(latent_boundary_counts > 0)),
        "mean_latents_abs_ge_3p8_divergent": float(np.mean(latent_boundary_counts[divergent])),
        "mean_latents_abs_ge_3p8_nondivergent": float(np.mean(latent_boundary_counts[~divergent])),
        "trace_file": trace_path.name,
        "trace_sha256": sha256(trace_path),
    }

    occupancy_rows: list[dict[str, Any]] = []
    for latent_name, values in (("z_logd", z_logd), ("z_logeta", z_logeta)):
        boundary = np.abs(values) >= 3.8
        for galaxy_index, galaxy in enumerate(galaxies):
            galaxy_values = values[:, :, galaxy_index]
            galaxy_boundary = boundary[:, :, galaxy_index]
            occupancy_rows.append(
                {
                    "model": model,
                    "latent": latent_name,
                    "galaxy": galaxy,
                    "fraction_abs_ge_3p8": float(np.mean(galaxy_boundary)),
                    "fraction_abs_ge_3p8_divergent": float(np.mean(galaxy_boundary[divergent])),
                    "fraction_abs_ge_3p8_nondivergent": float(np.mean(galaxy_boundary[~divergent])),
                    "posterior_mean": float(np.mean(galaxy_values)),
                    "posterior_sd": float(np.std(galaxy_values, ddof=1)),
                    "posterior_q01": float(np.quantile(galaxy_values, 0.01)),
                    "posterior_q99": float(np.quantile(galaxy_values, 0.99)),
                }
            )
    return rows, summary, occupancy_rows


def markdown_report(
    rows: list[dict[str, Any]],
    summaries: list[dict[str, Any]],
    occupancy_rows: list[dict[str, Any]],
) -> str:
    lines = [
        "# SPARC Full Non-centred NUTS Divergence Audit",
        "",
        "This post-run audit preserves the preregistered zero-divergence gate. It locates every divergent transition; it does not redefine convergence after seeing the outcome.",
        "",
        "## Counts",
        "",
        "| Model | Retained draws | Divergences | Rate | By chain | All draws near boundary |",
        "|---|---:|---:|---:|---|---:|",
    ]
    for item in summaries:
        lines.append(
            f"| {item['model']} | {item['retained_draws']:,} | {item['divergences']} | "
            f"{item['divergence_rate']:.6f} | {item['divergences_by_chain']} | "
            f"{item['all_draws_with_latent_abs_ge_3p8_fraction']:.3f} |"
        )

    lines.extend(
        [
            "",
            "A transition is labelled near the retained four-standard-deviation nuisance boundary when either maximum latent absolute value is at least 3.8.",
            "",
            "## Transition Locations",
            "",
            "| Model | Chain | Draw | Depth | Divergence energy error | Max abs(z_D) galaxy | Max abs(z_eta) galaxy | Latents >= 3.8 |",
            "|---|---:|---:|---:|---:|---|---|---:|",
        ]
    )
    for row in rows:
        lines.append(
            f"| {row['model']} | {row['chain']} | {row['draw']} | {row['depth']} | "
            f"{row['divergence_energy_error']:.3f} | "
            f"{row['max_abs_z_logd_galaxy']} ({row['z_logd_at_max']:.3f}) | "
            f"{row['max_abs_z_logeta_galaxy']} ({row['z_logeta_at_max']:.3f}) | "
            f"{row['n_latents_abs_ge_3p8']} |"
        )

    repeated = Counter(
        (row["model"], row[key])
        for row in rows
        for key in ("max_abs_z_logd_galaxy", "max_abs_z_logeta_galaxy")
    )
    repeated_text = ", ".join(
        f"{model}: {galaxy} ({count})"
        for (model, galaxy), count in sorted(repeated.items())
        if count > 1
    ) or "None"

    near_boundary = sum(int(item["divergent_draws_with_latent_abs_ge_3p8"]) for item in summaries)
    boundary_comparison = "; ".join(
        f"{item['model']}: divergent {item['mean_latents_abs_ge_3p8_divergent']:.2f}, "
        f"non-divergent {item['mean_latents_abs_ge_3p8_nondivergent']:.2f}"
        for item in summaries
    )
    top_occupancy = sorted(
        occupancy_rows,
        key=lambda item: float(item["fraction_abs_ge_3p8"]),
        reverse=True,
    )[:10]
    lines.extend(
        [
            "",
            "## Readout",
            "",
            f"All {near_boundary} divergent transitions lie near the explicit nuisance boundary, but so does every retained draw in both models. The mean numbers of latent variables with absolute value at least 3.8 are {boundary_comparison}. There is therefore no transition-level enrichment of boundary occupancy among the nine divergences.",
            f"Repeated nearest-boundary galaxies across divergent draws: {repeated_text}.",
            "",
            "The truncation is nevertheless active throughout the posterior. The ten highest marginal boundary occupancies are:",
            "",
            "| Model | Latent | Galaxy | Fraction with abs(z) >= 3.8 | Mean z |",
            "|---|---|---|---:|---:|",
        ]
    )
    for item in top_occupancy:
        lines.append(
            f"| {item['model']} | {item['latent']} | {item['galaxy']} | "
            f"{item['fraction_abs_ge_3p8']:.3f} | {item['posterior_mean']:.3f} |"
        )
    lines.extend(
        [
            "",
            "The divergence rates are small, and the rank-normalised R-hat, effective sample size, E-BFMI and tree-depth gates pass in the main report. Nevertheless, both models fail the separately preregistered zero-divergence requirement. The appropriate rescue is an unchanged-model rerun at higher target acceptance, registered before execution; these draws cannot be relabelled as a strict convergence pass.",
        ]
    )
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--outputs", type=Path, default=DEFAULT_OUTPUTS)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    inventory_path = args.run_dir / f"{PREFIX}_sample_inventory.csv"
    galaxies = pd.read_csv(inventory_path)["galaxy"].astype(str).tolist()
    if len(galaxies) != 157:
        raise ValueError(f"Expected 157 galaxies, found {len(galaxies)}")

    rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    occupancy_rows: list[dict[str, Any]] = []
    for model, config in MODELS.items():
        model_rows, model_summary, model_occupancy = audit_model(model, config, args.run_dir, galaxies)
        rows.extend(model_rows)
        summaries.append(model_summary)
        occupancy_rows.extend(model_occupancy)

    transitions_path = args.run_dir / f"{PREFIX}_divergence_transitions.csv"
    summary_path = args.run_dir / f"{PREFIX}_divergence_summary.csv"
    occupancy_path = args.run_dir / f"{PREFIX}_boundary_occupancy.csv"
    report_path = args.run_dir / f"{PREFIX}_divergence_audit.md"
    pd.DataFrame(rows).to_csv(transitions_path, index=False)
    pd.DataFrame(summaries).to_csv(summary_path, index=False)
    pd.DataFrame(occupancy_rows).to_csv(occupancy_path, index=False)
    report_path.write_text(markdown_report(rows, summaries, occupancy_rows), encoding="utf-8")
    manifest_path = refresh_manifest(args.run_dir)

    args.outputs.mkdir(parents=True, exist_ok=True)
    dated = {
        transitions_path: args.outputs / f"{PREFIX}_divergence_transitions_2026-07-16.csv",
        summary_path: args.outputs / f"{PREFIX}_divergence_summary_2026-07-16.csv",
        occupancy_path: args.outputs / f"{PREFIX}_boundary_occupancy_2026-07-16.csv",
        report_path: args.outputs / f"{PREFIX}_divergence_audit_2026-07-16.md",
    }
    for source, destination in dated.items():
        shutil.copy2(source, destination)
        print(f"Saved: {destination}")
    print(f"Refreshed: {manifest_path}")


if __name__ == "__main__":
    main()
