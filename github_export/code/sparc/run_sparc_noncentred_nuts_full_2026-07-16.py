#!/usr/bin/env python3
r"""Run the separately preregistered 157-galaxy SPARC NUTS posterior.

The physical model and non-centred parameterisation are imported from the
likelihood-equivalent pilot programme. This runner fixes the full filtered
all-Q2 sample, larger sampling budget, full-run filenames, and confirmatory
convergence boundary before sampling begins.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import os
import platform
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

os.environ.setdefault("PYTENSOR_FLAGS", "cxx=")

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm
import nutpie

import fit_sparc_hierarchical_map as fit
import sample_sparc_hierarchical_posterior as legacy


ROOT = Path(__file__).resolve().parents[3]
TASK_ROOT = Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu")
OUTPUTS = TASK_ROOT / "outputs"
PILOT_PROGRAMME = Path(__file__).with_name("run_sparc_noncentred_nuts_pilot_2026-07-16.py")
FILTERED_INPUT_DIR = TASK_ROOT / "work" / "sparc_filtered_persistent_negative"
DEFAULT_SAMPLE = FILTERED_INPUT_DIR / "sparc_galaxy_sample_without_persistent_negative_high_leverage.csv"
DEFAULT_POINTS = FILTERED_INPUT_DIR / "sparc_rotation_points_without_persistent_negative_high_leverage.csv"
DEFAULT_MAP_DIR = (
    ROOT
    / "plamb_runs"
    / "diagnostics"
    / "sparc_hierarchical_map"
    / "optical_depth_minus_persistent_negative_20260714"
)
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "sparc_noncentred_nuts_full_20260716"
SUPPORTED_MODELS = ("RAR", "PLAMB_OPTICAL_DEPTH_KAPPA_P")
RUN_DATE = "2026-07-16"


def load_pilot_module() -> Any:
    spec = importlib.util.spec_from_file_location("sparc_noncentred_nuts_pilot", PILOT_PROGRAMME)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load {PILOT_PROGRAMME}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


pilot = load_pilot_module()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def full_sample_inventory(prep: fit.Prepared) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for galaxy in prep.galaxies:
        mask = prep.data.galaxy == galaxy
        distance_method = int(round(float(np.median(prep.data.distance_method[mask]))))
        inclination = float(np.median(prep.data.inc_deg[mask]))
        rows.append(
            {
                "galaxy": galaxy,
                "distance_method": distance_method,
                "distance_group": "hubble" if distance_method == 1 else "non_hubble",
                "median_inclination_deg": inclination,
                "inclination_group": "high_inc" if inclination >= 60.0 else "low_inc",
                "n_points": int(np.sum(mask)),
                "included": True,
            }
        )
    return rows


def preregistration_markdown(config: dict[str, Any]) -> str:
    gates = config["gates"]
    sampling = config["sampling"]
    return "\n".join(
        [
            "# SPARC Full Filtered All-Q2 Non-centred NUTS Preregistration",
            "",
            f"- **Written before full sampling:** `{config['written_before_sampling_utc']}`",
            f"- **Galaxies:** `{config['full_filtered_galaxies']}`",
            f"- **Rotation points:** `{config['full_filtered_points']}`",
            f"- **Models:** `{', '.join(config['models'])}`",
            "",
            "## Fixed Scope",
            "",
            "Every galaxy remaining after the previously documented six-galaxy persistent-negative-family exclusion is included. There is no further selection, outcome ranking, stress-subset filtering or residual-based removal.",
            "",
            "The July 14 likelihood, fixed H0, distance treatment, velocity-error floor, stellar mass-to-light priors and physical parameter bounds are unchanged. The verified non-centred coordinates are",
            "",
            r"\[",
            r"\begin{aligned}",
            r"z_{D,g}        &\sim \mathcal{N}(0,1), & \log D_g        &= \mu_{D,g} + \sigma_{D,g}z_{D,g}, \\",
            r"z_{\eta,g}     &\sim \mathcal{N}(0,1), & \log\eta_g     &= \sigma_{\log\eta}z_{\eta,g}, \\",
            r"z_{\Upsilon d} &\sim \mathcal{N}(0,1), & \log\Upsilon_d &= \log(0.5)+\sigma_{\Upsilon}z_{\Upsilon d}, \\",
            r"z_{\Upsilon b} &\sim \mathcal{N}(0,1), & \log\Upsilon_b &= \log(0.7)+\sigma_{\Upsilon}z_{\Upsilon b}.",
            r"\end{aligned}",
            r"\]",
            "",
            "The same four-standard-deviation nuisance truncation is retained.",
            "",
            "## Sampling Budget",
            "",
            f"- four chains, `{sampling['tune_per_chain']}` warm-up and `{sampling['draws_per_chain']}` retained draws per chain;",
            f"- target acceptance `{sampling['target_accept']}`;",
            f"- maximum tree depth `{sampling['max_tree_depth']}`; and",
            f"- fixed seed `{sampling['seed']}`.",
            "",
            "## Locked Gates",
            "",
            f"- likelihood equivalence within `{gates['prediction_abs_tolerance_kms']}` km/s;",
            "- zero post-warm-up divergences;",
            f"- maximum rank-normalised split R-hat `<= {gates['rhat_target']}`;",
            f"- absolute R-hat failure ceiling `{gates['rhat_absolute_ceiling']}`;",
            f"- minimum bulk and tail ESS `>= {gates['minimum_bulk_and_tail_ess']}`; and",
            f"- minimum chain E-BFMI `>= {gates['minimum_bfmi']}`.",
            "",
            "A strict pass clears the convergence gate only. It is not a full-sample PLAMB win, a Bayes factor, or held-out predictive validation. Model promotion remains blocked until galaxy-held-out comparison is completed.",
            "",
        ]
    )


def write_report(
    path: Path,
    config: dict[str, Any],
    equivalence_rows: list[dict[str, Any]],
    gates: list[dict[str, Any]],
    global_rows: list[dict[str, Any]],
) -> None:
    lines = [
        "# SPARC Full Filtered All-Q2 Non-centred NUTS Report",
        "",
            f"- **Completed:** `{utc_now()}`",
            f"- **Galaxies:** `{config['full_filtered_galaxies']}`",
            f"- **Rotation points:** `{config['full_filtered_points']}`",
            f"- **Chains:** `{config['sampling']['chains']}`",
        "",
        "## Convergence Decision",
        "",
    ]
    if all(bool(gate["ready_for_full_run"]) for gate in gates):
        lines.append("Both full posteriors pass every strict preregistered convergence gate.")
    elif all(gate["status"] != "FAIL" for gate in gates):
        lines.append("Neither model fails an absolute gate, but at least one misses the strict R-hat target; extend sampling before interpretation.")
    else:
        lines.append("At least one full posterior fails a preregistered convergence gate and must not be interpreted as stable.")
    lines.extend(
        [
            "",
            "Passing this report clears only the full-sample sampler-convergence gate. It does not establish out-of-sample superiority over RAR.",
            "",
            "## Gate Results",
            "",
            "| Model | Status | Max R-hat | Min bulk ESS | Min tail ESS | Divergences | Min E-BFMI | Max depth hits |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for gate in gates:
        lines.append(
            f"| {gate['model']} | {gate['status']} | {gate['max_rhat']:.5f} | "
            f"{gate['min_bulk_ess']:.1f} | {gate['min_tail_ess']:.1f} | {gate['divergences']} | "
            f"{gate['minimum_bfmi']:.4f} | {gate['max_tree_depth_hits']} |"
        )
    lines.extend(
        [
            "",
            "## Likelihood Equivalence",
            "",
            "| Model | Maximum velocity difference (km/s) | Chi-squared difference | Pass |",
            "| --- | ---: | ---: | --- |",
        ]
    )
    for row in equivalence_rows:
        lines.append(
            f"| {row['model']} | {float(row['max_abs_prediction_difference_kms']):.3e} | "
            f"{float(row['chi2_data_absolute_difference']):.3e} | {row['pass']} |"
        )
    lines.extend(
        [
            "",
            "## Full Posterior Global Parameters",
            "",
            "These are well-mixed in-sample posterior estimates from chains that pass the R-hat, ESS, E-BFMI and tree-depth criteria but fail the locked zero-divergence criterion. They are not predictive model-comparison scores or a strict convergence pass.",
            "",
            "| Model | Parameter | Mean | SD | 3% HDI | 97% HDI | R-hat | Bulk ESS | Tail ESS |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in global_rows:
        lines.append(
            f"| {row['model']} | {row['parameter']} | {row['mean']:.6g} | {row['sd']:.6g} | "
            f"{row['hdi_3%']:.6g} | {row['hdi_97%']:.6g} | {row['r_hat']:.5f} | "
            f"{row['ess_bulk']:.1f} | {row['ess_tail']:.1f} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "- A strict pass replaces the old random-walk convergence uncertainty for the same filtered all-Q2 likelihood.",
            "- It does not validate the previous prior-redrawn stress scores.",
            "- It does not replace galaxy-held-out or grouped cross-validation.",
            "- The standing claim remains `subset wins, not a full-sample win` until predictive validation is passed.",
            "",
            "## Reproducibility",
            "",
            f"- PyMC: `{config['environment']['pymc']}`",
            f"- nutpie: `{config['environment']['nutpie']}`",
            f"- ArviZ: `{config['environment']['arviz']}`",
            f"- Python: `{config['environment']['python']}`",
            f"- sampling seed: `{config['sampling']['seed']}`",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Full 157-galaxy non-centred NUTS posterior for filtered all-Q2 SPARC.")
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--points", type=Path, default=DEFAULT_POINTS)
    parser.add_argument("--map-dir", type=Path, default=DEFAULT_MAP_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--models", default=",".join(SUPPORTED_MODELS))
    parser.add_argument("--seed", type=int, default=16072611)
    parser.add_argument("--chains", type=int, default=4)
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--tune", type=int, default=2000)
    parser.add_argument("--draws", type=int, default=2000)
    parser.add_argument("--target-accept", type=float, default=0.90)
    parser.add_argument("--max-tree-depth", type=int, default=12)
    parser.add_argument("--sigma-ln-ml", type=float, default=0.25)
    parser.add_argument("--sigma-global-ml", type=float, default=0.55)
    parser.add_argument("--distance-floor-frac", type=float, default=0.03)
    parser.add_argument("--hubble-prior-mode", choices=["model_h0_rescale", "published"], default="model_h0_rescale")
    parser.add_argument("--err-floor-kms", type=float, default=3.0)
    parser.add_argument("--rhat-target", type=float, default=1.01)
    parser.add_argument("--rhat-ceiling", type=float, default=1.05)
    parser.add_argument("--min-ess", type=float, default=400.0)
    parser.add_argument("--min-bfmi", type=float, default=0.30)
    parser.add_argument("--expected-galaxies", type=int, default=157)
    parser.add_argument("--expected-points", type=int, default=3001)
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--analyse-existing", action="store_true")
    parser.add_argument("--copy-to-outputs", action=argparse.BooleanOptionalAction, default=True)
    return parser.parse_args()


def read_equivalence(path: Path) -> list[dict[str, Any]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def analyse_existing(args: argparse.Namespace, models: list[str]) -> int:
    prereg_path = args.out_dir / "sparc_noncentred_nuts_full_preregistration.json"
    equivalence_path = args.out_dir / "sparc_noncentred_nuts_full_likelihood_equivalence.csv"
    if not prereg_path.exists() or not equivalence_path.exists():
        raise FileNotFoundError("Analysis-only mode requires existing full-run preregistration and equivalence files")
    prereg = json.loads(prereg_path.read_text(encoding="utf-8"))
    registered_sampling = prereg["sampling"]
    args.chains = int(registered_sampling["chains"])
    args.tune = int(registered_sampling["tune_per_chain"])
    args.draws = int(registered_sampling["draws_per_chain"])
    args.max_tree_depth = int(registered_sampling["max_tree_depth"])
    registered_gates = prereg["gates"]
    args.rhat_target = float(registered_gates["rhat_target"])
    args.rhat_ceiling = float(registered_gates["rhat_absolute_ceiling"])
    args.min_ess = float(registered_gates["minimum_bulk_and_tail_ess"])
    args.min_bfmi = float(registered_gates["minimum_bfmi"])
    equivalence_rows = read_equivalence(equivalence_path)
    equivalence_by_model = {row["model"]: row for row in equivalence_rows}
    gates: list[dict[str, Any]] = []
    global_rows: list[dict[str, Any]] = []
    for model_name in models:
        trace_path = args.out_dir / f"sparc_noncentred_nuts_full_{model_name}_trace.nc"
        idata = az.from_netcdf(trace_path)
        gate_path = args.out_dir / f"sparc_noncentred_nuts_full_{model_name}_gate.json"
        old_gate = json.loads(gate_path.read_text(encoding="utf-8")) if gate_path.exists() else {}
        summary, gate = pilot.diagnose(
            idata,
            model_name,
            args,
            float(old_gate.get("elapsed_seconds", 0.0)),
            str(equivalence_by_model[model_name]["pass"]).strip().lower() in {"true", "1", "yes"},
        )
        gate["compile_seconds"] = float(old_gate.get("compile_seconds", 0.0))
        summary.to_csv(args.out_dir / f"sparc_noncentred_nuts_full_{model_name}_diagnostics.csv", index=False)
        gate_path.write_text(json.dumps(gate, indent=2), encoding="utf-8")
        gates.append(gate)
        global_rows.extend(pilot.global_summary_rows(model_name, idata))
    gate_table = args.out_dir / "sparc_noncentred_nuts_full_gates.csv"
    global_table = args.out_dir / "sparc_noncentred_nuts_full_global_summary.csv"
    report_path = args.out_dir / "sparc_noncentred_nuts_full_report.md"
    pilot.write_csv(gate_table, gates)
    pilot.write_csv(global_table, global_rows)
    write_report(report_path, prereg, equivalence_rows, gates, global_rows)
    pilot.create_manifest(args.out_dir)
    if args.copy_to_outputs:
        copy_compact_outputs(args.out_dir)
    print(f"Reanalysed existing full traces: {report_path}", flush=True)
    return 0


def copy_compact_outputs(out_dir: Path) -> None:
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    names = [
        "sparc_noncentred_nuts_full_preregistration.md",
        "sparc_noncentred_nuts_full_sample_inventory.csv",
        "sparc_noncentred_nuts_full_likelihood_equivalence.csv",
        "sparc_noncentred_nuts_full_gates.csv",
        "sparc_noncentred_nuts_full_global_summary.csv",
        "sparc_noncentred_nuts_full_report.md",
    ]
    for name in names:
        path = out_dir / name
        shutil.copy2(path, OUTPUTS / f"{path.stem}_{RUN_DATE}{path.suffix}")


def main() -> int:
    args = parse_args()
    models = pilot.parse_models(args.models)
    if args.analyse_existing:
        return analyse_existing(args, models)
    if args.chains != 4 and not args.preflight_only:
        raise ValueError("The registered full run requires exactly four chains")
    required = [args.sample, args.points, args.map_dir / "sparc_hierarchical_map_summary.csv", PILOT_PROGRAMME]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required inputs: {missing}")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    map_rows = legacy.load_map_rows(args.map_dir)
    h0_by_model = {model: legacy.parse_float(map_rows[("all_Q2", model)]["H0_trial"]) for model in models}
    if len(set(h0_by_model.values())) != 1:
        raise ValueError(f"The full run requires a common fixed H0; found {h0_by_model}")
    h0 = next(iter(h0_by_model.values()))
    split_config = next(item for item in fit.split_configs() if item["split"] == "all_Q2")
    prep = fit.prepare(
        args.sample,
        args.points,
        int(split_config["quality_max"]),
        str(split_config["distance_method"]),
        args.err_floor_kms,
        [h0],
        args.sigma_ln_ml,
        args.distance_floor_frac,
        args.hubble_prior_mode,
    )
    if len(prep.galaxies) != args.expected_galaxies or len(prep.data.vobs) != args.expected_points:
        raise ValueError(
            f"Full-sample cardinality changed: got {len(prep.galaxies)} galaxies/{len(prep.data.vobs)} points, "
            f"expected {args.expected_galaxies}/{args.expected_points}"
        )
    inventory_path = args.out_dir / "sparc_noncentred_nuts_full_sample_inventory.csv"
    pilot.write_csv(inventory_path, full_sample_inventory(prep))

    environment = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "pymc": pm.__version__,
        "pytensor": __import__("pytensor").__version__,
        "nutpie": nutpie.__version__,
        "arviz": az.__version__,
    }
    prereg = {
        "date": RUN_DATE,
        "written_before_sampling_utc": utc_now(),
        "analysis": "full_filtered_all_Q2_noncentred_nuts_posterior",
        "scientific_use": "full posterior convergence only; no model promotion without held-out validation",
        "models": models,
        "sample": str(args.sample),
        "points": str(args.points),
        "map_dir": str(args.map_dir),
        "input_sha256": {
            "sample": pilot.sha256(args.sample),
            "points": pilot.sha256(args.points),
            "map_summary": pilot.sha256(args.map_dir / "sparc_hierarchical_map_summary.csv"),
            "full_runner": pilot.sha256(Path(__file__)),
            "verified_model_programme": pilot.sha256(PILOT_PROGRAMME),
        },
        "full_filtered_galaxies": len(prep.galaxies),
        "full_filtered_points": len(prep.data.vobs),
        "sample_rule": "all galaxies in the fixed filtered all_Q2 inputs; no additional selection",
        "fixed_h0": h0,
        "priors": {
            "log_ydisk": "Normal(log(0.5), 0.55), truncated to log(0.05)..log(1.50)",
            "log_ybul": "Normal(log(0.7), 0.55), truncated to log(0.05)..log(1.80)",
            "distance_nuisance": "Normal(mu_logd(H0), sigma_logd), truncated at +/-4 sigma",
            "stellar_ml_nuisance": "Normal(0, 0.25), truncated at +/-4 sigma",
            "RAR_log10_gdagger": "Uniform(-11.6, -9.2)",
            "PLAMB_log10_kappa": "Uniform(-0.50, 0.50)",
            "PLAMB_bridge_exponent": "Uniform(0.15, 1.20)",
        },
        "sampling": {
            "sampler": "nutpie NUTS compiled from the likelihood-equivalent PyMC model",
            "chains": args.chains,
            "cores": args.cores,
            "tune_per_chain": args.tune,
            "draws_per_chain": args.draws,
            "target_accept": args.target_accept,
            "max_tree_depth": args.max_tree_depth,
            "seed": args.seed,
        },
        "gates": {
            "prediction_abs_tolerance_kms": 1e-8,
            "rhat_target": args.rhat_target,
            "rhat_absolute_ceiling": args.rhat_ceiling,
            "minimum_bulk_and_tail_ess": args.min_ess,
            "maximum_divergences": 0,
            "minimum_bfmi": args.min_bfmi,
        },
        "environment": environment,
    }
    prereg_json = args.out_dir / "sparc_noncentred_nuts_full_preregistration.json"
    prereg_md = args.out_dir / "sparc_noncentred_nuts_full_preregistration.md"
    prereg_json.write_text(json.dumps(prereg, indent=2), encoding="utf-8")
    prereg_md.write_text(preregistration_markdown(prereg), encoding="utf-8")
    print(f"Preregistered full sample: galaxies={len(prep.galaxies)} points={len(prep.data.vobs)}", flush=True)

    equivalence_rows: list[dict[str, Any]] = []
    built_models: dict[str, pm.Model] = {}
    for model_name in models:
        pymc_model = pilot.build_model(prep, model_name, h0, args.sigma_global_ml)
        equivalence = pilot.equivalence_check(prep, model_name, h0, args.sigma_global_ml)
        equivalence_rows.append(equivalence)
        built_models[model_name] = pymc_model
        print(
            f"{model_name} equivalence max_abs={equivalence['max_abs_prediction_difference_kms']:.3e} "
            f"chi2_diff={equivalence['chi2_data_absolute_difference']:.3e} pass={equivalence['pass']}",
            flush=True,
        )
        if not equivalence["pass"]:
            raise RuntimeError(f"Full likelihood-equivalence check failed for {model_name}")
    equivalence_path = args.out_dir / "sparc_noncentred_nuts_full_likelihood_equivalence.csv"
    pilot.write_csv(equivalence_path, equivalence_rows)
    if args.preflight_only:
        pilot.create_manifest(args.out_dir)
        print(f"Full-run preflight complete: {args.out_dir}", flush=True)
        return 0

    gates: list[dict[str, Any]] = []
    global_rows: list[dict[str, Any]] = []
    for model_index, model_name in enumerate(models):
        print(
            f"=== Full sampling {model_name}: chains={args.chains} tune={args.tune} draws={args.draws} "
            f"galaxies={len(prep.galaxies)} points={len(prep.data.vobs)} ===",
            flush=True,
        )
        compile_start = time.perf_counter()
        compiled = nutpie.compile_pymc_model(built_models[model_name], backend="numba")
        compile_seconds = time.perf_counter() - compile_start
        sample_start = time.perf_counter()
        idata = nutpie.sample(
            compiled,
            draws=args.draws,
            tune=args.tune,
            chains=args.chains,
            cores=args.cores,
            seed=args.seed + model_index * 100003,
            target_accept=args.target_accept,
            maxdepth=args.max_tree_depth,
            save_warmup=True,
            progress_bar=True,
        )
        sample_seconds = time.perf_counter() - sample_start
        trace_path = args.out_dir / f"sparc_noncentred_nuts_full_{model_name}_trace.nc"
        pilot.sanitise_netcdf_attributes(idata)
        idata.to_netcdf(trace_path)
        summary, gate = pilot.diagnose(
            idata,
            model_name,
            args,
            sample_seconds,
            bool(equivalence_rows[model_index]["pass"]),
        )
        gate["compile_seconds"] = compile_seconds
        summary.to_csv(args.out_dir / f"sparc_noncentred_nuts_full_{model_name}_diagnostics.csv", index=False)
        (args.out_dir / f"sparc_noncentred_nuts_full_{model_name}_gate.json").write_text(
            json.dumps(gate, indent=2),
            encoding="utf-8",
        )
        gates.append(gate)
        global_rows.extend(pilot.global_summary_rows(model_name, idata))
        print(
            f"{model_name} status={gate['status']} max_rhat={gate['max_rhat']:.5f} "
            f"min_bulk_ess={gate['min_bulk_ess']:.1f} min_tail_ess={gate['min_tail_ess']:.1f} "
            f"divergences={gate['divergences']} min_bfmi={gate['minimum_bfmi']:.4f}",
            flush=True,
        )

    gate_table = args.out_dir / "sparc_noncentred_nuts_full_gates.csv"
    global_table = args.out_dir / "sparc_noncentred_nuts_full_global_summary.csv"
    report_path = args.out_dir / "sparc_noncentred_nuts_full_report.md"
    pilot.write_csv(gate_table, gates)
    pilot.write_csv(global_table, global_rows)
    write_report(report_path, prereg, equivalence_rows, gates, global_rows)
    pilot.create_manifest(args.out_dir)
    if args.copy_to_outputs:
        copy_compact_outputs(args.out_dir)
    print(f"Saved full report: {report_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
