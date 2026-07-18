#!/usr/bin/env python3
"""Generate and evaluate an exact-LINX design for hierarchical BBN inference."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from scipy.special import ndtr, ndtri
from scipy.stats import qmc


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
LINX_PATH = REPO_ROOT / "plamb_runs" / "deps" / "LINX"
OUTDIR = REPO_ROOT / "plamb_runs" / "diagnostics" / "bbn_lithium_hierarchical_surrogate_20260718"
ARCHIVE = REPO_ROOT / "github_export" / "results" / "2026-07-18" / "bbn_lithium"
REACTIONS = [
    "npdg", "dpHe3g", "ddHe3n", "ddtp", "tpag", "tdan",
    "taLi7g", "He3ntp", "He3dap", "He3aBe7g", "Be7nLi7p", "Li7paa",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["register", "design", "validate"], required=True)
    parser.add_argument("--outdir", type=Path, default=OUTDIR)
    parser.add_argument("--archive-dir", type=Path, default=ARCHIVE)
    parser.add_argument("--linx-path", type=Path, default=LINX_PATH)
    parser.add_argument("--n-design", type=int, default=4800)
    parser.add_argument("--n-train", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=20260718)
    parser.add_argument("--input-points", type=Path)
    parser.add_argument("--output-csv", type=Path)
    parser.add_argument("--max-points", type=int)
    parser.add_argument("--summary-every", type=int, default=25)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def now() -> str:
    return datetime.now(UTC).astimezone().isoformat(timespec="seconds")


def truncated_normal(u: np.ndarray, low: float, high: float) -> np.ndarray:
    return ndtri(ndtr(low) + u * (ndtr(high) - ndtr(low)))


def design_path(args: argparse.Namespace) -> Path:
    return args.outdir / "exact_linx_prior_design.csv"


def make_design(args: argparse.Namespace) -> Path:
    if not 0 < args.n_train < args.n_design:
        raise ValueError("Require 0 < n_train < n_design")
    args.outdir.mkdir(parents=True, exist_ok=True)
    path = design_path(args)
    sampler = qmc.LatinHypercube(d=14, scramble=True, seed=args.seed)
    unit = sampler.random(args.n_design)
    z = np.empty_like(unit)
    z[:, 0] = truncated_normal(unit[:, 0], -5.0, 5.0)
    z[:, 1] = truncated_normal(unit[:, 1], -5.0, 5.0)
    z[:, 2:] = truncated_normal(unit[:, 2:], -4.0, 4.0)
    order = np.random.default_rng(args.seed + 1).permutation(args.n_design)
    split = np.full(args.n_design, "train", dtype=object)
    split[order[args.n_train:]] = "holdout"
    fields = ["point_id", "split", "eta_fac", "tau_seconds"] + [f"q_{name}" for name in REACTIONS]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for i, row in enumerate(z):
            record: dict[str, object] = {
                "point_id": f"lhs_{i:05d}",
                "split": split[i],
                "eta_fac": 1.0 + row[0] * (0.04 / 6.12),
                "tau_seconds": 878.4 + row[1] * 0.5,
            }
            record.update({f"q_{name}": row[j + 2] for j, name in enumerate(REACTIONS)})
            writer.writerow(record)
    return path


def write_registration(args: argparse.Namespace) -> None:
    args.archive_dir.mkdir(parents=True, exist_ok=True)
    path = make_design(args)
    payload = {
        "registered_at": now(),
        "reason": "Direct LINX NUTS failed preposterior gradient smoke tests; use gradient-free surrogate inference with exact validation.",
        "programme": str(SCRIPT_PATH),
        "programme_sha256": hashlib.sha256(SCRIPT_PATH.read_bytes()).hexdigest(),
        "design": str(path),
        "design_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "n_design": args.n_design,
        "n_train": args.n_train,
        "n_holdout": args.n_design - args.n_train,
        "seed": args.seed,
        "network": "key_PRIMAT_2023",
        "features": ["eta_fac", "tau_seconds"] + REACTIONS,
        "priors": {
            "eta_fac": "Normal(1, 0.04/6.12), truncated at +/-5 sigma",
            "tau_seconds": "Normal(878.4, 0.5 s), truncated at +/-5 sigma",
            "rate_pulls": "12 independent Normal(0,1), truncated at +/-4",
            "lithium_depletion": "HalfNormal(0.30 dex) in depletion arm",
        },
        "numerics": {"rtol": 1e-5, "atol": 1e-9, "sampling_nTOp": 150},
        "posterior_arms": ["D+He+CMB predictive Li", "D+He+CMB+Li with depletion"],
        "surrogate_gate": {
            "holdout_rmse_max_observational_sigma": 0.25,
            "exact_validation_q95_abs_error_max_observational_sigma": 0.50,
            "importance_ess_min": 100,
        },
        "claim_boundary": "Homogeneous standard BBN; no matter-antimatter domains or annihilation physics.",
    }
    json_path = args.archive_dir / "bbn_lithium_hierarchical_surrogate_registration_2026-07-18.json"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md_path = args.archive_dir / "bbn_lithium_hierarchical_surrogate_preregistration_2026-07-18.md"
    md_path.write_text(
        """# Hierarchical LINX BBN Surrogate Posterior: Preregistration

Direct reverse-mode NUTS failed before retaining posterior draws because the
implicit Kvaerno solve became non-finite; forward-mode is unsupported by LINX's
`custom_vjp`. The registered fallback is therefore gradient-free.

An exact `key_PRIMAT_2023` Latin-hypercube design covers eta, neutron lifetime
and all 12 rate pulls. A fixed 4,000/800 train/holdout split is made before
evaluation. Surrogate posteriors are not promoted unless holdout and selected
posterior-draw exact-LINX gates pass. Both the D+He predictive-Li arm and the
explicit positive lithium-depletion arm are retained.

This remains homogeneous standard BBN. It is not a matter-antimatter domain or
annihilation calculation.
""",
        encoding="utf-8",
    )
    print(f"Saved design: {path}", flush=True)
    print(f"Saved registration: {md_path}", flush=True)


def load_done(path: Path) -> set[str]:
    if not path.exists():
        return set()
    with path.open("r", newline="", encoding="utf-8") as handle:
        return {row["point_id"] for row in csv.DictReader(handle) if row.get("status") == "ok"}


def evaluate(args: argparse.Namespace, input_path: Path, output_path: Path) -> None:
    sys.path.insert(0, str(args.linx_path.resolve()))
    import jax
    import jax.numpy as jnp
    from linx.abundances import AbundanceModel
    from linx.background import BackgroundModel
    from linx.nuclear import NuclearRates
    import linx.const as const

    jax.config.update("jax_enable_x64", True)
    rows = list(csv.DictReader(input_path.open("r", newline="", encoding="utf-8")))
    if args.max_points is not None:
        rows = rows[: args.max_points]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    done = set() if args.force else load_done(output_path)
    fields = list(rows[0]) + ["status", "message", "elapsed_s", "D_H_1e5", "Yp_mass", "Li7_H_1e10", "run_time"]
    write_header = not output_path.exists() or args.force
    mode = "w" if write_header else "a"

    rates = NuclearRates(nuclear_net="key_PRIMAT_2023")
    if list(rates.reactions_names) != REACTIONS:
        raise RuntimeError(f"Reaction ordering changed: {rates.reactions_names}")
    abundance = AbundanceModel(rates, throw=False)
    background = BackgroundModel()
    t, a, rho_g, rho_nu, rho_np, p_np, _ = background(jnp.asarray(0.0))
    arrays = tuple(jnp.asarray(x) for x in (t, a, rho_g, rho_nu, rho_np, p_np))
    t, a, rho_g, rho_nu, rho_np, p_np = arrays

    with output_path.open(mode, newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        completed = 0
        for row in rows:
            if row["point_id"] in done:
                continue
            started = time.perf_counter()
            record = dict(row)
            try:
                q = jnp.asarray([float(row[f"q_{name}"]) for name in REACTIONS])
                y = abundance(
                    rho_g, rho_nu, rho_np, p_np,
                    t_vec=t, a_vec=a,
                    eta_fac=jnp.asarray(float(row["eta_fac"])),
                    tau_n_fac=jnp.asarray(float(row["tau_seconds"]) / const.tau_n),
                    nuclear_rates_q=q,
                    sampling_nTOp=150,
                    rtol=1e-5,
                    atol=1e-9,
                )
                proton = float(y[1])
                record.update({
                    "status": "ok",
                    "message": "",
                    "D_H_1e5": float(y[2]) / proton * 1e5,
                    "Yp_mass": 4.0 * float(y[5]),
                    "Li7_H_1e10": (float(y[6]) + float(y[7])) / proton * 1e10,
                })
            except Exception as exc:
                record.update({"status": "failed", "message": f"{type(exc).__name__}: {exc}"})
            record["elapsed_s"] = time.perf_counter() - started
            record["run_time"] = now()
            writer.writerow(record)
            handle.flush()
            completed += 1
            if completed == 1 or completed % args.summary_every == 0:
                print(f"[{now()}] {output_path.name}: completed={completed}/{len(rows)} last={row['point_id']} status={record['status']} elapsed={record['elapsed_s']:.3f}s", flush=True)
    print(f"Saved exact LINX results: {output_path}", flush=True)


def main() -> int:
    args = parse_args()
    args.outdir = args.outdir.resolve()
    args.archive_dir = args.archive_dir.resolve()
    if args.mode == "register":
        write_registration(args)
        return 0
    if args.mode == "design":
        path = design_path(args)
        if not path.exists():
            raise FileNotFoundError(f"Run --mode register first: {path}")
        output = args.output_csv or (args.outdir / "exact_linx_design_results.csv")
        evaluate(args, path, output.resolve())
        return 0
    if args.input_points is None:
        raise ValueError("--input-points is required for validate mode")
    output = args.output_csv or (args.outdir / "exact_linx_posterior_validation_results.csv")
    evaluate(args, args.input_points.resolve(), output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
