#!/usr/bin/env python3
"""Preregistered hierarchical LINX/NumPyro BBN lithium posterior.

The baseline is deliberately conditional: homogeneous standard BBN with a
Gaussian CMB baryon-density summary, a laboratory neutron-lifetime prior, all
12 key-network LINX rate pulls, and current D, He and Li abundance summaries.
It does not model matter-antimatter domains, annihilation or entropy injection.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
DEFAULT_LINX_PATH = REPO_ROOT / "plamb_runs" / "deps" / "LINX"
DEFAULT_OUTDIR = (
    REPO_ROOT
    / "plamb_runs"
    / "diagnostics"
    / "bbn_lithium_hierarchical_numpyro_20260718"
)
DEFAULT_ARCHIVE = REPO_ROOT / "github_export" / "results" / "2026-07-18" / "bbn_lithium"

OBS = {
    "D_H_1e5": {"mean": 2.533, "sigma": 0.024},
    "Yp_mass": {"mean": 0.2458, "sigma": 0.0013},
    "Li7_H_1e10": {"mean": 1.45, "sigma": 0.25},
}
PRIORS = {
    "eta_fac": {"distribution": "TruncatedNormal", "mean": 1.0, "sigma": 0.04 / 6.12, "bounds_sigma": 5.0},
    "tau_seconds": {"distribution": "TruncatedNormal", "mean": 878.4, "sigma": 0.5, "bounds_sigma": 5.0},
    "rate_pull": {"distribution": "TruncatedNormal", "mean": 0.0, "sigma": 1.0, "low": -4.0, "high": 4.0},
    "delta_li_dex": {"distribution": "HalfNormal", "sigma": 0.30},
}
GATES = {"max_rhat": 1.01, "min_bulk_ess": 400.0, "max_divergences": 0}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["register", "smoke", "overnight", "analyse"], required=True)
    parser.add_argument("--arm", choices=["both", "predictive", "depletion"], default="both")
    parser.add_argument("--linx-path", type=Path, default=DEFAULT_LINX_PATH)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--archive-dir", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--chains", type=int)
    parser.add_argument("--warmup", type=int)
    parser.add_argument("--samples", type=int)
    parser.add_argument("--seed", type=int, default=20260718)
    parser.add_argument("--target-accept", type=float, default=0.90)
    parser.add_argument("--step-size", type=float, default=0.01)
    parser.add_argument("--max-tree-depth", type=int, default=4)
    parser.add_argument("--rtol", type=float)
    parser.add_argument("--atol", type=float, default=1e-9)
    parser.add_argument("--sampling-ntop", type=int)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def now() -> str:
    return datetime.now(UTC).astimezone().isoformat(timespec="seconds")


def script_sha256() -> str:
    return hashlib.sha256(SCRIPT_PATH.read_bytes()).hexdigest()


def settings(args: argparse.Namespace) -> dict[str, object]:
    smoke = args.mode == "smoke"
    return {
        "chains": args.chains if args.chains is not None else (1 if smoke else 4),
        "warmup": args.warmup if args.warmup is not None else (8 if smoke else 1000),
        "samples": args.samples if args.samples is not None else (8 if smoke else 2000),
        "rtol": args.rtol if args.rtol is not None else 1e-4,
        "atol": args.atol,
        "sampling_ntop": args.sampling_ntop if args.sampling_ntop is not None else (60 if smoke else 150),
        "target_accept": args.target_accept,
        "step_size": args.step_size,
        "max_tree_depth": args.max_tree_depth,
        "chain_method": "sequential",
        "gradient_mode": "forward",
        "max_steps": 8092 * 12,
        "network": "key_PRIMAT_2023",
    }


def registration_payload(args: argparse.Namespace) -> dict[str, object]:
    frozen = settings(argparse.Namespace(**{**vars(args), "mode": "overnight"}))
    return {
        "registered_at": now(),
        "script": str(SCRIPT_PATH),
        "script_sha256": script_sha256(),
        "linx_path": str(args.linx_path.resolve()),
        "outdir": str(args.outdir.resolve()),
        "archive_dir": str(args.archive_dir.resolve()),
        "jax_compilation_cache": str((args.outdir / "jax_compilation_cache").resolve()),
        "frozen_overnight_settings": frozen,
        "observations": OBS,
        "priors": PRIORS,
        "convergence_gates": GATES,
        "preposterior_numerical_amendment": (
            "The initial depth-8 smoke attempt was stopped during XLA compilation at "
            "0/16 transitions and produced no posterior draws. Maximum tree depth was "
            "set to the upstream LINX example value of 4, and chains were made sequential "
            "to control peak memory. A second smoke attempt failed at warm-up transition 2 "
            "after automatic initial step-size behaviour generated a non-finite implicit "
            "linear solve. Before any posterior was retained, the upstream LINX initial "
            "step size 0.01 was adopted; eta and tau were truncated at five prior standard "
            "deviations and rate pulls at +/-4 to keep trajectories inside the previously "
            "validated LINX matrix. That reverse-mode smoke still failed on its first "
            "trajectory, before retaining a draw. The final registered gradient trial uses "
            "NumPyro forward-mode differentiation and the upstream LINX Kvaerno3 maximum "
            "step budget at rtol 1e-4. A persistent run-local JAX compilation cache was "
            "enabled before the final smoke so compatible compiled modules can be reused "
            "by the overnight process. Likelihood anchors were unchanged."
        ),
        "arms": {
            "predictive": "D + He + Gaussian CMB eta summary; Li is posterior predictive only",
            "depletion": "D + He + Gaussian CMB eta summary + Li with positive dex depletion",
        },
        "claim_boundary": (
            "Homogeneous standard-BBN inference. It is not a matter-antimatter domain, "
            "annihilation, antimatter-galaxy, or self-consistent FR early-universe likelihood."
        ),
    }


def write_registration(args: argparse.Namespace) -> Path:
    args.archive_dir.mkdir(parents=True, exist_ok=True)
    payload = registration_payload(args)
    json_path = args.archive_dir / "bbn_lithium_hierarchical_numpyro_registration_2026-07-18.json"
    md_path = args.archive_dir / "bbn_lithium_hierarchical_numpyro_preregistration_2026-07-18.md"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md = """# Hierarchical LINX BBN Lithium Posterior: Preregistration

Registered: {payload['registered_at']}

## Registered Question

Can homogeneous standard BBN, after marginalising all 12 `key_PRIMAT_2023`
nuclear-rate pulls and the laboratory neutron-lifetime uncertainty, predict the
observed lithium abundance? If an explicit positive stellar-depletion nuisance
is admitted, what depletion is required?

## Likelihood And Priors

\[
\begin{aligned}
y_{\rm D}        &= 10^5({\rm D/H}), \\
Y_{\rm p}        &= 4Y_{{}^4{\rm He}}, \\
y_{{\rm Li},0}   &= 10^{10}[({}^7{\rm Li}+{}^7{\rm Be})/{\rm H}], \\
y_{{\rm Li},\star} &= y_{{\rm Li},0}\,10^{-\delta_{\rm Li}}.
\end{aligned}
\]

The abundance likelihoods are Normal summaries with D/H = 2.533 +/- 0.024,
Yp = 0.2458 +/- 0.0013 and Li/H = 1.45 +/- 0.25 in the scaled units above.
The CMB information is a five-sigma truncated Gaussian
`eta_fac = 1.0 +/- __ETA_SIGMA__`
summary, not a matched map-level CMB likelihood. Neutron lifetime is
`878.4 +/- 0.5 s`, truncated at five sigma; each LINX nuclear rate pull is a
standard Normal truncated at +/-4. The depletion
arm uses `delta_Li ~ HalfNormal(0.30 dex)`.

## Arms

1. `predictive`: condition on D, He and the CMB eta summary; do not condition on Li.
2. `depletion`: add the Li likelihood and positive stellar depletion.
3. Reweight the depletion posterior to 0.15 and 0.60 dex HalfNormal controls,
   reporting importance effective sample size before interpreting either control.

## Numerical Gates

The headline posterior requires maximum split Rhat <= 1.01, minimum bulk ESS >=
400 and zero divergent transitions. Failed gates make the output exploratory.

Maximum tree depth is 4, matching LINX's supplied NumPyro example. Four chains
run sequentially to control peak memory. An initial depth-8 smoke attempt was
stopped during XLA compilation at 0/16 transitions and supplied no posterior
draws. A second smoke attempt failed at warm-up transition 2 when automatic
initial step-size behaviour drove the implicit solver to a non-finite result.
Before retaining any posterior draw, the upstream LINX initial step size 0.01
and negligible-tail truncations described above were adopted. No likelihood
anchor was changed after observing results. The bounded reverse-mode smoke then
failed on its first trajectory, also before retaining a draw. The final gradient
trial uses forward-mode differentiation, Kvaerno3, `max_steps = 8092 x 12` and
`rtol = 1e-4`, following the upstream numerical envelope. A persistent run-local
JAX cache is enabled before the final smoke and reused by the overnight process.

## Matter-Antimatter Boundary

This run assumes homogeneous standard BBN with one net baryon-density field.
It contains no matter-antimatter domains, boundary annihilation, entropy
injection, spectral-distortion or antimatter-galaxy population parameter.
Consequently it can test standard-BBN and nuisance explanations of lithium, but
cannot validate or exclude an FR matter-antimatter formation mechanism.
""".replace("__ETA_SIGMA__", f"{PRIORS['eta_fac']['sigma']:.8f}")
    md_path.write_text(md, encoding="utf-8")
    print(f"Saved preregistration: {md_path}", flush=True)
    print(f"Saved frozen JSON: {json_path}", flush=True)
    return json_path


def ensure_registration(args: argparse.Namespace) -> None:
    path = args.archive_dir / "bbn_lithium_hierarchical_numpyro_registration_2026-07-18.json"
    if not path.exists():
        raise FileNotFoundError(f"Preregistration is required before overnight sampling: {path}")


def serialise_summary(value: object) -> object:
    array = np.asarray(value)
    if array.ndim == 0:
        return float(array)
    return array.tolist()


def quantiles(values: np.ndarray) -> dict[str, float]:
    finite = np.asarray(values, dtype=float).reshape(-1)
    finite = finite[np.isfinite(finite)]
    q = np.quantile(finite, [0.05, 0.50, 0.95])
    return {"q05": float(q[0]), "median": float(q[1]), "q95": float(q[2]), "mean": float(np.mean(finite)), "sd": float(np.std(finite, ddof=1))}


def selected_arms(name: str) -> list[str]:
    return ["predictive", "depletion"] if name == "both" else [name]


def run_arm(args: argparse.Namespace, arm: str, cfg: dict[str, object]) -> None:
    arm_dir = args.outdir / ("smoke" if args.mode == "smoke" else "overnight") / arm
    samples_path = arm_dir / "posterior_samples.npz"
    if samples_path.exists() and not args.force:
        print(f"Skipping completed arm: {arm} ({samples_path})", flush=True)
        return
    arm_dir.mkdir(parents=True, exist_ok=True)

    chains = int(cfg["chains"])
    cache_dir = args.outdir / "jax_compilation_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("XLA_FLAGS", "--xla_force_host_platform_device_count=1")
    sys.path.insert(0, str(args.linx_path.resolve()))

    import jax
    import jax.numpy as jnp
    import diffrax as dfx
    import numpyro
    import numpyro.distributions as dist
    from numpyro.diagnostics import summary
    from numpyro.infer import MCMC, NUTS, init_to_value
    from linx.abundances import AbundanceModel
    from linx.background import BackgroundModel
    from linx.nuclear import NuclearRates
    import linx.const as const

    jax.config.update("jax_enable_x64", True)
    jax.config.update("jax_compilation_cache_dir", str(cache_dir))
    jax.config.update("jax_persistent_cache_min_compile_time_secs", 1.0)
    numpyro.set_host_device_count(1)

    rates = NuclearRates(nuclear_net=str(cfg["network"]))
    reactions = list(rates.reactions_names)
    abundance = AbundanceModel(rates, throw=False)
    background = BackgroundModel()
    t_vec, a_vec, rho_g, rho_nu, rho_np, p_np, _ = background(jnp.asarray(0.0))
    t_vec, a_vec = jnp.asarray(t_vec), jnp.asarray(a_vec)
    rho_g, rho_nu = jnp.asarray(rho_g), jnp.asarray(rho_nu)
    rho_np, p_np = jnp.asarray(rho_np), jnp.asarray(p_np)

    def model() -> None:
        eta_mean = float(PRIORS["eta_fac"]["mean"])
        eta_sigma = float(PRIORS["eta_fac"]["sigma"])
        eta_bound = float(PRIORS["eta_fac"]["bounds_sigma"]) * eta_sigma
        tau_mean = float(PRIORS["tau_seconds"]["mean"])
        tau_sigma = float(PRIORS["tau_seconds"]["sigma"])
        tau_bound = float(PRIORS["tau_seconds"]["bounds_sigma"]) * tau_sigma
        eta_fac = numpyro.sample(
            "eta_fac",
            dist.TruncatedNormal(eta_mean, eta_sigma, low=eta_mean - eta_bound, high=eta_mean + eta_bound),
        )
        tau_seconds = numpyro.sample(
            "tau_seconds",
            dist.TruncatedNormal(tau_mean, tau_sigma, low=tau_mean - tau_bound, high=tau_mean + tau_bound),
        )
        rate_pull = numpyro.sample(
            "rate_pull",
            dist.TruncatedNormal(
                0.0,
                1.0,
                low=float(PRIORS["rate_pull"]["low"]),
                high=float(PRIORS["rate_pull"]["high"]),
            ).expand([len(reactions)]).to_event(1),
        )
        yields = abundance(
            rho_g,
            rho_nu,
            rho_np,
            p_np,
            t_vec=t_vec,
            a_vec=a_vec,
            eta_fac=eta_fac,
            tau_n_fac=tau_seconds / const.tau_n,
            nuclear_rates_q=rate_pull,
            sampling_nTOp=int(cfg["sampling_ntop"]),
            rtol=float(cfg["rtol"]),
            atol=float(cfg["atol"]),
            solver=dfx.Kvaerno3(),
            max_steps=int(cfg["max_steps"]),
        )
        proton = yields[1]
        d_h = yields[2] / proton * 1e5
        yp = 4.0 * yields[5]
        li_primordial = (yields[6] + yields[7]) / proton * 1e10
        numpyro.deterministic("omega_b_h2", eta_fac * const.Omegabh2)
        numpyro.deterministic("D_H_1e5", d_h)
        numpyro.deterministic("Yp_mass", yp)
        numpyro.deterministic("Li7_H_1e10", li_primordial)
        numpyro.sample("D_observed", dist.Normal(d_h, OBS["D_H_1e5"]["sigma"]), obs=OBS["D_H_1e5"]["mean"])
        numpyro.sample("He_observed", dist.Normal(yp, OBS["Yp_mass"]["sigma"]), obs=OBS["Yp_mass"]["mean"])
        if arm == "depletion":
            delta = numpyro.sample("delta_li_dex", dist.HalfNormal(PRIORS["delta_li_dex"]["sigma"]))
            li_surface = li_primordial * jnp.power(10.0, -delta)
            numpyro.deterministic("Li7_surface_H_1e10", li_surface)
            numpyro.sample("Li_observed", dist.Normal(li_surface, OBS["Li7_H_1e10"]["sigma"]), obs=OBS["Li7_H_1e10"]["mean"])

    initial = {
        "eta_fac": 1.0,
        "tau_seconds": PRIORS["tau_seconds"]["mean"],
        "rate_pull": jnp.zeros(len(reactions)),
    }
    if arm == "depletion":
        initial["delta_li_dex"] = 0.55
    kernel = NUTS(
        model,
        step_size=float(cfg["step_size"]),
        target_accept_prob=float(cfg["target_accept"]),
        max_tree_depth=int(cfg["max_tree_depth"]),
        dense_mass=True,
        forward_mode_differentiation=str(cfg["gradient_mode"]) == "forward",
        init_strategy=init_to_value(values=initial),
    )
    method = str(cfg["chain_method"])
    mcmc = MCMC(
        kernel,
        num_warmup=int(cfg["warmup"]),
        num_samples=int(cfg["samples"]),
        num_chains=chains,
        chain_method=method,
        progress_bar=True,
    )
    seed = args.seed + (0 if arm == "predictive" else 1000)
    started = time.time()
    print(f"[{now()}] Starting {arm}: {cfg}", flush=True)
    mcmc.run(jax.random.PRNGKey(seed))
    elapsed = time.time() - started
    samples = mcmc.get_samples(group_by_chain=True)
    extras = mcmc.get_extra_fields(group_by_chain=True)
    np.savez_compressed(samples_path, **{name: np.asarray(value) for name, value in samples.items()})

    latent_names = ["eta_fac", "tau_seconds", "rate_pull"]
    if arm == "depletion":
        latent_names.append("delta_li_dex")
    diagnostics_raw = summary({name: samples[name] for name in latent_names}, prob=0.90, group_by_chain=True)
    diagnostics = {
        name: {key: serialise_summary(value) for key, value in values.items()}
        for name, values in diagnostics_raw.items()
    }
    rhats = np.concatenate([np.asarray(values["r_hat"]).reshape(-1) for values in diagnostics_raw.values()])
    esses = np.concatenate([np.asarray(values["n_eff"]).reshape(-1) for values in diagnostics_raw.values()])
    divergences = int(np.asarray(extras.get("diverging", [])).sum())
    gate = {
        "max_rhat": float(np.nanmax(rhats)),
        "min_bulk_ess": float(np.nanmin(esses)),
        "divergences": divergences,
    }
    gate["passed"] = bool(
        gate["max_rhat"] <= GATES["max_rhat"]
        and gate["min_bulk_ess"] >= GATES["min_bulk_ess"]
        and divergences <= GATES["max_divergences"]
    )
    result = {
        "arm": arm,
        "completed_at": now(),
        "elapsed_seconds": elapsed,
        "settings": cfg,
        "reaction_names": reactions,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "jax": jax.__version__,
            "numpyro": numpyro.__version__,
        },
        "gate": gate,
        "posterior": {
            name: quantiles(np.asarray(samples[name]))
            for name in ["D_H_1e5", "Yp_mass", "Li7_H_1e10", "omega_b_h2"]
            + (["delta_li_dex", "Li7_surface_H_1e10"] if arm == "depletion" else [])
        },
        "diagnostics": diagnostics,
    }
    (arm_dir / "run_result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report = [
        f"# Hierarchical LINX Posterior: {arm}",
        "",
        f"Completed: {result['completed_at']}",
        f"Elapsed: {elapsed / 3600.0:.3f} h",
        "",
        "## Convergence Gate",
        "",
        f"- Passed: `{gate['passed']}`",
        f"- Maximum Rhat: `{gate['max_rhat']:.6g}` (required <= {GATES['max_rhat']})",
        f"- Minimum bulk ESS: `{gate['min_bulk_ess']:.6g}` (required >= {GATES['min_bulk_ess']})",
        f"- Divergences: `{divergences}` (required 0)",
        "",
        "## Posterior Summaries",
        "",
        "| quantity | 5% | median | 95% |",
        "|:--|--:|--:|--:|",
    ]
    for name, values in result["posterior"].items():
        report.append(f"| {name} | {values['q05']:.7g} | {values['median']:.7g} | {values['q95']:.7g} |")
    report.extend([
        "",
        "## Claim Boundary",
        "",
        "This is a summary-likelihood homogeneous standard-BBN posterior. It is not a matched map-level CMB likelihood and contains no matter-antimatter domain or annihilation physics.",
    ])
    (arm_dir / "posterior_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"[{now()}] Completed {arm}; gate={gate}; samples={samples_path}", flush=True)


def weighted_quantile(values: np.ndarray, weights: np.ndarray, probability: float) -> float:
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cdf = np.cumsum(weights) / np.sum(weights)
    return float(np.interp(probability, cdf, values))


def analyse(args: argparse.Namespace) -> None:
    root = args.outdir / "overnight"
    output_rows: list[dict[str, object]] = []
    reports: list[dict[str, object]] = []
    for arm in ["predictive", "depletion"]:
        result_path = root / arm / "run_result.json"
        if result_path.exists():
            reports.append(json.loads(result_path.read_text(encoding="utf-8")))
    depletion_path = root / "depletion" / "posterior_samples.npz"
    if depletion_path.exists():
        with np.load(depletion_path) as data:
            delta = np.asarray(data["delta_li_dex"]).reshape(-1)
            li0 = np.asarray(data["Li7_H_1e10"]).reshape(-1)
            lis = np.asarray(data["Li7_surface_H_1e10"]).reshape(-1)
        old_sigma = float(PRIORS["delta_li_dex"]["sigma"])
        for sigma in [0.15, 0.30, 0.60]:
            logw = -0.5 * delta**2 / sigma**2 - math.log(sigma)
            logw -= -0.5 * delta**2 / old_sigma**2 - math.log(old_sigma)
            logw -= np.max(logw)
            weights = np.exp(logw)
            weights /= np.sum(weights)
            output_rows.append({
                "depletion_prior_sigma_dex": sigma,
                "importance_ess": float(1.0 / np.sum(weights**2)),
                "delta_li_dex_median": weighted_quantile(delta, weights, 0.5),
                "depletion_factor_median": float(10.0 ** weighted_quantile(delta, weights, 0.5)),
                "primordial_li_median": weighted_quantile(li0, weights, 0.5),
                "surface_li_median": weighted_quantile(lis, weights, 0.5),
            })
    args.archive_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.archive_dir / "bbn_lithium_depletion_prior_sensitivity_2026-07-18.csv"
    if output_rows:
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
            writer.writeheader()
            writer.writerows(output_rows)
    report_path = args.archive_dir / "bbn_lithium_hierarchical_numpyro_readout_2026-07-18.md"
    lines = ["# Hierarchical LINX BBN Lithium Posterior Readout", "", f"Generated: {now()}", ""]
    for report in reports:
        gate = report["gate"]
        li = report["posterior"]["Li7_H_1e10"]
        lines.extend([
            f"## {report['arm'].title()} Arm",
            "",
            f"Gate passed: `{gate['passed']}`; max Rhat `{gate['max_rhat']:.5g}`; minimum ESS `{gate['min_bulk_ess']:.5g}`; divergences `{gate['divergences']}`.",
            f"Primordial Li/H x 1e10: median `{li['median']:.5g}`, 90% interval `[{li['q05']:.5g}, {li['q95']:.5g}]`.",
            "",
        ])
    if output_rows:
        lines.extend(["## Depletion-Prior Sensitivity", "", "| prior sigma (dex) | importance ESS | median depletion (dex) | median factor |", "|--:|--:|--:|--:|"])
        for row in output_rows:
            lines.append(f"| {row['depletion_prior_sigma_dex']:.2f} | {row['importance_ess']:.1f} | {row['delta_li_dex_median']:.4g} | {row['depletion_factor_median']:.4g} |")
        lines.append("")
    lines.extend([
        "## Interpretation Boundary",
        "",
        "These results remain conditional on homogeneous standard BBN. Matter-antimatter domains would require a separate transport and annihilation likelihood before any inference could be transferred to that sector of FR.",
    ])
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Saved analysis report: {report_path}", flush=True)
    if output_rows:
        print(f"Saved sensitivity table: {csv_path}", flush=True)


def main() -> int:
    args = parse_args()
    args.outdir = args.outdir.resolve()
    args.archive_dir = args.archive_dir.resolve()
    if args.mode == "register":
        write_registration(args)
        return 0
    if args.mode == "analyse":
        analyse(args)
        return 0
    if args.mode == "overnight":
        ensure_registration(args)
    cfg = settings(args)
    for arm in selected_arms(args.arm):
        run_arm(args, arm, cfg)
    if args.mode == "overnight" and args.arm == "both":
        analyse(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
