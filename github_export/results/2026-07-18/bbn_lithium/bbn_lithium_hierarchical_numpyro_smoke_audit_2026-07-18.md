# Hierarchical LINX NumPyro Smoke Audit

Date: 18 July 2026

## Outcome

The direct-gradient NumPyro route is **not validated** on the frozen Windows
LINX/JAX stack. No posterior draw was retained and no scientific result was
inferred from these attempts.

## Attempts

1. A depth-8 reverse-mode smoke was stopped during XLA compilation at 0/16
   transitions because its recursive graph and memory allocation were excessive.
2. The upstream depth-4 setting reached warm-up transition 2, then a large
   automatically initialised NUTS step caused the Kvaerno implicit linear solve
   to return NaN/Inf.
3. The upstream initial step size `0.01`, five-sigma eta/tau truncation and
   rate-pull bounds `[-4,4]` were applied before retaining any draw. Reverse mode
   still failed on its first trajectory with a non-finite implicit solve.
4. Forward-mode differentiation was tested with the upstream Kvaerno3 maximum
   step budget. JAX rejected it because LINX exposes a `custom_vjp` for which a
   forward-mode JVP is not defined.

## Decision

Longer wall time cannot repair these gradient failures. The registered fallback
uses a 4,800-point exact-LINX prior-space design over eta, neutron lifetime and
all 12 key-network rate pulls; a fixed 4,000/800 train/holdout split; a
gradient-free ensemble posterior; and selected posterior draws returned to exact
LINX for error and importance-weight gates.

The disabled NumPyro programme is retained as an audit artefact. The active
pipeline is `run_bbn_lithium_hierarchical_surrogate_pipeline_2026-07-18.ps1`.
