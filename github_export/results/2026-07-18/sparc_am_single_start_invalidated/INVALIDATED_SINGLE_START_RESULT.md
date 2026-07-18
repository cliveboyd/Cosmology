# Invalidated Single-Start Result

This result bundle is retained for audit history but must not be used for
scientific interpretation.

The post-run numerical audit found that the `stellar_ml_tight` profile for
`UGC06787` converged to different local minima for RAR and PLAMB. The resulting
single-start delta objective was `-430.4497`. Deterministic multi-start
profiling reduced it to approximately `-10.44`, while reproducing the baseline
and combined-conventional profiles to numerical precision.

The replacement analysis is preregistered and exported in the sibling
`sparc_am` directory after this bundle is archived as
`sparc_am_single_start_invalidated`.
