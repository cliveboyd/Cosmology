# Hierarchical LINX BBN Surrogate Posterior: Preregistration

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
