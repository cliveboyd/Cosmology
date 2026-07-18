# SPARC FR Environmental-Asymmetry Preregistration

Date: 2026-07-18
Written before outcome profiling: `2026-07-18T07:19:13.135626+00:00`

## Scientific Question

Does an independently measured, charge-blind galaxy environment predict where the fixed-global PLAMB bridge differs from RAR after conventional galaxy properties and nuisance structure are controlled?

This is the only rotation-curve consequence left open by the FR charge-conjugation derivation: a host-only change in a fixed asymmetric environment can alter the magnitude of the FR background asymmetry. The observable 2MRS environment is not a signed matter/antimatter field and cannot identify antimatter.

## Locked Samples

| sample | N_outcome_blind |
| --- | --- |
| primary_10_40_mpc_mk22 | 66 |
| near_5_40_mpc_mk22 | 97 |
| wide_10_80_mpc_mk23 | 80 |
| latitude_strict_10_40_mpc_mk22 | 61 |

The primary sample uses SPARC Q<=2 galaxies at 10-40 Mpc, |b|>=10 degrees and a 2MRS neighbour threshold M_K<=-22. The six previously outcome-selected development galaxies and all five reserved replication galaxies are excluded from primary inference.

## Environment Composite

$$
\begin{aligned}
E_1 &= \log(1+N_{2\,{\rm Mpc}}), \\
E_2 &= \log(1+N_{5\,{\rm Mpc}}), \\
E_3 &= -\log_{10} R_{\rm nearest}, \\
E_4 &= \log_{10}\!\left[10^{-6}+\sum_j\frac{L_{K,j}/L_{K,-23}}{(R_j^2+0.05^2)^{3/2}}\right], \\
E   &= \frac{1}{4}\sum_{k=1}^{4}\frac{E_k-\operatorname{median}(E_k)}{\operatorname{IQR}(E_k)}.
\end{aligned}
$$

Neighbours must lie within |Delta v|<=500 km/s. A 1000<=|Delta v|<2000 km/s sideband enters the host-control model as a projection/selection control. The negative-control environment repeats the calculation at RA+90 degrees.

## Locked Outcome

$$
\begin{aligned}
\Delta_g &= \frac{Q_{{\rm PLAMB},g}-Q_{{\rm RAR},g}}{N_g}, \\
Y_g      &= \operatorname{sign}(\Delta_g)\log(1+|\Delta_g|).
\end{aligned}
$$

The primary outcome uses the combined-conventional nuisance profile. The baseline profile is secondary. Positive values mean that PLAMB fits worse than RAR.

## Locked Test

Ridge models with fixed alpha=1 use ten deterministic folds. The primary association is a two-sided 20,000-permutation Spearman test of cross-fitted outcome and environment residuals after the locked host controls. Predictive value is the fractional reduction in cross-fitted RMSE when E is added.

## Development Gate

Every condition below is required before a separate replication preregistration may be written:

1. primary permutation p<=0.01;
2. cross-fitted RMSE improvement>=5 per cent;
3. environment coefficient has one sign in at least 8/10 folds;
4. the partial-correlation sign agrees in the wide-distance, strict-latitude and baseline-outcome controls; and
5. the shifted-position negative-control |rho| is smaller than the actual-environment |rho|.

Passing this gate would establish an environment residual worth replication. It would not establish a signed antimatter field.

## Reproducibility Lock

```json
{
  "date": "2026-07-18",
  "written_before_outcomes_utc": "2026-07-18T07:19:13.135626+00:00",
  "programme": "C:\\Users\\clive\\Documents\\Cosmology\\github_export\\code\\sparc\\run_sparc_fr_environment_asymmetry_2026-07-18.py",
  "programme_sha256": "e6264583ff4ad9e0eeb5fbb951aee71394abba8de7167dd2142a36344cf76f78",
  "outcome_opened_before_lock": false,
  "theory_boundary": "FR complete charge conjugation is charge-even. This study tests whether an independently observed charge-blind environment predicts model residuals; it does not classify antimatter.",
  "development_exclusions": {
    "previously_outcome_selected": [
      "UGC06787",
      "UGC03580",
      "UGC02487",
      "NGC5985",
      "UGC09133",
      "NGC2903"
    ],
    "reserved_replication": [
      "NGC2841",
      "UGC02953",
      "UGC07399",
      "UGC00128",
      "NGC2403"
    ]
  },
  "catalogues": {
    "central_positions_and_redshifts": {
      "source": "NASA/IPAC NED ObjectLookup",
      "url": "https://ned.ipac.caltech.edu/Documents/Guides/Interface/ObjectLookup",
      "cache": "C:\\Users\\clive\\Documents\\Cosmology\\external_datasets\\current_cosmology_sources\\SPARC\\sparc_ned_positions_environment_nonreserved_20260718.csv"
    },
    "environment_tracer": {
      "source": "2MASS Redshift Survey table 3",
      "url": "https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJS/199/26",
      "paper": "https://doi.org/10.1088/0067-0049/199/2/26",
      "catalogue": "C:\\Users\\clive\\Documents\\Cosmology\\external_datasets\\current_cosmology_sources\\2MRS\\2mrs_table3_environment_20260718.tsv",
      "catalogue_limit": "K_s <= 11.75 mag; 97.6 per cent redshift complete over 91 per cent of sky"
    }
  },
  "input_hashes": {
    "sample": "48da73dd7306536c1bdffb33d4b319e077dc53d5d11ba519f6f01aab07d20a69",
    "points": "7904b885e3244dec9500e40f343f7c30311d45b446137882580d8c2123f95dac",
    "map_summary": "b9d8d0b19252a433377ef19987ad50cb799a060c23be35d67e976f7c1842b732",
    "2mrs": "ff3d3f9ae22ac9e9ed2402843c762d76b25d6d6580495db55292ca6b4c18ee86",
    "positions": "c1e51431768fefa741bbc8cdecb57ae8e2c2b98b4b75c94e0056713c8263398f",
    "fr_rule": "baaca554cb19dac0f4c23a59b5a1054b0db07d9b8f531890683b1d02187504f0",
    "am_profile_programme": "c5dae3af6cc6b63ac5db29f2ea29db4f8e3ffb1bdfd5d8d0fbc43ccca08740a9",
    "environment_covariates": "72381709dc851aef18998d592dc4059be82a7c7c1323c6a6254cdfaec7c0b33e"
  },
  "samples": {
    "primary": {
      "name": "primary_10_40_mpc_mk22",
      "distance_min_mpc": 10.0,
      "distance_max_mpc": 40.0,
      "abs_galactic_latitude_min_deg": 10.0,
      "absolute_k_limit": -22.0
    },
    "sensitivities": [
      {
        "name": "near_5_40_mpc_mk22",
        "distance_min_mpc": 5.0,
        "distance_max_mpc": 40.0,
        "abs_galactic_latitude_min_deg": 10.0,
        "absolute_k_limit": -22.0
      },
      {
        "name": "wide_10_80_mpc_mk23",
        "distance_min_mpc": 10.0,
        "distance_max_mpc": 80.0,
        "abs_galactic_latitude_min_deg": 10.0,
        "absolute_k_limit": -23.0
      },
      {
        "name": "latitude_strict_10_40_mpc_mk22",
        "distance_min_mpc": 10.0,
        "distance_max_mpc": 40.0,
        "abs_galactic_latitude_min_deg": 20.0,
        "absolute_k_limit": -22.0
      }
    ],
    "outcome_blind_counts": {
      "primary_10_40_mpc_mk22": 66,
      "near_5_40_mpc_mk22": 97,
      "wide_10_80_mpc_mk23": 80,
      "latitude_strict_10_40_mpc_mk22": 61
    },
    "quality": "SPARC Q <= 2; finite NED position/redshift; stated distance and latitude cuts"
  },
  "environment": {
    "line_of_sight_window_kms": 500.0,
    "projected_radii_mpc": [
      2.0,
      5.0
    ],
    "sideband_kms": [
      1000.0,
      2000.0
    ],
    "self_exclusion": {
      "angle_deg": 0.008333333333333333,
      "delta_v_kms": 300.0
    },
    "tidal_softening_mpc": 0.05,
    "composite_components": [
      "log1p_count_2mpc",
      "log1p_count_5mpc",
      "minus_log10_nearest_mpc",
      "log10_tidal_strength"
    ],
    "composite_rule": "mean of robust-standardised log1p N(<2 Mpc), log1p N(<5 Mpc), -log10 nearest-neighbour distance and log10 softened K-luminosity tidal sum",
    "negative_control": "same proxy at RA shifted by +90 degrees with Dec and central redshift fixed"
  },
  "outcomes": {
    "primary": "signed_log1p((profile_objective_PLAMB - profile_objective_RAR) / N_points) under the locked combined_conventional nuisance scenario",
    "secondary": "same transformed contrast under the baseline nuisance scenario",
    "global_parameters": "fixed July 14 all-Q2 MAP values",
    "optimisation": "deterministic multi-start galaxy profile inherited unchanged from the July 18 systematics audit"
  },
  "host_controls": [
    "T",
    "log_L36",
    "log_SBeff",
    "log_MHI",
    "log_Vflat",
    "Inc_deg",
    "frac_distance_error",
    "bulge_proxy",
    "log_n_points",
    "D_Mpc",
    "abs_gal_b_deg",
    "log1p_sideband_count_5mpc",
    "f_D one-hot indicators"
  ],
  "statistics": {
    "ridge_alpha": 1.0,
    "outer_folds": 10,
    "fold_seed": 20260718,
    "partial_test": "Spearman correlation of cross-fitted host-control residuals",
    "permutations": 20000,
    "p_value": "two-sided deterministic permutation p=(1+exceedances)/(1+permutations)",
    "predictive_metric": "cross-fitted RMSE baseline versus baseline plus environment composite"
  },
  "development_gate_for_future_replication": {
    "all_required": true,
    "primary_partial_permutation_p_max": 0.01,
    "primary_cv_rmse_fractional_improvement_min": 0.05,
    "environment_coefficient_sign_fraction_min": 0.8,
    "same_partial_correlation_sign_in": [
      "wide_10_80_mpc_mk23",
      "latitude_strict_10_40_mpc_mk22",
      "baseline_outcome_primary_sample"
    ],
    "negative_control": "absolute shifted-position partial rho must be smaller than actual-environment rho",
    "replication_policy": "never unseal automatically; write a separate replication preregistration first"
  },
  "claim_boundary": [
    "A positive result would show ordinary environment dependence of model residuals, not antimatter.",
    "No fitted matter/antimatter galaxy sign is permitted.",
    "Candidate-galaxy environment comparisons are descriptive because those galaxies were outcome-selected.",
    "Reserved replication galaxies remain unopened in this programme."
  ]
}
```
