# SPARC Antimatter Identifiability and Conventional-Systematics Report

Date: 2026-07-18
Completed: `2026-07-18T06:00:25.460364+00:00`

## Bottom Line

The current RAR/PLAMB SPARC likelihood is matter/antimatter non-identifiable: it has no particle-sign or charge-conjugation branch. Rotation-curve persistence can therefore nominate an anomaly for further study, but cannot identify antimatter.

After the locked conventional-systematics profile, `3` of six development galaxies meet the preregistered persistent-tension rule: UGC02487, UGC09133, NGC2903. `3` are erased or reversed: UGC06787, UGC03580, NGC5985.

In the combined-conventional scenario, the mean candidate-minus-matched-control Delta objective is `0.374294`. Positive values mean that the candidate set remains relatively more difficult for PLAMB than its matched controls.

## Gate 1: FR Antimatter Identifiability

**NOT IDENTIFIABLE.** The supplied FR/SPARC implementation does not specify a matter-to-antimatter transformation for gbar, g0, kappa, p or an additional field. An arbitrary fitted sign is prohibited because it would be an empirical residual absorber rather than an FR prediction.

## Outcome-Blind Matched Controls

The maximum absolute robust-standardised mean difference across the ten locked matching variables is `0.225`. Matching used no RAR or PLAMB residual.

| feature | candidate_median_original | control_median_original | standardised_mean_difference |
| --- | --- | --- | --- |
| T | 2 | 2.5 | -0.183333 |
| log_L36 | 2.15597 | 2.04947 | 0.0107022 |
| log_SBeff | 2.85683 | 3.08483 | -0.0949277 |
| log_MHI | 0.882751 | 1.01968 | -0.102665 |
| log_Vflat | 2.37513 | 2.3298 | 0.0492577 |
| Inc_deg | 61.5 | 64 | -0.147436 |
| frac_distance_error | 0.249568 | 0.249806 | 0.126853 |
| bulge_proxy | 0.323914 | 0.31097 | -0.149505 |
| log_n_points | 1.60179 | 1.55897 | 0.22502 |
| Q | 1 | 1 | 0.0833333 |

## Conventional-Systematics Matrix

Positive Delta objective means the fixed-global PLAMB branch profiles worse than fixed-global RAR. The sign-test column is descriptive only because the six candidates were selected previously from their residuals.

| scenario | candidate_mean_delta | control_mean_delta | mean_matched_block_contrast | positive_candidate_blocks | descriptive_one_sided_sign_p | optimiser_failures |
| --- | --- | --- | --- | --- | --- | --- |
| baseline | 1.58906 | 0.369722 | 1.21933 | 3 | 0.65625 | 0 |
| combined_conventional | 0.388161 | 0.0138672 | 0.374294 | 3 | 0.65625 | 0 |
| disk_bulge_decoupled | 4.75957 | 0.129007 | 4.63057 | 4 | 0.34375 | 0 |
| distance_floor_10pct | 1.58906 | 0.369722 | 1.21933 | 3 | 0.65625 | 0 |
| distance_published_frame | 1.14471 | 0.387174 | 0.757541 | 3 | 0.65625 | 0 |
| distance_student_t3 | 1.59065 | 0.404854 | 1.1858 | 3 | 0.65625 | 0 |
| gas_scale_20pct | 2.52715 | 0.339051 | 2.1881 | 4 | 0.34375 | 0 |
| inclination_nuisance | 2.55995 | 0.386641 | 2.17331 | 4 | 0.34375 | 0 |
| radial_correlation_rho035 | 0.949249 | 0.146178 | 0.803071 | 3 | 0.65625 | 0 |
| stellar_ml_tight | -69.1816 | 0.391055 | -69.5726 | 3 | 0.65625 | 0 |
| stellar_ml_wide | 4.43005 | 0.365565 | 4.06449 | 4 | 0.34375 | 0 |
| velocity_jitter_5kms | 0.473844 | 0.157762 | 0.316082 | 3 | 0.65625 | 0 |

## Candidate Stability

| galaxy | fraction_scenarios_plamb_worse | baseline_delta_profile_objective | combined_delta_profile_objective | robustness_status |
| --- | --- | --- | --- | --- |
| UGC06787 | 0.25 | -3.0228 | -1.14102 | erased_or_reversed |
| UGC03580 | 0 | -2.69788 | -0.131142 | erased_or_reversed |
| UGC02487 | 0.916667 | 0.547865 | 0.0546327 | persistent_conventional_tension |
| NGC5985 | 0.25 | -3.66068 | 0.0702048 | erased_or_reversed |
| UGC09133 | 1 | 12.151 | 1.26055 | persistent_conventional_tension |
| NGC2903 | 1 | 6.2168 | 2.21574 | persistent_conventional_tension |

## Interpretation

A persistent flag means only that this particular nuisance matrix did not remove the candidate's relative PLAMB tension. It is not evidence for antimatter. A reversed or unstable flag supports an ordinary calibration, geometry or baryonic-decomposition explanation.

The candidate profiles hold the July 14 all-Q2 global parameters fixed. This provides a transparent and symmetric local challenge but is not a refitted hierarchical posterior. Publication-grade follow-up would require a newly derived FR antimatter transformation, a converged hierarchical nuisance model, and validation on the untouched replication galaxies.

## Locked Claim Boundary

- Do not label any SPARC galaxy as matter or antimatter from this analysis.
- Do not treat the descriptive sign p-value as a discovery probability.
- Do not inspect the reserved replication set until an FR transformation and decision rule are frozen.
- External gamma-ray analysis is an independent contact-annihilation test, not confirmation of a rotation-curve label by itself.

## Configuration

```json
{
  "date": "2026-07-18",
  "written_before_outcomes_utc": "2026-07-18T06:00:01.108268+00:00",
  "script": "C:\\Users\\clive\\Documents\\Cosmology\\github_export\\code\\sparc\\run_sparc_am_identifiability_systematics_2026-07-18.py",
  "script_sha256": "17c6279397f27ed39ecf5ac2bb95f344e4258749a8b6da74de7f8265ed6fbb50",
  "sample": "C:\\Users\\clive\\Documents\\Cosmology\\external_datasets\\current_cosmology_sources\\SPARC\\sparc_galaxy_sample.csv",
  "sample_sha256": "48da73dd7306536c1bdffb33d4b319e077dc53d5d11ba519f6f01aab07d20a69",
  "points": "C:\\Users\\clive\\Documents\\Cosmology\\external_datasets\\current_cosmology_sources\\SPARC\\sparc_rotation_points.csv",
  "points_sha256": "7904b885e3244dec9500e40f343f7c30311d45b446137882580d8c2123f95dac",
  "map_summary": "C:\\Users\\clive\\Documents\\Cosmology\\plamb_runs\\diagnostics\\sparc_hierarchical_map\\optical_depth_hierarchical_20260714\\sparc_hierarchical_map_summary.csv",
  "map_summary_sha256": "b9d8d0b19252a433377ef19987ad50cb799a060c23be35d67e976f7c1842b732",
  "candidates": [
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
  ],
  "models": [
    "RAR",
    "PLAMB_OPTICAL_DEPTH_KAPPA_P"
  ],
  "global_parameter_treatment": "fixed at July 14 all-Q2 hierarchical MAP values",
  "matching": {
    "eligible_control_population": 83,
    "matches_per_candidate": 2,
    "unique_controls": 12,
    "features": [
      "T",
      "log_L36",
      "log_SBeff",
      "log_MHI",
      "log_Vflat",
      "Inc_deg",
      "frac_distance_error",
      "bulge_proxy",
      "log_n_points",
      "Q"
    ],
    "weights": {
      "T": 0.75,
      "log_L36": 1.25,
      "log_SBeff": 1.0,
      "log_MHI": 1.0,
      "log_Vflat": 1.25,
      "Inc_deg": 0.75,
      "frac_distance_error": 0.75,
      "bulge_proxy": 1.25,
      "log_n_points": 0.75,
      "Q": 0.5
    },
    "robust_scaling": {
      "T": {
        "median": 7.0,
        "scale_iqr": 5.0
      },
      "log_L36": {
        "median": 0.6127838567197355,
        "scale_iqr": 1.8143028868699165
      },
      "log_SBeff": {
        "median": 1.7528164311882715,
        "scale_iqr": 1.3519347877810939
      },
      "log_MHI": {
        "median": 0.3963737275365065,
        "scale_iqr": 0.7735825162989017
      },
      "log_Vflat": {
        "median": 1.9585638832219674,
        "scale_iqr": 0.7386017958185545
      },
      "Inc_deg": {
        "median": 56.0,
        "scale_iqr": 26.0
      },
      "frac_distance_error": {
        "median": 0.2508710801393728,
        "scale_iqr": 0.10074487895716946
      },
      "bulge_proxy": {
        "median": 0.0,
        "scale_iqr": 0.19195353079916438
      },
      "log_n_points": {
        "median": 1.1760912590556813,
        "scale_iqr": 0.3979400086720377
      },
      "Q": {
        "median": 1.0,
        "scale_iqr": 1.0
      }
    },
    "outcome_blind": true,
    "exact_match": "presence or absence of a non-zero SPARC bulge component",
    "excluded_replication_galaxies": [
      "NGC2841",
      "UGC02953",
      "UGC07399",
      "UGC00128",
      "NGC2403"
    ]
  },
  "scenarios": [
    {
      "name": "baseline",
      "description": "Published coupled stellar scaling, Gaussian distance prior and 3 km/s error floor.",
      "sigma_ln_ml": 0.25,
      "ml_mode": "coupled",
      "distance_floor_frac": 0.03,
      "distance_prior": "gaussian",
      "hubble_prior_mode": "model_h0_rescale",
      "inclination_nuisance": false,
      "velocity_jitter_kms": 0.0,
      "radial_correlation": 0.0,
      "gas_scale_sigma": 0.0
    },
    {
      "name": "stellar_ml_tight",
      "description": "Tighter 15 per cent log stellar mass-to-light prior.",
      "sigma_ln_ml": 0.15,
      "ml_mode": "coupled",
      "distance_floor_frac": 0.03,
      "distance_prior": "gaussian",
      "hubble_prior_mode": "model_h0_rescale",
      "inclination_nuisance": false,
      "velocity_jitter_kms": 0.0,
      "radial_correlation": 0.0,
      "gas_scale_sigma": 0.0
    },
    {
      "name": "stellar_ml_wide",
      "description": "Wider 50 per cent log stellar mass-to-light prior.",
      "sigma_ln_ml": 0.5,
      "ml_mode": "coupled",
      "distance_floor_frac": 0.03,
      "distance_prior": "gaussian",
      "hubble_prior_mode": "model_h0_rescale",
      "inclination_nuisance": false,
      "velocity_jitter_kms": 0.0,
      "radial_correlation": 0.0,
      "gas_scale_sigma": 0.0
    },
    {
      "name": "disk_bulge_decoupled",
      "description": "Independent galaxy-level disc and bulge mass-to-light multipliers.",
      "sigma_ln_ml": 0.25,
      "ml_mode": "separate",
      "distance_floor_frac": 0.03,
      "distance_prior": "gaussian",
      "hubble_prior_mode": "model_h0_rescale",
      "inclination_nuisance": false,
      "velocity_jitter_kms": 0.0,
      "radial_correlation": 0.0,
      "gas_scale_sigma": 0.0
    },
    {
      "name": "distance_published_frame",
      "description": "No model-H0 rescaling of Hubble-flow distance-prior centres.",
      "sigma_ln_ml": 0.25,
      "ml_mode": "coupled",
      "distance_floor_frac": 0.03,
      "distance_prior": "gaussian",
      "hubble_prior_mode": "published",
      "inclination_nuisance": false,
      "velocity_jitter_kms": 0.0,
      "radial_correlation": 0.0,
      "gas_scale_sigma": 0.0
    },
    {
      "name": "distance_floor_10pct",
      "description": "Minimum fractional distance uncertainty increased to 10 per cent.",
      "sigma_ln_ml": 0.25,
      "ml_mode": "coupled",
      "distance_floor_frac": 0.1,
      "distance_prior": "gaussian",
      "hubble_prior_mode": "model_h0_rescale",
      "inclination_nuisance": false,
      "velocity_jitter_kms": 0.0,
      "radial_correlation": 0.0,
      "gas_scale_sigma": 0.0
    },
    {
      "name": "distance_student_t3",
      "description": "Robust Student-t(3) distance-pull penalty.",
      "sigma_ln_ml": 0.25,
      "ml_mode": "coupled",
      "distance_floor_frac": 0.03,
      "distance_prior": "student_t3",
      "hubble_prior_mode": "model_h0_rescale",
      "inclination_nuisance": false,
      "velocity_jitter_kms": 0.0,
      "radial_correlation": 0.0,
      "gas_scale_sigma": 0.0
    },
    {
      "name": "inclination_nuisance",
      "description": "First-order line-of-sight inclination nuisance with the SPARC inclination error.",
      "sigma_ln_ml": 0.25,
      "ml_mode": "coupled",
      "distance_floor_frac": 0.03,
      "distance_prior": "gaussian",
      "hubble_prior_mode": "model_h0_rescale",
      "inclination_nuisance": true,
      "velocity_jitter_kms": 0.0,
      "radial_correlation": 0.0,
      "gas_scale_sigma": 0.0
    },
    {
      "name": "velocity_jitter_5kms",
      "description": "Five km/s non-circular-motion jitter added in quadrature.",
      "sigma_ln_ml": 0.25,
      "ml_mode": "coupled",
      "distance_floor_frac": 0.03,
      "distance_prior": "gaussian",
      "hubble_prior_mode": "model_h0_rescale",
      "inclination_nuisance": false,
      "velocity_jitter_kms": 5.0,
      "radial_correlation": 0.0,
      "gas_scale_sigma": 0.0
    },
    {
      "name": "radial_correlation_rho035",
      "description": "Within-galaxy AR(1) radial residual correlation fixed at rho=0.35.",
      "sigma_ln_ml": 0.25,
      "ml_mode": "coupled",
      "distance_floor_frac": 0.03,
      "distance_prior": "gaussian",
      "hubble_prior_mode": "model_h0_rescale",
      "inclination_nuisance": false,
      "velocity_jitter_kms": 0.0,
      "radial_correlation": 0.35,
      "gas_scale_sigma": 0.0
    },
    {
      "name": "gas_scale_20pct",
      "description": "Galaxy-level gas contribution scaling with a 20 per cent log prior.",
      "sigma_ln_ml": 0.25,
      "ml_mode": "coupled",
      "distance_floor_frac": 0.03,
      "distance_prior": "gaussian",
      "hubble_prior_mode": "model_h0_rescale",
      "inclination_nuisance": false,
      "velocity_jitter_kms": 0.0,
      "radial_correlation": 0.0,
      "gas_scale_sigma": 0.2
    },
    {
      "name": "combined_conventional",
      "description": "Combined wide separate stellar, robust distance, inclination, jitter, radial-correlation and gas terms.",
      "sigma_ln_ml": 0.5,
      "ml_mode": "separate",
      "distance_floor_frac": 0.1,
      "distance_prior": "student_t3",
      "hubble_prior_mode": "model_h0_rescale",
      "inclination_nuisance": true,
      "velocity_jitter_kms": 5.0,
      "radial_correlation": 0.35,
      "gas_scale_sigma": 0.2
    }
  ],
  "persistence_gate": "PLAMB profile objective worse in >=80% scenarios and in combined_conventional",
  "identifiability_gate": "current likelihood must contain an independently derived FR charge-conjugation transformation; current result is NOT IDENTIFIABLE",
  "selection_inference": "development-set descriptive only; no discovery p-value"
}
```
