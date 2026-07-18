# SPARC Antimatter-Candidate Conventional-Systematics Preregistration

- Written before outcome profiling: `2026-07-18T06:05:17.656977+00:00`
- Candidate galaxies: `6`
- Unique matched controls: `12`
- Reserved replication galaxies: `5`

## Claim Boundary

The candidate set was selected previously from model residuals. Its profile scores are therefore descriptive and cannot supply an unbiased discovery p-value. Controls are selected here without using any RAR or PLAMB outcome. The reserved replication set is not examined by this programme.

Gate 1 has already determined that the current rotation-curve likelihood is matter/antimatter non-identifiable. Consequently, persistence in this matrix means only that a conventional-systematics challenge did not erase a model-relative anomaly. It does not identify antimatter.

## Locked Candidate Set

`UGC06787`, `UGC03580`, `UGC02487`, `NGC5985`, `UGC09133`, `NGC2903`

## Locked Control Matching

2 unique f_D=1, Q<=2 controls per candidate are assigned by minimum total weighted robust-standardised Euclidean distance, with exact matching on the presence or absence of a SPARC bulge component. Candidate residuals and model scores are excluded from matching. The matching variables are morphology, luminosity, effective surface brightness, H I mass, flat speed, inclination, fractional distance error, bulge proxy, point count and quality.

| candidate | match_slot | control | weighted_match_distance |
| --- | --- | --- | --- |
| NGC2903 | 1 | NGC3521 | 0.608086 |
| NGC2903 | 2 | NGC4559 | 0.82885 |
| NGC5985 | 1 | NGC5033 | 1.37724 |
| NGC5985 | 2 | UGC03205 | 1.0972 |
| UGC02487 | 1 | UGC02885 | 1.57361 |
| UGC02487 | 2 | NGC6674 | 1.35873 |
| UGC03580 | 1 | UGC08699 | 0.743579 |
| UGC03580 | 2 | UGC03546 | 1.18552 |
| UGC06787 | 1 | UGC05253 | 1.29888 |
| UGC06787 | 2 | UGC02916 | 1.68114 |
| UGC09133 | 1 | NGC6195 | 1.53939 |
| UGC09133 | 2 | UGC06786 | 1.37052 |

## Locked Sensitivity Matrix

| name | description |
| --- | --- |
| baseline | Published coupled stellar scaling, Gaussian distance prior and 3 km/s error floor. |
| stellar_ml_tight | Tighter 15 per cent log stellar mass-to-light prior. |
| stellar_ml_wide | Wider 50 per cent log stellar mass-to-light prior. |
| disk_bulge_decoupled | Independent galaxy-level disc and bulge mass-to-light multipliers. |
| distance_published_frame | No model-H0 rescaling of Hubble-flow distance-prior centres. |
| distance_floor_10pct | Minimum fractional distance uncertainty increased to 10 per cent. |
| distance_student_t3 | Robust Student-t(3) distance-pull penalty. |
| inclination_nuisance | First-order line-of-sight inclination nuisance with the SPARC inclination error. |
| velocity_jitter_5kms | Five km/s non-circular-motion jitter added in quadrature. |
| radial_correlation_rho035 | Within-galaxy AR(1) radial residual correlation fixed at rho=0.35. |
| gas_scale_20pct | Galaxy-level gas contribution scaling with a 20 per cent log prior. |
| combined_conventional | Combined wide separate stellar, robust distance, inclination, jitter, radial-correlation and gas terms. |

## Profile Method

The all-Q2 July 14 MAP global RAR and PLAMB parameters are held fixed. For every candidate and control, local nuisances are re-profiled independently in every scenario. Positive Delta objective = objective_PLAMB - objective_RAR means PLAMB fits worse than RAR. Both the data term and the data-plus-prior objective are retained.

A candidate is labelled persistent conventional tension only when PLAMB is worse in at least 80 per cent of locked scenarios and remains worse in the combined-conventional scenario. This label is a robustness flag, not an antimatter classification.

## Reproducibility Locks

```json
{
  "date": "2026-07-18",
  "written_before_outcomes_utc": "2026-07-18T06:05:17.656977+00:00",
  "script": "C:\\Users\\clive\\Documents\\Cosmology\\github_export\\code\\sparc\\run_sparc_am_identifiability_systematics_2026-07-18.py",
  "script_sha256": "c5dae3af6cc6b63ac5db29f2ea29db4f8e3ffb1bdfd5d8d0fbc43ccca08740a9",
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
  "selection_inference": "development-set descriptive only; no discovery p-value",
  "optimisation": "deterministic multi-start L-BFGS-B: centre, per-parameter quartile starts and six hash-seeded joint starts; minimum objective retained"
}
```
