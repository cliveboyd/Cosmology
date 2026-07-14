# SkyMap HHT Integration Report A4

- Source: `SkyMap_HHT_Integration_Report_A4.docx`
- Source size: 30908 bytes
- Source modified: 2025-11-01T10:05:04
- Extracted: 2026-07-14
- Word count estimate: 912

## Extracted Text
SkyMap Dipole Sweep

Results, Context, and HHT/Chirp Integration

Report date: 2025-10-31

1. Executive Summary

We ran SkyMap dipole sweeps on the BAO ratio field R_BAO = DH_over_rd / DM_over_rd at NSIDE=32 with FWHM=6° smoothing, for Galactic latitude masks |b| ≥ {0,10,20,30}°. The catalog does not include plain DH/DM, so both the “FR-safe” sweep and the “raw” sweep used the ratio inputs.

Across masks, the measured dipole amplitudes A_eff were ~1.0–1.2×10⁻³, with permutation p-values in the ~0.036–0.263 range for the specific examples below. The best axis-closest case (minimum Δ to the FR reference axis) occurred at |b|≥30° (Δ≈60.55°).

2. Data and Configuration

• Catalog: /Users/boyde/Downloads/unified_catalog_REAL_MIN.csv (missing DH/DM; contains DH_over_rd and DM_over_rd).

• Mode: R_BAO = DH_over_rd / DM_over_rd.

• Smoothing: NSIDE=32, FWHM=6°.

• Permutations: 1000 per map for p-value estimation (p ≡ Pr(A_null ≥ A_obs)).

• Reference axes for deltas: FR (l,b)=(264°, +48°), KdS (l,b)=(300°, +30°).

• Runs created both JSON summaries and logs; the reporter aggregates both (JSON preferred for Δ_FR, Δ_KdS).

3. Results: sweeps_real (FR mode)

|b| mask (°) | FWHM (°) | A_eff | l (°) | b (°) | p_perm | Δ_FR (°) | Δ_KdS (°) | axis Δ_min (°)

0 | 6 | 0.001155740 | 151.600934 | -62.963841 | 0.035964 | 38.937337 | 38.680565 | 75.836429

10 | 6 | 0.001122610 | 108.159611 | -73.925408 | 0.094905 | 27.977049 | 44.345204 | 88.681841

20 | 6 | 0.001162770 | 76.548440 | -69.867344 | 0.135864 | 22.164627 | 46.698044 | 81.941380

30 | 6 | 0.000973177 | 168.527768 | -31.537745 | 0.262737 | 63.698486 | 41.379430 | 60.548418

Notes: These match your recent FR sweep outputs (report_FR.csv).

4. Results: sweeps_raw (with --fr-compatible rerun)

After re-running the raw sweep with --fr-compatible, the numerical dipole estimates matched the FR set because both ultimately used the same ratio inputs. The reporter may list duplicates when both JSON and logs exist; JSON rows include ‘fr_compatible=True’, while log-derived rows default to False.

|b| mask (°) | FWHM (°) | A_eff | p_perm | rbao_mode | fr_compatible

0 | 6 | 0.001155740 | 0.035964 | DH_over_rd/DM_over_rd | True/False (both present)

10 | 6 | 0.001122610 | 0.094905 | DH_over_rd/DM_over_rd | True/False (both present)

20 | 6 | 0.001162770 | 0.135864 | DH_over_rd/DM_over_rd | True/False (both present)

30 | 6 | 0.000973177 | 0.262737 | DH_over_rd/DM_over_rd | True/False (both present)

Deduping tip: keep one row per (mask,fwhm); prefer JSON over logs to retain Δ_FR/Δ_KdS.

5. Interpretation & Implications

• R_BAO dipole amplitude ~10⁻³ across masks.

• Axis varies with mask; closest to FR axis at |b|≥30° (≈60.55°).

• p-values ~0.036–0.263 here; increase permutations and vary NSIDE/FWHM to refine significance.

• No numerical difference between ‘FR’ and ‘raw’ runs because both used the ratio inputs in this catalog.

6. Where HHT & Chirp Structures Fit

Context: HHT on SN residuals shows IMF2/IMF3 features and IMF2 ‘locking’ to BAO-relevant scales; tested via robustness grids and surrogate analysis.

Connection to R_BAO dipole: Anisotropic coupling from BAO‑linked IMF2 energy could project to low‑ℓ dipole. Chirp vs log(1+z) suggests scale drift; mask‑dependent axis shifts are compatible with footprint/selection effects. Combine map‑level dipole tests with HHT ridge/surrogate significance for a joint detection framework.

Operationally: use per‑mask CSV snapshots; recompute HHT ridge metrics within each mask’s footprint; sweep NSIDE/FWHM alongside HHT grids.

7. Relation to the PLamb Theory Document

The PLamb document (FR/KdS focus, accretion‑driven acceleration A_acc) connects parameter‑level inferences to sky‑level anisotropy tests. Current working fit (subject to refinement): H₀≈71.25, Ωₘ≈0.279, Ω_Λ≈0.739, A_acc≈0.030. We use FR/KdS reference axes to compute Δ_FR and Δ_KdS, tying the observed R_BAO dipole to theory.

8. Note on CMB ‘Ripples’ from Bubble Collisions (Hypothesis)

A lecture suggested possible large‑scale CMB ripples from multiverse bubble collisions. Interesting but speculative here; if pursued, cross‑check axis alignments and low‑ℓ patterns against the R_BAO dipole with robust null tests.

9. Recommended Next Steps

• Obtain a catalog with true DH and DM to compare against the ratio mode directly.

• Increase permutations (e.g., 10k) and sweep NSIDE∈{16,32,64}, FWHM∈{4°,6°,8°}.

• Add deduping (prefer JSON over logs) in the reporter to avoid duplicates.

• Build joint panels per mask: {A_eff, Δ_FR, Δ_KdS, p} ↔ {IMF2 energy, chirp slope, surrogate p}.

• If CMB hypothesis pursued, design a falsifiable cross‑correlation test.

Appendix A — Reproducibility

%run "/Users/boyde/.spyder-py3/sm_quickgrid.py" --catalog "/Users/boyde/Downloads/unified_catalog_REAL_MIN.csv" --rbao-mode DH_over_rd/DM_over_rd --fr-compatible --nside 32 --fwhm-smooth-deg 6 --gal-mask-list 0,10,20,30 --use-gal-mask-in-analysis --use-gal-mask-in-display --perm-n 1000 --outdir "/Users/boyde/.spyder-py3/plamb_runs/SkyMap/figs/sweeps_real"

%run "/Users/boyde/.spyder-py3/sm_sweep_report.py" --dir "/Users/boyde/.spyder-py3/plamb_runs/SkyMap/figs/sweeps_real" --out "/Users/boyde/.spyder-py3/plamb_runs/SkyMap/figs/sweeps_real/report_FR" --truth-A 0.0010 --truth-axis "230,-10"

%run "/Users/boyde/.spyder-py3/sm_quickgrid.py" --catalog "/Users/boyde/Downloads/unified_catalog_REAL_MIN.csv" --rbao-mode DH/DM --fr-compatible --nside 32 --fwhm-smooth-deg 6 --gal-mask-list 0,10,20,30 --use-gal-mask-in-analysis --use-gal-mask-in-display --perm-n 1000 --outdir "/Users/boyde/.spyder-py3/plamb_runs/SkyMap/figs/sweeps_raw"

%run "/Users/boyde/.spyder-py3/sm_sweep_report.py" --dir "/Users/boyde/.spyder-py3/plamb_runs/SkyMap/figs/sweeps_raw" --out "/Users/boyde/.spyder-py3/plamb_runs/SkyMap/figs/sweeps_raw/report_raw" --truth-A 0.0010 --truth-axis "230,-10"
