# FR + HHT Update Report (Alpha & Skymap/Chirp) — 1V0

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\FR_HHT_Update_2025-11-04_1V0.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2025-11-04T00:00:00Z`
- DOCX modified: `2025-11-04T00:00:00Z`

## Extracted Text

FR + HHT Update Report (Alpha & Skymap/Chirp) — 1V0

Date: 04 Nov 2025 (Europe/Stockholm)

Executive Summary

• IMF2 HHT ‘locking’ favors ln(1+z) with α = 1.0 over α = 0.5 based on paired comparisons (mean |r|: 0.909 vs 0.879; identical median p owing to the test floor).

• Skymap/locking heatmaps from the 2‑bin BAO pack (DH_over_rd, DM_over_rd) are consistent with α = 1.0.

• FR + OU(Λ) with amp = 1e‑120 is robust to corrz (0.05–0.60) at fixed seed; χ²/dof ≈ 1.166 (SN‑only).

• Piecewise Λ: discrete‑events Δz = 1e‑4 matches continuous (χ² ≈ 1975.21), validating the integrator.

• PBH with Λ↔A_acc coupling γ ∈ {0.00, 0.02, 0.05, 0.10} yields identical SN‑only best fits (χ² ≈ 1974.51).

• SN+BAO FR with c_of_z slightly improves χ² vs c0 (2014.55 vs 2017.30) under spline‑Λ fluctuations.

• Soft‑link EM↔c (flat FR, log1pz) with κ≈−1, σ=0.10 increases χ² (~2005.95) — currently tension‑adding.

HHT Locking & Alpha Comparison

• Tools: hht_locking_heatmap.py; hht_alpha_compare_quick.py; hht_alpha_compare_report.py (1V0).

• CEEMDAN, IMF2, seeds {3,7,11}, ns=200 → α=0.5 win‑rate 0% vs α=1.0 on |r| and p; mean |r| 0.879 vs 0.909.

• Interpretation: IMF2 locking favors ln(1+z) (α=1).

Skymap / Chirp Notes

• HHT locking heatmaps (‘skymaps’) show stronger IMF2 ridge coherence along α=1.0 axis.

• Next: extend to IMF3 for chirp energy tracking and BAO (DH/DM) stability vs seeds/method.

FR — OU(Λ) Robustness & Amplitude

• FR + OU(Λ), amp=1e‑120; corrz ∈ {0.05,0.15,0.30,0.60}; emcee: walkers=20, steps=700, burn=200.

• Best fit repeats across corrz: H0≈82.65, Ωm≈0.575, ΩΛ≈0.964, A_acc≈−0.250, n_acc≈−5.91, γ_c≈+0.053, ε_M≈+0.255.

• Amplitude sweep (no‑MCMC): unchanged up to 1e‑70; at 1e‑2 χ² degrades to ≈2089.88.

FR — Piecewise Λ: Continuous vs Discrete

• Piecewise Λ (knots=8): continuous χ²_SN≈1975.24; discrete Δz=0.015 χ²_SN≈1980.43; Δz=1e‑4 (at continuous best‑fit) χ²_SN≈1975.21.

• Conclusion: discrete formulation converges to the continuous limit as Δz→0.

PBH — Λ↔A_acc Coupling Sweep

• PBH + OU(Λ), amp=1e‑120, corrz=0.25, Λ↔A_acc coupling γ ∈ {0.00,0.02,0.05,0.10}.

• All runs converge to the same SN‑only best fit: χ²≈1974.51; suggests γ is unconstrained without BAO/Planck.

SN + BAO (FR, spline‑Λ)

• FR + BAO (spline‑Λ): c0 → χ²_total≈2017.30; c_of_z → χ²_total≈2014.55.

• Both above SN‑only best levels; c_of_z is modestly better.

EM↔c Soft Link Tests

• Soft link EM↔c (flat FR, dM=log1pz, c=log1pz): κ≈−1, σ=0.10 → χ²_total≈2005.95.

• Grid over (κ,σ) pending to quantify ΔAIC/ΔBIC and find a sweet‑spot prior.

Paths & Artefacts

• /Users/boyde/.spyder-py3/plamb_runs/tools/hht_imf2_compare/locking_alpha_compare.csv

• /Users/boyde/.spyder-py3/plamb_runs/tools/hht_imf2_compare/HHT_Locking_Alpha_Compare_1V0.docx

• /Users/boyde/.spyder-py3/plamb_runs/SNBAO_baocmode_c0/FR_emcee_bestfit.txt

• /Users/boyde/.spyder-py3/plamb_runs/SNBAO_baocmode_c_of_z/FR_emcee_bestfit.txt

• /Users/boyde/.spyder-py3/plamb_runs/SN_OU_amp1e-120_corr0.30/FR_emcee_bestfit.txt

• /Users/boyde/.spyder-py3/plamb_runs/SN_piecewise_cont/FR_emcee_bestfit.txt

• /Users/boyde/.spyder-py3/plamb_runs/SN_piecewise_disc/FR_emcee_bestfit.txt

• /Users/boyde/.spyder-py3/plamb_runs/SN_piecewise_disc_at_cont/FR_none_bestfit.txt

• /Users/boyde/.spyder-py3/plamb_runs/SN_piecewise_disc_at_cont_DZ1e-4/FR_none_bestfit.txt

• /Users/boyde/.spyder-py3/plamb_runs/PBH_couple_gamma_0.00/PBH_emcee_bestfit.txt

Next Steps

• Add Planck distance priors (planck_c_mode ∈ {none,R_only,R_and_rs}).

• Run IMF3 locking + chirp ridge detection; document surrogate significance.

• SN+BAO: vary rs_power and c‑model (log1pz, piecewise, saturating).

• Formalize (κ,σ) soft‑link grid; compare via ΔAIC/ΔBIC.

• Automate CSV→DOCX tables; embed corner/traces from each run.
