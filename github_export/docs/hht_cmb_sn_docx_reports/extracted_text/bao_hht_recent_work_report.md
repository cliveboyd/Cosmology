# BAO_HHT_recent_work_report

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\BAO_HHT_recent_work_report.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2025-10-09T04:13:00Z`

## Extracted Text

Recent BAO–HHT Integration: Files, Results, and Next Steps

This report summarizes the recent development of helper tooling, updates to covariance handling, and key BAO-only and SN+BAO results, plus implications for HHT and recommended next actions.

Executive Summary

• Implemented helper (bao_hht_helper.py) for constant-ρ BAO covariance, residuals/whitening, χ², per‑z pairs, and analytic ρ*.

• BAO‑only prefers ρ* ≈ −0.118 ± 0.129; SN+BAO prefers ρ* ≈ +0.174 ± 0.064.

• Using ρ = +0.143338 in the combined run already captures ~97% of the attainable χ² reduction vs ρ=0.

• Whitened BAO residuals are now exportable for HHT in the same style as SNe.

BAO‑Only Results (FR, flat; r_d = 147.090 Mpc)

Quadratic fit near the minimum gives ρ* = −0.118 ± 0.129 (χ²* ≈ 59.60).

Combined Pantheon+SH0ES + BAO Results (FR, flat; r_d = 147.090 Mpc)

Analytic optimization at the combined best-fit cosmology gives ρ* = +0.174 ± 0.064 with χ²_BAO(ρ*) ≈ 236.1908.

Δχ² vs ρ = +0.143338 is only ≈ 0.227; re‑running at exactly ρ* would change posteriors negligibly.

Relation to HHT Analysis

• Exported whitened BAO residuals (using the same covariance as the fit) enable HHT decomposition parallel to SN analysis.

Per‑z χ² highlights z≈0.845 and 0.700 as heavy contributors—matching suspect zones in SN residual structure.

• The ρ sign flip (BAO‑only negative → SN‑anchored positive) is expected because ρ de‑correlates DM and DH within each pair, and the optimal de‑correlation depends on cosmology.

Implications

1) Constant‑ρ block covariance adequately summarizes BAO DM–DH cross‑correlation for N=14.

• Parameter inferences are robust to small mis‑specification of ρ near the optimum.

• Whitened BAO residuals align BAO with SN HHT workflows for mode‑level comparisons.

• z≈0.7–0.9 remains the primary BAO lever arm—focus HHT inspection here.

Next Steps (Extended Redshift Focus)

• A) Extend BAO coverage: add/update high‑z BAO (Lyα at z≳2), ensure proper DM/DH pairing and uncertainties.

•    – Introduce ρ(z): two‑bin (z<0.9 vs z≥0.9) or linear slope; compare AIC/BIC to constant ρ.

•    – r_d sensitivity: re‑run fit‑rho and combined fits at r_d ± Δr_d.

• B) HHT integration: run HHT on exported whitened BAO residuals; compare mode energies/frequencies to SN modes.

•    – Perform leave‑one‑z‑pair‑out HHT to gauge influence of localized BAO pairs.

• C) Model checks: try c_model=log1pz vs linear_z; cautiously add CC or Planck R/rs to monitor tension migration.

• D) Reproducibility: export covariance.txt, perz_pair_report.csv, and whitened_residuals.csv for each major run.

Paths & Artifacts

• BAO CSV: /Users/boyde/Downloads/bao_MASTER_LONG.csv

• SN data/cov: /Users/boyde/Downloads/Pantheon+SH0ES.dat ; Pantheon+SH0ES_STAT+SYS.cov

• Helper: ~/.spyder-py3/bao_hht_helper.py

• Cov examples: cov_sweeps/bao_MASTER_cov_blockrho-0.12.txt, cov_sweeps/bao_MASTER_cov_blockrhop0.143338.txt

• Outputs: plamb_runs/COMBINED_BAO_SN_rho_p0.143338, plamb_runs/tools/bao_hht_at_combined_best/

## Extracted Tables

### Table 1

ρ | H0 (km/s/Mpc) | Ωm | χ²_BAO | Notes
−0.10 | 68.731 | 0.29276 | 59.624 | BAO-only min in grid around −0.1
−0.12 | 68.7337 | 0.292775 | 59.603 | Consistent with analytic minimum

### Table 2

Setting | H0 (km/s/Mpc) | Ωm | χ²_SN / dof_SN | χ²_BAO | χ²_total / dof_total
ρ=+0.143338 | 73.1400 | 0.235889 | 2098.521 / 1699 ≈ 1.235 | 236.4175 | 2334.938 / 1713 ≈ 1.363
