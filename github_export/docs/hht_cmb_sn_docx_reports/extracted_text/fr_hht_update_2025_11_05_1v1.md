# FR_HHT_Update_2025-11-05_1V1

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\FR_HHT_Update_2025-11-05_1V1.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2013-12-23T23:15:00Z`

## Extracted Text

FR–HHT Locking & Chirp — Interim Update

Generated: 2025-11-04 20:03  |  Version: 1V1

Executive summary

HHT IMF2 locking: baseline axis ln(1+z) (α=1.0) outperforms FR-style axis (α=0.5) on |r|; p-values tied at 0.00498 for the tested pairs.

FR SN+BAO (spline Λ): allowing c(z) in BAO (‘c_of_z’) gives a modest χ² improvement vs constant-c in BAO (‘c0’).

Piecewise-Λ ‘discrete events’: Δz→1e−4 reproduces the continuous result within rounding; coarse Δz=0.015 is slightly worse (as expected).

OU Λ(z) sensitivity: results stable from amp=1e−120 … 1e−70; degrade clearly by amp≈1e−2. Corr-length sweep (0.05–0.60) showed no material change.

PBH Λ↔A_acc coupling γ∈{0,0.02,0.05,0.10}: no change in best-fit or χ² with SN-only — suggests the tested γ range is too small or SN-only lacks constraining power.

Soft-link (energy–luminosity) with flat FR and log1pz scalings: initial κ={−1.0, −0.75}, σ=0.10 runs converge to identical solutions; likely penalty not binding at this σ.

HHT IMF2 locking (α=0.5 vs α=1.0)

Setup: ceemdan, IMF=2, packs=bao_hht_2bin_z0p7_longbest, seeds={3,7,11}, ns=200.

Summary table: mean|r| → α=0.5: 0.8790809752, α=1.0: 0.9093863675; median p both 0.0049751244; win rate (α=0.5) = 0/3.

FR (spline Λ) with SN+BAO — BAO c‑mode comparison

c0: χ²=2017.2954, χ²/dof=1.18247, AIC=2031.2954, BIC=2069.4174, H0=59.576, Ωm=0.3435, Ωk=0.00779, w_eff≈−0.912.

c_of_z: χ²=2014.5488, χ²/dof=1.18086, AIC=2028.5488, BIC=2066.6709, H0=60.579, Ωm=0.1333, Ωk=0.25117, w_eff≈−1.080.

‘c_of_z’ gains ~2.75 in χ² vs ‘c0’ with these settings; parameters shift toward higher |Ωk| and more negative w_eff.

Piecewise‑Λ: continuous vs discrete‑events

continuous: χ²=1975.2424, AIC=1989.2424 (baseline)

discrete Δz=0.015: χ²=1980.4284, AIC=1994.4284 (coarser events worse)

discrete Δz=1e−4: χ²=1975.2100, AIC=1989.2100 (matches continuous within rounding)

OU Λ(z) sensitivity

Correlation‑length sweep (amp=1e−120): corrz ∈ {0.05, 0.15, 0.30, 0.60} → identical best‑fit & χ²=1974.8730 (SN‑only).

Amplitude sweep at fixed corrz=0.30 and fixed parameters (‘sampler=none’): 1e−120…1e−70 unchanged; 1e−2 degrades (χ²≈2089.88).

Probe MCMC at amp=1e−2 (short chain): χ²≈2020.38.

PBH Λ↔A_acc coupling (γ) — SN‑only

γ ∈ {0.00, 0.02, 0.05, 0.10} → χ²=1974.5117 for all; no effect across this range for SN‑only.

Action: include BAO/Planck and/or widen γ grid.

Soft‑link experiments (flat FR, log1pz scalings)

κ=−1.0, σ=0.10 → χ²=2005.9505, AIC=2017.9505;

κ=−0.75, σ=0.10 → χ²=2005.9505, AIC=2017.9505;

Identical outcomes suggest the soft constraint isn’t binding at σ=0.10. Action: repeat for σ∈{0.05, 0.02, 0.01}.

Data & environment notes

Planck prior: ‘planck_distance_prior.json’ missing in /Users/boyde/Downloads/ — add file to enable --planck and --planck-c-mode tests.

emcee: use nwalkers ≥ 14 for 7 parameters; corner is optional.

HHT α‑compare DOCX needed python‑docx (installed).

Recommended next steps

HHT: α∈{0.6…1.4} (step 0.1) for IMF2/IMF3; keep surrogate trials; report |r| and p‑value distributions per pad/method.

BAO/Planck: once Planck JSON present, compare c0 vs c_of_z under Planck distance priors; inspect Ωk and w_eff.

Discrete events: quantify convergence vs Δz ladder (1e−3, 5e−4, 1e−4) with AIC/BIC.

OU Λ(z): place informative priors on corrz and amplitude; test rejection thresholds for large amplitudes (≥1e−4).

PBH coupling: extend γ grid and include BAO/Planck; track posteriors of A_acc, n_acc, and γ.

Soft‑link: sweep κ∈[−1.5,+1.5] with σ∈{0.05,0.02,0.01}; ensure nsteps≥2500.

Repro command snippets

# FR + SN+BAO (spline Λ): c_of_z%run '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \  --model FR --data '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \  --cov '/Users/boyde/Downloads/Pantheon+SH0ES_STAT+SYS.cov' \  --bao '/Users/boyde/Downloads/bao_long.csv' --rd 147.09 \  --sampler emcee --nwalkers 24 --nsteps 800 --nburn 200 --nint 220 --seed 73 \  --fluct-Lambda --Lambda-model spline --Lambda-amp 1.0e-120 --Lambda-knots 7 \  --bao-c-mode c_of_z --out 'plamb_runs/SNBAO_baocmode_c_of_z'

# Piecewise-Λ: Δz→1e−4 matches continuous%run /Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py \  --model FR --data /Users/boyde/Downloads/Pantheon+SH0ES.dat \  --cov  /Users/boyde/Downloads/Pantheon+SH0ES_STAT+SYS.cov \  --sampler none --nint 400 \  --H0 82.5742 --Om 0.68207 --Ol 0.96443 --A_acc -0.25624 --n_acc -5.84476 \  --gamma_c -0.02980 --epsilon_M 0.38324 \  --fluct-Lambda --Lambda-model piecewise --Lambda-amp 1e-120 --Lambda-knots 8 \  --discrete-events --event-size 0.0001 \  --out 'plamb_runs/SN_piecewise_disc_at_cont_DZ1e-4'
