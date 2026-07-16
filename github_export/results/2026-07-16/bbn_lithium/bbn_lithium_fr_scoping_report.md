# BBN Lithium FR Scoping Gate

Generated: 2026-07-16T21:00:27

## Claim Boundary

This is a first-pass emulator gate. It is suitable for ranking simple FR/no-expansion hypotheses and identifying vetoes, but any apparently viable region must be rerun with a full BBN nuclear network.

## Input Anchors

- D/H x 1e6 = 25.08 +/- 0.29
- Yp(He-4) = 0.245 +/- 0.003
- Li/H x 1e10 = 1.45 +/- 0.25
- eta10_CMB = 6.12 +/- 0.04
- Theory floors: D 3%, He 0.0005, Li 15%

## Executive Readout

- Standard BBN at CMB eta predicts Li/H x 1e10 = 4.72, which is 3.26 times the plateau value used here.
- The best clock-only FR-style shift at fixed CMB eta has S = 1, Delta N_eff = 0, and chi2 = 19.9.
- Allowing eta plus clock with the CMB prior gives eta10 = 6.113, S = 1, and chi2 = 19.8.
- A Li-only production suppression would need factor 0.305; a stellar-depletion interpretation needs surface depletion factor 3.25.
- If Li is treated only as a lower bound, SBBN+CMB is no longer penalised by Li; best chi2 = 0.957. That is an astrophysical/systematic resolution, not an FR signal.

## Model Ranking

| model                         |   best_chi2 |   best_eta10 |   best_clock_factor |   best_delta_neff |   best_li_suppression |   best_stellar_depletion |   pred_yD |   pred_Yp |   pred_yLi_surface |
|:------------------------------|------------:|-------------:|--------------------:|------------------:|----------------------:|-------------------------:|----------:|----------:|-------------------:|
| li_as_lower_bound_cmb         |    0.956505 |       6.12   |               1     |          0        |                 1     |                     1    |   2.45575 |  0.247228 |            4.72    |
| stellar_depletion_cmb         |    0.956553 |       6.12   |               1     |          0        |                 1     |                     3.25 |   2.45575 |  0.247228 |            1.45231 |
| li_production_suppression_cmb |    0.957496 |       6.12   |               1     |          0        |                 0.305 |                     1    |   2.45575 |  0.247228 |            1.4396  |
| eta_clock_free                |   17.9341   |       5.9    |               0.985 |         -0.182904 |                 1     |                     1    |   2.54158 |  0.244498 |            4.45392 |
| eta_only_free                 |   18.6695   |       5.9875 |               1     |          0        |                 1     |                     1    |   2.54328 |  0.247016 |            4.51783 |
| eta_only_with_cmb_prior       |   19.8273   |       6.1125 |               1     |          0        |                 1     |                     1    |   2.46057 |  0.247216 |            4.70844 |
| eta_clock_with_cmb_prior      |   19.8273   |       6.1125 |               1     |          0        |                 1     |                     1    |   2.46057 |  0.247216 |            4.70844 |
| sbbn_cmb                      |   19.9235   |       6.12   |               1     |          0        |                 1     |                     1    |   2.45575 |  0.247228 |            4.72    |
| clock_only_cmb_eta            |   19.9235   |       6.12   |               1     |          0        |                 1     |                     1    |   2.45575 |  0.247228 |            4.72    |

## Eta Inference Cross-Check

| probe   |   eta10_inferred |   eta10_sigma_obs_only | note                                                                            |
|:--------|-----------------:|-----------------------:|:--------------------------------------------------------------------------------|
| D/H     |          6.04    |              0.0436503 | D/H baryometer; matches PDG BBN concordance eta by construction.                |
| Li/H    |          3.39207 |              0.29242   | Lithium as direct measurement; lower-bound interpretation weakens this tension. |
| CMB     |          6.12    |              0.04      | Planck TT,TE,EE+lowE+lensing PDG value.                                         |

## Interpretation

The scoping gate favours the conventional reading: deuterium and CMB eta pin the baryon density near eta10 ~ 6, while lithium alone prefers a much lower eta if treated as an undepleted primordial measurement. Simple FR-style clock changes do not move Li independently; they also move D/H and He-4, so deuterium is the hard veto.

The useful next experiment is therefore a full-network BBN run with explicit FR parameters mapped to dimensionless nuclear quantities: eta_BBN/eta_CMB, H_BBN/H_SBBN or equivalent clock factor, neutron-proton mass difference, deuteron binding energy, and the 3He(alpha,gamma)7Be / 7Be(n,p)7Li / 7Li(p,alpha)4He reaction controls.

## Grid

- Smoke mode: False
- eta10: 321 values from 4 to 8
- clock: 241 values from 0.85 to 1.15
- li_suppression: 211 values from 0.15 to 1.2
- stellar_depletion: 241 values from 1 to 5

## Sources

- pdg_2025_bbn: https://ccwww.kek.jp/pdg/2025/reviews/rpp2025-rev-bbang-nucleosynthesis.pdf
- fields_lithium_review: https://ned.ipac.caltech.edu/level5/Sept15/Fields/Fields2.html
- steigman_fit_context: https://ned.ipac.caltech.edu/level5/Sept07/Steigman/Steigman2.html
- arxiv_recent_lithium_review: https://arxiv.org/abs/2304.08032

