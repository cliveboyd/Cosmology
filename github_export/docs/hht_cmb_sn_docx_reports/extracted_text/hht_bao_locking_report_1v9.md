# HHT_BAO_Locking_Report_1V9

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\HHT_BAO_Locking_Report_1V9.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2025-10-12T02:04:00Z`

## Extracted Text

Cosmology HHT–BAO Locking Report

(Revision 1V9)

Executive Summary

Application of  Hilbert–Huang Transform (HHT) to SN residuals and subsequently test ‘locking’ to BAO effective scales.

New analyses (energy sweep, locking heatmaps, sliding-τ windows, and frequency-gated tests) consolidate that IMF2 carries a robust, statistically significant BAO-linked signal, while IMF1/IMF3 exhibit high but broadly shared co-variation that does not alter the IMF2 conclusion.

Data & Method (brief)

Signals: SN residuals vs log(1+z); BAO observables DM/rd and DH/rd.Decomposition: EMD/EEMD/CEEMDAN; instantaneous frequency via Hilbert phase; uniform τ-grid with cubic spline.Statistics: Pearson r for locking, IMF2 energy ratio, phase-randomized surrogates (ns≈1000 per job), BH-FDR across cells/windows.Grids: pads ∈ {0.08, 0.10, 0.12, 0.15}, ngrid ∈ {192, 256, 320, 384}; seeds ∈ {3,7, 11,23}.

Key Findings (1V9)

• IMF2 is consistently strongest: median r ≈ 0.973 at (pad=0.15, ngrid=384) with q≤0.0015 across methods; IMF2 energy p≈0.0015 (median across seeds).

• Sliding-τ windows (IMF2, CEEMDAN): 17/18 windows significant at q≤0.10; peak τ-centers ≈ 0.615–0.815 (z≈0.85–1.26).

• Method-agnostic: EMD/EEMD/CEEMDAN agree on top cells; ‘strong’ ensembles (more trials/noise) are numerically identical to baseline.

• Adjacent modes: IMF1 (18/18) and IMF3 (16/18) windows significant, but frequency-gated tests (15–85% and 45–55% bands) change r by ≤0.006, indicating multi‑scale covariance rather than IMF2 leakage.

Recommended Defaults

Method: CEEMDAN (EEMD/EMD equivalent here)Grid: pad=0.15, ngrid=384IMF: IMF2 for primary metrics (locking r, energy ratio)Surrogates: ns=1000; BH-FDR across grid/windowsProgress: --progress auto (rich/tqdm)

Consensus best cells (energy + locking):

Robustness

• Seed stability: near-zero IQR at top cells; results invariant across seeds (3,7, 11,23).

• Ensemble strength: ‘strong’ EEMD/CEEMDAN equal to baseline (Δr=0, ΔIQR=0 across grid).

• Frequency-gating: IMF1/IMF3 Δr≈0 to −0.006 as gate widens; IMF2 unchanged.

Next Steps

• Export paper figures: energy heatmap, locking heatmap, sliding-τ curve, and specificity panel.

• Redshift-sliced scans to localize contributions and rule out edge artifacts.

• Optional: distance-based gating (|k_imf−k2|/|k2|<δ) for an even stricter specificity test.

Baseline context: Cosmology HHT–BAO Locking Report (Revision 1V8).

## Extracted Tables

### Table 1

Method | pad | ngrid | IMF2 energy p (med) | Locking r (med)
EMD | 0.15 | 384 | 0.0015 | 0.973
EEMD | 0.15 | 384 | 0.0015 | 0.973
CEEMDAN | 0.15 | 384 | 0.0015 | 0.973
