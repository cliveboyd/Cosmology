# HHT resonance diagnostic

## Purpose

This is a diagnostic side test for the SU2/chirality programme. It checks whether a Hilbert-Huang workflow can recover a controlled SU2 -> chi_dot -> U(1) resonance chain before applying the method to noisy cosmology residuals or parity observables.

It is not a SN/BAO/Planck likelihood and it is not a full U(2) model.

## Synthetic setup

- omega_Q: 1
- m_chi: 1.05
- k_mode: 0.5
- first-band mismatch (2k - omega_Q)/omega_Q: 0
- bias: 0.03
- damping: 0
- samples: 8192
- duration periods: 80

## Coupling result

| metric | value |
|---|---:|
| Q-chi 1:1 phase-locking value | 0.99823 |
| chi-U1 plus 1:2 phase-locking value | 0.52107 |
| chi-U1 minus 1:2 phase-locking value | 0.6513 |
| max chi-U1 1:2 phase-locking value | 0.6513 |
| U1 plus energy growth | 3.7198e+07 |
| U1 minus energy growth | 2.6429e+07 |
| raw median final helicity asymmetry | 0.25726 |
| intrinsic chirality flag | True |
| resonance score | 0.87743 |

Interpretation: strong HHT recovery of the synthetic resonance chain, with bias-supported helicity imbalance.

## Selected IMF components

| series | selected IMF | selected frequency | expected frequency |
|---|---:|---:|---:|
| Q | 1 | 0.15943 | 0.15915 |
| chi_dot | 1 | 0.15881 | 0.15915 |
| U1_plus | 1 | 0.075613 | 0.079577 |
| U1_minus | 1 | 0.078974 | 0.079577 |

## Notes for project use

- A high score here only proves that HHT can detect the synthetic resonance pattern.
- Applying this to cosmology residuals needs surrogate tests: shuffled redshift/order controls, LCDM-vs-SU2 comparisons, and repeated noise realizations.
- In the exact zero-bias periodic case, the SU2 -> chi -> U(1) mechanism is not intrinsically chiral; a raw finite-time helicity asymmetry should be treated as a phase-window effect unless a bias or other symmetry-breaking ingredient is present.
- The natural next data-facing targets are residual series ordered by log(1+z), lookback time, chiral magnetic-field diagnostics, or parity-sensitive 4PCF/GW observables.

## Output files

- hht_resonance_synthetic_signals.csv
- hht_resonance_imf_summary.csv
- hht_resonance_coupling_summary.csv
- hht_resonance_overview.png
- hht_resonance_hilbert_spectrum.png
