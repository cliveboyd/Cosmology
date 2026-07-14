# HHT surrogate control

Input: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\hht_resonance\extended_realdata_pulls_for_hht.csv`
x column: `x_hht`
y column: `joint_delta_pull_SU2R_baocov_4D_minus_LCDM_baocov`
surrogates: 100

## Original metrics

| metric | value |
|---|---:|
| n IMFs | 5 |
| max IMF energy fraction | 0.375538 |
| sum IMF energy fraction | 0.802933 |
| dominant IMF frequency | 24.7429 |

## Shuffled-control comparison

| metric | original | surrogate median | surrogate p90 | upper-tail p | percentile |
|---|---:|---:|---:|---:|---:|
| max IMF energy fraction | 0.375538 | 0.574293 | 0.702701 | 0.970297 | 0.03 |
| sum IMF energy fraction | 0.802933 | 0.969447 | 1.12829 | 0.940594 | 0.06 |

## Interpretation

The original IMF energy is not unusual relative to shuffled-order controls. Treat the HHT IMFs as likely sampling/probe-mixing structure rather than resonance evidence.
