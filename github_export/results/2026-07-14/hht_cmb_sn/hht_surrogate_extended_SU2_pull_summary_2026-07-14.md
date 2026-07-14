# HHT surrogate control

Input: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\hht_resonance\extended_realdata_pulls_for_hht.csv`
x column: `x_hht`
y column: `joint_pull_SU2_baocov_4D`
surrogates: 100

## Original metrics

| metric | value |
|---|---:|
| n IMFs | 4 |
| max IMF energy fraction | 0.543719 |
| sum IMF energy fraction | 0.901987 |
| dominant IMF frequency | 30.5909 |

## Shuffled-control comparison

| metric | original | surrogate median | surrogate p90 | upper-tail p | percentile |
|---|---:|---:|---:|---:|---:|
| max IMF energy fraction | 0.543719 | 0.551616 | 0.677 | 0.554455 | 0.45 |
| sum IMF energy fraction | 0.901987 | 0.977922 | 1.10776 | 0.80198 | 0.2 |

## Interpretation

The original IMF energy is not unusual relative to shuffled-order controls. Treat the HHT IMFs as likely sampling/probe-mixing structure rather than resonance evidence.
