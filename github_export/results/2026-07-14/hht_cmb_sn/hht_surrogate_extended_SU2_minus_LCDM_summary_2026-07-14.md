# HHT surrogate control

Input: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\hht_resonance\extended_realdata_pulls_for_hht.csv`
x column: `x_hht`
y column: `joint_delta_pull_SU2_baocov_4D_minus_LCDM_baocov`
surrogates: 100

## Original metrics

| metric | value |
|---|---:|
| n IMFs | 6 |
| max IMF energy fraction | 0.376515 |
| sum IMF energy fraction | 0.782199 |
| dominant IMF frequency | 24.7666 |

## Shuffled-control comparison

| metric | original | surrogate median | surrogate p90 | upper-tail p | percentile |
|---|---:|---:|---:|---:|---:|
| max IMF energy fraction | 0.376515 | 0.57553 | 0.703517 | 0.970297 | 0.03 |
| sum IMF energy fraction | 0.782199 | 0.977461 | 1.13933 | 0.960396 | 0.04 |

## Interpretation

The original IMF energy is not unusual relative to shuffled-order controls. Treat the HHT IMFs as likely sampling/probe-mixing structure rather than resonance evidence.
