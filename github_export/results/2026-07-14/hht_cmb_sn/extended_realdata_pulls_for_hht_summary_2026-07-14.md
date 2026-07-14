# Extended real-data HHT export

CSV: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\hht_resonance\extended_realdata_pulls_for_hht.csv`

## Scope

This file combines standardized pull series from SN bins, BAO rows, and cosmic chronometer H(z) rows. It is intended for HHT screening only, not as a replacement for the formal likelihood.

- z range: 0.014244 to 2.33
- HHT rows: 75
- BAO: 14 rows
- CC: 32 rows
- SN_bin: 29 rows

## Included models

| label | model | note |
|---|---|---|
| LCDM_baocov | LCDM | Matched LCDM SN+BAO covariance+Planck baseline. |
| SU2_baocov_4D | SU2 | Matched SU2 4D SN+BAO covariance+Planck branch. |
| SU2R_baocov_4D | SU2R | Matched reparameterized SU2R 4D branch. |
| SU2_wa_free_baocov | SU2 | Matched SU2 run with wa_chi free. |

## Skipped models

| label | reason |
|---|---|
| SU2R_long_baocov | missing C:\Users\clive\Documents\Cosmology\plamb_runs\next_su2_tests\runs\baocov_SU2R_long_SN_BAO_Planck\SU2R_emcee_bestfit.txt |

## Recommended HHT commands

```powershell
python C:\Users\clive\Documents\Cosmology\diagnose_hht_resonance.py --mode csv --input-csv C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\hht_resonance\extended_realdata_pulls_for_hht.csv --x-col x_hht --y-col joint_delta_pull_SU2_baocov_4D_minus_LCDM_baocov --label extended_SU2_minus_LCDM
python C:\Users\clive\Documents\Cosmology\diagnose_hht_resonance.py --mode csv --input-csv C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\hht_resonance\extended_realdata_pulls_for_hht.csv --x-col x_hht --y-col joint_delta_pull_SU2R_baocov_4D_minus_LCDM_baocov --label extended_SU2R_minus_LCDM
python C:\Users\clive\Documents\Cosmology\diagnose_hht_resonance.py --mode csv --input-csv C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\hht_resonance\extended_realdata_pulls_for_hht.csv --x-col x_hht --y-col joint_pull_SU2_baocov_4D --label extended_SU2_pull
```

## Caution

The combined series mixes different probes after pull-standardization. A detected IMF would be a feature candidate only; it must be checked probe-by-probe and against shuffled-order controls.
