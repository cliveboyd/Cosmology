# SN residual export for HHT

CSV: `C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\hht_resonance\sn_residuals_for_hht.csv`

## Included models

| label | model | residual column | RMS [mag] | weighted RMS [mag] | note |
|---|---|---|---:|---:|---|
| LCDM_baocov | LCDM | `residual_LCDM_baocov` | 0.19425 | 0.15159 | Matched LCDM SN+BAO covariance+Planck baseline. |
| SU2_baocov_4D | SU2 | `residual_SU2_baocov_4D` | 0.19844 | 0.15548 | Matched SU2 4D SN+BAO covariance+Planck branch. |
| SU2R_baocov_4D | SU2R | `residual_SU2R_baocov_4D` | 0.19839 | 0.15542 | Matched reparameterized SU2R 4D branch. |
| SU2_wa_free_baocov | SU2 | `residual_SU2_wa_free_baocov` | 0.19842 | 0.15544 | Matched SU2 run with wa_chi free. |

## Skipped models

| label | reason |
|---|---|
| SU2R_long_baocov | missing C:\Users\clive\Documents\Cosmology\plamb_runs\next_su2_tests\runs\baocov_SU2R_long_SN_BAO_Planck\SU2R_emcee_bestfit.txt |

## Example HHT commands

```powershell
python C:\Users\clive\Documents\Cosmology\diagnose_hht_resonance.py --mode csv --input-csv C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\hht_resonance\sn_residuals_for_hht.csv --x-col log1pz --y-col residual_SU2_baocov_4D --label SU2_residual
python C:\Users\clive\Documents\Cosmology\diagnose_hht_resonance.py --mode csv --input-csv C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\hht_resonance\sn_residuals_for_hht.csv --x-col log1pz --y-col delta_residual_SU2_baocov_4D_minus_LCDM_baocov --label SU2_minus_LCDM_residual
```

Use these as exploratory diagnostics only. Any apparent residual oscillation should be checked against shuffled-redshift controls and multiple model baselines.
