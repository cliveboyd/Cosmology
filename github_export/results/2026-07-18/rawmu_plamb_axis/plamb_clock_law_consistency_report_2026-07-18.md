# PLAMB Clock-law Consistency Audit

Date: 18 July 2026

## Decision

**PASS CLOCK-LAW RECONCILIATION.**

The raw-MU and historical no-expansion branches now use one shared clock-law
kernel while retaining distinct, explicit physical assumptions.

## Registered branches

| branch | gamma_c | p | consequence at gamma_c=1 |
| --- | ---: | ---: | --- |
| fixed raw-MU PLAMB | 1 | 0 | `D_path=(c0/H0) z(1+z/2)` |
| historical no-expansion pilot | parameter | 1 | `D_path=(c0/H0) z` when unclumped |
| constant-c linear-redshift control | 0 | 0 | `D_path=(c0/H0) z` |

The branches are no longer mathematically conflated. The raw-MU model uses a
linear redshift rate, whereas the earlier pilot uses a fractional redshift-rate
factor.

## Numerical gates

| gate | maximum absolute difference | status |
| --- | ---: | --- |
| Peter p=0 identity | `1.776e-15` | PASS |
| fractional-clock p=1 identity | `6.661e-16` | PASS |
| constant-c control | `6.661e-16` | PASS |
| analytic expression versus quadrature | `1.776e-15` | PASS |
| raw-MU p=0 wiring [mag] | `0.000e+00` | PASS |
| historical pilot p=1 wiring [Mpc] | `1.819e-12` | PASS |

## Likelihood regression

Re-running the three-release fixed-law comparison gives
`Delta BIC(PLAMB-Lambda-CDM)=94.344983139340` for `N=3422`.
The difference from the pre-refactor value is
`9.095e-13`. The refactor therefore changes no
scientific result.

## Claim boundary

This audit proves analytic and implementation consistency. It does not decide
whether `p=0` or `p=1` is physically correct. That question belongs to the
separately preregistered nested clock-law likelihood with external time-dilation
and flux constraints.

## Reproduction

```powershell
python github_export/code/shared/test_plamb_clock_distance.py
python github_export/code/rawmu/run_rawmu_release_grounded_holdouts_2026_07_18.py --self-test-only
python github_export/code/rawmu/run_plamb_axis_comparison_2026_07_18.py
python github_export/code/rawmu/audit_plamb_clock_law_consistency_2026_07_18.py
```
