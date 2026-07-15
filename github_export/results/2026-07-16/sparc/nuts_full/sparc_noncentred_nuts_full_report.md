# SPARC Full Filtered All-Q2 Non-centred NUTS Report

- **Completed:** `2026-07-15T21:30:55.132828+00:00`
- **Galaxies:** `157`
- **Rotation points:** `3001`
- **Chains:** `4`

## Convergence Decision

At least one full posterior fails a preregistered convergence gate and must not be interpreted as stable.

Passing this report clears only the full-sample sampler-convergence gate. It does not establish out-of-sample superiority over RAR.

## Gate Results

| Model | Status | Max R-hat | Min bulk ESS | Min tail ESS | Divergences | Min E-BFMI | Max depth hits |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| RAR | FAIL | 1.00259 | 1133.8 | 2182.4 | 3 | 0.8540 | 0 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | FAIL | 1.00283 | 1513.0 | 2304.5 | 6 | 0.8525 | 0 |

## Likelihood Equivalence

| Model | Maximum velocity difference (km/s) | Chi-squared difference | Pass |
| --- | ---: | ---: | --- |
| RAR | 8.527e-14 | 0.000e+00 | True |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | 1.421e-13 | 0.000e+00 | True |

## Full Posterior Global Parameters

These are well-mixed in-sample posterior estimates from chains that pass the R-hat, ESS, E-BFMI and tree-depth criteria but fail the locked zero-divergence criterion. They are not predictive model-comparison scores or a strict convergence pass.

| Model | Parameter | Mean | SD | 3% HDI | 97% HDI | R-hat | Bulk ESS | Tail ESS |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| RAR | ydisk | 0.556386 | 0.0142927 | 0.528881 | 0.582385 | 1.00142 | 1251.1 | 2182.4 |
| RAR | ybul | 0.602788 | 0.0172689 | 0.569717 | 0.634419 | 1.00118 | 1133.8 | 2511.0 |
| RAR | log10_gdagger | -10.0501 | 0.00764723 | -10.0641 | -10.035 | 1.00060 | 2940.4 | 4141.4 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | ydisk | 0.579576 | 0.0156236 | 0.55171 | 0.610073 | 1.00277 | 1643.8 | 3231.6 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | ybul | 0.598726 | 0.0168231 | 0.566716 | 0.629039 | 1.00249 | 1513.0 | 3139.8 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | log10_kappa | -0.101506 | 0.0140437 | -0.127683 | -0.0755627 | 1.00047 | 3270.6 | 4627.8 |
| PLAMB_OPTICAL_DEPTH_KAPPA_P | bridge_exponent | 0.51892 | 0.00494365 | 0.509945 | 0.528357 | 1.00050 | 4104.3 | 4999.1 |

## Interpretation Boundary

- A strict pass replaces the old random-walk convergence uncertainty for the same filtered all-Q2 likelihood.
- It does not validate the previous prior-redrawn stress scores.
- It does not replace galaxy-held-out or grouped cross-validation.
- The standing claim remains `subset wins, not a full-sample win` until predictive validation is passed.

## Reproducibility

- PyMC: `5.28.5`
- nutpie: `0.16.11`
- ArviZ: `0.23.4`
- Python: `3.13.9 | packaged by Anaconda, Inc. | (main, Oct 21 2025, 19:09:58) [MSC v.1929 64 bit (AMD64)]`
- sampling seed: `16072611`
