# Preregistration: leakage-free Quaia x DR16Q spectroscopic-redshift angular likelihood

**Frozen:** 18 July 2026, before calculating or inspecting cross-match residual outcomes  
**Status:** outcome-blind and immutable  
**Primary catalogue:** local Quaia $G<20$ catalogue positionally matched to official SDSS DR16Q v4

## Scientific question

Does the previously locked SU(2) redshift-dependent angular form predict an independently defined Quaia photometric-redshift error after joint survey-template control and genuinely held-out spatial validation?

## Frozen sample

1. Match the local `quaia_G20.0.fits` catalogue to official `DR16Q_v4.fits` on ICRS coordinates.
2. Use a maximum separation of **1.0 arcsec**.
3. Retain only reciprocal nearest-neighbour matches; this enforces a one-to-one unique match. No ambiguous duplicate is resolved using redshift.
4. Require DR16Q `IS_QSO_FINAL == 1`, finite coordinates and finite best spectroscopic `Z`.
5. Require $1.0 \le z_{\mathrm{spec}} < 1.5$ and $|b|\ge15$ degrees.
6. The 15-degree latitude cut is primary because the spectroscopic cross-match is expected to be much smaller than the full Quaia sample and the DR16Q footprint already strongly avoids the Galactic plane. It will not be changed after outcomes are seen.
7. Apply one joint finite-value mask to response, all retained nuisance templates and signal predictors. A candidate nuisance column with zero variance is omitted deterministically and recorded.

## Response and leakage prohibition

The response is

$$
y_i = z_{\mathrm{Quaia},i}-z_{\mathrm{spec},i}.
$$

The locked scalar curve is

$$
s_i = \operatorname{standardise}\!\left[\mu_{\mathrm{SU2R}}(z_{\mathrm{spec},i})-
                                            \mu_{\Lambda\mathrm{CDM}}(z_{\mathrm{spec},i})\right].
$$

It may use DR16Q spectroscopic redshift only. Quaia photometric redshift may enter only through $y_i$ and must never enter the design matrix, folds, matching, masks or injected signal.

Before any decision is accepted, Quaia photometric redshifts will be permuted with a fixed independent random permutation and the complete predictor matrix will be rebuilt. Exact equality, zero maximum absolute difference and equal SHA-256 hashes are required. Failure invalidates the replacement analysis.

## Frozen models

Both models contain an intercept and the same nuisance matrix $T_i$:

$$
\begin{aligned}
M_T:\quad y_i &= \alpha + T_i\beta + \epsilon_i,\\
M_S:\quad y_i &= \alpha + T_i\beta +
                   \boldsymbol d\mathbin{\cdot}(s_i\boldsymbol n_i)+\epsilon_i,
\end{aligned}
$$

where $\boldsymbol n_i$ is the Galactic Cartesian unit vector. The nuisance model includes the unmodulated components of $\boldsymbol n_i$, a cubic trend in $z_{\mathrm{spec}}$, and the scalar main effect $s_i$, so the candidate tests only the redshift-modulated angular interaction.

The joint nuisance groups are frozen as follows:

- **dust:** SFD $E(B-V)$, `log1p` dust and squared dust;
- **depth:** Quaia selection function, random-catalogue density, WISE W1/W2 coverage, and Gaia source density where locally available;
- **Gaia scanning geometry:** DR3 scan count and the cosine, sine and resultant summaries of scan angle;
- **WISE colour:** W1-W2, Gaia BP-RP and their interaction;
- **catalogue quality:** Quaia photometric-redshift uncertainty, Gaia $G$, W1, W2 and proper-motion uncertainties, plus available DR16Q confidence, warning, spectroscopic signal-to-noise, photometric-error and WISE-fit quality fields.

All continuous columns are centred by their sample median and scaled by their sample standard deviation after the frozen joint mask. Collinear columns are retained through a Moore-Penrose solution; a singular-value audit is reported.

## Spatial validation and calibration

- Fixed sky cells: 12-degree longitude by equal-width $\sin b$ cells.
- Fixed held-out folds: 12 Galactic-longitude sectors, assigned from cell-centre longitude.
- The nuisance-only and candidate models are fitted on eleven sectors and scored on the untouched twelfth sector in turn.
- Primary predictive statistic: pooled $\Delta\mathrm{CVLPD}=\mathrm{CVLPD}(M_S)-\mathrm{CVLPD}(M_T)$, with Gaussian variance estimated only from each training set.
- Fixed random seed: `18071826`.
- Null generator: cell-level Rademacher wild bootstrap of nuisance-only residuals, preserving the frozen positions, predictors, mask and within-cell residual structure.
- Locked null split: the first 1,000 generated null replicates calibrate the threshold; the next 1,000 independently validate false-positive rate.
- Frozen threshold: the higher empirical 99th percentile of calibration-null $\Delta\mathrm{CVLPD}$ values.
- Injection scales: 0.5, 1.0, 1.5 and 2.0 times the locked Quaia reference vector.
- Replicates: 1,000 per injection scale, using fresh cell signs after the 2,000 null replicates.
- Locked reference: amplitude `0.002548691334890594`, Galactic longitude `138.99062177880762` degrees and latitude `-36.60410263509791` degrees.

A replicate is a joint recovery only when its score exceeds the frozen threshold, its recovered cross-fit direction is within 30 degrees of the locked direction, its recovered amplitude is within `[0.5, 1.5]` times the injected amplitude, and its orientation has the locked sign.

## Conjunctive promotion rule

Promotion requires **every** gate below:

1. the predictor-permutation leakage audit passes exactly;
2. validation-null false-positive rate is at most 1%;
3. joint recovery power at the locked 1x injection is at least 80%;
4. observed pooled held-out gain is positive and exceeds the frozen empirical threshold;
5. observed cross-fit direction is within 30 degrees of the locked direction, has the locked orientation, and its amplitude ratio to the locked reference lies in `[0.5, 1.5]`.

The threshold and gates will not be changed if power is poor. In-sample BIC, ordinary least-squares significance and individual-fold values are diagnostics only and cannot promote the candidate.

## Interpretation boundary

A pass would establish a conditional cross-catalogue predictive association, not a discovery of SU(2) cosmology. Failure of calibration or power makes the test underpowered. A calibrated negative observed result rejects only this locked candidate in this sample and model. Survey-template inadequacy, cross-match selection and DR16Q footprint remain limitations.
