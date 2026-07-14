# IMF2_HHT_log1pz_interpretation_report

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\IMF2_HHT_log1pz_interpretation_report.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2025-10-09T10:09:00Z`

## Extracted Text

Interpreting IMF2 on a log(1+z) Axis (HHT)

Concepts, mappings to redshift space, and energy projections

1) Why use log(1+z)?

Redshift z relates to the scale factor a by a = 1/(1+z). Taking x = ln(1+z) turns multiplicative changes in (1+z) (or equivalently 1/a) into additive displacements. Instantaneous frequency measured in x has a natural interpretation as cycles per unit ln(1+z), i.e., a multiplicative scale in a.

x = ln(1+z)   ↔   a = e^{-x}

On the x-axis, equal intervals Δx correspond to equal factors in (1+z). Thus, an oscillation with frequency f(x) (cycles per unit x) captures features occurring at nearly constant multiplicative scale in cosmic expansion.

2) IMF / HHT essentials on x = ln(1+z)

EMD separates a signal s(x) into intrinsic mode functions (IMFs). For an IMF_k(x), the Hilbert transform produces the analytic signal A_k(x) e^{iφ_k(x)} with instantaneous amplitude A_k and phase φ_k. The instantaneous frequency is ω_k(x)=dφ_k/dx and f_k=ω_k/2π.

IMF_k(x)  →  z_k(x) = A_k(x) e^{i φ_k(x)}

ω_k(x) = dφ_k/dx      f_k(x) = ω_k(x) / (2π)    [cycles per ln(1+z)]

Note:Bandwidth is the range of instantaneous frequencies an IMF occupies. One cycle corresponds to Δx = 1/f,  mapping to a multiplicative change in (1+z) by e^{Δx}.

3) Mapping frequency and bandwidth to redshift space

If an IMF has local frequency f at redshift z (x = ln(1+z)), one oscillation spans Δx ≈ 1/f, which maps to:

Δx = 1/f    ⇒    Δz ≈ (1+z)(e^{1/f} − 1)

Higher f means finer structure in redshift (smaller Δz); lower  f means broader modulations  (larger Δz).

4) Energy projection and Hilbert spectrum on x

For IMF_k, the Hilbert spectrum H_k(x,f) assigns energy |A_k(x)|^2 at instantaneous frequency f_k(x). Summing over k gives H(x,f). Integrate over f to get energy vs x, or over x for the marginal spectrum vs f.

H_k(x, f) ≈ |A_k(x)|^2 · δ(f − f_k(x))

H(x, f)   = Σ_k H_k(x, f)E(x)      = ∫ H(x, f) df M(f)      = ∫ H(x, f) dx

5) IMF2 interpretation in log(1+z) space

IMF2 often captures the first non-trivial oscillatory mode after the trend; lower frequency than IMF1 but still localized.

IMF2’s f_2(x) maps to Δz ≈ (1+z)(e^{1/f_2}−1). Coherent f_2 with high amplitude A_2(x) indicates a repeatable modulation at a characteristic multiplicative scale in (1+z).

In your HHT screens, an IMF2 novelty around z≈0.38–0.70 would be consistent with BAO residual leverage in that range (seen via LOOCV).

6) Mapping back to cosmic evolution

Since x = ln(1+z) = −ln a, uniform frequency in x corresponds to uniform frequency in ln a. An IMF with nearly constant f describes processes repeating per e-fold in the scale factor. For time interpretation, one may map z→t via a background H(z).

t(z) = ∫_z^∞ dz' / [(1+z') H(z')] H(z) = H0 √(Ω_m(1+z)^3 + Ω_Λ)

7) What can be learned from IMF2

Localization of structured residuals in specific redshift bands.

Scale diagnosis via f-range → Δz-range (fine vs broad structure).

Cross‑dataset checks: compare IMF2 energy in BAO vs SN residuals.

Model stress tests where IMF2 energy and χ² leverage coincide (e.g., z ≈ 0.38–0.70).

Appendix: Key equations and definitions

x ≡ ln(1+z) a = e^{-x}

Analytic signal:  z_k(x) = IMF_k(x) + i · H[IMF_k(x)] = A_k(x) e^{i φ_k(x)}

Instantaneous frequency:  f_k(x) = (1/2π) dφ_k/dx   [cycles per ln(1+z)]

Redshift span per cycle:  Δz ≈ (1+z)(e^{1/f_k}−1)

Hilbert spectrum:  H(x,f)                = Σ_k |A_k(x)|^2 δ(f−f_k(x))Energy   vs x:  E(x) = ∫ H(x,f) df Marginal vs f:  M(f) = ∫ H(x,f) dx
