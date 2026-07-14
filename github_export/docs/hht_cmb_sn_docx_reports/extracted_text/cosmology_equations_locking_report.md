# cosmology_equations_locking_report

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\cosmology_equations_locking_report.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2013-12-23T23:15:00Z`
- DOCX modified: `2013-12-23T23:15:00Z`

## Extracted Text

Cosmology Equations for Locking Report

1. Distances & Distance Modulus

χ(z) = ∫₀ᶻ [ c(z′) / H(z′) ] dz′

Sₖ(χ) = { sin(√|Ωₖ| χ)/√|Ωₖ|   if Ωₖ < 0;  χ if Ωₖ = 0; sinh(√|Ωₖ| χ)/√|Ωₖ|  if Ωₖ > 0 }

D_M(z) = Sₖ(χ)

D_H(z) = c(z) / H(z)

D_A(z) = D_M(z) / (1 + z)

D_L(z) = (1 + z) · D_M(z)

μ(z) = 5 log₁₀[ D_L(z) / 10 pc ]

Δμ(z) = μ_obs(z) − μ_model(z)

2. BAO Core Definitions

r_d = ∫_{z_d}^∞ [ c_s(z) / H(z) ] dz

c_s(z) = c(z) / √[ 3 ( 1 + R_b(z) ) ]

R_b(z) = 3 ρ_b(z) / (4 ρ_γ(z))

BAO observables: D_M/r_d, D_H/r_d

D_V(z) = [ z · D_H(z) · D_M²(z) ]^{1/3}

F_AP(z) = D_M(z) · H(z) / c(z)

3. Expansion History

E(z) = H(z) / H₀ = √[ Ω_r (1+z)^4 + Ω_m(z)(1+z)^3 + Ωₖ(1+z)^2 + Ω_DE(z) + … ]

Ω_m(z) = Ω_m0 · ℳ(z)

c(z) = c₀ · ℭ(z)

4. Growth of Structure

f(z) = d ln D₊(z) / d ln a

fσ₈(z) = f(z) · σ₈(z)

f(z) ≈ [Ω_m(z)]^γ  (γ ≈ 0.55 in GR+ΛCDM)

5. CMB Distance Priors

R = √[Ω_m H₀²] · D_M(z_*) / c(z_*)

ℓ_a = π · D_M(z_*) / r_s(z_*)

6. Locking & Statistics

r(ω̄, 1/k_eff) = Cov(ω̄, 1/k_eff) / (σ_ω̄ · σ_{1/k})

t = r √[(N−2)/(1−r²)]

χ² = (d − m)^T C⁻¹ (d − m)

AIC = 2k + χ²_min

BIC = k ln N + χ²_min

7. BAO/k-Space Link

k_BAO(z) ≈ π / r_d

λ_BAO(z) ≈ 2π / k_BAO ≈ 2 r_d

ω̄(log(1+z)) ↔ 1/k_eff
