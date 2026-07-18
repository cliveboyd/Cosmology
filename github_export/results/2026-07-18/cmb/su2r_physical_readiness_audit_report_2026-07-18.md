# SU2R Physical Perturbation Readiness Audit

Date: 18 July 2026

## Decision

**BLOCK PHYSICAL SU2R CMB/LENSING INTERPRETATION.**

The adapter refusal gate passed: `True`. The existing growth closure and CAMB PPF proxy remain useful diagnostics, but neither supplies the field-level equations required for a physical SU2R spectrum.

## Existing evidence

- The 66,150-row growth scan found only a `Delta chi2=-0.0971` raw improvement over its smooth closure; parameter penalties retained the smooth branch.
- The historical ACT PPF proxy strongly penalised its effective-w row, but it was an unmatched fluid proxy and is not a physical SU2R exclusion.

## Missing lensing requirements

- `action_or_field_equations`
- `background_field_equations`
- `scalar_dynamical_equations`
- `scalar_constraint_equations`
- `perturbed_stress_energy`
- `gauge_and_variable_dictionary`
- `initial_conditions`
- `weyl_potential_and_lensing_relation`
- `implementation.module_path`
- `implementation.sha256`

## Required next input from the model

Provide the action or complete field equations; define the background SU2R field; derive scalar dynamical and constraint equations in a stated gauge; give delta-rho, delta-p, velocity and anisotropic-stress terms; specify adiabatic/isocurvature initial conditions; and derive the Weyl-potential relation used for lensing. The implementation module must then be hash-locked in the registry before ACT or Planck evaluation.
