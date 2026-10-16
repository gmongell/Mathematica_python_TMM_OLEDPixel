# OLED dipole-outcoupling solver: first implementation specification

Prepared 2026-10-16. **Proposed scientific scope, not implemented or validated.**
The plane-wave preview remains separate and must not be relabeled as this solver.
Seller verification is a commercial gate, not a technical prerequisite to
research implementation. A future solver must pass its own scientific gates.

**First reference milestone implemented:** `dipole_reference.py` contains
homogeneous angular/power kernels and single-interface total source-work
integration, with [documented equations and benchmark scope](DIPOLE_FORMULATION.md).
This partially addresses O1/O2/O4/O6/O8/O9 only for the stated reference geometry;
it does not satisfy the complete finite-multilayer acceptance table below.

## Intended first problem

An electric dipole in a finite, homogeneous, lossless, isotropic emitting layer
of a planar multilayer stack, with passive isotropic surrounding layers.
Use a frequency-domain dyadic Green-function/angular-spectrum formulation and
outgoing-wave boundary conditions. Resolve dipoles parallel and perpendicular
to the interfaces and s/p polarization explicitly. A merely transmitted incident
plane wave is not an emission source and cannot establish outcoupling efficiency.

Use an explicitly documented time convention, complex-index convention and
square-root branch choice. For example, with `exp(-i omega t)`, passive media
have nonnegative imaginary permittivity and outgoing/decaying longitudinal waves.
Verify that the implementation and any dependency use compatible conventions.

### Bounded first release

- Flat laterally infinite layers; no scattering, grating or finite pixel.
- Lossless, positive-index, nonmagnetic emitter host and exterior half-spaces.
- Passive absorbing non-emitting films are allowed with documented
  complex permittivity conventions.
- Point dipoles strictly inside the emitter layer, not on an interface.
- No absorbing emitter-host singularity regularization in this first model.
- No exciton transport, electrical injection, thermal effects, stimulated
  emission, gain, anisotropy, microcavity lasing or coherence between emitters.
- No claim of full OLED external quantum efficiency from optics alone.

## Inputs and reproducibility contract

| Input | Required interpretation |
|---|---|
| Wavelengths | Vacuum wavelength in nm, finite positive ordered samples |
| Stack | Explicit top-to-bottom layer order; first/last are half-spaces |
| Thicknesses | Finite interior thicknesses in nm; positive emitter thickness |
| Material data | Complex n or epsilon convention, units, source, license, interpolation domain |
| Emitter location | Layer index and depth in nm measured from a documented interface |
| Orientation | Vertical fraction in [0,1]; horizontal fraction split equally over in-plane axes |
| Isotropic ensemble | Vertical fraction 1/3, horizontal fraction 2/3 |
| Spectrum | Optional nonnegative intrinsic spectral weights with explicit normalization |
| Intrinsic quantum yield | Optional value in [0,1], with stated reference-medium assumptions |
| Numerical controls | Relative/absolute tolerance, integration strategy and tail-convergence controls |

Reject nonfinite inputs, invalid orientation/yield, interfaces as dipole
locations, unlicensed/missing material metadata, unsupported host absorption,
and wavelengths outside the measured domain unless explicit, documented
extrapolation is implemented. No silent substitution of placeholder data.
Test numerical interpolation separately from the underlying measurement's validity.

## Outputs and definitions

1. Orientation-resolved spectral electromagnetic power relative to a clearly
   specified homogeneous-host radiated power for the same dipole.
2. Top and bottom far-field hemispherical powers and angular distributions
   with declared angular variables, solid-angle measures and polarization.
3. Total optical absorption per absorbing layer from a physical flux/volume-loss
   calculation, not an unexplained residual labeled as absorption.
4. Bound/guided contributions where resolvable, with pole/integration treatment.
5. Numerical balance error, estimated quadrature/tail error and convergence flag.
6. Normalized air-escape fraction and, separately, substrate escape.
   In a glass half-space, substrate power is not automatically air outcoupling.

Let `P_EM` be total source work transferred into electromagnetic channels
(far field, guided and lossy channels) under the specified source convention.
For an accepted passive case, define the optical top-escape fraction as
`P_top_air / P_EM`, not `T_s + T_p`, and not normalized independently per angle.
Do not claim a complete channel decomposition if guided contributions or
high-wavevector lossy tails are unresolved.

Report raw powers and denominators so users can audit normalization. If a
denominator is zero, integration fails, or errors exceed limits, return an
explicit error/non-converged result; do not clip to [0,1] and report success.

Only after deriving the rate model may intrinsic yield be combined with the
electromagnetic rate enhancement. For fixed intrinsic nonradiative decay and
`F = P_EM / P_host`, a candidate rate model gives an escape probability per
excitation `q0 * (P_top_air / P_host) / ((1-q0) + q0*F)`.
This requires reference-medium assumptions and special-case tests. It is not
EQE without separately supplied and justified electrical/exciton efficiencies.

## Numerical architecture

1. Keep existing preview entry point and outputs backward compatible.
2. Introduce an explicit stack/material/source model, separate from plotting.
3. Implement stable multilayer reflection/response using scattering matrices or
   a comparably stable recursion. Do not use exponentially growing transfer
   products across thick absorbing/evanescent layers.
4. Derive dipole source terms and Green-function self-interaction for each
   orientation. Record equations, references, conventions and units before coding.
5. Integrate in-plane wavevector, separating propagating, evanescent and
   guided/pole regions. Split branch points; resolve resonances and near-interface
   tails with explicit convergence checks.
6. Use outgoing field flux for far-field channels and physically derived
   absorption. Do not count reflected powers as separate emitted populations.
7. Add orientation and emission-depth/spectral ensemble averaging only after
   single-dipole normalization is demonstrated.
8. Expose CLI/CSV/JSON only with explicit model version, material provenance,
   reference power, tolerances, quadrature diagnostics and limitations.

## Required validation before commercial claims

The following are proposed acceptance criteria, not current results:

| ID | Physical or numerical requirement | Proposed criterion |
|---|---|---|
| O1 | Homogeneous lossless medium, all orientations | Total power / reference = 1 within 1e-6; equal hemispherical total power |
| O2 | Dipole angular patterns | Vertical sin-squared and in-plane azimuthal patterns; integrals within 1e-5 |
| O3 | Symmetric lossless stack, centered emitter | Equal top/bottom power within 1e-6 relative; orientation treated independently |
| O4 | Single-interface dipole reference | Independent analytic/literature comparison across propagating/evanescent regions, agreed tolerance <=1e-4 where reference is well-conditioned |
| O5 | Passive absorbing multilayer | Nonnegative physical powers within error budget; full energy balance residual <=1e-4 of reference power |
| O6 | Integration convergence | Tightened tolerances and enlarged tail domain change integrated outputs <=1e-4 relative with stated absolute floor |
| O7 | Stable thick-film limit | No overflow/NaN; limiting behavior matches independently evaluated boundary cases |
| O8 | Orientation ensemble | Isotropic result equals weighted 1/3 vertical + 2/3 horizontal powers before taking ratios |
| O9 | Near-interface regime | Flag lack of convergence; never report finite-cutoff artifacts as validated efficiency |
| O10 | Material input | Reject out-of-domain wavelengths and unsupported host absorption explicitly |
| O11 | Emission/rate normalization | q0=0 gives zero external photon probability; q0=1 and F=1 reduce to the correct optical limit |
| O12 | Reproducible output | Coordinates, units, layer order, source position, references and error diagnostics round-trip unambiguously |

Any tolerance adjustment must be justified by numerical/reference error, not
changed merely to obtain passing tests. References and data sets must be
versioned and reviewed for redistribution permissions.

## Research references and unresolved evidence

- [S. J. Byrnes, Multilayer optical calculations, arXiv:1603.02720](https://arxiv.org/abs/1603.02720)
  supports passive multilayer conventions and plane-wave flux accounting;
  it does not by itself supply a complete emitting-dipole model.
- Use the Chance-Prock-Silbey / planar Green-function literature to derive the
  source-work and decay-rate expressions; verify the exact source and equations
  before implementation. This document does not reproduce a validated derivation.
- Select an independently implemented dipole reference or reproducible published
  case with stated stack, constants, orientation, normalization and uncertainties.
  The existing `tmm` wrapper cannot serve as an independent dipole reference.
- Obtain measured, licensed optical constants for the intended OLED materials
  and measured-stack evidence before any device-specific accuracy claim.

## Implementation order and sale gates

**Next technical milestone:** equation/convention derivation plus homogeneous
and single-interface reference kernels/tests (O1, O2, O4), not a full paid solver.
Then stable multilayer source physics, channel balance/convergence, then UI/export.
Keep each stage reviewable and fail closed on unsupported cases.

Seller onboarding can proceed independently. No publication of an "OLED solver"
offer, pre-order with delivery promises, or reuse of preview test counts as solver
validation until the relevant solver milestones and commercial gates pass.
