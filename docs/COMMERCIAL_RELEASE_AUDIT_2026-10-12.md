# OLED TMM commercial release readiness audit (2026-10-12)

Status: **BLOCKED — research prototype; not a validated commercial solver.**

## Source reviewed
- `tmm_OLED9.py` on `main`, blob `320829789f7d4dab662b9dfad2b202eddcd7889e`.
- Existing README identifies earlier numerical and unit issues. This audit independently confirms the issues below from source inspection; it does not claim runtime testing.

## Reproducible static findings
1. `coh_tmm` is called but never imported or defined; running the sweep reaches a `NameError` once the first permitted angle is evaluated.
2. `n_list` contains air, one refractive index per wavelength (1000 values), and silver, while `d_list` contains four layer thicknesses. The TMM API expects one refractive index per physical layer at the *current wavelength*.
3. The call `coh_tmm('s', n_list, d_list, wavelength * 1e-9, theta)` reverses angle and vacuum-wavelength arguments relative to the standard `tmm.coh_tmm(pol, n_list, d_list, th_0, lam_vac)` signature.
4. The air-to-organic critical-angle condition is applied to all angles. A critical-angle cutoff is not appropriate for incidence from lower-index air into higher-index organic material.
5. The wavelength dispersion parameters have no stated dimensional calibration or provenance. The output is described as OLED outcoupling efficiency although summed TE/TM transmission is not a calibrated dipole-emission/outcoupling calculation.
6. Script performs a large plotting sweep at module import time; it lacks a guarded entry point and configurable resolution, making automated unit testing difficult.

## Acceptance tests required before sale
- Install/pin compatible `numpy`, `matplotlib`, `tmm`; document Python version and dependencies.
- Add deterministic Fresnel normal-incidence and zero-thickness checks, and lossless R+T=1 checks.
- Validate absorbing-film power accounting and TE/TM polarization normalization.
- Compare a representative multilayer stack to an independent solver with stated tolerances.
- Distinguish plane-wave transmittance from physically defined dipole-emission outcoupling.
- Add headless smoke test, parameter input validation, benchmark fixture, and CI.
- Review Apache-2.0/CC BY notices and all third-party material before packaging.
- Do not publish numerical efficiency claims or paid solver releases until tests pass.

This file is an audit artifact only. It does not change executable code or certify performance.
