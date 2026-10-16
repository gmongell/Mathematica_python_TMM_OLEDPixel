# Multilayer OLED transfer-matrix and Green-function explorations

Python scripts and Mathematica notebook revisions exploring multilayer optics, polarization, wavelength dependence, and OLED outcoupling concepts.

## Installable plane-wave preview (0.1.0a1)

**Release preparation only; not approved for paid publication.** The packaged
program is a headless plane-wave transmittance demonstration, not an OLED
dipole-emission, outcoupling-efficiency, chromaticity, or device-design solver.
The inherited organic and silver index formulas are illustrative placeholders,
not measured material data. Historical scripts and notebooks are not included
in the wheel and are not certified by these changes.

Python 3.11-3.14 is the declared compatibility range; CI checks the matrix below,
and only completed test runs establish support. From the repository, on Windows:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install ".[dev]"
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\mongelli-plane-wave.exe --wavelength-count 31 --angle-count 19 --output transmission.png
```

The command creates:

- `transmission.png`: a headless wavelength/angle contour plot.
- `transmission.csv`: columns `wavelength_nm`, `angle_deg`,
  `unpolarized_transmittance`; wavelength-major order, one row per sample.
- `transmission.json`: software version, coordinates, layer thicknesses,
  polarization normalization, and model limitations.

Lengths use nm, angles use degrees, and transmittance is dimensionless.
The stack is air / illustrative organic / 50 nm illustrative silver / air,
with semi-infinite air on both sides. The default grid is 400-700 nm and
0-89 degrees; unpolarized transmittance is `(T_s + T_p) / 2`, not their sum.
Use `--thickness-nm` to change the organic thickness. Counts must be at least
two. Output must end in `.png` and its parent directory must already exist.
Existing output files are protected unless `--overwrite` is explicitly passed.
Errors exit nonzero; a write failure may leave partial outputs.

For a built wheel, install with `python -m pip install <wheel-path>` and use
`mongelli-plane-wave` from the same environment. An internet connection is
needed to resolve dependencies unless wheels were separately supplied. This
is **not** a standalone Windows installer, and does not bundle Mathematica,
NumPy, Matplotlib, or tmm.

### Validation and distribution

```powershell
.\.venv\Scripts\python.exe -m build
.\.venv\Scripts\python.exe -m twine check dist\*
```

Tests cover analytic Fresnel interfaces, Airy thin-film transmittance, a separate
normal-incidence characteristic-matrix calculation for the absorbing stack,
lossless power balance, invalid inputs, and CLI export. These checks validate
bounded plane-wave behavior; they do not validate emitter physics or measured
material properties. See `tests/test_tmm_oled9.py` for cases and tolerances.

The GitHub workflow runs numerical tests, builds distributions, smoke-tests an
installed wheel outside the checkout, and uploads preview artifacts. Its
runtime-dependency audit emits a CycloneDX inventory and fails on known
vulnerabilities (or unavailable advisory data); it is not a full security review.
The inventory also includes environment tooling such as pip. Its
optional manual draft-release action accepts only a tag matching the package
version and creates a **draft prerelease**; it never publishes a paid listing.
No automatic tag push or public release is performed.

### GitHub + Gumroad commercial approval gates

Owner-review drafts are available in [commercial policies and storefront copy](docs/COMMERCIAL_DRAFTS.md),
[Gumroad seller verification guidance](docs/GUMROAD_SELLER_CHECKLIST.md), and
[the proposed OLED dipole-solver specification](docs/OLED_SOLVER_SPEC.md).
These documents are not active policies, completed verification, or implemented
solver functionality.

The five enrollment packets are supporting checklists, not enrollment approvals.
Use GitHub for technical documentation, demonstrations and versioned artifacts;
Gumroad is the proposed separate checkout and buyer-delivery channel.

Before charging customers:

1. Decide whether the product is this educational plane-wave tool or a future
   validated OLED solver. Do not use OLED efficiency claims for this preview.
2. Confirm code/contributor and material-data provenance and preserve the
   existing Apache-2.0 code and CC BY 4.0 documentation terms. A purchase cannot
   remove existing open-source rights; exclusive licensing is not implied.
3. Complete independent scientific review for the intended claims, dependency
   vulnerability scanning, license review, supported-platform testing, release
   checksums, a dependency SBOM, and signed-tag provenance.
4. Approve price, customer license wording, support contact and response policy,
   refund terms, updates policy, privacy notice, and truthful screenshots.
5. Complete seller identity, tax and payout verification directly in Gumroad.
   Never store tax identifiers, identity documents, or banking data here.
6. Test checkout and buyer delivery using the platform's supported test process.
   Upload paid deliverables to authorized buyer delivery, not a public GitHub
   release if they are intended to be access-restricted.

No price, exclusive license, seller approval, support SLA, paid listing, or
production readiness is asserted by this repository. Microsoft/AWS/GCP billing
integrations are intentionally outside this first downloadable-preview scope.

### Proposed storefront copy (draft, not a live offer)

**Title:** Mongelli Plane-Wave Optics Preview

**Audience:** Researchers, educators and scientific programmers exploring
planar thin-film transmission with Python.

**Description:** Generate wavelength/angle transmission maps for an illustrative
air/organic/silver/air stack. Export a PNG plot, numerical CSV and model metadata.
Includes source code and a Python command-line interface. Requires Python and
separately installed dependencies; no standalone Windows executable is supplied.

**Limitations:** Alpha research/educational preview. No measured material
database, emitting dipoles, calibrated OLED efficiency, GUI, or device-performance
guarantee. The independent checks cover specified plane-wave cases only.

**Rights:** Existing Apache-2.0 code and CC BY 4.0 documentation terms remain.
Any eventual charge must accurately describe the added value (for example,
packaging or a separately agreed support service), not imply exclusive rights.

**Not yet approved:** Price, support commitments, refund/updates terms, customer
license presentation, final screenshots, seller verification and publication.

## Functions and engineering applications

### Dispersive data, spectral EQE and ensemble planning

See [native optical data and spectral EQE derivation](docs/SPECTRAL_EQE_AND_DATA.md)
and [orientation/resource planning](docs/ENSEMBLE_RESOURCE_ESTIMATES.md).
The new APIs provide pinned native refractiveindex.info ingestion, local
measurement CSV ingestion, checksum-checked offline catalog records,
non-extrapolating n+ik evaluation, wavelength-resolved layer stacks, explicit
EQE rate composition, incoherent orientation comparisons, and CPU/GPU/RAM
workload estimates. No cloud account is required or provisioned.

**Scope boundary:** stack-dependent plane-wave R/T/A is wired; finite-layer
dipole emission/collected-flux physics is not yet implemented. EQE composition
requires externally derived physical source-work/escape-rate budgets, not
plane-wave T. The resource estimator is not an implemented GPU/FDTD solver.

### Research dipole reference API

`dipole_reference` provides a separate homogeneous and single-interface
benchmark milestone. It does **not** replace the preview CLI or implement
a finite-stack OLED solver. See [equations, conventions and reference limits](docs/DIPOLE_FORMULATION.md)
and [the larger solver specification](docs/OLED_SOLVER_SPEC.md).

```python
from dipole_reference import SingleInterface, homogeneous_pattern

pattern = homogeneous_pattern([0.0, 1.5707963267948966], 0.0, "vertical")
print(pattern.s, pattern.p)  # dP / (P_host dOmega), radians, per steradian

# Synthetic permittivity, not a measured silver data set.
interface = SingleInterface(
    wavelength_nm=550.0, distance_nm=20.0,
    substrate_epsilon=-5.0 + 0.5j, host_index=1.0,
)
result = interface.decay_rates(atol=1e-8, rtol=1e-7)
print(result.vertical, result.parallel, result.isotropic)
print(result.estimated_error, result.refinement_delta)
```

`homogeneous_power` returns macroscopic SI watts for a peak harmonic dipole
amplitude in C m; local-field corrections are not included.
`SingleInterface` requires a lossless real positive host index and positive
dipole distance. Permittivity is relative, with nonnegative imaginary part
under `exp(-i omega t)`; `"pec"` selects the ideal perfect-mirror benchmark.
Lossless nonpositive-permittivity substrates, absorbing hosts, and dimensionless
separations `2*k_host*d < 1e-6` are unsupported and rejected.

The result is total electromagnetic source work relative to homogeneous-host
power, **not** external efficiency, radiated-only power, or a channel
decomposition. Its propagating/evanescent corrections are integration regions;
in particular host-evanescent fields can transmit into a lossless high-index
substrate. Integration failures raise `DipoleConvergenceError` instead of
returning clipped or success-shaped values. Diagnostic errors are numerical
estimates, not guaranteed bounds or experimental validation.

The module uses SciPy and is included in the wheel alongside the unchanged
preview API. Tests are in `tests/test_dipole_reference.py`; no measured optical
constants or production OLED claims are supplied.

| Source / interface | Observed function | Inferred application |
|---|---|---|
| `python_tmm_oled_Green'sfcnApproach.py`: `OLEDOutcouplingTMM` | Layer definitions, aluminum Drude model, interface/propagation matrices, angle/wavelength calculations | Exploring material and thickness sensitivity in OLED stacks [1] |
| `oled_outcoupling_simulation.py`: `transfer_matrix_method(...)` | Normal-incidence characteristic-matrix calculation | Teaching and checking thin-film interference concepts [1] |
| `tmm_OLED9.py` | Installable plane-wave transmittance sweep, PNG/CSV/JSON export | Illustrative optics only; not OLED dipole outcoupling |
| `tmm_OLED8.py` | Historical wavelength/angle/thickness sweep | Unvalidated prototype, excluded from the package |
| `TMMandGreensFunctions*.nb`, `mathematica_tmm_oled.txt` | Symbolic/numerical multilayer and Green-function formulations | Comparing analytical formulations with numerical solvers [1,2] |

These are research prototypes. The class filename does not establish a complete dyadic Green-function emission calculation: its efficiency method explicitly uses a simplified transmission-amplitude approach.

## Example function calls

Python examples require NumPy, SciPy, and Matplotlib; the separate sweep scripts may also require `tmm`. The following uses the actual class API and avoids executing its full plotting sweep:

```python
import importlib.util
from pathlib import Path
path = Path("python_tmm_oled_Green'sfcnApproach.py")
spec = importlib.util.spec_from_file_location("oled_model", path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
oled = module.OLEDOutcouplingTMM()
oled.define_oled_structure()
value = oled.outcoupling_efficiency_angle_wavelength(550.0, 0.0)
print(value)  # Model output; not a validated device efficiency.
```

A smoke run of this class example completed during the review and returned approximately `7.22e-54`; that demonstrates execution only, not physical accuracy. Wavelength arguments are nm and angle arguments are degrees. The structure must be initialized before calling the angle/wavelength method. To launch its built-in plots instead:

```bash
python "python_tmm_oled_Green'sfcnApproach.py"
```

For notebook exploration, open `TMMandGreensFunctions.nb` in Mathematica and inspect/evaluate initialization and parameter cells before the dependent cells. There is no uniform package-level Wolfram API.

## Known issues and interpretation

`oled_outcoupling_simulation.py` mixes a wavelength in nm with layer thickness converted to meters inside the phase calculation. `tmm_OLED9.py` calls `coh_tmm` without importing it, constructs incompatible layer/index lists, and passes wavelength/angle arguments in an order that needs correction against the dependency API. These examples must be repaired before numerical interpretation. The class uses approximate material constants, simplifies emission to transmission, and catches calculation exceptions by returning zero for a polarization; a zero is therefore not conclusive evidence of no emitted power. None of these scripts was validated against experimental OLED data in this review. Check lossless energy balance, single-interface Fresnel limits, absorbing-media conventions, and agreement with an independent solver before design optimization.

## Review scope and software citation

Documentation reviewed on 2026-10-08 against source commit [`b4adb4a1d758`](https://github.com/gmongell/Mathematica_python_TMM_OLEDPixel/tree/b4adb4a1d758cf3a7a19148929c4e9a8bfdcdc05). “Observed” means supported by source inspection; engineering applications are reasoned possibilities unless explicitly demonstrated. Scholarly references provide methodological context and do not certify these implementations. Runtime validation is stated separately above.

For software attribution, cite Guy Francis Mongelli, *Mathematica_python_TMM_OLEDPixel*, the [repository](https://github.com/gmongell/Mathematica_python_TMM_OLEDPixel), the exact commit used, and your access date. Also cite the relevant method publications and any original third-party contributors. No unverified software DOI or release version is assigned by this documentation.

## Scholarly references

1. S. J. Byrnes (2016; revised 2020). “Multilayer optical calculations.” *arXiv:1603.02720*, preprint. [DOI: 10.48550/arXiv.1603.02720](https://doi.org/10.48550/arXiv.1603.02720). Reference for planar transfer matrices, polarization, absorption, and power normalization; it is not a validation of the repository’s emitter models.

2. A. F. Oskooi et al. (2010). “Meep: A flexible free-software package for electromagnetic simulations by the FDTD method.” *Computer Physics Communications* 181, 687–702. [DOI: 10.1016/j.cpc.2009.11.008](https://doi.org/10.1016/j.cpc.2009.11.008). Method/software reference for FDTD, dispersive media, sources, and boundary conditions.

## Ownership and existing license notices

Copyright (c) 2025 Guy Francis Mongelli

The existing project notice declares Apache License 2.0 for project code. Documentation, prose, and figures are declared CC BY 4.0; notebook code cells are Apache-2.0 and narrative/figures CC BY 4.0. Preserve all file-level and third-party notices. This README update does not change ownership or licensing terms.


## Reproducibility, source verification, and contribution policy (2026-10-08)

- **Observed versus proposed:** The function inventory above describes inspected source where indicated. Engineering applications identified as *inferred* are potential uses, not verified features or validated performance claims.
- **Usage examples:** Treat the documented commands and calls as illustrative until the referenced source file, runtime version, dependencies, required input data, and working directory have been checked. Do not execute notebook fragments or batch-scheduler directives as standalone programs without adapting their context.
- **Scientific citations:** References above identify relevant governing methods and computational background; citing a publication does not imply that its algorithm is implemented in this repository or that the publication was authored by this repository owner.
- **Missing or null artifacts:** Empty, placeholder, missing, or non-executable source files must not be represented as functional implementations. Candidate restorations from personal archives require provenance, content comparison, license review, and explicit verification before committing code.
- **Access model:** This repository is publicly readable. Public visibility does not grant anonymous push rights; write access is controlled separately through repository collaborators, credentials, apps, deploy keys, and branch rules. This README is descriptive and does not itself enforce permissions.
