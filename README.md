# Multilayer OLED transfer-matrix and Green-function explorations

Python scripts and Mathematica notebook revisions exploring multilayer optics, polarization, wavelength dependence, and OLED outcoupling concepts.

## Functions and engineering applications

| Source / interface | Observed function | Inferred application |
|---|---|---|
| `python_tmm_oled_Green'sfcnApproach.py`: `OLEDOutcouplingTMM` | Layer definitions, aluminum Drude model, interface/propagation matrices, angle/wavelength calculations | Exploring material and thickness sensitivity in OLED stacks [1] |
| `oled_outcoupling_simulation.py`: `transfer_matrix_method(...)` | Normal-incidence characteristic-matrix calculation | Teaching and checking thin-film interference concepts [1] |
| `tmm_OLED8.py`, `tmm_OLED9.py` | Wavelength/angle/thickness sweeps and image output | Exploratory design maps, after correcting implementation issues |
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


## Optical-constant databases and integration roadmap (2026-10)

The resources below are **major, useful optical-constant collections**, not a
verified ranking by record count: catalog sizes, material coverage, licensing,
and access terms are not directly comparable. Always retain the original
dataset, publication citation, measurement conditions, and source version.

| Resource | Official link | Strength / relevance to OLED stacks | Access / integration note |
|---|---|---|---|
| **RefractiveIndex.INFO** | [Website](https://refractiveindex.info/) · [GitHub database](https://github.com/polyanskiy/refractiveindex.info-database) · [releases](https://github.com/polyanskiy/refractiveindex.info-database/releases) | Broad, literature-sourced optical constants including metals, dielectrics, semiconductors and some organics; **first choice** for automated research integration | Database is **CC0**; YAML datasets may contain tabulated n,k or dispersion formulas. Pin a release/commit, record the selected *material/book/page* and primary literature reference. |
| **KLA / Filmetrics Refractive Index Database** | [Database](https://www.kla.com/products/instruments/refractive-index-database/Schott%2BN-SF10) | Large practical thin-film material catalog, including ITO, Ag, Al, SiO2, Si3N4 | Some datasets are proprietary or require a request; check reuse terms before downloading or redistributing. |
| **SCHOTT Optical Glass Catalog** | [Official downloads](https://www.schott.com/en-gb/products/optical-glass-p1000267/downloads) | Manufacturer dispersion and transmission for substrate and encapsulation glass | XLSX, optical-design formats and datasheets; formulas may need sampling onto the model's wavelength grid. Not an OLED organic-emitter n,k source. |
| **OHARA Optical Glass Catalog** | [Official downloads](https://www.ohara-inc.co.jp/en/product/catalog/) | Extensive glass dispersion data for substrates and optical components | CSV/XLSX and optical-design catalogs; check material-specific spectral validity and license. |
| **NASA Optical Constants Database (OCdb)** | [Official database](https://ocdb.smce.nasa.gov/) | Published n,k datasets, particularly ices, tholins and planetary/astro materials | Useful specialized reference, **not** a primary source for OLED device films. |
| **NIST optical-constants research** | [NIST program](https://www.nist.gov/programs-projects/calculation-optical-constants-and-dielectric-response) | Computational dielectric-response methods and cross-checking | Research portal rather than a general downloadable OLED n,k catalog. |

### Recommended data pipeline

1. Prefer **RefractiveIndex.INFO** for initial reproducible n,k ingestion, then
   compare at least one independent experimental dataset for each critical
   material. Use original peer-reviewed measurements for Alq3, emissive
   organics, ITO, Ag, and other stack layers whenever available.
2. Record **material name, source URL, publication DOI, dataset ID/page, source
   commit or release, film preparation, substrate, temperature, spectral
   interval, and license**. Different deposition conditions can produce
   substantially different n,k, especially for ITO and thin metals.
3. Convert selected **tabulated** data to a CSV with header
   `wavelength_nm,n,k`, ordered by increasing wavelength. The repository's
   `OpticalConstants.from_csv(path)` accepts this format; its
   `evaluate(wavelength_nm)` returns `n + 1j*k`, with no extrapolation.
   Formula-only data requires evaluation with its documented dispersion
   equation before CSV export. No direct YAML importer is implemented yet.
4. Preserve the passive-medium sign convention `n+ik` for
   `exp(-i omega t)`. Check wavelength units (nm vs micrometers), coverage,
   interpolation artifacts, and energy balance `R + T + A = 1` for
   non-scattering stacks. Avoid using bulk metal n,k as an unqualified model
   for ultrathin or discontinuous films.
5. Run wavelength/angle TMM benchmarks and sensitivity analysis across
   **alternative measured datasets** before any optimization. Plane-wave
   transmittance is **not** OLED dipole outcoupling or external quantum
   efficiency.

**Current implementation status:** `optical_constants.py` provides
validated CSV ingestion and bounded interpolation; automatic database
download, YAML/dispersion-formula parsing, per-layer assignment in the
OLED solver, and measurement-based model validation remain future work.
Do not treat the illustrative Alq3/silver formulas in `tmm_OLED9.py`
as experimental optical constants.

Useful scholarly attribution: M. N. Polyanskiy, *Scientific Data* **11**,
94 (2024), [doi:10.1038/s41597-023-02898-2](https://doi.org/10.1038/s41597-023-02898-2).
