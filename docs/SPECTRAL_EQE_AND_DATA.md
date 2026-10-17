# Dispersive stacks, native optical data and spectral EQE

See the [consolidated mathematical record](MATHEMATICAL_DERIVATIONS.md)
for rendered LaTeX equations and their dipole/ensemble/resource rationale.

Research milestone, 2026-10-16. Implemented foundations are in `spectral_eqe.py`.
**A planar finite emitting-layer source solver now exists in `multilayer_dipole.py`**
(one wavelength per solve; see section 12 of the derivations). It reports source work and
top/bottom exit flux; it does not separate guided from absorbed power, handle anisotropic
layers, model charge transport or optimize stacks. The data catalog, stack evaluator and
EQE composition are working components, not a complete device-performance prediction.

## What wavelength-specific EQE means

Let a stack be defined by thicknesses d_j and relative complex indices
$N_j(\lambda)=n_j(\lambda)+ik_j(\lambda)$, under $e^{-i\omega t}$.
For a specified emitter position, dipole orientation and collection geometry,
the emitting-source solver computes (or must still compute, where noted), each as a function of $(\lambda;\{d_j,N_j\},z,\text{orientation})$:

$$
\begin{aligned}
F_{\operatorname{total}}&=P_{\operatorname{EM}}/P_{\operatorname{host}}\\
F_{\operatorname{escape}}&=P_{\operatorname{collected}}/P_{\operatorname{host}}
\end{aligned}
$$

P_EM is source work into all electromagnetic channels, including loss;
P_collected is emitted flux into a specified external collection region.
Both use the SAME lossless homogeneous-host reference for the same dipole.
They are not the R/T coefficients of an externally incident plane wave.
Air escape, substrate emission and detector-NA collection are different regions.

Let q0 be intrinsic photon yield in that reference environment, eta_charge
the charge-balance factor and eta_exciton the fraction of injected pairs
forming the modeled emissive state. For an independently evaluated
monochromatic transition, with the intrinsic nonradiative rate held fixed:

$$
\operatorname{EQE}_{\operatorname{mono}}(\lambda)=\frac{\eta_{\operatorname{charge}}\,\eta_{\operatorname{exciton}}\,q_0\,F_{\operatorname{escape}}(\lambda)}{1-q_0+q_0F_{\operatorname{total}}(\lambda)}
$$

This follows from Gamma_rad,host/(Gamma_rad,host+Gamma_nr)=q0,
Gamma_EM=F_total*Gamma_rad,host and
Gamma_escape=F_escape*Gamma_rad,host. q0=0 yields no photons; q0=1
reduces to the escape/total optical-rate fraction (unless total rate vanishes,
in which case the steady-state yield is undefined).

Thickness and n+ik affect both optical rates. They do not determine
charge balance, exciton formation, intrinsic emitter loss or spatial excitation.
No universal 25% spin factor, unit charge balance or unit q0 is silently assumed.

### A broadband emitter has a different denominator

For one emitter population with intrinsic PHOTON spectral density S0(lambda)
normalized to integral S0 d(lambda)=1:

$$
\begin{aligned}
D&=1-q_0+q_0\int S_0(\lambda)F_{\operatorname{total}}(\lambda)\,d\lambda\\
\frac{d\operatorname{EQE}}{d\lambda}&=\frac{\eta_{\operatorname{charge}}\,\eta_{\operatorname{exciton}}\,q_0\,S_0(\lambda)F_{\operatorname{escape}}(\lambda)}{D}
\end{aligned}
$$

Integrate the latter to get dimensionless broadband EQE. It is incorrect
to apply an independent lifetime denominator to every spectral bin of this
same emitter. A spectrum in optical power per nm must be converted to photon
weights (multiply by lambda/(h*c)) and explicitly normalized first.
The supplied wavelength range must cover the intended emission spectrum;
truncated data normalized to one describe a conditional/truncated model.
The implementation requires a nonnegative already-normalized photon density;
it never silently renormalizes, extrapolates or mixes power/photon units.

## Native material-data catalog

`OpticalDataLake(root)` uses an explicit local or already-mounted filesystem
root. It stores immutable dataset manifests plus SHA-256-addressed original
source bytes. Manifest/source checksums are checked on load. These detect
accidental mutation, not malicious alteration of both object and checksum;
they are not cryptographic signatures or proof of material accuracy.

Implemented providers:

1. Direct, explicit HTTPS fetch of individual native
   [refractiveindex.info database](https://github.com/polyanskiy/refractiveindex.info-database)
   YAML datasets pinned to a complete Git commit.
2. Native YAML import from an explicitly identified copy of that database.
3. Local measurements or other licensed sources via exact CSV headers
   `wavelength_nm,n,k`, with required source, license, citation and conditions.

This is a local-first curated catalog, not a provisioned Azure/GCP/S3 service,
full database mirror, or integration with every refractive-index vendor.
New providers can import that CSV schema without losing source provenance.
An existing mounted data lake can be the root; cloud authorization, object-store
adapters, access governance and multi-user transactional writes are future work.
Do not point it at untrusted manifests or overwrite scientific source objects.

### Supported native formats and units

- `tabulated nk`, or separate `tabulated n` and `tabulated k`.
- Sellmeier `formula 1` and `formula 2` for n, optionally with tabulated k.
- Unsupported analytic formulas, anisotropic/tensor models, nonlinear n2
  and duplicate components raise errors rather than being approximated.
- Database wavelengths in micrometers are converted to nm. n and k tables
  are linearly interpolated separately, only inside their overlapping domains.
- Missing k is NOT evidence of zero absorption. An n-only dataset needs
  explicit `assume_lossless=True`, which is preserved in metadata.
- No automatic material selection, zero-loss fallback, extrapolation, or
  network access occurs during stack evaluation.

The upstream database declares CC0-1.0; retain database/original publication
citations and measurement conditions. That does not make its hosted articles,
figures or unrelated commercial databases freely redistributable. Import
only data you have rights to use. Sample preparation, film thickness,
temperature, deposition method, anisotropy and wavelength range affect suitability.

### Example: pinned measured silver data

```python
from pathlib import Path
from spectral_eqe import OpticalDataLake

lake = OpticalDataLake(Path.home() / "optical-material-catalog")
commit = "c5c2f188e848453def5970e347399d653df2ffc2"
silver = lake.fetch_rii(
    "Ag-Johnson-1972",
    "database/data/main/Ag/nk/Johnson.yml",
    commit,
)
print(silver.index([450.0, 550.0, 650.0]))  # n + i*k
print(silver.metadata["citation"])
```

The first fetch downloads only this public dataset. Subsequent same-source
calls use verified local files without network access. Reusing an ID for
another revision or lossless assumption is rejected; use a new ID.
Select your own root rather than committing a private material lake.

Fused silica is also supported through
`database/data/main/SiO2/nk/Malitson.yml`, using formula 1; treating its
missing k as zero requires an explicit lossless assumption. Its room-temperature
optical formula is not a substitute for an OLED organic emitter measurement.

Database citation:
M. N. Polyanskiy, Scientific Data 11, 94 (2024),
[doi:10.1038/s41597-023-02898-2](https://doi.org/10.1038/s41597-023-02898-2).
Upstream [database notice](https://refractiveindex.info/about),
checked 2026-10-16.

## Layer stack wiring

```python
from spectral_eqe import Layer, DispersiveStack

# Each ID must already be imported from appropriate measurement data.
stack = DispersiveStack(
    layers=(
        Layer("ambient", None),       # semi-infinite top half-space
        Layer("transport-top", 35.0),
        Layer("emitter-host", 25.0),
        Layer("transport-bottom", 40.0),
        Layer("metal-electrode", 80.0),
        Layer("substrate", None),     # semi-infinite bottom half-space
    ),
    emitter_layer=2,
    emitter_depth_nm=12.5,            # from the top of that layer
)
evaluated = stack.evaluate(lake, [450.0, 500.0, 550.0, 600.0])
print(evaluated["indices"])          # [wavelength, layer], complex
print(evaluated["thicknesses_nm"])

diagnostic = stack.plane_wave_spectrum(lake, [450.0, 500.0, 550.0, 600.0])
# diagnostic["R"], ["T"], ["A"] are NOT an EQE result.
```

Interior thicknesses must be nonnegative; the emitter position must be strictly
inside a positive-thickness emitting layer. For the current model contract,
the emitting host and exterior half-spaces must be lossless positive-index
media. Absorbing interior non-emitting films are allowed.
Every material's provenance and explicit lossless assumptions accompany
the evaluated stack. Sweeps construct new stacks, not mutate catalog datasets.

## EQE composition API

```python
from spectral_eqe import OpticalRateBudget

# ANALYTIC homogeneous example, not output from a finite OLED device solver.
# Total source work is 1; one homogeneous hemisphere collects 1/2.
budget = OpticalRateBudget(
    wavelengths_nm=[500.0, 550.0, 600.0],
    total=[1.0, 1.0, 1.0],
    escape=[0.5, 0.5, 0.5],
    reference="Analytic homogeneous lossless host, one hemisphere, same dipole",
)
print(budget.monochromatic_eqe(
    intrinsic_yield=0.8, charge_balance=0.9, exciton_yield=1.0,
))  # 0.36 at each wavelength for this stipulated model
```

The constructor checks rate bounds/grid/reference, not whether a caller's
source calculation is scientifically valid. Raw rates MUST be derived and
validated independently; no input label alone proves this.

`orientation_average` combines total and escape powers before calculating a
rapidly reorienting population's yield. Independently excited, fixed-orientation
dipoles have different lifetimes: compute each one's yield then average
populations. `dipole_ensemble_eqe` returns both results explicitly rather than
silently assuming they are equivalent. It currently treats monochromatic
transition arrays, planar isotropic media and azimuth-integrated collection.
Broadband, spatially distributed fixed-emitter populations require per-population
integrated rate denominators before averaging; they are not implemented here.

## Stack-derived optical rates and remaining physics

`emitter_rate_budget` solves the emitter in a `DispersiveStack` at each wavelength and returns an
`OpticalRateBudget` whose escape rate is the flux into the chosen exit half-space (`"exit_medium"`)
or a first-pass lower bound into air through a thick exit medium (`"air_first_pass"`):

```python
from multilayer_dipole import emitter_rate_budget

budget = emitter_rate_budget(
    stack, lake, wavelengths_nm, exit="bottom", collection="air_first_pass",
    vertical_fraction=0.2,            # orientation order parameter, not a re-orienting ensemble
)
eqe = budget.monochromatic_eqe(intrinsic_yield=0.8, charge_balance=0.9, exciton_yield=1.0)
```

The solver refuses lossless guiding layers (their guided-mode power is a delta function the
quadrature cannot see); give a guiding layer a small passive `k`. Still to be built and validated:
guided/plasmon/absorption separation, anisotropic layers, emission-zone weighting from charge
transport, and substrate-to-air recycling. Plane-wave transmittance must not replace the
collected-emission calculation, and measured-device comparison is required before any
performance claim.
