# Orientation ensembles and CPU/GPU/RAM planning

Prepared 2026-10-16. APIs: `spectral_eqe.py` and `simulation_resources.py`.
Estimates are explicit workload models, NOT a GPU implementation, hardware
performance guarantee, measured peak memory, or full-pixel field simulation.

## Three distinct source problems

1. **Incoherent dipoles in one plane:** common emitting depth; unrelated
   phases; sum emitted powers, not amplitudes.
2. **Incoherent volume ensemble:** distribution of depths and orientations.
   In an infinite planar isotropic stack, lateral translation does not change
   the optical response; depth and orientation distributions still matter.
3. **Coherent plane-wave illumination or phased emitter plane:** this is not
   the same ensemble. It has different source correlations, energy normalization,
   interference and computational requirements. The existing plane-wave preview
   computes external illumination, not emission.

For azimuth-integrated power in planar isotropic media, compute vertical and
parallel basis responses once per wavelength and depth. A tilted dipole has
vertical power weight cos(theta)^2. Random isotropic orientations have expected
vertical weight 1/3 and parallel weight 2/3. Horizontal-only, vertical-only and
tilted distributions can reuse these bases. Direction-resolved collection
in a restricted azimuth, anisotropic media or coherent emitters need more than
this two-basis integrated-power model.

```python
from spectral_eqe import random_dipoles, dipole_ensemble_eqe

directions = random_dipoles(10000, seed=42)  # uniform solid angle, reproducible
# vertical_budget and horizontal_budget must share wavelength/reference.
comparison = dipole_ensemble_eqe(
    vertical_budget, horizontal_budget, directions,
    intrinsic_yield=0.8, charge_balance=0.9, exciton_yield=1.0,
)
```

The result compares equal-excitation **independent fixed dipoles** (average
their yields) with **rapid reorientation** (average rates then form the yield).
These can differ; an orientation-averaged Purcell factor alone does not imply
the same EQE for both ensembles. Correlations, exciton transfer, saturation,
orientation-dependent intrinsic emission and superradiance are not included.

For a uniform laterally infinite stack, area increases emitter count, photons
and current together; normalized EQE does not inherently scale with area.
Finite pixel edges, contacts, scattering, current crowding and local heating
break that simplification. No finite-pixel EQE is claimed here.

## Planar workload model

Let L be wavelength count, Z sampled depths, C thickness configurations,
J finite layers, Q quadrature evaluations per orientation including refinements.

    basis integrals = 2 * L * Z * C
    layer-response evaluations = 2 * L * Z * C * Q * J

The estimator reports an adjustable arithmetic count, a streamed workspace,
retained rate-array storage, and optional N-dipole direction storage/averaging
cost. The workspace omits Python/SciPy interpreter/allocator overhead and any
large field/angle output. It is not process peak RAM. A streaming analytic
orientation distribution needs no N-by-wavelength matrix.

Increasing physical area or active-emitter density does not multiply these
planar basis integrals. A Cartesian sweep of thicknesses CAN explode:
C is the product of samples per swept layer, not the sum. Adaptive quadrature
cost depends on optical constants, separations and resonances; Q is not fixed
by geometry alone. Array count/FLOP factors are assumptions, not measured code.

The implemented single-interface kernels use scalar SciPy adaptive quadrature
on the CPU. No CUDA/ROCm/GPU kernel exists. Independent wavelengths/depths can
eventually be parallelized, but GPU performance cannot be inferred from the
device model name or CPU core count.

## Full 3-D pixel FDTD planning model

The alternative estimator sizes a uniform 3-D Yee grid with PML added to every
face. Defaults for comparison:

- Pixel x/y side: 1, 10, 50, 100 um; simulated z: 1 um.
- Active emitter region: 0.1 um; illustrative active density 1000 per um^3.
  This is not a measured exciton density or total molecular density.
- Grid spacing: 10 nm; 10 PML cells per face; single-precision real fields.
- 24 real arrays per cell, plus 30% memory headroom.
  Six are electromagnetic field components; auxiliary dispersive/PML storage
  is budgeted, not derived for a specific material implementation.
- Physical duration: 200 fs; Courant factor 0.99.
- 20 SAMPLE source positions, 3 independent dipole axes each, run sequentially.
  This is NOT every emitter in the volume and NOT a demonstrated converged
  sampling density. Finite pixels require converged lateral AND depth sampling.
- 300 assumed FLOPs and 192 bytes of traffic per cell-step. Change these for
  the implemented update scheme; double precision also needs a correspondingly
  revised traffic assumption, not merely doubled allocated field storage.

For spacing delta and minimum nondispersive index n_min, the reference
Courant timestep is `0.99*n_min*delta/(c*sqrt(3))`. The fastest medium controls
this bound; a maximum-index CFL formula would be wrong. Dispersive models can
require smaller timesteps. Spatial resolution must also resolve wavelength
inside the largest-index medium, skin depths, thin films and geometric edges.
Ten nm/200 fs is a scenario, not a convergence certification.

| Pixel side (um) | Active volume (um^3) | Planar sheet basis integrals | Planar 20-depth integrals | FDTD provisioned grid GiB |
|---|---:|---:|---:|---:|
| 1 | 0.1 | 62 | 1,240 | 0.201 |
| 10 | 10 | 62 | 1,240 | 14.511 |
| 50 | 250 | 62 | 1,240 | 351.482 |
| 100 | 1,000 | 62 | 1,240 | 1,400.333 |

The planar rows assume 31 wavelengths and one thickness configuration.
FDTD RAM is for ONE sequential run, not 60 resident copies. Work multiplies
by independent source runs; memory multiplies only if those runs are simultaneous.
Summing many phase-coherent sources in one run is not an incoherent ensemble.
Stochastic-source/reciprocity methods can reduce work but require separate
implementation and statistical validation.

The grid estimate excludes full-field snapshots, frequency-domain monitors,
geometry metadata, MPI/halo copies and CPU/GPU transfers. Volume snapshots
can easily dominate storage. A full-spectrum FDTD pulse can cover several
wavelengths in one run if duration/sampling and dispersion are adequate;
do not multiply its cost by planar wavelength count automatically.

## This workstation and honest runtime limits

Read-only inspection found an Intel Core i9-13900HX (24 cores, 32 logical
processors), approximately 31.69 GiB system RAM and an RTX 4070 Laptop GPU
reporting 8188 MiB VRAM. These are observed capacity facts, not measured
sustained bandwidth or floating-point throughput.

Under the table's model, the 10 um pixel exceeds that GPU's memory but fits
the nominal system RAM budget; other application allocations/OS overhead can
still make it impractical. The 50/100 um cases exceed both capacities by large
factors. No multi-GPU capability is implied.

`Hardware` optionally accepts USER-SUPPLIED sustained CPU/GPU GFLOP/s and GB/s.
The estimator then reports only the idealized lower bound:

    seconds >= max(assumed_FLOPs / sustained_FLOPs_per_second,
                   assumed_memory_traffic / sustained_bytes_per_second)

This roofline lower bound is not elapsed-time prediction. If a grid does not
fit, the single-device bound does not mean the run is feasible. Without supplied
rates, no CPU/GPU runtime is invented. The planning tool itself does not execute
FDTD on either processor.

## Reproduce the scenario report and CPU reference benchmark

```powershell
python -m simulation_resources --output pixel-resource-study.json --benchmark-cpu
```

Or use installed `mongelli-resource-study` with the same flags. Output creation
does not overwrite existing files. `--side-um`, `--emitter-depth-um`, and
`--density-per-um3` customize selected dimensions; use Python APIs for all
quadrature, geometry, precision and hardware assumptions.

The optional CPU benchmark measures only the implemented single-interface
reference at 550 nm/20 nm for matched and synthetic lossy substrates, with
warmup and five sequential repeats. It reports median/min/max wall time and
quadrature counts. It does NOT benchmark multilayer OLED emission, full
volumetric FDTD, GPU execution, concurrency scaling or peak memory.
It also composes horizontal, vertical, 45-degree and isotropic source-work
responses from these two computed bases, including a 10,000-dipole sample
with seed 42. These are electromagnetic decay/source-work comparisons,
not collected photon fractions or EQE. The lossy reference has relative
source work approximately 2.07018 horizontally, 7.75800 vertically and
3.96612 for an analytic isotropic ensemble.

For actual procurement/runtime estimates, next benchmark the completed source
solver on representative material stacks and convergence targets, then implement
and measure any GPU backend. Retain workload definitions and uncertainty ranges.
