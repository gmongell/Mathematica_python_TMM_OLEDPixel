"""Transparent workload/RAM estimates, not a GPU solver or measured speed prediction."""

from dataclasses import dataclass
import argparse
import json
import math
from numbers import Real
from pathlib import Path
import statistics
import time


def _positive(name: str, value: float) -> None:
    if not isinstance(value, Real) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")


def _count(name: str, value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")


@dataclass(frozen=True)
class Hardware:
    """User-supplied SUSTAINED rates, not vendor peak/nameplate specifications."""

    cpu_gflops: float
    cpu_bandwidth_gb_s: float
    gpu_gflops: float
    gpu_bandwidth_gb_s: float
    system_ram_gib: float
    gpu_vram_gib: float

    def __post_init__(self):
        for name, value in vars(self).items():
            _positive(name, value)


def _roofline(flops: float, traffic_bytes: float, gflops: float, gb_s: float) -> float:
    return max(flops / (gflops * 1e9), traffic_bytes / (gb_s * 1e9))


def planar_estimate(
    *, wavelengths: int, finite_layers: int, emitter_depths: int,
    quadrature_evaluations: int, thickness_configurations: int = 1,
    ensemble_emitters: int = 1, bytes_per_real: int = 8,
    flops_per_layer_evaluation: int = 200,
) -> dict:
    """Two orientation bases per depth; uniform lateral area needs no new solves.

    Quadrature count is per orientation and includes convergence refinements.
    Assumes streaming wavelength/depth integration, no full-field spatial output.
    """
    arguments = locals().copy()
    for name, value in arguments.items():
        _count(name, value)
    if bytes_per_real not in (4, 8):
        raise ValueError("Real precision must be 4 or 8 bytes")
    solves = 2 * wavelengths * emitter_depths * thickness_configurations
    layer_evaluations = solves * quadrature_evaluations * finite_layers
    flops = layer_evaluations * flops_per_layer_evaluation
    # Workspace, index tables and a reusable pair of rate spectra; not per emitter.
    working_bytes = bytes_per_real * (
        16 * (finite_layers + 2) + 2 * wavelengths * (finite_layers + 2) + 8 * wavelengths
    )
    retained_bytes = bytes_per_real * 4 * wavelengths * emitter_depths * thickness_configurations
    return {
        "method": "planar_incoherent_two_orientation_basis",
        "basis_integrals": solves,
        "layer_evaluations": layer_evaluations,
        "assumed_flops": flops,
        "streamed_workspace_bytes": working_bytes,
        "retained_rate_output_bytes": retained_bytes,
        "sampled_orientation_storage_bytes_if_retained": 3 * bytes_per_real * ensemble_emitters,
        "orientation_composition_operations": 20 * wavelengths * ensemble_emitters,
        "assumptions": arguments,
        "runtime_status": "No seconds predicted: complex recursion/quadrature throughput must be benchmarked",
        "gpu_status": "No GPU implementation; scalar adaptive SciPy quadrature currently executes on CPU",
    }


def emitter_population(*, area_um2: float, depth_um: float, density_per_um3: float) -> dict:
    for name, value in locals().copy().items():
        _positive(name, value)
    volume = area_um2 * depth_um
    expected = volume * density_per_um3
    if not math.isfinite(expected):
        raise ValueError("Emitter population overflowed")
    return {"area_um2": area_um2, "volume_um3": volume, "expected_emitters": expected}


def fdtd_estimate(
    *, size_um: tuple[float, float, float], cell_nm: float, duration_fs: float,
    minimum_index: float = 1.0, pml_cells: int = 10, bytes_per_real: int = 4,
    allocated_real_arrays_per_cell: int = 24, flops_per_cell_step: int = 300,
    bytes_transferred_per_cell_step: int = 192, independent_source_runs: int = 1,
    hardware: Hardware | None = None,
) -> dict:
    """Uniform 3-D Yee-grid planning model with configurable ADE/PML array budget.

    CFL uses the minimum nondispersive refractive index (fastest phase velocity),
    NOT the largest index. Dispersive media can impose stricter stability bounds.
    """
    if len(size_um) != 3:
        raise ValueError("Physical domain must have three dimensions in micrometers")
    for size in size_um:
        _positive("Domain dimension (um)", size)
    for name, value in (("Cell size (nm)", cell_nm), ("Duration (fs)", duration_fs), ("Minimum index", minimum_index)):
        _positive(name, value)
    for name, value in (
        ("PML cells", pml_cells), ("Array count", allocated_real_arrays_per_cell),
        ("FLOPs per cell step", flops_per_cell_step),
        ("Bytes transferred per cell step", bytes_transferred_per_cell_step),
        ("Independent source runs", independent_source_runs),
    ):
        _count(name, value)
    if allocated_real_arrays_per_cell < 6 or bytes_per_real not in (4, 8):
        raise ValueError("Need at least six field arrays and 4- or 8-byte real precision")
    shape = tuple(math.ceil(size * 1000 / cell_nm) + 2 * pml_cells for size in size_um)
    cells = math.prod(shape)
    dt_fs = 0.99 * minimum_index * cell_nm * 1e-9 / (299792458 * math.sqrt(3)) * 1e15
    steps = math.ceil(duration_fs / dt_fs)
    updates = cells * steps * independent_source_runs
    raw_bytes = cells * allocated_real_arrays_per_cell * bytes_per_real
    provision_bytes = math.ceil(raw_bytes * 1.3)
    result = {
        "method": "uniform_3d_fdtd_planning_only",
        "grid_shape_including_pml": shape, "cells": cells,
        "courant_step_fs": dt_fs, "time_steps_per_run": steps,
        "cell_updates_all_runs": updates,
        "raw_grid_bytes": raw_bytes, "provisioned_grid_bytes_30pct_headroom": provision_bytes,
        "assumed_flops": updates * flops_per_cell_step,
        "assumed_memory_traffic_bytes": updates * bytes_transferred_per_cell_step,
        "independent_source_runs": independent_source_runs,
        "limitations": (
            "No FDTD implementation. Excludes full-field snapshots, DFT monitors, geometry, "
            "MPI copies, device transfers, dispersive stability analysis and convergence runs. "
            "Array/FLOP/traffic counts are user-adjustable planning assumptions. "
            "One coherent source run does not represent an incoherent emitter ensemble."
        ),
    }
    if hardware is not None:
        result.update({
            "cpu_ideal_roofline_seconds_lower_bound": _roofline(result["assumed_flops"], result["assumed_memory_traffic_bytes"], hardware.cpu_gflops, hardware.cpu_bandwidth_gb_s),
            "gpu_ideal_roofline_seconds_lower_bound": _roofline(result["assumed_flops"], result["assumed_memory_traffic_bytes"], hardware.gpu_gflops, hardware.gpu_bandwidth_gb_s),
            "fits_cpu_ram_budget": provision_bytes <= hardware.system_ram_gib * 2**30,
            "fits_gpu_vram_budget": provision_bytes <= hardware.gpu_vram_gib * 2**30,
            "hardware_assumptions": vars(hardware),
        })
    return result


def pixel_study(
    *, side_lengths_um: tuple[float, ...] = (1, 10, 50, 100),
    emitter_depth_um: float = .1, active_density_per_um3: float = 1000,
    simulation_depth_um: float = 1, volume_depth_samples: int = 20,
    fdtd_source_positions: int = 20, hardware: Hardware | None = None,
) -> dict:
    """Compare sheet/volume planar ensembles with a finite-pixel reference workload."""
    _count("Volume depth samples", volume_depth_samples)
    _count("FDTD sampled source positions", fdtd_source_positions)
    _positive("Simulation depth (um)", simulation_depth_um)
    if not side_lengths_um:
        raise ValueError("At least one pixel side length is required")
    if emitter_depth_um > simulation_depth_um:
        raise ValueError("Emission region cannot exceed the simulated domain depth")
    rows = []
    for side in side_lengths_um:
        _positive("Pixel side (um)", side)
        population = emitter_population(
            area_um2=side**2, depth_um=emitter_depth_um, density_per_um3=active_density_per_um3,
        )
        common = dict(wavelengths=31, finite_layers=6, quadrature_evaluations=2000,
                      ensemble_emitters=max(1, math.ceil(population["expected_emitters"])))
        rows.append({
            "pixel_side_um": side, **population,
            "emitter_sheet": planar_estimate(**common, emitter_depths=1),
            "emitter_volume": planar_estimate(**common, emitter_depths=volume_depth_samples),
            "finite_pixel_fdtd": fdtd_estimate(
                size_um=(side, side, simulation_depth_um), cell_nm=10, duration_fs=200,
                independent_source_runs=3 * fdtd_source_positions, hardware=hardware,
            ),
        })
    return {
        "status": "Planning estimates, not an implemented finite-stack or full-wave simulation",
        "assumptions": {
            "planar_wavelengths": 31, "finite_layers": 6,
            "quadrature_evaluations_per_orientation": 2000,
            "volume_depth_samples": volume_depth_samples,
            "active_density_per_um3": active_density_per_um3,
            "fdtd_source_positions": fdtd_source_positions,
            "fdtd_source_runs_per_position": 3,
            "fdtd_cell_nm": 10, "fdtd_duration_fs": 200,
            "fdtd_simulation_depth_um": simulation_depth_um,
            "fdtd_real_arrays_per_cell": 24, "fdtd_precision_bytes": 4,
            "fdtd_pml_cells_per_face": 10,
        },
        "cases": rows,
    }


def benchmark_cpu_reference(*, repeats: int = 5) -> dict:
    """Measure only the implemented single-interface CPU kernels on this machine."""
    from dipole_reference import SingleInterface
    from spectral_eqe import random_dipoles

    _count("Benchmark repeats", repeats)
    directions = random_dipoles(10000, seed=42)
    sampled_vertical_fraction = float((directions[:, 2] ** 2).mean())
    records = []
    for name, epsilon in [("matched_homogeneous", 1.0), ("synthetic_lossy_interface", -5 + .5j)]:
        interface = SingleInterface(550, 20, epsilon)
        interface.decay_rates()
        samples = []
        for _ in range(repeats):
            start = time.perf_counter()
            result = interface.decay_rates()
            samples.append(time.perf_counter() - start)
        records.append({
            "case": name, "wavelength_nm": 550, "distance_nm": 20,
            "median_wall_seconds": statistics.median(samples),
            "min_wall_seconds": min(samples), "max_wall_seconds": max(samples),
            "repeats": repeats, "vertical_source_work": result.vertical,
            "parallel_source_work": result.parallel,
            "quadrature_evaluations_per_call": result.evaluations,
            "orientation_source_work_not_eqe": {
                "vertical": result.vertical,
                "horizontal_plane": result.parallel,
                "tilted_45_degrees": .5 * result.vertical + .5 * result.parallel,
                "isotropic_analytic": result.vertical / 3 + 2 * result.parallel / 3,
                "isotropic_sampled": sampled_vertical_fraction * result.vertical + (1 - sampled_vertical_fraction) * result.parallel,
            },
        })
    return {
        "method": "Measured sequential CPU single-interface reference; no finite-stack extrapolation",
        "cases": records,
        "orientation_sampling": {
            "count": len(directions), "seed": 42,
            "mean_vertical_fraction": sampled_vertical_fraction,
            "quantity": "Total electromagnetic source work / homogeneous-host power, not collected flux or EQE",
        },
        "gpu": "Not measured: no GPU kernel is implemented",
        "memory": "No process peak-RAM measurement; use the separately labeled planning model",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--side-um", nargs="+", type=float, default=[1, 10, 50, 100])
    parser.add_argument("--emitter-depth-um", type=float, default=.1)
    parser.add_argument("--density-per-um3", type=float, default=1000)
    parser.add_argument("--benchmark-cpu", action="store_true")
    args = parser.parse_args()
    report = pixel_study(side_lengths_um=tuple(args.side_um), emitter_depth_um=args.emitter_depth_um,
                         active_density_per_um3=args.density_per_um3)
    if args.benchmark_cpu:
        report["measured_reference_cpu"] = benchmark_cpu_reference()
    with args.output.open("x", encoding="utf-8") as output:
        json.dump(report, output, indent=2, allow_nan=False)
        output.write("\n")
    print(f"Created planning report: {args.output}; no GPU execution or runtime measurement claimed")


if __name__ == "__main__":
    main()
