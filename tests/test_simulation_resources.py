import math
import json
import sys

import pytest

from simulation_resources import Hardware, emitter_population, fdtd_estimate, planar_estimate, pixel_study, benchmark_cpu_reference, main


def test_planar_two_bases_reuse_random_ensemble():
    kwargs = dict(wavelengths=31, finite_layers=6, emitter_depths=10, quadrature_evaluations=1000)
    one = planar_estimate(**kwargs)
    many = planar_estimate(**kwargs, ensemble_emitters=1000000)
    assert one["basis_integrals"] == 620
    assert one["layer_evaluations"] == 620 * 1000 * 6
    assert many["assumed_flops"] == one["assumed_flops"]
    assert many["streamed_workspace_bytes"] == one["streamed_workspace_bytes"]
    assert many["sampled_orientation_storage_bytes_if_retained"] == 24_000_000
    assert many["orientation_composition_operations"] == 20 * 31 * 1000000
    sweep = planar_estimate(**kwargs, thickness_configurations=5)
    assert sweep["assumed_flops"] == 5 * one["assumed_flops"]
    assert sweep["retained_rate_output_bytes"] == 5 * one["retained_rate_output_bytes"]


def test_area_volume_population_scaling():
    small = emitter_population(area_um2=100, depth_um=.1, density_per_um3=1000)
    large = emitter_population(area_um2=400, depth_um=.2, density_per_um3=1000)
    assert small == {"area_um2": 100, "volume_um3": 10, "expected_emitters": 10000}
    assert large["expected_emitters"] == 8 * small["expected_emitters"]


def test_pixel_study_compares_area_volume_and_independent_source_work():
    result = pixel_study(side_lengths_um=(1, 10), volume_depth_samples=5, fdtd_source_positions=7)
    small, large = result["cases"]
    assert large["area_um2"] == 100 * small["area_um2"]
    assert large["volume_um3"] == 100 * small["volume_um3"]
    assert large["expected_emitters"] == 100 * small["expected_emitters"]
    assert small["emitter_volume"]["basis_integrals"] == 5 * small["emitter_sheet"]["basis_integrals"]
    assert large["emitter_volume"]["basis_integrals"] == small["emitter_volume"]["basis_integrals"]
    assert small["finite_pixel_fdtd"]["independent_source_runs"] == 21
    assert large["finite_pixel_fdtd"]["cells"] > small["finite_pixel_fdtd"]["cells"]
    assert "cpu_ideal_roofline_seconds_lower_bound" not in small["finite_pixel_fdtd"]


def test_cpu_benchmark_measures_reference_not_full_solver():
    report = benchmark_cpu_reference(repeats=1)
    homogeneous, lossy = report["cases"]
    assert homogeneous["vertical_source_work"] == homogeneous["parallel_source_work"] == 1
    assert lossy["vertical_source_work"] > 1
    assert lossy["parallel_source_work"] > 0
    assert homogeneous["repeats"] == 1
    assert homogeneous["median_wall_seconds"] > 0
    assert set(homogeneous["orientation_source_work_not_eqe"].values()) == {1.0}
    orientations = lossy["orientation_source_work_not_eqe"]
    assert orientations["isotropic_analytic"] == pytest.approx(3.966116501134272)
    assert orientations["tilted_45_degrees"] == pytest.approx(4.9140864628373)
    assert orientations["isotropic_sampled"] == pytest.approx(orientations["isotropic_analytic"], abs=.03)
    assert report["orientation_sampling"]["count"] == 10000
    assert report["orientation_sampling"]["seed"] == 42
    assert "no GPU kernel" in report["gpu"]
    with pytest.raises(ValueError):
        benchmark_cpu_reference(repeats=0)


def test_resource_cli_exports_requested_sizes_without_overwrite(tmp_path, monkeypatch):
    output = tmp_path / "report.json"
    monkeypatch.setattr(sys, "argv", ["resource-study", "--output", str(output), "--side-um", "2", "4"])
    main()
    report = json.loads(output.read_text())
    assert [case["pixel_side_um"] for case in report["cases"]] == [2, 4]
    assert report["cases"][0]["volume_um3"] == pytest.approx(.4)
    previous = output.read_bytes()
    with pytest.raises(FileExistsError):
        main()
    assert output.read_bytes() == previous


def test_fdtd_grid_cfl_memory_and_roofline_arithmetic():
    hardware = Hardware(100, 20, 1000, 100, 16, 8)
    result = fdtd_estimate(size_um=(1, 2, 3), cell_nm=100, duration_fs=10,
                           pml_cells=2, hardware=hardware)
    assert result["grid_shape_including_pml"] == (14, 24, 34)
    cells = 14 * 24 * 34
    assert result["cells"] == cells
    assert result["raw_grid_bytes"] == cells * 24 * 4
    assert result["provisioned_grid_bytes_30pct_headroom"] == math.ceil(cells * 24 * 4 * 1.3)
    dt = .99 * 100e-9 / (299792458 * math.sqrt(3)) * 1e15
    assert result["courant_step_fs"] == pytest.approx(dt)
    assert result["time_steps_per_run"] == math.ceil(10 / dt)
    work = result["cell_updates_all_runs"]
    assert result["cpu_ideal_roofline_seconds_lower_bound"] == pytest.approx(max(work * 300 / 100e9, work * 192 / 20e9))
    assert result["gpu_ideal_roofline_seconds_lower_bound"] == pytest.approx(max(work * 300 / 1000e9, work * 192 / 100e9))
    assert result["fits_cpu_ram_budget"] and result["fits_gpu_vram_budget"]


def test_fdtd_independent_emitters_scale_work_not_serial_memory():
    kwargs = dict(size_um=(10, 10, 1), cell_nm=10, duration_fs=100)
    one = fdtd_estimate(**kwargs)
    ensemble = fdtd_estimate(**kwargs, independent_source_runs=300)
    assert ensemble["cell_updates_all_runs"] == 300 * one["cell_updates_all_runs"]
    assert ensemble["raw_grid_bytes"] == one["raw_grid_bytes"]
    double = fdtd_estimate(**kwargs, bytes_per_real=8)
    assert double["raw_grid_bytes"] == 2 * one["raw_grid_bytes"]
    assert "gpu_ideal_roofline_seconds_lower_bound" not in one


def test_fdtd_memory_fit_can_fail_and_cfl_uses_fastest_index():
    result = fdtd_estimate(size_um=(100, 100, 1), cell_nm=10, duration_fs=100,
                           hardware=Hardware(100, 20, 1000, 100, 16, 8))
    assert not result["fits_cpu_ram_budget"]
    assert not result["fits_gpu_vram_budget"]
    slower = fdtd_estimate(size_um=(1, 1, 1), cell_nm=10, duration_fs=100, minimum_index=2)
    faster = fdtd_estimate(size_um=(1, 1, 1), cell_nm=10, duration_fs=100, minimum_index=1)
    assert slower["courant_step_fs"] == 2 * faster["courant_step_fs"]
    assert slower["time_steps_per_run"] <= faster["time_steps_per_run"]


@pytest.mark.parametrize("kwargs", [
    {"cell_nm": 0}, {"duration_fs": -1}, {"minimum_index": float("nan")},
    {"pml_cells": 0}, {"bytes_per_real": 3}, {"allocated_real_arrays_per_cell": 5},
    {"independent_source_runs": 1.5}, {"size_um": (1, 2)}, {"size_um": (1, 0, 1)},
])
def test_invalid_fdtd_assumptions(kwargs):
    args = dict(size_um=(1, 1, 1), cell_nm=10, duration_fs=100)
    args.update(kwargs)
    with pytest.raises(ValueError):
        fdtd_estimate(**args)


@pytest.mark.parametrize("value", [0, -1, 1.5, True])
def test_invalid_planar_counts(value):
    with pytest.raises(ValueError):
        planar_estimate(wavelengths=value, finite_layers=2, emitter_depths=1, quadrature_evaluations=100)


def test_invalid_hardware_and_population():
    with pytest.raises(ValueError):
        Hardware(0, 20, 1000, 100, 16, 8)
    with pytest.raises(ValueError):
        emitter_population(area_um2=1, depth_um=0, density_per_um3=1000)
