import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

import tmm_OLED9 as model


@pytest.mark.parametrize("polarization", ["s", "p"])
@pytest.mark.parametrize("angle", [0.0, 30.0, 65.0, 89.0])
def test_zero_thickness_is_air(polarization, angle):
    result = model.power_coefficients(550.0, angle, 0.0, polarization, silver_nm=0.0)
    assert result["R"] == pytest.approx(0.0, abs=1e-12)
    assert result["T"] == pytest.approx(1.0, abs=1e-12)
    assert result["A"] == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("polarization", ["s", "p"])
@pytest.mark.parametrize("angle", [0.0, 30.0, 65.0])
def test_interface_matches_fresnel(monkeypatch, polarization, angle):
    monkeypatch.setattr(model, "stack_at", lambda *args: ([1.0, 1.5], [np.inf, np.inf]))
    theta = np.deg2rad(angle)
    incident_cos = np.cos(theta)
    transmitted_cos = np.sqrt(1 - (np.sin(theta) / 1.5) ** 2)
    if polarization == "s":
        r = (incident_cos - 1.5 * transmitted_cos) / (incident_cos + 1.5 * transmitted_cos)
    else:
        r = (1.5 * incident_cos - transmitted_cos) / (1.5 * incident_cos + transmitted_cos)
    result = model.power_coefficients(550.0, angle, 0.0, polarization)
    assert result["R"] == pytest.approx(r**2, abs=1e-12)
    assert result["T"] == pytest.approx(1 - r**2, abs=1e-12)


@pytest.mark.parametrize("polarization", ["s", "p"])
def test_lossless_film_matches_airy_formula(polarization):
    wavelength, thickness = 550.0, 75.0
    index = 1.65 + 0.5 / 550.0**2
    interface_reflectance = ((index - 1) / (index + 1)) ** 2
    phase = 2 * np.pi * index * thickness / wavelength
    expected = 1 / (
        1 + 4 * interface_reflectance * np.sin(phase) ** 2 / (1 - interface_reflectance) ** 2
    )
    result = model.power_coefficients(wavelength, 0.0, thickness, polarization, silver_nm=0)
    assert result["T"] == pytest.approx(expected, abs=1e-12)
    assert result["R"] == pytest.approx(1 - expected, abs=1e-12)
    assert result["A"] == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("polarization", ["s", "p"])
def test_absorbing_stack_matches_characteristic_matrix(polarization):
    # Independent normal-incidence field-matrix formulation, not another coh_tmm call.
    matrix = np.eye(2, dtype=complex)
    for index, thickness in [(1.65 + 0.5 / 550**2, 75), (0.16 - 1.5 / 550**2 + 3.16j, 50)]:
        phase = 2 * np.pi * index * thickness / 550
        matrix = matrix @ np.array([
            [np.cos(phase), -1j * np.sin(phase) / index],
            [-1j * index * np.sin(phase), np.cos(phase)],
        ])
    b = matrix[0, 0] + matrix[0, 1]
    c = matrix[1, 0] + matrix[1, 1]
    expected_r = abs((b - c) / (b + c)) ** 2
    expected_t = abs(2 / (b + c)) ** 2
    result = model.power_coefficients(550.0, 0.0, 75.0, polarization)
    assert result["R"] == pytest.approx(expected_r, abs=1e-12)
    assert result["T"] == pytest.approx(expected_t, abs=1e-12)
    assert result["A"] == pytest.approx(1 - expected_r - expected_t, abs=1e-12)
    assert 0 < result["A"] < 1


@pytest.mark.parametrize("polarization", ["s", "p"])
def test_lossless_oblique_power_balance(polarization):
    for wavelength in [400.0, 550.0, 700.0]:
        for angle in [0.0, 45.0, 89.0]:
            result = model.power_coefficients(wavelength, angle, 100.0, polarization, silver_nm=0)
            assert result["R"] + result["T"] == pytest.approx(1.0, abs=1e-12)
            assert result["A"] == pytest.approx(0.0, abs=1e-12)


def test_sweep_shape_weighting_and_high_angle_transmission():
    values = model.sweep([450.0, 550.0], [0.0, 65.0], 75.0)
    assert values.shape == (2, 2)
    for i, wavelength in enumerate([450.0, 550.0]):
        for j, angle in enumerate([0.0, 65.0]):
            expected = sum(model.transmission(wavelength, angle, 75.0, p) for p in ("s", "p")) / 2
            assert values[i, j] == pytest.approx(expected, abs=1e-12)
    assert np.all(values > 0)
    assert np.all(values <= 1)


@pytest.mark.parametrize("wavelength,angle,organic,silver,polarization", [
    (0, 0, 75, 50, "s"), (-1, 0, 75, 50, "s"), (np.nan, 0, 75, 50, "s"),
    (np.inf, 0, 75, 50, "s"), (550, -1, 75, 50, "s"), (550, 90, 75, 50, "s"),
    (550, np.nan, 75, 50, "s"), (550, np.inf, 75, 50, "s"),
    (550, 0, -1, 50, "s"), (550, 0, np.nan, 50, "s"),
    (550, 0, np.inf, 50, "s"), (550, 0, 75, -1, "s"),
    (550, 0, 75, np.nan, "s"), (550, 0, 75, np.inf, "s"), (550, 0, 75, 50, "x"),
])
def test_invalid_scalar_inputs(wavelength, angle, organic, silver, polarization):
    with pytest.raises(ValueError):
        model.power_coefficients(wavelength, angle, organic, polarization, silver)


@pytest.mark.parametrize("wavelengths,angles,thickness", [
    ([], [0], 75), ([550], [], 75), ([[550]], [0], 75), ([550], [[0]], 75),
    ([0], [0], 75), ([np.nan], [0], 75), ([550], [90], 75),
    ([550], [np.inf], 75), ([550], [0], -1),
])
def test_invalid_sweep_inputs(wavelengths, angles, thickness):
    with pytest.raises(ValueError):
        model.sweep(wavelengths, angles, thickness)


def test_cli_exports_plot_csv_and_metadata(tmp_path, capsys):
    output = tmp_path / "preview.png"
    model.main(["--wavelength-count", "2", "--angle-count", "3", "--output", str(output)])
    assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    with output.with_suffix(".csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 6
    assert [(float(row["wavelength_nm"]), float(row["angle_deg"])) for row in rows] == [
        (400, 0), (400, 44.5), (400, 89), (700, 0), (700, 44.5), (700, 89),
    ]
    expected = model.sweep([400, 700], [0, 44.5, 89], 75)
    np.testing.assert_allclose(
        [float(row["unpolarized_transmittance"]) for row in rows], expected.ravel(), atol=1e-12
    )
    metadata = json.loads(output.with_suffix(".json").read_text(encoding="utf-8"))
    assert metadata["quantity"] == "unpolarized_plane_wave_transmittance"
    assert metadata["organic_thickness_nm"] == 75
    assert metadata["silver_thickness_nm"] == 50
    assert metadata["wavelengths_nm"] == [400, 700]
    assert metadata["angles_deg"] == [0, 44.5, 89]
    assert metadata["software_version"] == model.__version__
    assert "NOT a validated OLED solver" in capsys.readouterr().err


@pytest.mark.parametrize("arguments,message", [
    (["--wavelength-count", "1"], "sample counts"),
    (["--angle-count", "1"], "sample counts"),
    (["--thickness-nm", "nan"], "finite and nonnegative"),
    (["--thickness-nm", "-1"], "finite and nonnegative"),
])
def test_cli_rejects_invalid_inputs(tmp_path, capsys, arguments, message):
    with pytest.raises(SystemExit) as error:
        model.main(arguments + ["--output", str(tmp_path / "preview.png")])
    assert error.value.code != 0
    assert message in capsys.readouterr().err
    assert list(tmp_path.iterdir()) == []


def test_cli_does_not_overwrite_existing_sidecar(tmp_path, capsys):
    output = tmp_path / "preview.png"
    output.with_suffix(".csv").write_text("keep me", encoding="utf-8")
    with pytest.raises(SystemExit) as error:
        model.main(["--output", str(output)])
    assert error.value.code == 2
    assert "output already exists" in capsys.readouterr().err
    assert output.with_suffix(".csv").read_text(encoding="utf-8") == "keep me"
    assert not output.exists()


def test_cli_explicit_overwrite_replaces_outputs(tmp_path):
    output = tmp_path / "preview.png"
    for suffix in [".png", ".csv", ".json"]:
        output.with_suffix(suffix).write_text("old content", encoding="utf-8")
    model.main([
        "--wavelength-count", "2", "--angle-count", "2",
        "--output", str(output), "--overwrite",
    ])
    assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert output.with_suffix(".csv").read_text(encoding="utf-8").startswith("wavelength_nm,")
    assert json.loads(output.with_suffix(".json").read_text(encoding="utf-8"))["schema_version"] == 1


@pytest.mark.parametrize("relative_path,message", [
    ("bad.svg", ".png extension"),
    ("missing/preview.png", "parent directory"),
])
def test_cli_rejects_invalid_output_paths(tmp_path, capsys, relative_path, message):
    with pytest.raises(SystemExit) as error:
        model.main(["--output", str(tmp_path / relative_path)])
    assert error.value.code != 0
    assert message in capsys.readouterr().err
    assert list(tmp_path.iterdir()) == []


def test_cli_reports_write_failure(tmp_path, monkeypatch, capsys):
    def fail_open(*args, **kwargs):
        raise PermissionError("test write denied")

    monkeypatch.setattr(Path, "open", fail_open)
    with pytest.raises(SystemExit) as error:
        model.main([
            "--wavelength-count", "2", "--angle-count", "2",
            "--output", str(tmp_path / "preview.png"),
        ])
    assert error.value.code == 1
    assert "test write denied" in capsys.readouterr().err


def test_import_has_no_plot_or_file_side_effects(tmp_path):
    env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]))
    result = subprocess.run(
        [sys.executable, "-c", "import sys; import tmm_OLED9; assert 'matplotlib.pyplot' not in sys.modules"],
        cwd=tmp_path, env=env, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    assert list(tmp_path.iterdir()) == []
