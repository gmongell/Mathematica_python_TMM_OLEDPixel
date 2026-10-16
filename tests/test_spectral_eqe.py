import hashlib
import json
from unittest.mock import patch
import io

import numpy as np
import pytest

from spectral_eqe import (
    DispersiveStack, Layer, OpticalDataLake, OpticalRateBudget,
    dipole_ensemble_eqe, random_dipoles,
)


TABLE = b"""REFERENCES: Synthetic test data, not a measured material
COMMENTS: Unit test only
DATA:
  - type: tabulated nk
    data: |
      0.4 1.5 0.01
      0.6 1.7 0.03
"""


def import_table(lake, name="film", raw=TABLE):
    return lake.import_rii(name, raw, source="synthetic fixture", revision="test-v1")


def constant(lake, name, n, k=0):
    raw = f"REFERENCES: Synthetic test\nDATA:\n  - type: tabulated nk\n    data: |\n      0.4 {n} {k}\n      0.6 {n} {k}\n".encode()
    return import_table(lake, name, raw)


def test_catalog_csv_interpolation_and_integrity(tmp_path):
    source = tmp_path / "input.csv"
    source.write_text("wavelength_nm,n,k\n400,1.5,0.01\n600,1.7,0.03\n")
    lake = OpticalDataLake(tmp_path / "lake")
    material = lake.import_csv("film", source, source="measurement.csv", license="test-only",
                               citation="Synthetic", conditions="test")
    np.testing.assert_allclose(material.index([400, 500, 600]), [1.5 + .01j, 1.6 + .02j, 1.7 + .03j])
    assert lake.list_datasets() == ["film"]
    assert lake.load("film").metadata["conditions"] == "test"
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    (lake.root / f"{digest}.source").write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="Raw material checksum"):
        lake.load("film")


def test_catalog_record_checksum_and_immutable_ids(tmp_path):
    lake = OpticalDataLake(tmp_path)
    import_table(lake)
    import_table(lake)
    with pytest.raises(FileExistsError):
        constant(lake, "film", 2)
    record = json.loads((tmp_path / "film.json").read_text())
    record["payload"]["metadata"]["citation"] = "tampered"
    (tmp_path / "film.json").write_text(json.dumps(record))
    with pytest.raises(ValueError, match="Catalog checksum"):
        lake.load("film")


@pytest.mark.parametrize("kind,pole", [("formula 1", ".1"), ("formula 2", ".01")])
def test_rii_native_tables_and_formula(tmp_path, kind, pole):
    lake = OpticalDataLake(tmp_path)
    formula = f"REFERENCES: Synthetic formula\nDATA:\n - type: {kind}\n   wavelength_range: 0.4 0.6\n   coefficients: 0 1 {pole}\n".encode()
    material = lake.import_rii("formula", formula, source="synthetic", revision="v1", assume_lossless=True)
    expected = np.sqrt(1 + .5**2 / (.5**2 - .1**2))
    assert material.index([500])[0] == pytest.approx(expected + 0j)
    assert material.metadata["license"] == "CC0-1.0"
    assert material.assume_lossless
    nk = import_table(lake)
    assert nk.index([500])[0] == pytest.approx(1.6 + .02j)


def test_separate_n_and_k_tables_use_intersection(tmp_path):
    raw = b"""REFERENCES: Synthetic separate grids
DATA:
 - type: tabulated n
   data: |
     0.4 1.5
     0.6 1.7
 - type: tabulated k
   data: |
     0.45 0.01
     0.55 0.03
"""
    material = import_table(OpticalDataLake(tmp_path), raw=raw)
    assert material.index([500])[0] == pytest.approx(1.6 + .02j)
    with pytest.raises(ValueError, match="domain"):
        material.index([400])


def test_missing_k_requires_explicit_lossless_assumption(tmp_path):
    raw = b"REFERENCES: Synthetic\nDATA:\n - type: tabulated n\n   data: |\n     .4 1.5\n     .6 1.5\n"
    material = import_table(OpticalDataLake(tmp_path), raw=raw)
    with pytest.raises(ValueError, match="Extinction coefficient"):
        material.index([500])


def test_fetch_is_pinned_and_offline(tmp_path):
    lake = OpticalDataLake(tmp_path)
    commit = "a" * 40
    path = "database/data/main/Ag/nk/Johnson.yml"
    with patch("spectral_eqe.urlopen", return_value=io.BytesIO(TABLE)) as remote:
        material = lake.fetch_rii("ag", path, commit)
        assert remote.call_args.args[0] == f"https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/{commit}/{path}"
        assert material.metadata["revision"] == commit
    with patch("spectral_eqe.urlopen", side_effect=AssertionError("must be offline")):
        assert lake.fetch_rii("ag", path, commit).index([500])[0] == pytest.approx(1.6 + .02j)
    with pytest.raises(FileExistsError):
        lake.fetch_rii("ag", path, "b" * 40)


@pytest.mark.parametrize("path,commit", [
    ("database/data/main/Ag/nk/Johnson.yml", "main"),
    ("database/data/../private.yml", "a" * 40),
    ("https://example.com/file.yml", "a" * 40),
    ("database/data/a.yml?x=1", "a" * 40),
    ("database\\data\\a.yml", "a" * 40),
])
def test_fetch_rejects_unpinned_and_noncanonical_paths(tmp_path, path, commit):
    with patch("spectral_eqe.urlopen", side_effect=AssertionError("no network")):
        with pytest.raises(ValueError):
            OpticalDataLake(tmp_path).fetch_rii("sample", path, commit)


def test_fetch_network_errors_are_not_hidden(tmp_path):
    with patch("spectral_eqe.urlopen", side_effect=OSError("offline")):
        with pytest.raises(OSError, match="offline"):
            OpticalDataLake(tmp_path).fetch_rii("sample", "database/data/a.yml", "a" * 40)


@pytest.mark.parametrize("size", [4_000_000, 4_000_001])
def test_fetch_dataset_size_boundary(tmp_path, size):
    raw = TABLE + b"#" + b" " * (size - len(TABLE) - 1)
    lake = OpticalDataLake(tmp_path)
    with patch("spectral_eqe.urlopen", return_value=io.BytesIO(raw)):
        if size > 4_000_000:
            with pytest.raises(ValueError, match="4 MB"):
                lake.fetch_rii("sample", "database/data/a.yml", "a" * 40)
            assert lake.list_datasets() == []
        else:
            material = lake.fetch_rii("sample", "database/data/a.yml", "a" * 40)
            assert material.index([500])[0] == pytest.approx(1.6 + .02j)


def test_catalog_rejects_path_ids_and_missing_csv_provenance(tmp_path):
    lake = OpticalDataLake(tmp_path / "lake")
    with pytest.raises(ValueError, match="identifier"):
        lake.load("../outside")
    source = tmp_path / "data.csv"
    source.write_text("wavelength_nm,n,k\n400,1,0\n600,1,0\n")
    with pytest.raises(ValueError, match="license"):
        lake.import_csv("sample", source, source="source", license="", citation="Test", conditions="Test")
    assert lake.list_datasets() == []
    source.write_text("wavelength_um,n,k\n.4,1,0\n.6,1,0\n")
    with pytest.raises(ValueError, match="columns"):
        lake.import_csv("sample", source, source="source", license="test", citation="Test", conditions="Test")


def test_formula_domain_poles_and_malformed_coefficients_are_rejected(tmp_path):
    lake = OpticalDataLake(tmp_path)
    raw = b"REFERENCES: Test\nDATA:\n - type: formula 1\n   wavelength_range: 0.4 0.6\n   coefficients: 0 1 0.5\n"
    material = lake.import_rii("pole", raw, source="fixture", revision="1", assume_lossless=True)
    with pytest.raises(ValueError, match="pole"):
        material.index([500])
    with pytest.raises(ValueError, match="domain"):
        material.index([700])
    with pytest.raises(ValueError):
        lake.import_rii("malformed", raw.replace(b"0 1 0.5", b"0 1 nonsense"), source="fixture", revision="1")


def test_cached_lossless_assumption_cannot_change(tmp_path):
    lake = OpticalDataLake(tmp_path)
    with patch("spectral_eqe.urlopen", return_value=io.BytesIO(TABLE)):
        lake.fetch_rii("sample", "database/data/a.yml", "a" * 40)
    with pytest.raises(FileExistsError, match="assumption"):
        lake.fetch_rii("sample", "database/data/a.yml", "a" * 40, assume_lossless=True)


@pytest.mark.parametrize("raw", [
    b"REFERENCES: Test\nDATA: []",
    b"REFERENCES: Test\nDATA:\n - type: formula 9\n",
    b"DATA:\n - type: tabulated nk\n   data: |\n     .4 1.5 0\n     .6 1.5 0\n",
    b"REFERENCES: Test\nDATA:\n - type: tabulated nk\n   data: |\n     .4 1.5 -1\n     .6 1.5 0\n",
    b"REFERENCES: Test\nDATA:\n - type: tabulated nk\n   data: |\n     .6 1.5 0\n     .4 1.5 0\n",
])
def test_invalid_imports_are_explicit(tmp_path, raw):
    with pytest.raises(ValueError):
        import_table(OpticalDataLake(tmp_path), raw=raw)


@pytest.mark.parametrize("grid", [[399], [601], [500, 400], [500, 500], [], [np.nan], [[500]]])
def test_material_rejects_invalid_or_out_of_domain_wavelengths(tmp_path, grid):
    material = import_table(OpticalDataLake(tmp_path))
    with pytest.raises(ValueError):
        material.index(grid)


def make_stack(lake, thickness=80):
    constant(lake, "air", 1)
    constant(lake, "emitter", 1.7)
    import_table(lake)
    return DispersiveStack((Layer("air", None), Layer("emitter", 40),
                            Layer("film", thickness), Layer("air", None)), 1, 20)


def test_stack_dispersion_and_thickness(tmp_path):
    lake = OpticalDataLake(tmp_path)
    stack = make_stack(lake)
    data = stack.evaluate(lake, [400, 500, 600])
    assert data["indices"].shape == (3, 4)
    np.testing.assert_allclose(data["indices"][:, 2], [1.5 + .01j, 1.6 + .02j, 1.7 + .03j])
    np.testing.assert_equal(data["thicknesses_nm"], [np.inf, 40, 80, np.inf])
    assert data["emitter_depth_nm"] == 20
    assert data["materials"][2]["dataset_id"] == "film"


def test_plane_wave_thickness_spectrum_is_not_eqe(tmp_path):
    lake = OpticalDataLake(tmp_path)
    thin = make_stack(lake, 10).plane_wave_spectrum(lake, [400, 500, 600])
    thick = make_stack(lake, 100).plane_wave_spectrum(lake, [400, 500, 600])
    np.testing.assert_allclose(thin["R"] + thin["T"] + thin["A"], 1, atol=1e-12)
    assert thin["quantity"] == "plane_wave_power_not_EQE"
    assert np.max(np.abs(thin["T"] - thick["T"])) > .01
    assert np.all(thin["A"] >= 0)
    assert "EQE" not in thin


@pytest.mark.parametrize("layers,emitter,depth", [
    ((Layer("a", None), Layer("b", None)), 1, 1),
    ((Layer("a", 1), Layer("b", 10), Layer("c", None)), 1, 1),
    ((Layer("a", None), Layer("b", -1), Layer("c", None)), 1, 1),
    ((Layer("a", None), Layer("b", 10), Layer("c", None)), 0, 1),
    ((Layer("a", None), Layer("b", 10), Layer("c", None)), 1, 0),
    ((Layer("a", None), Layer("b", 10), Layer("c", None)), 1, 10),
    ((Layer("a", None), Layer("b", 10), Layer("c", None)), 1, np.nan),
])
def test_invalid_stack_geometry(layers, emitter, depth):
    with pytest.raises(ValueError):
        DispersiveStack(layers, emitter, depth)


def test_absorbing_host_and_halfspace_are_rejected(tmp_path):
    lake = OpticalDataLake(tmp_path)
    constant(lake, "air", 1)
    constant(lake, "lossy", 1.7, .1)
    for layers in [
        (Layer("air", None), Layer("lossy", 10), Layer("air", None)),
        (Layer("lossy", None), Layer("air", 10), Layer("air", None)),
    ]:
        with pytest.raises(ValueError, match="lossless"):
            DispersiveStack(layers, 1, 5).evaluate(lake, [500])


def test_monochromatic_rate_equation_and_yield_boundaries():
    budget = OpticalRateBudget([400, 500], [1, 2], [.5, .4], "P_host")
    kwargs = dict(charge_balance=.8, exciton_yield=.25)
    np.testing.assert_allclose(budget.monochromatic_eqe(intrinsic_yield=.5, **kwargs), [.05, .2 * .2 / 1.5])
    np.testing.assert_array_equal(budget.monochromatic_eqe(intrinsic_yield=0, **kwargs), [0, 0])
    np.testing.assert_allclose(budget.monochromatic_eqe(intrinsic_yield=1, **kwargs), [.1, .04])
    np.testing.assert_array_equal(
        budget.monochromatic_eqe(intrinsic_yield=1, charge_balance=0, exciton_yield=1), [0, 0])


def test_broadband_uses_global_rate_denominator():
    budget = OpticalRateBudget([400, 500, 600], [1, 2, 3], [.5, 1, 1.5], "P_host")
    result = budget.broadband_eqe([.005, .005, .005], intrinsic_yield=.5, charge_balance=.8, exciton_yield=.25)
    assert result["rate_denominator"] == pytest.approx(1.5)
    np.testing.assert_allclose(result["eqe_density_per_nm"], .2 * .5 * .005 * np.array([.5, 1, 1.5]) / 1.5)
    assert result["integrated_eqe"] == pytest.approx(.2 * .5 / 1.5)
    zero = budget.broadband_eqe([.005] * 3, intrinsic_yield=0, charge_balance=1, exciton_yield=1)
    assert zero["integrated_eqe"] == 0
    np.testing.assert_array_equal(zero["eqe_density_per_nm"], [0, 0, 0])
    perfect = budget.broadband_eqe([.005] * 3, intrinsic_yield=1, charge_balance=1, exciton_yield=1)
    assert perfect["integrated_eqe"] == pytest.approx(.5)


def test_orientation_averages_powers_before_ratios():
    vertical = OpticalRateBudget([500], [4], [1], "P_host")
    horizontal = OpticalRateBudget([500], [1], [.5], "P_host")
    average = OpticalRateBudget.orientation_average(vertical, horizontal, 1 / 3)
    assert average.total[0] == pytest.approx(2)
    assert average.escape[0] == pytest.approx(2 / 3)
    assert average.monochromatic_eqe(intrinsic_yield=1, charge_balance=1, exciton_yield=1)[0] == pytest.approx(1 / 3)
    compare = dipole_ensemble_eqe(vertical, horizontal, [[0, 0, 1], [1, 0, 0]],
                                intrinsic_yield=1, charge_balance=1, exciton_yield=1)
    assert compare["independent_fixed_dipoles_eqe"][0] == pytest.approx((.25 + .5) / 2)
    assert compare["rapid_reorientation_eqe"][0] == pytest.approx(.75 / 2.5)
    assert compare["mean_vertical_fraction"] == .5


def test_random_orientations_are_reproducible_and_isotropic():
    a = random_dipoles(10000, seed=42)
    np.testing.assert_array_equal(a, random_dipoles(10000, seed=42))
    np.testing.assert_allclose(np.linalg.norm(a, axis=1), 1, atol=1e-14)
    np.testing.assert_allclose(np.mean(a**2, axis=0), 1 / 3, atol=.01)
    assert not np.array_equal(a, random_dipoles(10000, seed=43))


@pytest.mark.parametrize("total,escape", [([-1], [0]), ([1], [-1]), ([1], [2]), ([np.nan], [0]), ([1, 2], [.5])])
def test_invalid_optical_budgets(total, escape):
    with pytest.raises(ValueError):
        OpticalRateBudget([500], total, escape, "P_host")


@pytest.mark.parametrize("yield_value", [-1, 1.1, np.nan, np.inf])
def test_invalid_eqe_factors(yield_value):
    budget = OpticalRateBudget([500], [1], [.5], "P_host")
    for parameter in ["intrinsic_yield", "charge_balance", "exciton_yield"]:
        kwargs = dict(intrinsic_yield=1, charge_balance=1, exciton_yield=1)
        kwargs[parameter] = yield_value
        with pytest.raises(ValueError):
            budget.monochromatic_eqe(**kwargs)


def test_dark_emitter_and_invalid_spectrum_fail_explicitly():
    budget = OpticalRateBudget([400, 600], [0, 0], [0, 0], "P_host")
    with pytest.raises(ValueError, match="Zero total"):
        budget.monochromatic_eqe(intrinsic_yield=1, charge_balance=1, exciton_yield=1)
    for spectrum in [[1, 1], [-.005, .015], [np.nan, .005], [.005]]:
        with pytest.raises(ValueError):
            budget.broadband_eqe(spectrum, intrinsic_yield=.5, charge_balance=1, exciton_yield=1)
    with pytest.raises(ValueError, match="Zero integrated"):
        budget.broadband_eqe([.005, .005], intrinsic_yield=1, charge_balance=1, exciton_yield=1)


def test_spectral_overflow_is_not_reported_as_zero_efficiency():
    # A narrow normalized peak exposes intermediate spectrum-times-rate overflow.
    budget = OpticalRateBudget([500, 500.0000001], [1.7e308] * 2, [1.7e308] * 2, "P_host")
    density = 1 / (500.0000001 - 500)
    with pytest.raises(ValueError, match="overflowed"):
        budget.broadband_eqe([density] * 2, intrinsic_yield=1, charge_balance=1, exciton_yield=1)


def test_budget_copies_inputs_and_requires_matching_references():
    total = np.array([1.])
    budget = OpticalRateBudget([500], total, [.5], "host A")
    total[0] = 2
    assert budget.total[0] == 1
    with pytest.raises(ValueError):
        budget.total[0] = 2
    for other in [OpticalRateBudget([500], [1], [.5], "host B"),
                  OpticalRateBudget([600], [1], [.5], "host A")]:
        with pytest.raises(ValueError):
            OpticalRateBudget.orientation_average(budget, other, .5)


@pytest.mark.parametrize("directions", [[], [[0, 0, 0]], [[np.nan, 0, 1]], [[1, 0]]])
def test_invalid_ensemble_directions(directions):
    budget = OpticalRateBudget([500], [1], [.5], "P_host")
    with pytest.raises(ValueError):
        dipole_ensemble_eqe(budget, budget, directions, intrinsic_yield=1, charge_balance=1, exciton_yield=1)


@pytest.mark.parametrize("count", [0, -1, True, 1.5])
def test_random_orientation_count_validation(count):
    with pytest.raises(ValueError):
        random_dipoles(count, seed=42)
