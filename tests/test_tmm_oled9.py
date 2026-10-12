import importlib.util
from pathlib import Path
import numpy as np
import pytest
from tmm import coh_tmm

spec = importlib.util.spec_from_file_location("oled9", Path(__file__).resolve().parents[1] / "tmm_OLED9.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_stack_has_matching_lengths():
    n, d = module.stack_at(550.0, 75.0)
    assert len(n) == len(d) == 4
    assert np.isinf(d[0]) and np.isinf(d[-1])


def test_transmission_matches_reference_tmm():
    n, d = module.stack_at(550.0, 75.0)
    expected = coh_tmm("s", n, d, np.deg2rad(30.0), 550.0)["T"]
    assert module.transmission(550.0, 30.0, 75.0, "s") == pytest.approx(expected)


def test_sweep_shape_and_finite_values():
    result = module.sweep([450.0, 550.0], [0.0, 30.0], 75.0)
    assert result.shape == (2, 2)
    assert np.isfinite(result).all()
    assert np.all(result >= 0)


@pytest.mark.parametrize("wavelength,angle,thickness", [(0, 0, 50), (550, 90, 50), (550, 0, -1)])
def test_invalid_inputs(wavelength, angle, thickness):
    with pytest.raises(ValueError):
        module.transmission(wavelength, angle, thickness, "s")


def test_no_plot_on_import():
    assert not hasattr(module, "plt")
