import numpy as np
import pytest
from tmm import coh_tmm

from dipole_reference import DipoleConvergenceError, SingleInterface
from multilayer_dipole import EmitterStack, _kz, _outward, emitter_rate_budget
from spectral_eqe import DispersiveStack, Layer, OpticalDataLake

LAM = 550.0
FILM = b"""REFERENCES: Synthetic test data
COMMENTS: Unit test only
DATA:
  - type: tabulated nk
    data: |
      0.4 {n} {k}
      0.8 {n} {k}
"""


def lossless_cavity(**kwargs):
    base = dict(wavelength_nm=LAM, indices=(1.8, 1.2, 1.0, 1.1, 1.5),
                thicknesses_nm=(None, 90.0, 120.0, 70.0, None), emitter_layer=2, emitter_depth_nm=40.0)
    base.update(kwargs)
    return EmitterStack(**base)


@pytest.mark.parametrize("epsilon", [-5 + 0.5j, 2.25 + 0j, -12 + 1.2j])
def test_matched_top_reduces_to_single_interface_kernel(epsilon):
    index = np.sqrt(epsilon)
    index = index if index.imag >= 0 else -index
    result = EmitterStack(LAM, (1.0, 1.0, index), (None, 300.0, None), 1, 150.0).rates()
    reference = SingleInterface(LAM, 150.0, epsilon, 1.0).decay_rates()
    assert result.total == pytest.approx((reference.vertical, reference.parallel), rel=1e-6)


@pytest.mark.parametrize("indices,thicknesses,layer,depth", [
    ((1.8, 1.0, 1.5), (None, 200.0, None), 1, 80.0),
    ((1.8, 1.2, 1.0, 1.1, 1.5), (None, 90.0, 120.0, 70.0, None), 2, 40.0),
    ((1.0, 1.4, 1.0, 1.6), (None, 60.0, 250.0, None), 2, 30.0),
])
def test_lossless_source_work_equals_exit_flux(indices, thicknesses, layer, depth):
    # Independent energy balance: source work (field at the dipole) vs radiated flux (far-field amplitudes).
    result = EmitterStack(LAM, indices, thicknesses, layer, depth).rates()
    for total, top, bottom in zip(result.total, result.top, result.bottom):
        assert top + bottom == pytest.approx(total, abs=1e-6)
    assert result.other == pytest.approx((0.0, 0.0), abs=1e-6)


def test_flipping_stack_swaps_exit_channels():
    stack = lossless_cavity()
    flipped = lossless_cavity(indices=(1.5, 1.1, 1.0, 1.2, 1.8), thicknesses_nm=(None, 70.0, 120.0, 90.0, None),
                              emitter_depth_nm=80.0)
    a, b = stack.rates(), flipped.rates()
    assert a.top == pytest.approx(b.bottom, rel=1e-6)
    assert a.bottom == pytest.approx(b.top, rel=1e-6)
    assert a.total == pytest.approx(b.total, rel=1e-6)


def test_homogeneous_limit_has_unit_rates_and_half_space_split():
    result = EmitterStack(LAM, (1.0, 1.0, 1.0), (None, 400.0, None), 1, 200.0).rates()
    assert result.total == pytest.approx((1.0, 1.0), abs=1e-9)
    assert result.top == pytest.approx((0.5, 0.5), abs=1e-9)
    assert result.bottom == pytest.approx((0.5, 0.5), abs=1e-9)


def test_recursion_matches_tmm_package_power_coefficients():
    # Plane-wave R and T from the recursion versus the independent `tmm` package at propagating kappa.
    eps = [1.8**2, (1.5 + 0.1j) ** 2, 1.75**2, 1.0]
    d = [0.0, 90.0, 120.0, 0.0]
    n_list = [np.sqrt(e) for e in eps]
    angle = np.deg2rad(20.0)
    kappa = np.array([1.75 * np.sin(angle)])
    k0 = 2 * np.pi / LAM
    kz = [_kz(e, kappa) for e in eps]
    for pol in ("s", "p"):
        # emitter layer index 2; outward list runs from exit (layer 0) to the emitter layer
        r, t = _outward(pol, eps[:3], kz[:3], d[:3], k0)
        reference = coh_tmm(pol, [n_list[2], n_list[1], n_list[0]][::-1][::-1] and
                            [1.75, 1.5 + 0.1j, 1.8], [np.inf, 90.0, np.inf], angle, LAM)
        assert abs(r[0]) ** 2 == pytest.approx(reference["R"], rel=1e-9)
        y = (lambda k, e: np.real(k / e) if pol == "p" else np.real(k))
        transmitted = y(kz[0], eps[0])[0] * abs(t[0]) ** 2 / y(kz[2], eps[2])[0]
        assert transmitted == pytest.approx(reference["T"], rel=1e-9)


def test_lossy_stack_balances_into_other_and_orientation_ordering():
    cathode = np.sqrt(-30 + 10j)
    cathode = cathode if cathode.imag > 0 else -cathode
    result = EmitterStack(LAM, (1.5, 1.9 + 0.01j, 1.75, cathode), (None, 100.0, 150.0, None), 2, 75.0).rates()
    assert result.other[0] > 0.5 and result.other[1] > 0.1
    assert result.bottom == (0.0, 0.0)
    assert result.top[0] < result.top[1]  # vertical dipoles couple poorly to the escape cone
    assert result.oriented("total", 0.0) == result.total[1]
    assert result.oriented("total", 1 / 3) == pytest.approx(result.total[0] / 3 + 2 * result.total[1] / 3)


def test_weak_loss_limit_is_continuous_for_guided_power():
    guide = lambda k: EmitterStack(LAM, (1.0, 2.0 + 1j * k, 1.7, 1.5), (None, 300.0, 100.0, None), 2, 50.0).rates()
    strong, weak = guide(3e-3), guide(3e-4)
    assert strong.other == pytest.approx(weak.other, rel=0.05)
    assert strong.other[0] > 0.5


def test_lossless_guided_layer_is_refused():
    with pytest.raises(DipoleConvergenceError, match="guid"):
        EmitterStack(LAM, (1.0, 2.0, 1.0, 1.5), (None, 300.0, 100.0, None), 2, 50.0).rates()


def test_air_first_pass_is_lower_than_exit_flux():
    result = EmitterStack(LAM, (1.5, 1.0, 1.0), (None, 400.0, None), 1, 200.0).rates(air_index=1.0)
    assert result.air_top[0] <= result.top[0] and result.air_top[1] <= result.top[1]
    assert result.air_bottom == pytest.approx(result.bottom)


@pytest.mark.parametrize("kwargs,message", [
    (dict(wavelength_nm=0.0), "wavelength"),
    (dict(indices=(1.0, 1.0)), "half-spaces"),
    (dict(emitter_layer=0), "interior"),
    (dict(emitter_depth_nm=120.0), "depth"),
    (dict(indices=(1.8, 1.2, 1.0 + 0.1j, 1.1, 1.5)), "lossless"),
    (dict(indices=(1.8, 1.2, 1.0, -1.1, 1.5)), "passive"),
    (dict(thicknesses_nm=(None, 90.0, -1.0, 70.0, None)), "positive"),
])
def test_invalid_inputs_are_rejected(kwargs, message):
    with pytest.raises(ValueError, match=message):
        lossless_cavity(**kwargs)


def test_invalid_tolerances_are_rejected():
    with pytest.raises(ValueError):
        lossless_cavity().rates(atol=1.0)
    with pytest.raises(ValueError):
        lossless_cavity().rates(air_index=0)
    with pytest.raises(ValueError):
        lossless_cavity().rates(max_intervals=5)


def test_orientation_helper_validates_inputs():
    result = lossless_cavity().rates()
    with pytest.raises(ValueError):
        result.oriented("nonsense", 0.5)
    with pytest.raises(ValueError):
        result.oriented("total", 1.5)


def test_dispersive_stack_budget_feeds_eqe_without_exceeding_total(tmp_path):
    lake = OpticalDataLake(tmp_path)
    for name, n, k in (("glass", 1.5, 0), ("organic", 1.75, 0), ("ito", 1.9, 0.02)):
        raw = FILM.replace(b"{n}", str(n).encode()).replace(b"{k}", str(k).encode())
        lake.import_rii(name, raw, source="fixture", revision="v1")
    stack = DispersiveStack((Layer("glass", None), Layer("organic", 150.0), Layer("ito", 100.0), Layer("glass", None)),
                            emitter_layer=1, emitter_depth_nm=75.0)
    wavelengths = [500.0, 550.0, 600.0]
    budget = emitter_rate_budget(stack, lake, wavelengths, exit="bottom", collection="exit_medium",
                                 vertical_fraction=0.2)
    assert np.all(budget.escape <= budget.total + 1e-9)
    assert np.all(budget.escape > 0)
    eqe = budget.monochromatic_eqe(intrinsic_yield=0.9, charge_balance=1.0, exciton_yield=1.0)
    assert np.all((eqe > 0) & (eqe < 1))
    with pytest.raises(ValueError):
        emitter_rate_budget(stack, lake, wavelengths, exit="side", collection="exit_medium", vertical_fraction=0.2)
