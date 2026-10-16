import numpy as np
import pytest
from scipy.constants import c, epsilon_0
from scipy.integrate import quad

from dipole_reference import (
    DipoleConvergenceError,
    SingleInterface,
    homogeneous_pattern,
    homogeneous_power,
)


@pytest.mark.parametrize("orientation", [
    "vertical", "horizontal_x", "horizontal_y", "horizontal", "isotropic",
])
def test_homogeneous_sphere_and_hemisphere_normalization(orientation):
    nodes, weights = np.polynomial.legendre.leggauss(40)
    phi = np.linspace(0, 2 * np.pi, 80, endpoint=False)
    for lo, hi, expected in [(-1, 1, 1), (0, 1, 0.5), (-1, 0, 0.5)]:
        cosine = lo + (nodes + 1) * (hi - lo) / 2
        density = homogeneous_pattern(np.arccos(cosine[:, None]), phi, orientation)
        integral = np.sum(weights[:, None] * density.total) * (hi - lo) / 2 * 2 * np.pi / len(phi)
        assert integral == pytest.approx(expected, abs=1e-12)
        assert np.all(density.s >= 0) and np.all(density.p >= 0)


def test_homogeneous_directional_patterns_and_polarization():
    normalization = 3 / (8 * np.pi)
    vertical = homogeneous_pattern([0, np.pi / 2, np.pi], 0, "vertical")
    np.testing.assert_allclose(vertical.total, [0, normalization, 0], atol=1e-15)
    np.testing.assert_array_equal(vertical.s, [0, 0, 0])
    x = homogeneous_pattern(np.pi / 2, [0, np.pi / 2], "horizontal_x")
    y = homogeneous_pattern(np.pi / 2, [0, np.pi / 2], "horizontal_y")
    np.testing.assert_allclose(x.total, [0, normalization], atol=1e-15)
    np.testing.assert_allclose(y.total, [normalization, 0], atol=1e-15)
    np.testing.assert_allclose(x.p, [0, 0], atol=1e-15)
    isotropic = homogeneous_pattern([0, 0.7, np.pi], [1, 2, 3], "isotropic")
    np.testing.assert_allclose(isotropic.s, 1 / (8 * np.pi), atol=1e-15)
    np.testing.assert_allclose(isotropic.p, 1 / (8 * np.pi), atol=1e-15)


def test_isotropic_pattern_is_orientation_power_average():
    theta, phi = np.linspace(0, np.pi, 17)[:, None], np.linspace(0, 2 * np.pi, 19)
    vertical = homogeneous_pattern(theta, phi, "vertical")
    horizontal = homogeneous_pattern(theta, phi, "horizontal")
    x = homogeneous_pattern(theta, phi, "horizontal_x")
    y = homogeneous_pattern(theta, phi, "horizontal_y")
    isotropic = homogeneous_pattern(theta, phi, "isotropic")
    for part in ("s", "p"):
        np.testing.assert_allclose(getattr(horizontal, part), (getattr(x, part) + getattr(y, part)) / 2, atol=1e-15)
        np.testing.assert_allclose(
            getattr(isotropic, part), (getattr(vertical, part) + 2 * getattr(horizontal, part)) / 3,
            atol=1e-15,
        )


def test_homogeneous_power_si_normalization_and_scaling():
    wavelength_m, amplitude = 550e-9, 1e-29
    # Equivalent k0^4*c expression checks wavelength conversion and harmonic amplitude.
    expected = (2 * np.pi / wavelength_m) ** 4 * c * amplitude**2 / (12 * np.pi * epsilon_0)
    vacuum = homogeneous_power(550, 1, amplitude)
    assert vacuum == pytest.approx(expected, rel=1e-13)
    assert homogeneous_power(550, 1.5, amplitude) == pytest.approx(1.5 * vacuum)
    assert homogeneous_power(1100, 1, 2j * amplitude) == pytest.approx(vacuum / 4)
    assert homogeneous_power(550, 1, 0) == 0


@pytest.mark.parametrize("host_index", [1.0, 1.7])
def test_matched_interface_has_unit_rate_and_zero_spectrum(host_index):
    interface = SingleInterface(550, 37, host_index**2, host_index)
    for q in [0, 0.5, 0.999, 1.001, 2, 10]:
        assert interface.spectral_correction(q) == (0, 0)
    result = interface.decay_rates()
    assert result.vertical == result.parallel == result.isotropic == 1
    assert result.propagating_correction == result.evanescent_correction == (0, 0)
    assert result.estimated_error == result.tail_magnitude == (0, 0)


@pytest.mark.parametrize("a", [0.001, 0.1, 1.0, 4.0, 20.0, 100.0])
def test_perfect_mirror_rates_match_image_dipole(a):
    interface = SingleInterface(600, a * 600 / (4 * np.pi * 1.5), "pec", 1.5)
    result = interface.decay_rates(atol=1e-10, rtol=1e-9)
    if a < 0.01:
        expected_vertical = 2 - a**2 / 10 + a**4 / 280
        expected_parallel = a**2 / 5 - 3 * a**4 / 280
    else:
        expected_vertical = 1 + 3 * (np.sin(a) - a * np.cos(a)) / a**3
        expected_parallel = 1 - 1.5 * (np.sin(a) / a + np.cos(a) / a**2 - np.sin(a) / a**3)
    assert result.vertical == pytest.approx(expected_vertical, abs=1e-10)
    assert result.parallel == pytest.approx(expected_parallel, abs=1e-10)
    assert result.evanescent_correction == (0, 0)
    assert result.evaluations > 0


def test_mirror_raw_spectrum_has_correct_parallel_sign():
    interface = SingleInterface(600, 50, "pec")
    q = 0.6
    w = np.sqrt(1 - q**2)
    phase_real = np.cos(4 * np.pi * 50 / 600 * w)
    vertical, parallel = interface.spectral_correction(q)
    assert vertical == pytest.approx(1.5 * q**3 / w * phase_real)
    assert parallel == pytest.approx(-0.75 * q / w * (1 + w**2) * phase_real)


def _alternative_angular_integral(epsilon, distance, wavelength, host):
    """Different theta/cosh coordinates and absolute wavevectors for reference."""
    k0, k1 = 2 * np.pi / wavelength, 2 * np.pi * host / wavelength

    def reflection(q):
        kz1 = k1 * np.sqrt(1 - q * q + 0j)
        kz2 = np.sqrt(k0 * k0 * epsilon - (q * k1) ** 2 + 0j)
        rs = (kz1 - kz2) / (kz1 + kz2)
        rp = (epsilon * kz1 - host**2 * kz2) / (epsilon * kz1 + host**2 * kz2)
        return rs, rp

    values = []
    for vertical in [True, False]:
        def angular(theta):
            rs, rp = reflection(np.sin(theta))
            cosine = np.cos(theta)
            amplitude = np.sin(theta) ** 3 * rp if vertical else np.sin(theta) * (rs - cosine**2 * rp)
            return np.real(amplitude * np.exp(2j * k1 * distance * cosine))

        def hyperbolic(t):
            rs, rp = reflection(np.cosh(t))
            sinh = np.sinh(t)
            amplitude = np.cosh(t) ** 3 * rp if vertical else np.cosh(t) * (rs + sinh**2 * rp)
            return np.imag(amplitude) * np.exp(-2 * k1 * distance * sinh)

        prop_points, evan_points = [], []
        if complex(epsilon).imag == 0:
            ratio = float(np.real(epsilon)) / host**2
            if ratio < 1:
                prop_points = [np.arcsin(np.sqrt(ratio))]
            elif ratio > 1:
                evan_points = [np.arccosh(np.sqrt(ratio))]
        prop = quad(angular, 0, np.pi / 2, points=prop_points, epsabs=1e-11, epsrel=1e-11, limit=300)[0]
        upper = np.arcsinh(80 / (2 * k1 * distance))
        evan = quad(hyperbolic, 0, upper, points=[p for p in evan_points if p < upper],
                    epsabs=1e-11, epsrel=1e-11, limit=300)[0]
        values.append(1 + (1.5 if vertical else 0.75) * (prop + evan))
    return values


@pytest.mark.parametrize("epsilon,host", [
    (2.25, 1.0), (1.0, 1.5), (2.25 + 0.4j, 1.0), (-5 + 0.5j, 1.0),
])
def test_interface_matches_alternative_angular_integration(epsilon, host):
    result = SingleInterface(550, 45, epsilon, host).decay_rates(atol=1e-10, rtol=1e-9)
    expected = _alternative_angular_integral(epsilon, 45, 550, host)
    np.testing.assert_allclose([result.vertical, result.parallel], expected, rtol=1e-8, atol=1e-9)


def test_high_index_lossless_substrate_has_host_evanescent_source_work():
    result = SingleInterface(550, 40, 2.25).decay_rates()
    assert result.evanescent_correction[0] > 0.01
    assert result.evanescent_correction[1] > 0.01
    # This contribution is transmission into a dielectric, not substrate absorption.
    assert result.vertical > 0 and result.parallel > 0


def test_lossy_quasistatic_limit_and_orientation_ratio():
    epsilon = -5 + 0.5j
    interface = SingleInterface(550, 0.05, epsilon)
    result = interface.decay_rates(atol=1e-8, rtol=1e-7)
    a = 4 * np.pi * 0.05 / 550
    reflection_imag = np.imag((epsilon - 1) / (epsilon + 1))
    assert (result.vertical - 1) * a**3 / 3 == pytest.approx(reflection_imag, rel=2e-5)
    assert (result.parallel - 1) * a**3 / 1.5 == pytest.approx(reflection_imag, rel=2e-5)


def test_refinement_and_diagnostics_for_lossy_interface():
    interface = SingleInterface(550, 12, -3 + 0.2j)
    loose = interface.decay_rates()
    tight = interface.decay_rates(atol=1e-10, rtol=1e-9)
    np.testing.assert_allclose([loose.vertical, loose.parallel], [tight.vertical, tight.parallel], rtol=1e-7)
    for i, rate in enumerate([tight.vertical, tight.parallel]):
        assert tight.estimated_error[i] <= max(1e-10, 1e-9 * rate)
        assert tight.refinement_delta[i] <= max(1e-10, 1e-9 * rate)
        assert 0 <= tight.tail_magnitude[i] < 1e-12
        assert rate == pytest.approx(1 + tight.propagating_correction[i] + tight.evanescent_correction[i])
    assert tight.oriented(0) == tight.parallel
    assert tight.oriented(1) == tight.vertical
    assert tight.isotropic == pytest.approx((tight.vertical + 2 * tight.parallel) / 3)


def test_dimensionless_interface_scaling_preserves_source_work():
    baseline = SingleInterface(550, 20, -5 + 0.5j).decay_rates()
    scaled_lengths = SingleInterface(1100, 40, -5 + 0.5j).decay_rates()
    scaled_medium = SingleInterface(1100, 20, 4 * (-5 + 0.5j), 2).decay_rates()
    for other in [scaled_lengths, scaled_medium]:
        np.testing.assert_allclose(
            [other.vertical, other.parallel], [baseline.vertical, baseline.parallel], rtol=1e-12,
        )


@pytest.mark.parametrize("epsilon,host", [(2.25, 1.0), (1.0, 1.5)])
def test_weak_dielectric_loss_approaches_lossless_branch_limit(epsilon, host):
    lossless = SingleInterface(550, 20, epsilon, host).decay_rates()
    weak_loss = SingleInterface(550, 20, epsilon + 1e-7j, host).decay_rates()
    np.testing.assert_allclose(
        [weak_loss.vertical, weak_loss.parallel], [lossless.vertical, lossless.parallel],
        rtol=1e-5, atol=1e-7,
    )


def test_weak_loss_surface_resonance_is_resolved():
    interface = SingleInterface(550, 20, -5 + 1e-5j)
    result = interface.decay_rates()
    refined = interface.decay_rates(atol=1e-10, rtol=1e-9)
    np.testing.assert_allclose(
        [result.vertical, result.parallel], [refined.vertical, refined.parallel], rtol=1e-7,
    )
    assert result.evanescent_correction[0] > 0.1
    assert result.evanescent_correction[1] > 0.01


def test_far_mirror_tends_to_homogeneous_rate():
    result = SingleInterface(550, 100 * 550 / (4 * np.pi), "pec").decay_rates()
    assert result.vertical == pytest.approx(1, abs=4e-4)
    assert result.parallel == pytest.approx(1, abs=0.016)


def test_lossy_raw_evanescent_spectrum_is_nonzero_and_positive():
    interface = SingleInterface(550, 20, -5 + 0.5j)
    vertical, parallel = interface.spectral_correction(2)
    assert vertical > 0.01
    assert parallel > 0.01
    # At high q the p-dipole near-field terms approach a 2:1 ratio.
    vertical, parallel = interface.spectral_correction(20)
    assert vertical / parallel == pytest.approx(2, rel=0.01)


def test_underresolved_quadrature_fails_explicitly():
    interface = SingleInterface(550, 2000, "pec")
    with pytest.raises(DipoleConvergenceError, match="integration failed"):
        interface.decay_rates(max_subdivisions=1)


@pytest.mark.parametrize("wavelength,distance,epsilon,host", [
    (0, 20, 2.25, 1), (np.nan, 20, 2.25, 1), (np.inf, 20, 2.25, 1),
    (550, 0, 2.25, 1), (550, -1, 2.25, 1), (550, np.nan, 2.25, 1),
    (550, 1e-12, 2.25, 1), (550, 20, 2.25, 0), (550, 20, 2.25, np.inf),
    (550, 20, -1, 1), (550, 20, 0, 1), (550, 20, 2.25 - 0.1j, 1),
    (550, 20, complex(np.nan, 1), 1), (550, 20, complex(1, np.inf), 1),
    (550, 20, "mirror", 1),
    (550, 20, 2.25, 1 + 0.1j),
])
def test_invalid_interface_inputs(wavelength, distance, epsilon, host):
    with pytest.raises(ValueError):
        SingleInterface(wavelength, distance, epsilon, host)


@pytest.mark.parametrize("q", [-1, 1, np.nan, np.inf])
def test_invalid_spectral_coordinates(q):
    with pytest.raises(ValueError, match="q="):
        SingleInterface(550, 20, 2.25).spectral_correction(q)


@pytest.mark.parametrize("kwargs", [
    {"atol": 0}, {"rtol": np.nan}, {"rtol": np.inf}, {"atol": 1e-15},
    {"max_subdivisions": 0}, {"max_subdivisions": 2.5}, {"max_subdivisions": True},
])
def test_invalid_integration_controls(kwargs):
    with pytest.raises(ValueError):
        SingleInterface(550, 20, 2.25).decay_rates(**kwargs)


@pytest.mark.parametrize("theta,phi,orientation", [
    (-1, 0, "vertical"), (np.pi + 0.01, 0, "vertical"), (np.nan, 0, "vertical"),
    (0, np.inf, "vertical"), (0, 0, "unknown"),
])
def test_invalid_pattern_inputs(theta, phi, orientation):
    with pytest.raises(ValueError):
        homogeneous_pattern(theta, phi, orientation)


@pytest.mark.parametrize("wavelength,host,amplitude", [
    (0, 1, 1e-29), (550, 0, 1e-29), (550, 1, np.nan), (550, 1, np.inf),
    (550, 1 + 0.1j, 1e-29), (1e-100, 1, 1e100),
])
def test_invalid_homogeneous_power_inputs(wavelength, host, amplitude):
    with pytest.raises(ValueError):
        homogeneous_power(wavelength, host, amplitude)


@pytest.mark.parametrize("fraction", [-0.01, 1.01, np.nan, np.inf])
def test_invalid_orientation_fraction(fraction):
    with pytest.raises(ValueError):
        SingleInterface(550, 20, "pec").decay_rates().oriented(fraction)
