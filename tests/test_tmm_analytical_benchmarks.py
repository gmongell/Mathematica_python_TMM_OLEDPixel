"""Independent analytical benchmarks for planar thin-film transmission.

A zero-thickness, lossless film must reduce to the Fresnel single-interface
power transmission, independent of its nominal refractive index.
"""
import unittest
import numpy as np
from tmm import coh_tmm
from tmm_OLED9 import calculate_transmission


class ThinFilmBenchmarks(unittest.TestCase):
    def test_zero_thickness_stack_is_transparent(self):
        wavelengths = np.array([450., 550., 650.])
        angles = np.array([0., 30., 60.])
        t = calculate_transmission(wavelengths, angles, 0., 0.)
        np.testing.assert_allclose(t, 1., atol=1e-10)

    def test_fresnel_interface_normal_incidence(self):
        n0, n1 = 1.0, 1.5
        expected = 4 * n0 * n1 / (n0 + n1)**2
        for pol in ("s", "p"):
            result = coh_tmm(pol, [n0, n1], [np.inf, np.inf], 0., 550.)
            self.assertAlmostEqual(result["T"], expected, places=12)
            self.assertAlmostEqual(result["R"] + result["T"], 1., places=12)

    def test_lossless_multilayer_energy_conservation(self):
        layers = [1., 1.7, 1.45, 1.]
        thicknesses = [np.inf, 90., 130., np.inf]
        for wavelength in (430., 550., 680.):
            for angle in (0., np.pi / 6, np.pi / 3):
                for pol in ("s", "p"):
                    result = coh_tmm(pol, layers, thicknesses, angle, wavelength)
                    self.assertAlmostEqual(result["R"] + result["T"], 1., places=10)

    def test_passive_absorber_energy_budget(self):
        layers = [1., 1.7, 0.15 + 3.1j, 1.]
        thicknesses = [np.inf, 80., 20., np.inf]
        for pol in ("s", "p"):
            result = coh_tmm(pol, layers, thicknesses, np.pi / 4, 550.)
            self.assertGreaterEqual(result["R"], 0.)
            self.assertGreaterEqual(result["T"], 0.)
            self.assertLessEqual(result["R"] + result["T"], 1. + 1e-10)


if __name__ == "__main__":
    unittest.main()
