"""Tests for the white-spectrum plane-wave optimization proxy."""
import unittest
import numpy as np
from white_oled_proxy_optimizer import reference_spectrum, evaluate_proxy, optimize_proxy


class WhiteProxyTests(unittest.TestCase):
    def setUp(self):
        self.wavelengths = np.linspace(420, 680, 9)
        self.angles = np.linspace(0, 80, 7)

    def test_spectrum_normalization(self):
        s = reference_spectrum(self.wavelengths)
        self.assertTrue(np.all(s >= 0))
        self.assertAlmostEqual(float(np.trapz(s, self.wavelengths)), 1.0, places=10)

    def test_metrics(self):
        m = evaluate_proxy(60, 20, self.wavelengths, self.angles)
        self.assertGreaterEqual(m["hemispherical_transmission_proxy"], 0)
        self.assertLessEqual(m["hemispherical_transmission_proxy"], 1.000001)
        self.assertGreaterEqual(m["angular_coefficient_of_variation"], 0)

    def test_reject_invalid_spectrum(self):
        with self.assertRaises(ValueError):
            evaluate_proxy(60, 20, self.wavelengths, self.angles,
                           spectrum=np.zeros_like(self.wavelengths))

    def test_small_reproducible_optimization(self):
        a = optimize_proxy(self.wavelengths, self.angles, maxiter=1, seed=17)
        b = optimize_proxy(self.wavelengths, self.angles, maxiter=1, seed=17)
        self.assertAlmostEqual(a["alq3_nm"], b["alq3_nm"], places=8)
        self.assertAlmostEqual(a["silver_nm"], b["silver_nm"], places=8)


if __name__ == "__main__":
    unittest.main()
