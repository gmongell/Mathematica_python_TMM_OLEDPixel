"""Numerical colorimetry tests with synthetic CMFs (not CIE reference data)."""
import unittest
import numpy as np
from white_oled_colorimetry import calculate_angular_color


class AngularColorTests(unittest.TestCase):
    def setUp(self):
        self.w = np.array([450., 550., 650.])
        self.a = np.array([0., 30., 60.])
        self.cmf = np.array([[0.3, 0.1, 0.8], [0.4, 1.0, 0.2], [0.8, 0.3, 0.1]])
        self.spd = np.array([[1., 2., 3.], [2., 4., 6.], [3., 6., 9.]])

    def test_angle_independent_color_when_only_intensity_changes(self):
        result = calculate_angular_color(self.w, self.a, self.spd, self.cmf)
        self.assertAlmostEqual(result["max_delta_uv"], 0., places=12)
        self.assertEqual(len(result["xy"]), 3)

    def test_spectral_shape_change_shifts_color(self):
        altered = self.spd.copy()
        altered[0, 2] *= 4
        result = calculate_angular_color(self.w, self.a, altered, self.cmf)
        self.assertGreater(result["max_delta_uv"], 0.)

    def test_zero_spectrum_is_rejected(self):
        with self.assertRaises(ValueError):
            calculate_angular_color(self.w, self.a, np.zeros_like(self.spd), self.cmf)

    def test_shape_mismatch_is_rejected(self):
        with self.assertRaises(ValueError):
            calculate_angular_color(self.w, self.a, self.spd.T[:2], self.cmf)


if __name__ == "__main__":
    unittest.main()
