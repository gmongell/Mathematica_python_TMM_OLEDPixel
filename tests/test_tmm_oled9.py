"""Lightweight regression checks for the corrected plane-wave TMM example.

Run: python -m unittest discover -s tests -v
Requires numpy, matplotlib, and tmm.
"""
import unittest
import numpy as np
from tmm_OLED9 import calculate_transmission


class TransmissionTests(unittest.TestCase):
    def test_shape_and_physical_bounds(self):
        t = calculate_transmission([450.0, 550.0], [0.0, 30.0], 50.0)
        self.assertEqual(t.shape, (2, 2))
        self.assertTrue(np.all(np.isfinite(t)))
        self.assertTrue(np.all(t >= -1e-10))
        self.assertTrue(np.all(t <= 1.0 + 1e-10))

    def test_invalid_input(self):
        for wavelengths, angles, thickness in [
            ([0.0], [0.0], 50.0),
            ([550.0], [90.0], 50.0),
            ([550.0], [0.0], -1.0),
            ([float("nan")], [0.0], 50.0),
        ]:
            with self.subTest(wavelengths=wavelengths, angles=angles, thickness=thickness):
                with self.assertRaises(ValueError):
                    calculate_transmission(wavelengths, angles, thickness)


if __name__ == "__main__":
    unittest.main()
