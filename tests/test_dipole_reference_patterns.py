"""Analytical limits for homogeneous-medium electric-dipole radiation."""
import unittest
import numpy as np
from dipole_reference_patterns import (
    dipole_angular_pattern, integrated_fraction, orientation_mixture_fraction
)


class DipoleReferenceTests(unittest.TestCase):
    def test_vertical_dipole_has_axis_null(self):
        self.assertAlmostEqual(float(dipole_angular_pattern(0., 0., "vertical")), 0.)
        self.assertAlmostEqual(float(dipole_angular_pattern(np.pi / 2, 0., "vertical")), 1.)

    def test_horizontal_dipole_has_axis_maximum(self):
        self.assertAlmostEqual(float(dipole_angular_pattern(0., 0., "horizontal")), 1.)
        self.assertAlmostEqual(float(dipole_angular_pattern(np.pi / 2, 0., "horizontal")), 0.)

    def test_full_sphere_integrates_to_unity(self):
        for orientation in ("vertical", "horizontal"):
            self.assertAlmostEqual(integrated_fraction(np.pi, orientation), 1., delta=3e-5)

    def test_upper_hemisphere_is_half(self):
        for orientation in ("vertical", "horizontal"):
            self.assertAlmostEqual(integrated_fraction(np.pi / 2, orientation), .5, delta=2e-5)

    def test_isotropic_orientation_average(self):
        # 2/3 horizontal and 1/3 vertical makes the ensemble isotropic.
        theta = np.pi / 3
        fraction = orientation_mixture_fraction(theta)
        self.assertAlmostEqual(fraction, (1 - np.cos(theta)) / 2, delta=2e-5)

    def test_invalid_orientation_rejected(self):
        with self.assertRaises(ValueError):
            dipole_angular_pattern(.3, .5, "diagonal")


if __name__ == "__main__":
    unittest.main()
