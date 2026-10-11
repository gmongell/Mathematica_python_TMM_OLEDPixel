"""Tests for the idealized geometric single-interface dipole escape cone."""
import unittest
import numpy as np
from dipole_escape_cone import critical_angle, escape_cone_fraction, orientation_escape_cone_fractions


class EscapeConeTests(unittest.TestCase):
    def test_equal_indices_have_hemisphere(self):
        self.assertAlmostEqual(critical_angle(1.5, 1.5), np.pi/2)
        self.assertAlmostEqual(escape_cone_fraction(1.5, 1.5), .5)

    def test_total_internal_reflection_critical_angle(self):
        self.assertAlmostEqual(critical_angle(1.5, 1.0), np.arcsin(2/3))

    def test_higher_index_reduces_cone_fraction(self):
        self.assertGreater(escape_cone_fraction(1.5, 1.0), escape_cone_fraction(2.0, 1.0))

    def test_horizontal_and_vertical_distinct(self):
        r = orientation_escape_cone_fractions(1.7, 1.0)
        self.assertGreater(r["horizontal"], r["vertical"])

    def test_invalid_indices_rejected(self):
        for values in ((0., 1.), (1., -1.), (float("nan"), 1.)):
            with self.assertRaises(ValueError):
                critical_angle(*values)

    def test_orientation_fraction_validation(self):
        with self.assertRaises(ValueError):
            escape_cone_fraction(1.7, 1., 1.1)


if __name__ == "__main__":
    unittest.main()
