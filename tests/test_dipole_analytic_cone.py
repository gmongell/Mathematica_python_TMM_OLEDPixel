"""Compare numerical dipole angular quadrature against analytic cone formulas."""
import unittest
import numpy as np
from dipole_analytic_cone import analytic_collection_fraction, analytic_mixture_fraction
from dipole_reference_patterns import integrated_fraction


class AnalyticConeTests(unittest.TestCase):
    def test_endpoints(self):
        for orientation in ("vertical", "horizontal"):
            self.assertAlmostEqual(analytic_collection_fraction(0, orientation), 0.)
            self.assertAlmostEqual(analytic_collection_fraction(np.pi, orientation), 1.)

    def test_numerical_quadrature_agrees_with_analytic(self):
        for orientation in ("vertical", "horizontal"):
            for theta in (np.pi / 12, np.pi / 4, np.pi / 2, 3*np.pi / 4):
                numeric = integrated_fraction(theta, orientation, 301, 361)
                exact = analytic_collection_fraction(theta, orientation)
                self.assertAlmostEqual(numeric, exact, delta=4e-5)

    def test_isotropic_mixture(self):
        for theta in (0., np.pi / 6, np.pi / 2, np.pi):
            self.assertAlmostEqual(
                analytic_mixture_fraction(theta, 2/3),
                (1 - np.cos(theta))/2,
                places=12,
            )

    def test_invalid_angles(self):
        with self.assertRaises(ValueError):
            analytic_collection_fraction(-.1)


if __name__ == "__main__":
    unittest.main()
