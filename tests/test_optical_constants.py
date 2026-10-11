"""Regression tests for measured optical-constant CSV ingestion."""
import os
import tempfile
import unittest
import numpy as np
from optical_constants import OpticalConstants


class OpticalConstantsTests(unittest.TestCase):
    def test_complex_interpolation(self):
        table = OpticalConstants([400., 500., 600.], [1.4, 1.6, 1.8], [0., .1, .2])
        self.assertAlmostEqual(table.evaluate(450.).real, 1.5)
        self.assertAlmostEqual(table.evaluate(450.).imag, .05)
        np.testing.assert_allclose(table.evaluate([400., 600.]), [1.4+0j, 1.8+.2j])

    def test_csv_roundtrip(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "measured.csv")
            with open(path, "w", encoding="utf-8") as f:
                f.write("wavelength_nm,n,k\n450,1.7,0.01\n550,1.8,0.03\n")
            self.assertAlmostEqual(OpticalConstants.from_csv(path).evaluate(500.).imag, .02)

    def test_extrapolation_rejected(self):
        table = OpticalConstants([400., 600.], [1.5, 1.6], [0., 0.])
        for wavelength in (399., 601., np.nan):
            with self.assertRaises(ValueError):
                table.evaluate(wavelength)

    def test_invalid_physics_rejected(self):
        cases = [
            ([500., 400.], [1., 1.], [0., 0.]),
            ([400., 500.], [1., 1.], [-.1, 0.]),
            ([400., 500.], [0., 1.], [0., 0.]),
            ([400.], [1.], [0.]),
        ]
        for args in cases:
            with self.assertRaises(ValueError):
                OpticalConstants(*args)


if __name__ == "__main__":
    unittest.main()
