"""Regression checks for coupled white-light transmission/color optimization.

Uses synthetic color matching functions solely to test mathematical behavior.
"""
import unittest
import numpy as np
from white_oled_color_optimizer import evaluate_design, optimize_color_uniformity


class CoupledOptimizationTests(unittest.TestCase):
    def setUp(self):
        self.w = np.linspace(430, 670, 9)
        self.a = np.linspace(0, 75, 6)
        self.cmf = np.column_stack((
            np.exp(-((self.w - 620) / 80)**2),
            np.exp(-((self.w - 540) / 65)**2),
            np.exp(-((self.w - 460) / 60)**2),
        ))

    def test_metrics_are_finite(self):
        result = evaluate_design(75, 20, self.w, self.a, self.cmf)
        self.assertTrue(np.isfinite(result["proxy"]["hemispherical_transmission_proxy"]))
        self.assertTrue(np.isfinite(result["color"]["max_delta_uv"]))
        self.assertGreaterEqual(result["color"]["max_delta_uv"], 0)

    def test_optimizer_reproducible_with_fixed_seed(self):
        first = optimize_color_uniformity(self.w, self.a, self.cmf, maxiter=1, seed=3)
        second = optimize_color_uniformity(self.w, self.a, self.cmf, maxiter=1, seed=3)
        self.assertAlmostEqual(first["alq3_nm"], second["alq3_nm"], places=8)
        self.assertAlmostEqual(first["silver_nm"], second["silver_nm"], places=8)

    def test_negative_penalties_rejected(self):
        with self.assertRaises(ValueError):
            optimize_color_uniformity(self.w, self.a, self.cmf, color_penalty=-1)


if __name__ == "__main__":
    unittest.main()
