"""Independent oblique-incidence Fresnel benchmarks for coherent TMM.

Both media must be real, positive, nonabsorbing. Angles are incident angles
in radians; only propagating transmitted waves are tested here.
"""
import unittest
import numpy as np
from tmm import coh_tmm


def fresnel_interface_power(n0, n1, theta0, polarization):
    """Analytical reflectance and transmittance for one lossless interface."""
    if not np.isfinite(n0) or not np.isfinite(n1) or min(n0, n1) <= 0:
        raise ValueError("Indices must be positive and finite")
    if not np.isfinite(theta0) or not 0 <= theta0 < np.pi/2:
        raise ValueError("Incident angle must be in [0, pi/2)")
    if polarization not in ("s", "p"):
        raise ValueError("Polarization must be s or p")
    sin_theta1 = n0 * np.sin(theta0) / n1
    if sin_theta1 >= 1:
        raise ValueError("Reference function requires propagating transmission")
    c0, c1 = np.cos(theta0), np.sqrt(1 - sin_theta1**2)
    if polarization == "s":
        r = (n0*c0 - n1*c1) / (n0*c0 + n1*c1)
        t = 2*n0*c0 / (n0*c0 + n1*c1)
    else:
        r = (n1*c0 - n0*c1) / (n1*c0 + n0*c1)
        t = 2*n0*c0 / (n1*c0 + n0*c1)
    return float(r*r), float((n1*c1)/(n0*c0) * t*t)


class FresnelObliqueBenchmarks(unittest.TestCase):
    def test_te_tm_against_independent_fresnel(self):
        for n0, n1 in ((1., 1.5), (1.5, 1.), (1.33, 1.8)):
            for angle_deg in (0., 10., 25., 40.):
                theta = np.deg2rad(angle_deg)
                if n0*np.sin(theta) >= n1:
                    continue
                for pol in ("s", "p"):
                    expected_r, expected_t = fresnel_interface_power(n0, n1, theta, pol)
                    result = coh_tmm(pol, [n0, n1], [np.inf, np.inf], theta, 550.)
                    self.assertAlmostEqual(result["R"], expected_r, places=11)
                    self.assertAlmostEqual(result["T"], expected_t, places=11)

    def test_brewster_angle_zero_tm_reflectance(self):
        n0, n1 = 1., 1.7
        theta = np.arctan(n1/n0)
        p = coh_tmm("p", [n0, n1], [np.inf, np.inf], theta, 550.)
        s = coh_tmm("s", [n0, n1], [np.inf, np.inf], theta, 550.)
        self.assertLess(p["R"], 1e-12)
        self.assertGreater(s["R"], 0.01)

    def test_total_internal_reflection(self):
        theta = np.deg2rad(60.)
        for pol in ("s", "p"):
            result = coh_tmm(pol, [1.5, 1.], [np.inf, np.inf], theta, 550.)
            self.assertAlmostEqual(result["R"], 1., places=10)
            self.assertAlmostEqual(result["T"], 0., places=10)

    def test_reference_validation(self):
        for args in ((0., 1., .2, "s"), (1., 1.5, -.1, "s"),
                     (1., 1.5, .2, "x"), (1.5, 1., 1.0, "p")):
            with self.assertRaises(ValueError):
                fresnel_interface_power(*args)


if __name__ == "__main__":
    unittest.main()
