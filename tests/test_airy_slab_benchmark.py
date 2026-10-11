"""Independent Airy/Fresnel analytical verification of coherent TMM.

The reference equations below are evaluated directly from interface Fresnel
coefficients rather than by calling the transfer-matrix implementation.
Assumptions: normal incidence, real positive indices, coherent planar film,
lossless incident and exit media, no scattering.
"""
import unittest
import numpy as np
from tmm import coh_tmm


def airy_slab_transmission(n_incident, n_film, n_exit, thickness_nm, wavelength_nm):
    """Normal-incidence power transmittance of one lossless film."""
    indices = (n_incident, n_film, n_exit)
    if not all(np.isfinite(n) and n > 0 for n in indices):
        raise ValueError("Indices must be real, finite and positive")
    if not np.isfinite(thickness_nm) or thickness_nm < 0:
        raise ValueError("Film thickness must be finite and nonnegative")
    if not np.isfinite(wavelength_nm) or wavelength_nm <= 0:
        raise ValueError("Wavelength must be finite and positive")
    r01 = (n_incident - n_film) / (n_incident + n_film)
    r12 = (n_film - n_exit) / (n_film + n_exit)
    t01 = 2 * n_incident / (n_incident + n_film)
    t12 = 2 * n_film / (n_film + n_exit)
    phase = 2 * np.pi * n_film * thickness_nm / wavelength_nm
    amplitude = t01 * t12 * np.exp(1j * phase) / (
        1 + r01 * r12 * np.exp(2j * phase)
    )
    return float((n_exit / n_incident) * abs(amplitude)**2)


class AirySlabBenchmarks(unittest.TestCase):
    def test_air_spacer_air_against_closed_form(self):
        for n_film in (1.2, 1.5, 2.1):
            for wavelength in (405., 532., 650.):
                for thickness in (0., 20., wavelength/(4*n_film),
                                  wavelength/(2*n_film), 240.):
                    expected = airy_slab_transmission(
                        1., n_film, 1., thickness, wavelength
                    )
                    for polarization in ("s", "p"):
                        actual = coh_tmm(
                            polarization, [1., n_film, 1.],
                            [np.inf, thickness, np.inf], 0., wavelength
                        )["T"]
                        self.assertAlmostEqual(actual, expected, places=11)

    def test_asymmetric_interfaces(self):
        for n0, n1, n2 in ((1., 1.6, 1.3), (1.5, 2.2, 1.)):
            for thickness in (0., 77., 180.):
                expected = airy_slab_transmission(n0, n1, n2, thickness, 550.)
                actual = coh_tmm(
                    "s", [n0, n1, n2],
                    [np.inf, thickness, np.inf], 0., 550.
                )["T"]
                self.assertAlmostEqual(actual, expected, places=11)

    def test_half_wave_resonance_in_symmetric_medium(self):
        wavelength, n_film = 600., 1.8
        thickness = wavelength / (2 * n_film)
        self.assertAlmostEqual(
            airy_slab_transmission(1., n_film, 1., thickness, wavelength),
            1., places=12
        )

    def test_invalid_reference_inputs(self):
        for args in ((1., -1., 1., 80., 550.),
                     (1., 1.5, 1., -2., 550.),
                     (1., 1.5, 1., 80., 0.)):
            with self.assertRaises(ValueError):
                airy_slab_transmission(*args)


if __name__ == "__main__":
    unittest.main()
