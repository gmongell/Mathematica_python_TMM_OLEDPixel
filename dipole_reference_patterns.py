"""Radiation-pattern integration for orientation-resolved electric dipoles.

This module handles the *homogeneous, lossless* reference medium only.
It is not a multilayer Green-tensor or LDOS solver. The orientation
convention is a dipole aligned to z (vertical) or x (horizontal), with z
as the optical axis and theta measured from z.
"""
import numpy as np


def dipole_angular_pattern(theta_rad, phi_rad, orientation="vertical"):
    """Return dimensionless far-field dP/dOmega proportionality."""
    theta = np.asarray(theta_rad, dtype=float)
    phi = np.asarray(phi_rad, dtype=float)
    if not np.all(np.isfinite(theta)) or not np.all(np.isfinite(phi)):
        raise ValueError("Angles must be finite")
    if np.any(theta < 0) or np.any(theta > np.pi):
        raise ValueError("Polar angles must lie in [0, pi]")
    if orientation == "vertical":
        return np.sin(theta)**2 + np.zeros(np.broadcast_shapes(theta.shape, phi.shape))
    if orientation == "horizontal":
        return 1. - (np.sin(theta) * np.cos(phi))**2
    raise ValueError("orientation must be vertical or horizontal")


def integrated_fraction(theta_max_rad, orientation="vertical", theta_samples=401,
                        phi_samples=721):
    """Fraction of total homogeneous-medium radiation in 0 <= theta <= theta_max.

    Integrates dP/dOmega sin(theta) dtheta dphi. Both orientations have
    total radiation integral 8*pi/3 over the full sphere.
    """
    if not np.isfinite(theta_max_rad) or not 0 <= theta_max_rad <= np.pi:
        raise ValueError("theta_max_rad must lie in [0, pi]")
    if theta_samples < 3 or phi_samples < 3:
        raise ValueError("At least three samples are required per angle")
    theta = np.linspace(0., theta_max_rad, theta_samples)
    phi = np.linspace(0., 2 * np.pi, phi_samples)
    pattern = dipole_angular_pattern(theta[:, None], phi[None, :], orientation)
    azimuthal = np.trapz(pattern, phi, axis=1)
    power = np.trapz(azimuthal * np.sin(theta), theta)
    return float(power / (8 * np.pi / 3))


def orientation_mixture_fraction(theta_max_rad, horizontal_fraction=2./3.,
                                 theta_samples=401, phi_samples=721):
    """Mixture of x/y-equivalent horizontal and z-oriented dipoles."""
    if not np.isfinite(horizontal_fraction) or not 0 <= horizontal_fraction <= 1:
        raise ValueError("horizontal_fraction must lie in [0, 1]")
    h = integrated_fraction(theta_max_rad, "horizontal", theta_samples, phi_samples)
    v = integrated_fraction(theta_max_rad, "vertical", theta_samples, phi_samples)
    return horizontal_fraction * h + (1 - horizontal_fraction) * v
