"""Closed-form collection-cone fractions for homogeneous-medium dipoles.

These equations apply to free-space angular radiation patterns, NOT to a
stratified optical cavity or a dipole near an absorbing interface.
"""
import numpy as np


def analytic_collection_fraction(theta_max_rad, orientation="vertical"):
    """Fraction of total radiated power collected inside a polar cone.

    Vertical z dipole: (2 - 3 cos(theta) + cos(theta)^3)/4.
    Horizontal x/y dipole (azimuth integrated):
        (4 - 3 cos(theta) - cos(theta)^3)/8.
    Both equal 1/2 at pi/2 and 1 at pi.
    """
    theta = np.asarray(theta_max_rad, dtype=float)
    if not np.all(np.isfinite(theta)) or np.any(theta < 0) or np.any(theta > np.pi):
        raise ValueError("Cone angle must lie within [0, pi]")
    c = np.cos(theta)
    if orientation == "vertical":
        result = (2 - 3*c + c**3) / 4
    elif orientation == "horizontal":
        result = (4 - 3*c - c**3) / 8
    else:
        raise ValueError("Unknown dipole orientation")
    return float(result) if result.ndim == 0 else result


def analytic_mixture_fraction(theta_max_rad, horizontal_fraction=2/3):
    """Orientation-weighted homogeneous-medium collection fraction."""
    if not np.isfinite(horizontal_fraction) or not 0 <= horizontal_fraction <= 1:
        raise ValueError("horizontal_fraction must lie within [0, 1]")
    return (horizontal_fraction * analytic_collection_fraction(theta_max_rad, "horizontal")
            + (1-horizontal_fraction) * analytic_collection_fraction(theta_max_rad, "vertical"))
