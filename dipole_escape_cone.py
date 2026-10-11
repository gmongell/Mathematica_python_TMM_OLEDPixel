"""Single-interface dipole radiation escape-cone estimate.

A deliberately limited geometric-optics model: emission is homogeneous and
isotropic in azimuth inside a lossless high-index emitter, and rays above
the critical angle are totally internally reflected. This does NOT
include Fresnel transmission, multiple reflections, LDOS, cavity effects,
absorption, or frustrated total internal reflection.
"""
import numpy as np
from dipole_analytic_cone import analytic_collection_fraction, analytic_mixture_fraction


def critical_angle(n_emitter, n_exit):
    """Critical angle in radians for lossless real positive indices.

    Returns pi/2 when there is no total-internal-reflection cutoff.
    """
    if not np.isfinite(n_emitter) or not np.isfinite(n_exit):
        raise ValueError("Refractive indices must be finite")
    if n_emitter <= 0 or n_exit <= 0:
        raise ValueError("Refractive indices must be positive")
    if n_emitter <= n_exit:
        return np.pi / 2
    return float(np.arcsin(n_exit / n_emitter))


def escape_cone_fraction(n_emitter, n_exit, horizontal_fraction=2/3):
    """Fraction of total dipole radiation inside the forward escape cone.

    This is an upper bound only for the specified geometric escape-cone
    selection model; it is NOT an upper bound for arbitrary OLED stacks.
    """
    theta = critical_angle(n_emitter, n_exit)
    return float(analytic_mixture_fraction(theta, horizontal_fraction))


def orientation_escape_cone_fractions(n_emitter, n_exit):
    """Return separate horizontal and vertical dipole cone fractions."""
    theta = critical_angle(n_emitter, n_exit)
    return {
        "horizontal": float(analytic_collection_fraction(theta, "horizontal")),
        "vertical": float(analytic_collection_fraction(theta, "vertical")),
        "critical_angle_deg": float(np.rad2deg(theta)),
    }
