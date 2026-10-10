"""Multiobjective *plane-wave proxy* optimization for white OLED optical stacks.

This does not model dipole LDOS, orientation, Purcell factor, internal quantum
efficiency, color-rendering index, or absolute external quantum efficiency.
Use measured wavelength-dependent n,k and emitter spectra before design use.

Run: python white_oled_proxy_optimizer.py --samples 30 --seed 7
Requires: numpy, scipy, matplotlib, tmm.
"""
import argparse
import json
import numpy as np
from scipy.optimize import differential_evolution
from tmm_OLED9 import calculate_transmission

# Illustrative normalized RGB-like white spectral mixture; not a calibrated white SPD.
CENTERS_NM = (460.0, 530.0, 620.0)
WIDTHS_NM = (22.0, 38.0, 42.0)
WEIGHTS = (0.25, 0.45, 0.30)


def reference_spectrum(wavelengths_nm):
    w = np.asarray(wavelengths_nm, dtype=float)
    s = sum(a * np.exp(-0.5 * ((w - mu) / sigma)**2)
            for a, mu, sigma in zip(WEIGHTS, CENTERS_NM, WIDTHS_NM))
    return s / np.trapz(s, w)


def evaluate_proxy(thickness_nm, silver_nm, wavelengths_nm, angles_deg, spectrum=None):
    """Return hemispherical weighted transmission and angular nonuniformity.

    Angular transmission is spectrum-weighted. The hemispherical measure
    integrates T(theta) sin(theta)cos(theta) and normalizes by its integral.
    Angular coefficient of variation uses solid-angle weights sin(theta),
    which avoids overemphasizing uniform samples near normal incidence.
    """
    w = np.asarray(wavelengths_nm, float)
    a = np.asarray(angles_deg, float)
    if w.ndim != 1 or a.ndim != 1 or len(w) < 2 or len(a) < 2:
        raise ValueError("Need 1D wavelength and angle grids with at least two points")
    if np.any(np.diff(w) <= 0) or np.any(np.diff(a) <= 0):
        raise ValueError("Grids must be strictly increasing")
    if a[0] < 0 or a[-1] >= 90:
        raise ValueError("Angle grid must lie in [0, 90) degrees")
    s = reference_spectrum(w) if spectrum is None else np.asarray(spectrum, float)
    if s.shape != w.shape or not np.all(np.isfinite(s)) or np.any(s < 0):
        raise ValueError("Spectrum must be finite, nonnegative, and match wavelengths")
    norm = np.trapz(s, w)
    if norm <= 0:
        raise ValueError("Spectrum integral must be positive")
    s = s / norm
    t = calculate_transmission(w, a, thickness_nm, silver_nm)
    spectral_t = np.trapz(t * s[:, None], w, axis=0)
    theta = np.deg2rad(a)
    solid_weight = np.sin(theta)
    projection_weight = solid_weight * np.cos(theta)
    flux = np.trapz(spectral_t * projection_weight, theta) / np.trapz(projection_weight, theta)
    mean = np.trapz(spectral_t * solid_weight, theta) / np.trapz(solid_weight, theta)
    variance = np.trapz((spectral_t - mean)**2 * solid_weight, theta) / np.trapz(solid_weight, theta)
    cv = float(np.sqrt(max(variance, 0.0)) / max(mean, 1e-12))
    return {"hemispherical_transmission_proxy": float(flux),
            "angular_coefficient_of_variation": cv,
            "normal_incidence_transmission": float(spectral_t[0]) if a[0] == 0 else None}


def optimize_proxy(wavelengths_nm, angles_deg, uniformity_penalty=0.3,
                   bounds=((25.0, 150.0), (5.0, 60.0)), maxiter=30, seed=7):
    """Optimize Alq3 and Ag thickness in nm, with a reproducible search seed."""
    if not np.isfinite(uniformity_penalty) or uniformity_penalty < 0:
        raise ValueError("uniformity_penalty must be nonnegative")
    def objective(x):
        m = evaluate_proxy(x[0], x[1], wavelengths_nm, angles_deg)
        return -m["hemispherical_transmission_proxy"] + uniformity_penalty * m["angular_coefficient_of_variation"]
    result = differential_evolution(objective, bounds, seed=seed, maxiter=maxiter,
                                    popsize=6, polish=True, workers=1)
    metrics = evaluate_proxy(*result.x, wavelengths_nm, angles_deg)
    return {"alq3_nm": float(result.x[0]), "silver_nm": float(result.x[1]),
            "objective": float(result.fun), "optimizer_success": bool(result.success),
            "optimizer_message": str(result.message), "metrics": metrics,
            "note": "Plane-wave proxy only; NOT OLED EQE, dipole outcoupling, or white-point/color uniformity."}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--samples", type=int, default=30, help="Wavelength grid points")
    p.add_argument("--angles", type=int, default=20, help="Angular grid points")
    p.add_argument("--iterations", type=int, default=30)
    p.add_argument("--penalty", type=float, default=0.3)
    p.add_argument("--seed", type=int, default=7)
    args = p.parse_args()
    if args.samples < 2 or args.angles < 2 or args.iterations < 1:
        p.error("samples and angles must be >= 2; iterations >= 1")
    wavelengths = np.linspace(400, 700, args.samples)
    angles = np.linspace(0, 85, args.angles)
    print(json.dumps(optimize_proxy(wavelengths, angles, args.penalty,
                                    maxiter=args.iterations, seed=args.seed), indent=2))


if __name__ == "__main__":
    main()
