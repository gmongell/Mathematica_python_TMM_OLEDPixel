"""Illustrative plane-wave transmission through an Alq3/silver stack.

Requires: numpy, matplotlib, tmm (pip install tmm).
All lengths, including wavelength and layer thickness, are in nanometers.
The toy refractive-index formulas are not experimentally validated dispersion.
Transmission is NOT OLED dipole outcoupling efficiency.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tmm import coh_tmm


def n_alq3(wavelength_nm):
    return 1.65 + 0.5 / wavelength_nm**2


def n_silver(wavelength_nm):
    return 0.16 - 1.5 / wavelength_nm**2 + 3.16j


def calculate_transmission(wavelengths_nm, angles_deg, thickness_nm, silver_nm=50.0):
    """Return unpolarized plane-wave power transmission [wavelength, angle]."""
    wavelengths_nm = np.asarray(wavelengths_nm, dtype=float)
    angles_deg = np.asarray(angles_deg, dtype=float)
    if wavelengths_nm.ndim != 1 or angles_deg.ndim != 1:
        raise ValueError("Wavelengths and angles must be one-dimensional")
    if not np.all(np.isfinite(wavelengths_nm)) or np.any(wavelengths_nm <= 0):
        raise ValueError("Wavelengths must be finite and positive")
    if not np.all(np.isfinite(angles_deg)) or np.any((angles_deg < 0) | (angles_deg >= 90)):
        raise ValueError("Incidence angles must be in [0, 90) degrees")
    if not np.isfinite(thickness_nm) or not np.isfinite(silver_nm) or min(thickness_nm, silver_nm) < 0:
        raise ValueError("Thicknesses must be finite and nonnegative")

    result = np.empty((len(wavelengths_nm), len(angles_deg)))
    d_list = [np.inf, float(thickness_nm), float(silver_nm), np.inf]
    for i, wavelength_nm in enumerate(wavelengths_nm):
        # coh_tmm requires one refractive index per layer, evaluated at this wavelength.
        n_list = [1.0, n_alq3(wavelength_nm), n_silver(wavelength_nm), 1.0]
        for j, theta_deg in enumerate(angles_deg):
            theta = np.deg2rad(theta_deg)
            ts = coh_tmm("s", n_list, d_list, wavelength_nm, theta)["T"]
            tp = coh_tmm("p", n_list, d_list, wavelength_nm, theta)["T"]
            result[i, j] = 0.5 * (ts + tp)
    return result


def main():
    wavelengths_nm = np.linspace(400, 700, 100)
    angles_deg = np.linspace(0, 89, 90)
    for thickness_nm in np.linspace(25, 150, 10):
        transmission = calculate_transmission(wavelengths_nm, angles_deg, thickness_nm)
        fig, ax = plt.subplots(figsize=(10, 6))
        image = ax.contourf(angles_deg, wavelengths_nm, transmission, levels=50)
        fig.colorbar(image, ax=ax, label="Unpolarized plane-wave transmittance")
        ax.set_xlabel("Angle of incidence (degrees)")
        ax.set_ylabel("Wavelength (nm)")
        ax.set_title(f"Alq3 thickness = {thickness_nm:.1f} nm; illustrative indices")
        fig.savefig(f"transmittance_thickness_{thickness_nm:.1f}nm.png")
        plt.close(fig)


if __name__ == "__main__":
    main()
