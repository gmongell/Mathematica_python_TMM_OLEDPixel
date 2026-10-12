"""Plane-wave multilayer transmission sweep (not OLED dipole outcoupling)."""
from pathlib import Path
import argparse
import numpy as np
from tmm import coh_tmm


def n_alq3(wavelength_nm):
    """Illustrative dispersion only; not calibrated material data."""
    return 1.65 + 0.5 / wavelength_nm**2


def n_silver(wavelength_nm):
    """Illustrative complex refractive index; replace with measured data."""
    return 0.16 - 1.5 / wavelength_nm**2 + 3.16j


def stack_at(wavelength_nm, organic_nm, silver_nm=50.0):
    if wavelength_nm <= 0 or organic_nm < 0 or silver_nm < 0:
        raise ValueError("Wavelength must be positive and thicknesses nonnegative")
    return ([1.0, n_alq3(wavelength_nm), n_silver(wavelength_nm), 1.0],
            [np.inf, organic_nm, silver_nm, np.inf])


def transmission(wavelength_nm, angle_deg, organic_nm, polarization):
    if not 0 <= angle_deg < 90:
        raise ValueError("Angle must be in [0, 90) degrees")
    if polarization not in ("s", "p"):
        raise ValueError("Polarization must be s or p")
    indices, thicknesses = stack_at(wavelength_nm, organic_nm)
    result = coh_tmm(polarization, indices, thicknesses,
                     np.deg2rad(angle_deg), wavelength_nm)
    return float(result["T"])


def sweep(wavelengths_nm, angles_deg, organic_nm):
    wavelengths_nm = np.asarray(wavelengths_nm, dtype=float)
    angles_deg = np.asarray(angles_deg, dtype=float)
    if wavelengths_nm.ndim != 1 or angles_deg.ndim != 1:
        raise ValueError("Sweep coordinates must be 1-D")
    if not np.all(np.isfinite(wavelengths_nm)) or not np.all(np.isfinite(angles_deg)):
        raise ValueError("Sweep coordinates must be finite")
    values = np.empty((len(wavelengths_nm), len(angles_deg)))
    for i, wavelength in enumerate(wavelengths_nm):
        for j, angle in enumerate(angles_deg):
            values[i, j] = 0.5 * (
                transmission(wavelength, angle, organic_nm, "s")
                + transmission(wavelength, angle, organic_nm, "p"))
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wavelength-count", type=int, default=31)
    parser.add_argument("--angle-count", type=int, default=19)
    parser.add_argument("--thickness-nm", type=float, default=75.0)
    parser.add_argument("--output", type=Path, default=Path("transmission.png"))
    args = parser.parse_args()
    if args.wavelength_count < 2 or args.angle_count < 2:
        parser.error("sample counts must be >= 2")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    wavelengths = np.linspace(400.0, 700.0, args.wavelength_count)
    angles = np.linspace(0.0, 89.0, args.angle_count)
    values = sweep(wavelengths, angles, args.thickness_nm)
    fig, ax = plt.subplots()
    plot = ax.contourf(angles, wavelengths, values, levels=30)
    fig.colorbar(plot, ax=ax, label="Unpolarized plane-wave transmittance")
    ax.set(xlabel="Incidence angle (degrees)", ylabel="Vacuum wavelength (nm)",
           title=f"Illustrative stack: organic {args.thickness_nm:g} nm")
    fig.savefig(args.output, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
