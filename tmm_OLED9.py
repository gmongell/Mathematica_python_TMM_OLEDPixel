"""Illustrative plane-wave multilayer optics; NOT OLED dipole outcoupling."""

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from tmm import coh_tmm

__version__ = "0.1.0a1"


def _positive_wavelength(wavelength_nm: float) -> None:
    if not np.isfinite(wavelength_nm) or wavelength_nm <= 0:
        raise ValueError("Wavelength must be finite and positive (nm)")


def _nonnegative_thickness(thickness_nm: float) -> None:
    if not np.isfinite(thickness_nm) or thickness_nm < 0:
        raise ValueError("Thickness must be finite and nonnegative (nm)")


def n_alq3(wavelength_nm: float) -> float:
    """Illustrative dispersion only, not measured or calibrated Alq3 data."""
    _positive_wavelength(wavelength_nm)
    return 1.65 + 0.5 / wavelength_nm**2


def n_silver(wavelength_nm: float) -> complex:
    """Illustrative passive index, not measured or calibrated silver data."""
    _positive_wavelength(wavelength_nm)
    return 0.16 - 1.5 / wavelength_nm**2 + 3.16j


def stack_at(
    wavelength_nm: float, organic_nm: float, silver_nm: float = 50.0
) -> tuple[list[complex], list[float]]:
    _positive_wavelength(wavelength_nm)
    _nonnegative_thickness(organic_nm)
    _nonnegative_thickness(silver_nm)
    return (
        [1.0, n_alq3(wavelength_nm), n_silver(wavelength_nm), 1.0],
        [np.inf, organic_nm, silver_nm, np.inf],
    )


def power_coefficients(
    wavelength_nm: float,
    angle_deg: float,
    organic_nm: float,
    polarization: str,
    silver_nm: float = 50.0,
) -> dict[str, float]:
    """Return R, T and total A for air/organic/silver/air plane-wave incidence.

    A = 1 - R - T applies here because both external media are lossless air.
    All lengths are nm; angle is in degrees. No emitting dipoles are modeled.
    """
    if not np.isfinite(angle_deg) or not 0 <= angle_deg < 90:
        raise ValueError("Angle must be finite and in [0, 90) degrees")
    if polarization not in ("s", "p"):
        raise ValueError("Polarization must be s or p")
    indices, thicknesses = stack_at(wavelength_nm, organic_nm, silver_nm)
    result = coh_tmm(
        polarization, indices, thicknesses, np.deg2rad(angle_deg), wavelength_nm
    )
    reflection, transmittance = float(result["R"]), float(result["T"])
    absorption = 1.0 - reflection - transmittance
    values = np.array([reflection, transmittance, absorption])
    if not np.all(np.isfinite(values)) or np.any(values < -1e-9) or np.any(values > 1 + 1e-9):
        raise ArithmeticError("Solver returned nonphysical passive-stack power coefficients")
    return {"R": reflection, "T": transmittance, "A": absorption}


def transmission(
    wavelength_nm: float, angle_deg: float, organic_nm: float, polarization: str
) -> float:
    return power_coefficients(wavelength_nm, angle_deg, organic_nm, polarization)["T"]


def sweep(
    wavelengths_nm: ArrayLike, angles_deg: ArrayLike, organic_nm: float
) -> NDArray[np.float64]:
    """Return equally weighted s/p transmittance, indexed [wavelength, angle]."""
    wavelengths = np.asarray(wavelengths_nm, dtype=float)
    angles = np.asarray(angles_deg, dtype=float)
    _nonnegative_thickness(organic_nm)
    if wavelengths.ndim != 1 or angles.ndim != 1 or not wavelengths.size or not angles.size:
        raise ValueError("Sweep coordinates must be nonempty 1-D arrays")
    if not np.all(np.isfinite(wavelengths)) or np.any(wavelengths <= 0):
        raise ValueError("Wavelengths must be finite and positive (nm)")
    if not np.all(np.isfinite(angles)) or np.any(angles < 0) or np.any(angles >= 90):
        raise ValueError("Angles must be finite and in [0, 90) degrees")
    values = np.empty((len(wavelengths), len(angles)))
    for i, wavelength in enumerate(wavelengths):
        for j, angle in enumerate(angles):
            values[i, j] = 0.5 * (
                transmission(wavelength, angle, organic_nm, "s")
                + transmission(wavelength, angle, organic_nm, "p")
            )
    return values


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("--wavelength-count", type=int, default=31)
    parser.add_argument("--angle-count", type=int, default=19)
    parser.add_argument("--thickness-nm", type=float, default=75.0)
    parser.add_argument("--output", type=Path, default=Path("transmission.png"))
    parser.add_argument("--overwrite", action="store_true", help="Replace existing output files")
    args = parser.parse_args(argv)
    if args.wavelength_count < 2 or args.angle_count < 2:
        parser.error("sample counts must be >= 2")
    if args.output.suffix.lower() != ".png":
        parser.error("--output must have a .png extension")
    paths = [args.output, args.output.with_suffix(".csv"), args.output.with_suffix(".json")]
    if not args.overwrite and any(path.exists() for path in paths):
        parser.error("output already exists; choose another path or use --overwrite")
    try:
        _nonnegative_thickness(args.thickness_nm)
        if not args.output.parent.is_dir():
            raise ValueError("Output parent directory must already exist")
        wavelengths = np.linspace(400.0, 700.0, args.wavelength_count)
        angles = np.linspace(0.0, 89.0, args.angle_count)
        values = sweep(wavelengths, angles, args.thickness_nm)
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        try:
            plot = ax.contourf(angles, wavelengths, values, levels=30)
            fig.colorbar(plot, ax=ax, label="Unpolarized plane-wave transmittance")
            ax.set(
                xlabel="Incidence angle (degrees)",
                ylabel="Vacuum wavelength (nm)",
                title=f"Illustrative stack: organic {args.thickness_nm:g} nm",
            )
            fig.savefig(paths[0], dpi=150)
        finally:
            plt.close(fig)
        with paths[1].open("w", newline="", encoding="utf-8") as output:
            writer = csv.writer(output)
            writer.writerow(["wavelength_nm", "angle_deg", "unpolarized_transmittance"])
            for i, wavelength in enumerate(wavelengths):
                for j, angle in enumerate(angles):
                    writer.writerow([wavelength, angle, values[i, j]])
        metadata = {
            "schema_version": 1,
            "software_version": __version__,
            "quantity": "unpolarized_plane_wave_transmittance",
            "polarization_weighting": "0.5 * (T_s + T_p)",
            "stack": ["air", "illustrative organic", "illustrative silver", "air"],
            "organic_thickness_nm": args.thickness_nm,
            "silver_thickness_nm": 50.0,
            "wavelengths_nm": wavelengths.tolist(),
            "angles_deg": angles.tolist(),
            "limitations": "Illustrative indices; not calibrated materials or OLED dipole outcoupling.",
        }
        paths[2].write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    except (ValueError, ArithmeticError, OSError) as error:
        parser.exit(1, f"error: {error}; outputs may be incomplete\n")
    print("WARNING: illustrative plane-wave preview; NOT a validated OLED solver.", file=sys.stderr)
    print("Created: " + ", ".join(str(path) for path in paths))


if __name__ == "__main__":
    main()
