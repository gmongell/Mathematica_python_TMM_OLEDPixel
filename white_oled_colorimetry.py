"""Angle-resolved white-light chromaticity and uniformity utilities.

Supply measured spectral radiant intensity per steradian on a wavelength/angle
grid and tabulated CIE 1931 2-degree xbar,ybar,zbar values. No CMF data is
fabricated or embedded. This module does not solve dipole electrodynamics.

CSV format: wavelength_nm,xbar,ybar,zbar (strictly increasing wavelengths).
All wavelengths must match the intensity wavelength grid. No interpolation is
performed silently.
"""
import csv
import numpy as np


def read_cmf_csv(path):
    with open(path, newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) < 2:
        raise ValueError("At least two CMF rows required")
    columns = ("wavelength_nm", "xbar", "ybar", "zbar")
    values = np.array([[float(row[k]) for k in columns] for row in rows])
    if not np.all(np.isfinite(values)) or np.any(np.diff(values[:, 0]) <= 0):
        raise ValueError("CMF wavelengths must increase and values be finite")
    if np.any(values[:, 1:] < 0):
        raise ValueError("CMF values must be nonnegative")
    return values[:, 0], values[:, 1:]


def calculate_angular_color(wavelengths_nm, angles_deg, radiant_intensity, cmf,
                            reference_angle_index=0):
    """Return XYZ, xy, u'v' and max Δu'v' relative to a reference angle.

    radiant_intensity shape: (n_wavelength, n_angle). XYZ are relative
    tristimulus values in arbitrary units; no absolute photometric calibration.
    """
    w = np.asarray(wavelengths_nm, dtype=float)
    a = np.asarray(angles_deg, dtype=float)
    s = np.asarray(radiant_intensity, dtype=float)
    c = np.asarray(cmf, dtype=float)
    if w.ndim != 1 or a.ndim != 1 or len(w) < 2 or len(a) < 1:
        raise ValueError("Expected 1D wavelength and angle grids")
    if s.shape != (len(w), len(a)) or c.shape != (len(w), 3):
        raise ValueError("Intensity/CMF array shapes do not match grids")
    if not np.all(np.isfinite(w)) or not np.all(np.isfinite(a)):
        raise ValueError("Grids must be finite")
    if np.any(np.diff(w) <= 0) or np.any(np.diff(a) <= 0):
        raise ValueError("Grids must be strictly increasing")
    if not np.all(np.isfinite(s)) or np.any(s < 0):
        raise ValueError("Radiant intensity must be finite and nonnegative")
    if not np.all(np.isfinite(c)) or np.any(c < 0):
        raise ValueError("Color matching functions must be finite and nonnegative")
    if not 0 <= reference_angle_index < len(a):
        raise ValueError("Reference angle index is out of bounds")
    xyz = np.trapz(s[:, :, None] * c[:, None, :], w, axis=0)
    xyz_sum = xyz.sum(axis=1)
    uv_den = xyz[:, 0] + 15 * xyz[:, 1] + 3 * xyz[:, 2]
    if np.any(xyz_sum <= 0) or np.any(uv_den <= 0):
        raise ValueError("Zero/undefined tristimulus at one or more angles")
    xy = xyz[:, :2] / xyz_sum[:, None]
    uv = np.column_stack((4 * xyz[:, 0] / uv_den, 9 * xyz[:, 1] / uv_den))
    delta_uv = np.linalg.norm(uv - uv[reference_angle_index], axis=1)
    return {"XYZ": xyz.tolist(), "xy": xy.tolist(), "uv_prime": uv.tolist(),
            "delta_uv": delta_uv.tolist(),
            "max_delta_uv": float(np.max(delta_uv)),
            "reference_angle_deg": float(a[reference_angle_index])}
