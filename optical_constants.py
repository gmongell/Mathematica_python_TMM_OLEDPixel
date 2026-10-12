"""Load tabulated wavelength-dependent complex optical constants.

CSV schema: wavelength_nm,n,k with strictly increasing wavelengths.
The convention is N=n+ik, k>=0 for passive materials with exp(-i omega t).
Linear interpolation is confined to the measured range: extrapolation fails.
No material values are bundled or fabricated by this module.
"""
import csv
import numpy as np


class OpticalConstants:
    def __init__(self, wavelength_nm, n, k):
        w = np.asarray(wavelength_nm, dtype=float)
        nr = np.asarray(n, dtype=float)
        ki = np.asarray(k, dtype=float)
        if any(a.ndim != 1 for a in (w, nr, ki)):
            raise ValueError("Optical-constant arrays must be one-dimensional")
        if len(w) < 2 or len(nr) != len(w) or len(ki) != len(w):
            raise ValueError("At least two matched wavelength,n,k records required")
        if not all(np.all(np.isfinite(a)) for a in (w, nr, ki)):
            raise ValueError("Optical constants must be finite")
        if np.any(w <= 0) or np.any(np.diff(w) <= 0):
            raise ValueError("Wavelengths must be positive and strictly increasing")
        if np.any(nr <= 0) or np.any(ki < 0):
            raise ValueError("Require n>0 and passive extinction coefficient k>=0")
        self.wavelength_nm = w.copy()
        self.n = nr.copy()
        self.k = ki.copy()

    @classmethod
    def from_csv(cls, path):
        with open(path, newline="", encoding="utf-8-sig") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None or not {"wavelength_nm", "n", "k"}.issubset(reader.fieldnames):
                raise ValueError("CSV must contain wavelength_nm,n,k columns")
            rows = list(reader)
        try:
            return cls([float(r["wavelength_nm"]) for r in rows],
                       [float(r["n"]) for r in rows],
                       [float(r["k"]) for r in rows])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("Invalid optical-constant CSV records") from exc

    def evaluate(self, wavelength_nm):
        w = np.asarray(wavelength_nm, dtype=float)
        if not np.all(np.isfinite(w)) or np.any(w < self.wavelength_nm[0]) or np.any(w > self.wavelength_nm[-1]):
            raise ValueError("Requested wavelengths must lie within tabulated range")
        value = np.interp(w, self.wavelength_nm, self.n) + 1j*np.interp(w, self.wavelength_nm, self.k)
        return complex(value) if w.ndim == 0 else value
