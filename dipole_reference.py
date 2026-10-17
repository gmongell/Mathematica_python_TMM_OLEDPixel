"""Research dipole reference kernels; not a multilayer OLED outcoupling solver."""

from dataclasses import dataclass
from numbers import Real
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.constants import c, epsilon_0
from scipy.integrate import quad

Orientation = Literal["vertical", "horizontal_x", "horizontal_y", "horizontal", "isotropic"]
Substrate = complex | Literal["pec"]


class DipoleConvergenceError(ArithmeticError):
    """The requested source-work accuracy could not be established."""


def _positive(name: str, value: float) -> None:
    if not isinstance(value, Real) or not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")


@dataclass(frozen=True)
class AngularPattern:
    """Power per steradian divided by homogeneous-host total power."""

    s: NDArray[np.float64]
    p: NDArray[np.float64]

    @property
    def total(self) -> NDArray[np.float64]:
        return self.s + self.p


def homogeneous_pattern(
    theta_rad: ArrayLike, phi_rad: ArrayLike, orientation: Orientation
) -> AngularPattern:
    """Broadcast polar/azimuth angles; theta is in [0, pi], phi is periodic."""
    theta, phi = np.broadcast_arrays(
        np.asarray(theta_rad, dtype=float), np.asarray(phi_rad, dtype=float)
    )
    if not np.all(np.isfinite(theta)) or np.any(theta < 0) or np.any(theta > np.pi):
        raise ValueError("Polar angles must be finite and in [0, pi] radians")
    if not np.all(np.isfinite(phi)):
        raise ValueError("Azimuth angles must be finite radians")
    if orientation == "vertical":
        s, p = np.zeros_like(theta), np.sin(theta) ** 2
    elif orientation == "horizontal_x":
        s, p = np.sin(phi) ** 2, np.cos(theta) ** 2 * np.cos(phi) ** 2
    elif orientation == "horizontal_y":
        s, p = np.cos(phi) ** 2, np.cos(theta) ** 2 * np.sin(phi) ** 2
    elif orientation == "horizontal":
        s, p = np.full_like(theta, 0.5), 0.5 * np.cos(theta) ** 2
    elif orientation == "isotropic":
        s, p = np.full_like(theta, 1 / 3), np.full_like(theta, 1 / 3)
    else:
        raise ValueError("Unknown dipole orientation")
    return AngularPattern(3 / (8 * np.pi) * s, 3 / (8 * np.pi) * p)


def homogeneous_power(wavelength_nm: float, host_index: float, dipole_c_m: complex) -> float:
    """Macroscopic SI watts for peak harmonic dipole amplitude; no local-field factor."""
    _positive("Vacuum wavelength (nm)", wavelength_nm)
    _positive("Host refractive index", host_index)
    if not np.isfinite(dipole_c_m):
        raise ValueError("Dipole amplitude must be finite (C m)")
    with np.errstate(over="raise", divide="raise", invalid="raise"):
        try:
            omega = 2 * np.pi * c / (np.float64(wavelength_nm) * 1e-9)
            power = host_index * omega**4 * np.abs(np.complex128(dipole_c_m)) ** 2 / (12 * np.pi * epsilon_0 * c**3)
        except FloatingPointError as error:
            raise ValueError("Homogeneous power overflowed; check input units") from error
    if not np.isfinite(power):
        raise ValueError("Homogeneous power overflowed; check input units")
    return float(power)


@dataclass(frozen=True)
class DecayResult:
    """Total electromagnetic source-work ratios, not escape efficiencies.

    Tuple diagnostics are ordered (vertical, parallel). Corrections can be
    negative and are integration regions, not individually positive channels.
    """

    vertical: float
    parallel: float
    propagating_correction: tuple[float, float]
    evanescent_correction: tuple[float, float]
    estimated_error: tuple[float, float]
    refinement_delta: tuple[float, float]
    tail_magnitude: tuple[float, float]
    evaluations: int

    def oriented(self, vertical_fraction: float) -> float:
        if not np.isfinite(vertical_fraction) or not 0 <= vertical_fraction <= 1:
            raise ValueError("Vertical dipole fraction must be finite and in [0, 1]")
        return vertical_fraction * self.vertical + (1 - vertical_fraction) * self.parallel

    @property
    def isotropic(self) -> float:
        return self.oriented(1 / 3)


@dataclass(frozen=True)
class SingleInterface:
    """Dipole at z=d>0 above a substrate at z<0; host is lossless."""

    wavelength_nm: float
    distance_nm: float
    substrate_epsilon: Substrate
    host_index: float = 1.0

    def __post_init__(self) -> None:
        _positive("Vacuum wavelength (nm)", self.wavelength_nm)
        _positive("Dipole-interface distance (nm)", self.distance_nm)
        _positive("Lossless host refractive index", self.host_index)
        if not np.isfinite(self.phase_distance) or self.phase_distance < 1e-6:
            raise ValueError("Dimensionless separation 2*k_host*d must be finite and >= 1e-6")
        if isinstance(self.substrate_epsilon, str):
            if self.substrate_epsilon != "pec":
                raise ValueError("The only ideal-boundary sentinel is 'pec'")
        else:
            epsilon = complex(self.substrate_epsilon)
            if not np.isfinite(epsilon) or epsilon.imag < 0:
                raise ValueError("Substrate permittivity must be finite and passive (Im >= 0)")
            if epsilon.imag == 0 and epsilon.real <= 0:
                raise ValueError("Lossless nonpositive-epsilon substrates require unsupported pole/ENZ treatment")

    @property
    def phase_distance(self) -> float:
        return 4 * np.pi * self.host_index * self.distance_nm / self.wavelength_nm

    def _reflection(self, w: complex) -> tuple[complex, complex]:
        if self.substrate_epsilon == "pec":
            return -1 + 0j, 1 + 0j
        relative = complex(self.substrate_epsilon) / self.host_index**2
        if relative == 1:
            return 0j, 0j
        v = np.sqrt(relative - 1 + w * w + 0j)
        return (w - v) / (w + v), (relative * w - v) / (relative * w + v)

    def spectral_correction(self, q: float) -> tuple[float, float]:
        """Vertical/parallel d(F-1)/dq; rejects the singular host light line q=1."""
        if not np.isfinite(q) or q < 0 or q == 1:
            raise ValueError("q=k_parallel/k_host must be finite, nonnegative, and not 1")
        w = np.sqrt(1 - q * q + 0j)
        rs, rp = self._reflection(w)
        phase = np.exp(1j * self.phase_distance * w)
        values = (
            1.5 * np.real(q**3 / w * rp * phase),
            0.75 * np.real(q / w * (rs - w * w * rp) * phase),
        )
        if not np.all(np.isfinite(values)):
            raise DipoleConvergenceError("Nonfinite raw spectral kernel; use a supported q range")
        return float(values[0]), float(values[1])

    def _breakpoints(self) -> tuple[list[float], list[float]]:
        prop: list[float] = []
        evan: list[float] = []
        if self.substrate_epsilon == "pec":
            return prop, evan
        relative = complex(self.substrate_epsilon) / self.host_index**2
        if relative.imag == 0:
            if 0 < relative.real < 1:
                prop.append(float(np.sqrt(1 - relative.real)))
            elif relative.real > 1:
                evan.append(float(self.phase_distance * np.sqrt(relative.real - 1)))
        elif relative.real < -1:
            pole = self.phase_distance * np.sqrt(-1 / (relative + 1))
            center, width = float(pole.real), abs(float(pole.imag))
            offsets = (-32768, -4096, -512, -64, -8, -1, 0, 1, 8, 64, 512, 4096, 32768)
            evan.extend(center + offset * width for offset in offsets)
        elif relative.imag > 0:
            # Resolve a weakly lossy dielectric's smoothed substrate light line.
            if relative.real < 1:
                branch = np.sqrt(1 - relative)
                prop.extend(
                    point for offset in (-8, -1, 0, 1, 8)
                    if 0 < (point := float(branch.real + offset * abs(branch.imag))) < 1
                )
            elif relative.real > 1:
                branch = self.phase_distance * np.sqrt(relative - 1)
                evan.extend(float(branch.real + offset * abs(branch.imag)) for offset in (-8, -1, 0, 1, 8))
        return sorted(set(prop)), sorted(set(point for point in evan if point > 0))

    def decay_rates(
        self, *, atol: float = 1e-8, rtol: float = 1e-7, max_subdivisions: int = 300
    ) -> DecayResult:
        """Integrate source work, repeating transformed infinite-tail quadrature.

        Raises DipoleConvergenceError on failed error estimates/refinement.
        No quadrature estimate is a rigorous bound on all possible resonances.
        """
        if not 1e-12 <= atol <= 1e-2 or not 1e-12 <= rtol <= 1e-2:
            raise ValueError("atol and rtol must be finite and in [1e-12, 1e-2]")
        if isinstance(max_subdivisions, bool) or not isinstance(max_subdivisions, int) or not 1 <= max_subdivisions <= 10000:
            raise ValueError("max_subdivisions must be an integer in [1, 10000]")
        prop_points, evan_points = self._breakpoints()
        a = self.phase_distance
        evaluations = 0

        def integrate(function, boundaries, tolerance):
            nonlocal evaluations
            total = error = 0.0
            budget = atol * tolerance / (len(boundaries) - 1)
            for left, right in zip(boundaries[:-1], boundaries[1:]):
                outcome = quad(
                    function, left, right, epsabs=budget,
                    epsrel=rtol * tolerance, limit=max_subdivisions, full_output=1,
                )
                value, estimate, info = outcome[:3]
                evaluations += info["neval"]
                if len(outcome) != 3 or not np.isfinite(value) or not np.isfinite(estimate):
                    message = outcome[3] if len(outcome) > 3 else "nonfinite quadrature"
                    raise DipoleConvergenceError(f"Single-interface integration failed: {message}")
                total += value
                error += estimate
            return total, error

        corrections, errors, deltas, tails = [], [], [], []
        for vertical in (True, False):
            factor = 1.5 if vertical else 0.75

            def prop(u):
                rs, rp = self._reflection(complex(u))
                amplitude = (1 - u * u) * rp if vertical else rs - u * u * rp
                return float(factor * np.real(amplitude * np.exp(1j * a * u)))

            def evan(x):
                if x > 700:
                    return 0.0  # Beyond floating-point exponential support at the requested accuracy.
                t = x / a
                rs, rp = self._reflection(1j * t)
                amplitude = (1 + t * t) * rp if vertical else rs + t * t * rp
                return float(factor / a * np.exp(-x) * np.imag(amplitude))

            passes = []
            for split, tolerance in ((32.0, 0.05), (64.0, 0.0125)):
                p, ep = integrate(prop, [0.0, *prop_points, 1.0], tolerance)
                points = sorted(set([0.0, split, *evan_points]))
                e, ee = integrate(evan, [*points, np.inf], tolerance)
                tail, _ = integrate(evan, [points[-1], np.inf], tolerance)
                passes.append((p, e, ep + ee, abs(tail)))
            p, e, error, tail = passes[-1]
            rate = 1 + p + e
            delta = abs(p - passes[0][0]) + abs(e - passes[0][1])
            budget = max(atol, rtol * abs(rate))
            if not np.isfinite(rate) or rate < 0 or max(error, delta) > budget:
                raise DipoleConvergenceError(
                    f"Source-work accuracy not established: rate={rate}, error={error}, refinement={delta}"
                )
            corrections.append((p, e))
            errors.append(error)
            deltas.append(delta)
            tails.append(tail)
        return DecayResult(
            vertical=1 + sum(corrections[0]),
            parallel=1 + sum(corrections[1]),
            propagating_correction=(corrections[0][0], corrections[1][0]),
            evanescent_correction=(corrections[0][1], corrections[1][1]),
            estimated_error=(errors[0], errors[1]),
            refinement_delta=(deltas[0], deltas[1]),
            tail_magnitude=(tails[0], tails[1]),
            evaluations=evaluations,
        )
