"""Dipole source work and exit-channel power for an emitter inside a planar multilayer.

Convention: time dependence exp(-i omega t); p reflection/transmission are
magnetic-field amplitude coefficients (a perfect electric conductor has
r_p = +1), s coefficients are electric-field amplitudes. All rates are
normalized to the power the same dipole radiates in a lossless homogeneous
medium of the emitter-layer index.

Channels reported: source work into all electromagnetic channels (total),
flux into the top and bottom semi-infinite media, and the remainder
(guided modes + absorption + plasmons, not separated here).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Sequence

import numpy as np
from numpy.typing import NDArray

from dipole_reference import DipoleConvergenceError

_CHANNELS = ("total", "top", "bottom", "air_top", "air_bottom", "other")
_OUTPUTS = 10  # (total, top, bottom, air_top, air_bottom) x (vertical, parallel)


def _kz(epsilon: complex, kappa: NDArray) -> NDArray:
    """Normalized normal wavenumber with Im >= 0 (decaying/outgoing branch)."""
    root = np.sqrt(epsilon - kappa**2 + 0j)
    flip = (root.imag < 0) | ((root.imag == 0) & (root.real < 0))
    return np.where(flip, -root, root)


def _interface(polarization: str, eps_i: complex, k_i: NDArray, eps_t: complex, k_t: NDArray):
    """Reflection r and transmission t for a wave going from medium i into t."""
    if polarization == "s":
        return (k_i - k_t) / (k_i + k_t), 2 * k_i / (k_i + k_t)
    denominator = eps_t * k_i + eps_i * k_t
    return (eps_t * k_i - eps_i * k_t) / denominator, 2 * eps_t * k_i / denominator


def _outward(polarization: str, eps: Sequence[complex], kz: Sequence[NDArray],
             thickness_nm: Sequence[float], k0: float) -> tuple[NDArray, NDArray]:
    """Effective reflection/transmission seen from the last layer toward layer 0.

    Layers are ordered from the exit half-space to the emitter layer. The
    reflection is referenced to the emitter layer's interface; the
    transmission is the amplitude entering the exit half-space per unit
    up-going amplitude at that interface.
    """
    reflection, transmission = _interface(polarization, eps[1], kz[1], eps[0], kz[0])
    for k in range(2, len(eps)):
        phase = np.exp(2j * kz[k - 1] * k0 * thickness_nm[k - 1])
        r, t = _interface(polarization, eps[k], kz[k], eps[k - 1], kz[k - 1])
        denominator = 1 + r * reflection * phase
        transmission = t * np.exp(1j * kz[k - 1] * k0 * thickness_nm[k - 1]) * transmission / denominator
        reflection = (r + reflection * phase) / denominator
    return reflection, transmission


def _flux_weight(polarization: str, eps: complex, kz: NDArray) -> NDArray:
    return np.real(kz / eps) if polarization == "p" else np.real(kz)


@dataclass(frozen=True)
class EmitterRates:
    """Rates / P_host for vertical and in-plane (azimuth-averaged) dipoles.

    Tuple entries are (vertical, parallel). `air_*` is a first-pass lower
    bound for emission into air through a thick exit medium: it applies the
    exit-medium/air Fresnel power transmission once and ignores light
    returned into the stack. `other` = total - top - bottom is guided +
    absorbed power; it is not split further.
    """

    total: tuple[float, float]
    top: tuple[float, float]
    bottom: tuple[float, float]
    air_top: tuple[float, float]
    air_bottom: tuple[float, float]
    estimated_error: float
    evaluations: int

    @property
    def other(self) -> tuple[float, float]:
        return (self.total[0] - self.top[0] - self.bottom[0], self.total[1] - self.top[1] - self.bottom[1])

    def oriented(self, channel: str, vertical_fraction: float) -> float:
        if channel not in _CHANNELS:
            raise ValueError(f"channel must be one of {_CHANNELS}")
        if not np.isfinite(vertical_fraction) or not 0 <= vertical_fraction <= 1:
            raise ValueError("Vertical dipole fraction must be finite and in [0, 1]")
        vertical, parallel = getattr(self, channel)
        return vertical_fraction * vertical + (1 - vertical_fraction) * parallel


@dataclass(frozen=True)
class EmitterStack:
    """Top-to-bottom stack; interior layers are finite and the emitter sits at `emitter_depth_nm`
    below the top interface of `emitter_layer`. The emitter host must be lossless; other layers
    and half-spaces may be passive. A lossy half-space reports no exit flux (its power is "other")."""

    wavelength_nm: float
    indices: tuple[complex, ...]
    thicknesses_nm: tuple[float | None, ...]
    emitter_layer: int
    emitter_depth_nm: float

    def __post_init__(self) -> None:
        if not np.isfinite(self.wavelength_nm) or self.wavelength_nm <= 0:
            raise ValueError("Vacuum wavelength (nm) must be finite and positive")
        count = len(self.indices)
        if count < 3 or len(self.thicknesses_nm) != count:
            raise ValueError("Need two half-spaces, at least one interior layer, and matching thicknesses")
        if self.thicknesses_nm[0] is not None or self.thicknesses_nm[-1] is not None:
            raise ValueError("Half-space thicknesses must be None")
        for n, d in zip(self.indices, self.thicknesses_nm):
            value = complex(n)
            if not np.isfinite(value) or value.imag < 0 or value.real < 0 or value == 0:
                raise ValueError("Refractive indices must be finite and passive (Re, Im >= 0)")
            if d is not None and (not np.isfinite(d) or d <= 0):
                raise ValueError("Interior thicknesses must be finite and positive (nm)")
            if value.imag == 0 and value.real == 0:
                raise ValueError("Lossless zero index is unsupported")
        if isinstance(self.emitter_layer, bool) or not isinstance(self.emitter_layer, int) \
                or not 0 < self.emitter_layer < count - 1:
            raise ValueError("Emitter must occupy an interior layer")
        host = complex(self.indices[self.emitter_layer])
        if host.imag != 0 or host.real <= 0:
            raise ValueError("Emitting host must be lossless with positive index")
        for edge in (self.indices[0], self.indices[-1]):
            if complex(edge).real <= 0:
                raise ValueError("Exterior half-space index must have positive real part")
        for n in self.indices:
            if complex(n) ** 2 == (complex(n) ** 2).real and (complex(n) ** 2).real <= 0:
                raise ValueError("Unsupported lossless nonpositive permittivity")
        thickness = self.thicknesses_nm[self.emitter_layer]
        if not np.isfinite(self.emitter_depth_nm) or not 0 < self.emitter_depth_nm < thickness:
            raise ValueError("Emitter depth must be strictly inside its layer (nm)")

    @classmethod
    def from_dispersive(cls, stack, lake, wavelength_nm: float) -> "EmitterStack":
        """Evaluate a DispersiveStack and data-lake materials at one wavelength."""
        evaluated = stack.evaluate(lake, [wavelength_nm])
        return cls(
            wavelength_nm=float(evaluated["wavelengths_nm"][0]),
            indices=tuple(complex(n) for n in evaluated["indices"][0]),
            thicknesses_nm=tuple(None if np.isinf(d) else float(d) for d in evaluated["thicknesses_nm"]),
            emitter_layer=evaluated["emitter_layer"],
            emitter_depth_nm=evaluated["emitter_depth_nm"],
        )

    @property
    def host_index(self) -> float:
        return complex(self.indices[self.emitter_layer]).real

    @property
    def distances_nm(self) -> tuple[float, float]:
        """(distance to the top interface, distance to the bottom interface)."""
        thickness = self.thicknesses_nm[self.emitter_layer]
        return self.emitter_depth_nm, thickness - self.emitter_depth_nm

    def _breakpoints(self) -> list[float]:
        n_e = self.host_index
        points = {0.0, 1.0}
        for n in self.indices:
            points.add(complex(n).real / n_e)
        top = max(points)
        k_e = 2 * np.pi * n_e / self.wavelength_nm
        scale = 1 / (2 * k_e * min(self.distances_nm))
        points.update(top + scale * f for f in (0.25, 0.5, 1, 2, 3, 5, 7, 10, 14, 20))
        return sorted(p for p in points if p >= 0)

    def _resonances(self) -> tuple[list[float], float]:
        """Locate cavity-denominator minima |1 - rho| on a dense scan.

        Returns refined minimum positions (used as quadrature breakpoints)
        and the smallest minimum value found.
        """
        from scipy.optimize import minimize_scalar

        points = self._breakpoints()
        top = max([1.0] + [complex(n).real / self.host_index for n in self.indices]) * 1.5
        grid = np.unique(np.concatenate([np.linspace(1e-6, top, 20001), np.geomspace(top, points[-1], 4000)]))

        def denominator(u: NDArray) -> NDArray:
            return np.minimum(*(np.abs(1 - self._rho(pol, u)) for pol in ("s", "p")))

        values = denominator(grid)
        found, smallest = [], float("inf")
        for i in np.nonzero((values[1:-1] < values[:-2]) & (values[1:-1] <= values[2:]) & (values[1:-1] < 0.5))[0] + 1:
            result = minimize_scalar(lambda x: float(denominator(np.array([x]))[0]),
                                     bounds=(grid[i - 1], grid[i + 1]), method="bounded",
                                     options={"xatol": 1e-13})
            found.append(float(result.x))
            smallest = min(smallest, float(result.fun))
        return found, smallest

    def _rho(self, polarization: str, u: NDArray) -> NDArray:
        n_e = self.host_index
        k0 = 2 * np.pi / self.wavelength_nm
        kappa = n_e * u
        eps = [complex(n) ** 2 for n in self.indices]
        kz = [_kz(e, kappa) for e in eps]
        d1, d2 = self.distances_nm
        e = self.emitter_layer
        thicknesses = [0.0 if d is None else d for d in self.thicknesses_nm]
        up, down = list(range(0, e + 1)), list(range(len(eps) - 1, e - 1, -1))
        r_t, _ = _outward(polarization, [eps[i] for i in up], [kz[i] for i in up], [thicknesses[i] for i in up], k0)
        r_b, _ = _outward(polarization, [eps[i] for i in down], [kz[i] for i in down],
                          [thicknesses[i] for i in down], k0)
        return r_t * r_b * np.exp(2j * kz[e] * k0 * (d1 + d2))

    def _integrand(self, u: NDArray, air_index: float) -> NDArray:
        n_e = self.host_index
        k0 = 2 * np.pi / self.wavelength_nm
        kappa = n_e * u
        eps = [complex(n) ** 2 for n in self.indices]
        kz = [_kz(e, kappa) for e in eps]
        d1, d2 = self.distances_nm
        e_index = self.emitter_layer
        w = kz[e_index] / n_e
        e1 = np.exp(2j * kz[e_index] * k0 * d1)
        e2 = np.exp(2j * kz[e_index] * k0 * d2)
        phase1, phase2 = np.exp(1j * kz[e_index] * k0 * d1), np.exp(1j * kz[e_index] * k0 * d2)
        up = list(range(0, e_index + 1))
        down = list(range(len(eps) - 1, e_index - 1, -1))
        thicknesses = [0.0 if d is None else d for d in self.thicknesses_nm]
        out = np.zeros((_OUTPUTS, u.size))
        air_eps = air_index**2
        k_air = _kz(air_eps, kappa)
        weights = {"vp": 0.75 * n_e * u**3 / np.abs(w) ** 2, "ss": 0.375 * u / (n_e * np.abs(w) ** 2),
                   "sp": 0.375 * n_e * u}
        parts = {}
        for pol in ("s", "p"):
            r_t, t_t = _outward(pol, [eps[i] for i in up], [kz[i] for i in up],
                                [thicknesses[i] for i in up], k0)
            r_b, t_b = _outward(pol, [eps[i] for i in down], [kz[i] for i in down],
                                [thicknesses[i] for i in down], k0)
            parts[pol] = (r_t, t_t, r_b, t_b, r_t * r_b * e1 * e2)
        z = [(0, "p", "vp", +1), (1, "s", "ss", +1), (1, "p", "sp", -1)]
        for column, pol, key, alpha in z:
            r_t, t_t, r_b, t_b, rho = parts[pol]
            cavity = (r_t * e1 + r_b * e2 + 2 * alpha * rho) / (1 - rho)
            scattered = cavity if alpha > 0 else -cavity
            if key == "vp":
                total = 1.5 * u**3 / w * (1 + scattered)
            elif key == "ss":
                total = 0.75 * u / w * (1 + scattered)
            else:
                total = 0.75 * u / w * w**2 * (1 + scattered)
            out[column] += np.real(total)
            tau_top = t_t * phase1 * (1 + alpha * r_b * e2) / (1 - rho)
            tau_bottom = t_b * phase2 * (alpha + r_t * e1) / (1 - rho)
            for exit_index, tau, base in ((0, tau_top, 2), (-1, tau_bottom, 4)):
                lossy = eps[exit_index].imag != 0
                y_exit = np.zeros_like(u) if lossy else _flux_weight(pol, eps[exit_index], kz[exit_index])
                flux = weights[key] * y_exit * np.abs(tau) ** 2
                out[base + column] += flux
                r_air, t_air = _interface(pol, eps[exit_index] + 0j, kz[exit_index], air_eps, k_air)
                y_air = _flux_weight(pol, air_eps, k_air)
                with np.errstate(divide="ignore", invalid="ignore"):
                    transmittance = np.where(y_exit > 0, y_air * np.abs(t_air) ** 2 / y_exit, 0.0)
                out[base + 4 + column] += flux * transmittance
        return out

    def _panel_integrals(self, spans: NDArray, air_index: float, order: int = 10):
        """Gauss-Legendre integrals over (panel_lo, panel_hi, theta_lo, theta_hi) rows.

        The map u = a + (b-a)(1 - cos(theta))/2 removes inverse-square-root
        endpoint singularities at light lines and interface branch points.
        """
        x, wx = np.polynomial.legendre.leggauss(order)
        lo, hi = spans[:, 2:3], spans[:, 3:4]
        theta = (lo + hi) / 2 + (hi - lo) / 2 * x
        a, b = spans[:, 0:1], spans[:, 1:2]
        nodes = a + (b - a) * (1 - np.cos(theta)) / 2
        weights = (hi - lo) / 2 * wx * (b - a) / 2 * np.sin(theta)
        values = self._integrand(nodes.ravel(), air_index).reshape(_OUTPUTS, *nodes.shape)
        return np.einsum("cij,ij->ic", values, weights), nodes.size

    def rates(self, *, air_index: float = 1.0, atol: float = 1e-7, rtol: float = 1e-6,
              max_intervals: int = 6000) -> EmitterRates:
        """Adaptively integrate over in-plane wavenumber.

        Raises DipoleConvergenceError when the summed bisection error estimate
        does not meet the tolerance (for example a lossless guided mode, whose
        pole lies on the integration axis: give guiding layers a small passive
        Im(n)). The estimate is not a rigorous bound for arbitrarily narrow resonances.
        """
        if not np.isfinite(air_index) or air_index <= 0:
            raise ValueError("Air index must be finite and positive")
        if not 1e-12 <= atol <= 1e-2 or not 1e-12 <= rtol <= 1e-2:
            raise ValueError("atol and rtol must be finite and in [1e-12, 1e-2]")
        if isinstance(max_intervals, bool) or not isinstance(max_intervals, int) or not 100 <= max_intervals <= 100000:
            raise ValueError("max_intervals must be an integer in [100, 100000]")
        if all(complex(n).imag == 0 for n in self.indices) and \
                max(complex(n).real for n in self.indices[1:-1]) > max(complex(self.indices[0]).real, complex(self.indices[-1]).real):
            raise DipoleConvergenceError("Lossless guiding layer: guided-mode delta power is not integrable; add a small passive Im(n)")
        resonances, smallest = self._resonances()
        if smallest < 1e-6:
            raise DipoleConvergenceError(
                f"Cavity denominator reaches {smallest:.2g}: lossless guided pole; add a small passive Im(n)")
        points = sorted(set(self._breakpoints()) | set(resonances))
        edges = np.linspace(0, np.pi, 5)
        spans = np.array([[a, b, lo, hi] for a, b in zip(points[:-1], points[1:])
                          for lo, hi in zip(edges[:-1], edges[1:])])
        evaluations = 0
        while True:
            mid = (spans[:, 2] + spans[:, 3]) / 2
            left, right = spans.copy(), spans.copy()
            left[:, 3], right[:, 2] = mid, mid
            whole, n1 = self._panel_integrals(spans, air_index)
            half_l, n2 = self._panel_integrals(left, air_index)
            half_r, n3 = self._panel_integrals(right, air_index)
            evaluations += n1 + n2 + n3
            refined = half_l + half_r
            error = np.abs(refined - whole)
            value = refined.sum(axis=0)
            if not np.isfinite(value).all():
                raise DipoleConvergenceError("Non-finite multilayer integral")
            tolerance = np.maximum(atol, rtol * np.abs(value))
            if np.all(error.sum(axis=0) <= tolerance):
                break
            if 2 * len(spans) > max_intervals:
                raise DipoleConvergenceError(
                    f"Adaptive error {np.max(error.sum(axis=0) / tolerance):.3g} x tolerance with {len(spans)} "
                    "intervals; possible lossless guided pole or unresolved resonance")
            score = np.max(error / tolerance, axis=1)
            split = score > 0.5 / len(spans)
            spans = np.concatenate([left[split], right[split], spans[~split]])
        pairs = [(float(value[i]), float(value[i + 1])) for i in range(0, _OUTPUTS, 2)]
        return EmitterRates(*pairs, estimated_error=float(error.sum(axis=0).max()), evaluations=evaluations)

def emitter_rate_budget(
    stack, lake, wavelengths_nm, *, exit: Literal["top", "bottom"], collection: Literal["exit_medium", "air_first_pass"],
    vertical_fraction: float, air_index: float = 1.0,
):
    """Build an OpticalRateBudget by solving the emitter at every wavelength.

    `vertical_fraction` weights vertical against in-plane dipole rates
    (a fixed orientation order parameter, not a re-orienting ensemble).
    """
    from spectral_eqe import OpticalRateBudget

    if exit not in ("top", "bottom") or collection not in ("exit_medium", "air_first_pass"):
        raise ValueError("exit must be top/bottom and collection exit_medium/air_first_pass")
    channel = exit if collection == "exit_medium" else f"air_{exit}"
    grid = np.asarray(stack.evaluate(lake, wavelengths_nm)["wavelengths_nm"], dtype=float)
    total, escape = [], []
    for wavelength in grid:
        rates = EmitterStack.from_dispersive(stack, lake, float(wavelength)).rates(air_index=air_index)
        total.append(rates.oriented("total", vertical_fraction))
        escape.append(rates.oriented(channel, vertical_fraction))
    total_array, escape_array = np.array(total), np.array(escape)
    escape_array = np.where((escape_array > total_array) & (escape_array <= total_array * (1 + 1e-6) + 1e-9),
                            total_array, escape_array)
    return OpticalRateBudget(
        wavelengths_nm=grid, total=total_array, escape=escape_array,
        reference=f"multilayer dipole solver; exit={exit}; collection={collection}; vertical_fraction={vertical_fraction}",
    )
