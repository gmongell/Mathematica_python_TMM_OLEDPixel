"""Dispersive optical-data catalog and rate-model EQE; no finite-stack emitter solver."""

import csv
from dataclasses import dataclass
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
from typing import Mapping
from urllib.request import urlopen

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import trapezoid
from tmm import coh_tmm
import yaml


def _grid(values: ArrayLike) -> NDArray[np.float64]:
    grid = np.asarray(values, dtype=float)
    if grid.ndim != 1 or grid.size == 0 or not np.all(np.isfinite(grid)) or np.any(grid <= 0):
        raise ValueError("Wavelengths must be a nonempty finite positive 1-D array (nm)")
    if np.any(np.diff(grid) <= 0):
        raise ValueError("Wavelengths must be strictly increasing")
    return grid


def _probability(name: str, value: float) -> None:
    if not np.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"{name} must be finite and in [0, 1]")


def _table(values: ArrayLike, columns: int) -> NDArray[np.float64]:
    data = np.asarray(values, dtype=float)
    if data.ndim != 2 or data.shape[1] != columns or len(data) < 2 or not np.all(np.isfinite(data)):
        raise ValueError("Optical table must have at least two finite rows and the required columns")
    _grid(data[:, 0])
    if np.any(data[:, 1:] < 0):
        raise ValueError("Passive optical tables require n >= 0 and k >= 0")
    return data


def _within(grid: NDArray[np.float64], lower: float, upper: float) -> None:
    if np.any(grid < lower) or np.any(grid > upper):
        raise ValueError(f"Wavelength outside dataset domain [{lower:g}, {upper:g}] nm; extrapolation forbidden")


@dataclass(frozen=True)
class Material:
    """A verified catalog record. Data are interpreted only through index()."""

    dataset_id: str
    metadata: Mapping
    blocks: tuple[dict, ...]
    assume_lossless: bool

    def index(self, wavelengths_nm: ArrayLike) -> NDArray[np.complex128]:
        grid = _grid(wavelengths_nm)
        n = k = None
        for block in self.blocks:
            kind = block["type"]
            if kind in ("tabulated nk", "tabulated n", "tabulated k"):
                table = _table(block["rows_nm"], 3 if kind == "tabulated nk" else 2)
                _within(grid, table[0, 0], table[-1, 0])
                if kind != "tabulated k":
                    if n is not None:
                        raise ValueError("Multiple refractive-index components in dataset")
                    n = np.interp(grid, table[:, 0], table[:, 1])
                if kind != "tabulated n":
                    if k is not None:
                        raise ValueError("Multiple extinction components in dataset")
                    k = np.interp(grid, table[:, 0], table[:, -1])
            elif kind in ("formula 1", "formula 2"):
                if n is not None:
                    raise ValueError("Multiple refractive-index components in dataset")
                lower, upper = np.asarray(block["range_nm"], dtype=float)
                if not np.isfinite([lower, upper]).all() or not 0 < lower < upper:
                    raise ValueError("Formula range must be finite, positive and increasing")
                _within(grid, lower, upper)
                coefficients = np.asarray(block["coefficients"], dtype=float)
                if coefficients.ndim != 1 or len(coefficients) < 3 or len(coefficients) % 2 != 1 or not np.isfinite(coefficients).all():
                    raise ValueError("Sellmeier coefficients must be a finite constant plus coefficient pairs")
                wavelength_squared = (grid / 1000) ** 2
                n_squared = np.full_like(grid, 1 + coefficients[0])
                for strength, pole in zip(coefficients[1::2], coefficients[2::2]):
                    denominator = wavelength_squared - (pole**2 if kind == "formula 1" else pole)
                    if np.any(denominator == 0):
                        raise ValueError("Wavelength coincides with a dispersion-formula pole")
                    n_squared += strength * wavelength_squared / denominator
                if np.any(n_squared <= 0) or not np.isfinite(n_squared).all():
                    raise ValueError("Dispersion formula produced invalid dielectric refractive index")
                n = np.sqrt(n_squared)
            else:
                raise ValueError(f"Unsupported optical-data block: {kind}")
        if n is None:
            raise ValueError("Dataset has no refractive-index component")
        if k is None:
            if not self.assume_lossless:
                raise ValueError("Extinction coefficient is missing; explicit assume_lossless=True is required")
            k = np.zeros_like(n)
        result = n + 1j * k
        if not np.isfinite(result).all():
            raise ValueError("Material evaluation produced nonfinite optical constants")
        return result


class OpticalDataLake:
    """Immutable content-addressed local catalog; simulations never fetch implicitly.

    The root may be a local directory or an already-mounted shared filesystem.
    This is not a cloud provisioning client or a database-wide mirror.
    """

    def __init__(self, root: str | Path):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _identifier(dataset_id: str) -> str:
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", dataset_id):
            raise ValueError("Dataset ID must be a simple catalog identifier, not a path")
        return dataset_id

    def list_datasets(self) -> list[str]:
        return sorted(path.stem for path in self.root.glob("*.json"))

    def _store(self, dataset_id: str, raw: bytes, blocks: list[dict], metadata: dict, assume_lossless: bool) -> Material:
        dataset_id = self._identifier(dataset_id)
        for field in ("source", "license", "citation"):
            if not isinstance(metadata.get(field), str) or not metadata[field].strip():
                raise ValueError(f"Material provenance requires {field}")
        if not isinstance(assume_lossless, bool):
            raise ValueError("assume_lossless must be explicitly boolean")
        if not blocks:
            raise ValueError("Dataset contains no supported optical blocks")
        for block in blocks:
            if block["type"].startswith("tabulated "):
                _table(block["rows_nm"], 3 if block["type"] == "tabulated nk" else 2)
            elif block["type"] not in ("formula 1", "formula 2"):
                raise ValueError(f"Unsupported optical-data block: {block['type']}")
        raw_sha = hashlib.sha256(raw).hexdigest()
        record = {
            "schema_version": 1, "dataset_id": dataset_id,
            "metadata": metadata, "blocks": blocks,
            "assume_lossless": assume_lossless, "raw_sha256": raw_sha,
        }
        payload = json.dumps(record, sort_keys=True, allow_nan=False).encode("utf-8")
        envelope = {"payload": record, "sha256": hashlib.sha256(payload).hexdigest()}
        path = self.root / f"{dataset_id}.json"
        encoded = json.dumps(envelope, indent=2, sort_keys=True, allow_nan=False).encode("utf-8")
        if path.exists():
            if path.read_bytes() != encoded:
                raise FileExistsError(f"Dataset ID already exists with different content: {dataset_id}")
            return self.load(dataset_id)
        raw_path = self.root / f"{raw_sha}.source"
        if raw_path.exists():
            if hashlib.sha256(raw_path.read_bytes()).hexdigest() != raw_sha:
                raise ValueError("Existing raw material object failed checksum validation")
        else:
            with raw_path.open("xb") as output:
                output.write(raw)
        with path.open("xb") as output:
            output.write(encoded)
        return self.load(dataset_id)

    def load(self, dataset_id: str) -> Material:
        dataset_id = self._identifier(dataset_id)
        envelope = json.loads((self.root / f"{dataset_id}.json").read_text(encoding="utf-8"))
        record = envelope["payload"]
        payload = json.dumps(record, sort_keys=True, allow_nan=False).encode("utf-8")
        if hashlib.sha256(payload).hexdigest() != envelope["sha256"]:
            raise ValueError(f"Catalog checksum mismatch: {dataset_id}")
        if record["schema_version"] != 1 or record["dataset_id"] != dataset_id:
            raise ValueError("Unsupported or mismatched material record")
        digest = record["raw_sha256"]
        if not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValueError("Invalid raw material digest")
        if hashlib.sha256((self.root / f"{digest}.source").read_bytes()).hexdigest() != digest:
            raise ValueError(f"Raw material checksum mismatch: {dataset_id}")
        return Material(dataset_id, record["metadata"], tuple(record["blocks"]), record["assume_lossless"])

    def import_csv(
        self, dataset_id: str, path: str | Path, *, source: str, license: str,
        citation: str, conditions: str,
    ) -> Material:
        """Import exact wavelength_nm,n,k headers, preserving original bytes."""
        raw = Path(path).read_bytes()
        reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
        if reader.fieldnames != ["wavelength_nm", "n", "k"]:
            raise ValueError("CSV columns must be exactly wavelength_nm,n,k")
        rows = [[float(row[key]) for key in reader.fieldnames] for row in reader]
        return self._store(dataset_id, raw, [{"type": "tabulated nk", "rows_nm": rows}], {
            "source": source, "license": license, "citation": citation, "conditions": conditions,
            "format": "CSV wavelength_nm,n,k",
        }, False)

    def import_rii(
        self, dataset_id: str, raw: bytes, *, source: str, revision: str,
        assume_lossless: bool = False,
    ) -> Material:
        """Import native refractiveindex.info YAML from an explicit source/revision."""
        if not source.strip() or not revision.strip():
            raise ValueError("Source and revision are required")
        document = yaml.safe_load(raw)
        if not isinstance(document, dict) or not isinstance(document.get("DATA"), list):
            raise ValueError("Expected native refractiveindex.info DATA blocks")
        blocks = []
        for item in document["DATA"]:
            kind = item["type"]
            if kind in ("tabulated nk", "tabulated n", "tabulated k"):
                columns = 3 if kind == "tabulated nk" else 2
                rows = np.loadtxt(io.StringIO(item["data"]), ndmin=2)
                rows = _table(rows, columns).copy()
                rows[:, 0] *= 1000  # Database wavelengths are micrometers, API wavelengths are nm.
                blocks.append({"type": kind, "rows_nm": rows.tolist()})
            elif kind in ("formula 1", "formula 2"):
                blocks.append({
                    "type": kind,
                    "range_nm": [1000 * float(value) for value in str(item["wavelength_range"]).split()],
                    "coefficients": [float(value) for value in str(item["coefficients"]).split()],
                })
            else:
                raise ValueError(f"Unsupported refractiveindex.info block: {kind}")
        return self._store(dataset_id, raw, blocks, {
            "source": source, "revision": revision, "license": "CC0-1.0",
            "citation": document.get("REFERENCES", ""),
            "conditions": document.get("CONDITIONS", {}),
            "comments": document.get("COMMENTS", ""),
            "database_citation": "M. N. Polyanskiy, Scientific Data 11, 94 (2024), doi:10.1038/s41597-023-02898-2",
        }, assume_lossless)

    def fetch_rii(
        self, dataset_id: str, database_path: str, commit: str, *,
        assume_lossless: bool = False,
    ) -> Material:
        """Explicit HTTPS fetch from the official database, pinned to a full Git commit."""
        self._identifier(dataset_id)
        if not re.fullmatch(r"[0-9a-f]{40}", commit):
            raise ValueError("Use an immutable full lowercase 40-character Git commit")
        path = PurePosixPath(database_path)
        if (not database_path.startswith("database/data/") or path.suffix != ".yml"
                or ".." in path.parts or "\\" in database_path
                or not re.fullmatch(r"[A-Za-z0-9_./-]+", database_path)):
            raise ValueError("Expected a relative database/data/.../*.yml source path")
        url = f"https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/{commit}/{database_path}"
        if (self.root / f"{dataset_id}.json").exists():
            material = self.load(dataset_id)
            if material.metadata["source"] != url or material.assume_lossless != assume_lossless:
                raise FileExistsError("Dataset ID is bound to another source or lossless assumption")
            return material
        with urlopen(url, timeout=30) as response:
            raw = response.read(4_000_001)
        if len(raw) > 4_000_000:
            raise ValueError("Material YAML exceeds the 4 MB per-dataset limit")
        return self.import_rii(dataset_id, raw, source=url, revision=commit, assume_lossless=assume_lossless)


@dataclass(frozen=True)
class Layer:
    material_id: str
    thickness_nm: float | None


@dataclass(frozen=True)
class DispersiveStack:
    """Ordered top-to-bottom stack with an explicit finite lossless emitter layer."""

    layers: tuple[Layer, ...]
    emitter_layer: int
    emitter_depth_nm: float

    def __post_init__(self):
        if len(self.layers) < 3 or self.layers[0].thickness_nm is not None or self.layers[-1].thickness_nm is not None:
            raise ValueError("Stack needs two explicit half-spaces and at least one finite layer")
        for layer in self.layers[1:-1]:
            if layer.thickness_nm is None or not np.isfinite(layer.thickness_nm) or layer.thickness_nm < 0:
                raise ValueError("Interior layer thicknesses must be finite and nonnegative (nm)")
        if isinstance(self.emitter_layer, bool) or not isinstance(self.emitter_layer, int) or not 0 < self.emitter_layer < len(self.layers) - 1:
            raise ValueError("Emitter must occupy an interior layer")
        thickness = self.layers[self.emitter_layer].thickness_nm
        if not np.isfinite(self.emitter_depth_nm) or not 0 < self.emitter_depth_nm < thickness:
            raise ValueError("Emitter depth must be strictly inside its layer (nm)")

    def evaluate(self, lake: OpticalDataLake, wavelengths_nm: ArrayLike) -> dict:
        grid = _grid(wavelengths_nm)
        materials = [lake.load(layer.material_id) for layer in self.layers]
        indices = np.column_stack([material.index(grid) for material in materials])
        if np.any(indices[:, self.emitter_layer].imag != 0) or np.any(indices[:, self.emitter_layer].real <= 0):
            raise ValueError("Emitting host must be lossless and positive index")
        if np.any(indices[:, [0, -1]].imag != 0) or np.any(indices[:, [0, -1]].real <= 0):
            raise ValueError("Exterior half-spaces must be lossless and positive index")
        return {
            "wavelengths_nm": grid,
            "indices": indices,
            "thicknesses_nm": np.array([np.inf if layer.thickness_nm is None else layer.thickness_nm for layer in self.layers]),
            "emitter_layer": self.emitter_layer,
            "emitter_depth_nm": self.emitter_depth_nm,
            "materials": [{"dataset_id": material.dataset_id, **dict(material.metadata),
                           "assume_lossless": material.assume_lossless} for material in materials],
        }

    def plane_wave_spectrum(
        self, lake: OpticalDataLake, wavelengths_nm: ArrayLike, *, angle_deg: float = 0,
    ) -> dict:
        """Bounded dispersive R/T/A integration check; NEVER an emitter EQE result."""
        if not np.isfinite(angle_deg) or not 0 <= angle_deg < 90:
            raise ValueError("Incidence angle must be finite and in [0, 90) degrees")
        evaluated = self.evaluate(lake, wavelengths_nm)
        reflection, transmission = [], []
        for wavelength, indices in zip(evaluated["wavelengths_nm"], evaluated["indices"]):
            results = [coh_tmm(p, indices, evaluated["thicknesses_nm"], np.deg2rad(angle_deg), wavelength)
                       for p in ("s", "p")]
            reflection.append(sum(result["R"] for result in results) / 2)
            transmission.append(sum(result["T"] for result in results) / 2)
        r, t = np.array(reflection), np.array(transmission)
        a = 1 - r - t
        if not np.isfinite([r, t, a]).all() or np.any(np.array([r, t, a]) < -1e-9):
            raise ArithmeticError("Dispersive plane-wave power balance is invalid")
        return {**evaluated, "quantity": "plane_wave_power_not_EQE", "R": r, "T": t, "A": a}


@dataclass(frozen=True)
class OpticalRateBudget:
    """Source-work F_total and collected emitted power F_escape, both / P_host.

    Arrays must come from a physically validated emitting-source calculation,
    not from incident-plane-wave transmission. No such finite-stack solver
    is manufactured by this data container.
    """

    wavelengths_nm: ArrayLike
    total: ArrayLike
    escape: ArrayLike
    reference: str

    def __post_init__(self):
        grid = _grid(self.wavelengths_nm).copy()
        total, escape = np.asarray(self.total, dtype=float), np.asarray(self.escape, dtype=float)
        if total.shape != grid.shape or escape.shape != grid.shape or not np.isfinite([total, escape]).all():
            raise ValueError("Optical rates must be finite arrays on the wavelength grid")
        if np.any(total < 0) or np.any(escape < 0) or np.any(escape > total):
            raise ValueError("Require 0 <= escape source-work rate <= total rate")
        if not isinstance(self.reference, str) or not self.reference.strip():
            raise ValueError("Optical rate normalization/source reference is required")
        total, escape = total.copy(), escape.copy()
        for array in (grid, total, escape):
            array.setflags(write=False)
        object.__setattr__(self, "wavelengths_nm", grid)
        object.__setattr__(self, "total", total)
        object.__setattr__(self, "escape", escape)

    @classmethod
    def orientation_average(
        cls, vertical: "OpticalRateBudget", horizontal: "OpticalRateBudget", vertical_fraction: float,
    ) -> "OpticalRateBudget":
        # This averages rates, not the quantum yields of independently decaying fixed dipoles.
        _probability("Vertical fraction", vertical_fraction)
        if not np.array_equal(vertical.wavelengths_nm, horizontal.wavelengths_nm):
            raise ValueError("Orientation budgets must use exactly the same wavelength grid")
        if vertical.reference != horizontal.reference:
            raise ValueError("Orientation budgets require the same homogeneous reference normalization")
        weight = vertical_fraction
        return cls(
            vertical.wavelengths_nm,
            weight * vertical.total + (1 - weight) * horizontal.total,
            weight * vertical.escape + (1 - weight) * horizontal.escape,
            vertical.reference,
        )

    def monochromatic_eqe(
        self, *, intrinsic_yield: float, charge_balance: float, exciton_yield: float,
    ) -> NDArray[np.float64]:
        """Independent monochromatic transitions; NOT a broadband spectral density."""
        for name, value in (("Intrinsic yield", intrinsic_yield), ("Charge balance", charge_balance), ("Exciton yield", exciton_yield)):
            _probability(name, value)
        if intrinsic_yield == 0:
            return np.zeros_like(self.total)
        denominator = (1 - intrinsic_yield) + intrinsic_yield * self.total
        if np.any(denominator <= 0):
            raise ValueError("Zero total decay rate: steady-state photon yield is undefined")
        return charge_balance * exciton_yield * intrinsic_yield * self.escape / denominator

    def broadband_eqe(
        self, intrinsic_photon_density_per_nm: ArrayLike, *, intrinsic_yield: float,
        charge_balance: float, exciton_yield: float,
    ) -> dict:
        """One emitter population: integrate rates BEFORE constructing spectral EQE."""
        for name, value in (("Intrinsic yield", intrinsic_yield), ("Charge balance", charge_balance), ("Exciton yield", exciton_yield)):
            _probability(name, value)
        spectrum = np.asarray(intrinsic_photon_density_per_nm, dtype=float)
        grid = self.wavelengths_nm
        if len(grid) < 2 or spectrum.shape != grid.shape or not np.isfinite(spectrum).all() or np.any(spectrum < 0):
            raise ValueError("Photon spectrum must be nonnegative and finite on at least two wavelengths")
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            try:
                normalization = trapezoid(spectrum, grid)
            except FloatingPointError as error:
                raise ValueError("Photon-spectrum integral overflowed") from error
        if not np.isclose(normalization, 1, rtol=1e-6, atol=1e-9):
            raise ValueError("Intrinsic PHOTON density must integrate to 1 per nm; normalization is explicit")
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            try:
                if intrinsic_yield == 0:
                    denominator, density = 1.0, np.zeros_like(spectrum)
                else:
                    denominator = (1 - intrinsic_yield) + intrinsic_yield * trapezoid(spectrum * self.total, grid)
                    if denominator <= 0:
                        raise ValueError("Zero integrated decay rate: spectral EQE is undefined")
                    density = charge_balance * exciton_yield * intrinsic_yield * spectrum * self.escape / denominator
                integrated_eqe = float(trapezoid(density, grid))
            except FloatingPointError as error:
                raise ValueError("Spectral optical-rate integration overflowed") from error
        if not np.isfinite(denominator) or not np.isfinite(density).all() or not np.isfinite(integrated_eqe):
            raise ValueError("Nonfinite spectral EQE calculation")
        return {
            "wavelengths_nm": grid.copy(), "eqe_density_per_nm": density,
            "integrated_eqe": integrated_eqe,
            "rate_denominator": float(denominator), "reference": self.reference,
        }


def random_dipoles(count: int, *, seed: int) -> NDArray[np.float64]:
    """Reproducible uniform orientations on the unit sphere, not coherent phases."""
    if isinstance(count, bool) or not isinstance(count, int) or count < 1:
        raise ValueError("Dipole count must be a positive integer")
    rng = np.random.default_rng(seed)
    cosine = rng.uniform(-1, 1, count)
    azimuth = rng.uniform(0, 2 * np.pi, count)
    sine = np.sqrt(1 - cosine**2)
    return np.column_stack((sine * np.cos(azimuth), sine * np.sin(azimuth), cosine))


def dipole_ensemble_eqe(
    vertical: OpticalRateBudget, horizontal: OpticalRateBudget, directions: ArrayLike, *,
    intrinsic_yield: float, charge_balance: float, exciton_yield: float,
) -> dict:
    """Compare independent fixed dipoles with a rapidly reorienting population.

    Planar isotropic media and azimuth-integrated collection only. Incoherent,
    equal excitation probability per emitter; no superradiance or phase array.
    Streams emitters to avoid allocating emitter-by-wavelength arrays.
    """
    vectors = np.asarray(directions, dtype=float)
    if vectors.ndim != 2 or vectors.shape[1] != 3 or len(vectors) == 0 or not np.isfinite(vectors).all():
        raise ValueError("Directions must be a nonempty finite N-by-3 array")
    norms = np.linalg.norm(vectors, axis=1)
    if np.any(norms == 0) or not np.isfinite(norms).all():
        raise ValueError("Dipole orientation vectors must have finite nonzero norm")
    weights = (vectors[:, 2] / norms) ** 2
    averaged = OpticalRateBudget.orientation_average(vertical, horizontal, float(np.mean(weights)))
    kwargs = dict(intrinsic_yield=intrinsic_yield, charge_balance=charge_balance, exciton_yield=exciton_yield)
    independent = np.zeros_like(averaged.total)
    for weight in weights:
        budget = OpticalRateBudget.orientation_average(vertical, horizontal, float(weight))
        independent += budget.monochromatic_eqe(**kwargs) / len(weights)
    return {
        "independent_fixed_dipoles_eqe": independent,
        "rapid_reorientation_eqe": averaged.monochromatic_eqe(**kwargs),
        "mean_vertical_fraction": float(np.mean(weights)),
        "dipole_count": len(weights),
        "assumptions": "Incoherent equal-excitation ensemble; planar isotropic media; azimuth-integrated collection",
    }
