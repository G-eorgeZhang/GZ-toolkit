"""Alloy EAM, FS and signed-density HE tables and a small-cell evaluator.

FS/HE densities depend on both receiver and donor; alloy EAM only on donor.
This evaluator includes periodic images (including self images), not just
minimum-image pairs. It is for fitting small cells, not production dynamics.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from itertools import product
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline

from .data import Configuration


def pair_key(a: str, b: str) -> str:
    return "-".join(sorted((a, b)))


def density_key(receiver: str, donor: str) -> str:
    return f"{receiver}<-{donor}"


@dataclass
class Element:
    symbol: str
    atomic_number: int
    mass: float
    lattice_constant: float
    lattice: str = "BCC"


@dataclass
class Prediction:
    energy: float
    forces: np.ndarray
    stress: np.ndarray
    densities: np.ndarray


@dataclass
class EAMModel:
    elements: list[Element]
    radial_grid: list[float]
    density_grid: list[float]
    embedding: dict[str, list[float]]
    density: dict[str, list[float]]
    pair: dict[str, list[float]]
    metadata: dict
    style: str = "eam/alloy"
    cutoff_radius: float | None = None

    @property
    def symbols(self):
        return [e.symbol for e in self.elements]

    @property
    def cutoff(self):
        return float(self.radial_grid[-1] if self.cutoff_radius is None else self.cutoff_radius)

    def validate(self):
        if not self.elements or len(set(self.symbols)) != len(self.symbols):
            raise ValueError("Model elements must be unique and nonempty.")
        for e in self.elements:
            if not e.symbol.isalpha() or e.atomic_number < 1 or not np.isfinite(e.mass) or e.mass <= 0:
                raise ValueError("Invalid element metadata.")
            if not e.lattice or any(ch.isspace() for ch in e.lattice) or not np.isfinite(e.lattice_constant) or e.lattice_constant < 0:
                raise ValueError("Invalid lattice metadata.")
        if self.style not in {"eam/alloy", "eam/fs", "eam/he"}:
            raise ValueError("Supported styles: eam/alloy, eam/fs, eam/he.")
        for index, grid in enumerate((self.radial_grid, self.density_grid)):
            g = np.asarray(grid)
            if len(g) < 4 or not np.isfinite(g).all() or np.any(np.diff(g) <= 0):
                raise ValueError("Grids need >=4 finite increasing values.")
            if (index == 0 or self.style != "eam/he") and abs(g[0]) > 1e-10:
                raise ValueError("Only eam/he embedding grids may start below/above zero.")
            if not np.allclose(np.diff(g), g[1] - g[0], rtol=1e-8, atol=1e-12):
                raise ValueError("setfl grids must be uniformly spaced.")
        if not np.isfinite(self.cutoff) or self.cutoff <= 0 or self.cutoff > self.radial_grid[-1] + (self.radial_grid[1] - self.radial_grid[0]) + 1e-8:
            raise ValueError("Cutoff must lie within the radial table or at most one interval beyond it.")
        pairs = {pair_key(a, b) for a in self.symbols for b in self.symbols}
        density_keys = set(self.symbols) if self.style == "eam/alloy" else {density_key(a, b) for a in self.symbols for b in self.symbols}
        for tables, keys, n in ((self.embedding, set(self.symbols), len(self.density_grid)),
                                (self.density, density_keys, len(self.radial_grid)),
                                (self.pair, pairs, len(self.radial_grid))):
            if set(tables) != keys:
                raise ValueError("Missing or extra EAM functions.")
            for v in tables.values():
                if len(v) != n or not np.isfinite(v).all():
                    raise ValueError("Invalid table values.")
        if self.style != "eam/he" and any(np.any(np.asarray(v) < 0) for v in self.density.values()):
            raise ValueError("Signed densities require eam/he.")

    def save(self, path):
        self.validate()
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2, allow_nan=False), encoding="utf-8")
        return path

    @classmethod
    def load(cls, path):
        d = json.loads(Path(path).read_text(encoding="utf-8"))
        d["elements"] = [Element(**e) for e in d["elements"]]
        model = cls(**d)
        model.validate()
        return model

    def _splines(self):
        embedding = {s: CubicSpline(self.density_grid, v, bc_type="natural") for s, v in self.embedding.items()}
        density = {s: CubicSpline(self.radial_grid, v, bc_type="natural") for s, v in self.density.items()}
        pair = {}
        for s, v in self.pair.items():
            rp = np.asarray(self.radial_grid) * v
            rp[0] = self.metadata.get("rphi_at_zero", {}).get(s, 0.0)
            pair[s] = CubicSpline(self.radial_grid, rp, bc_type="natural")
        return embedding, density, pair

    def neighbors(self, c: Configuration):
        c.validate()
        if not set(c.species).issubset(self.symbols):
            raise ValueError("Configuration contains elements outside the model.")
        cell = np.asarray(c.cell, float)
        inv = np.linalg.inv(cell)
        frac = np.asarray(c.positions) @ inv
        for axis, periodic in enumerate(c.pbc):
            if periodic:
                frac[:, axis] %= 1.0
        positions = frac @ cell
        # Reciprocal-vector bounds guarantee coverage even for skewed cells.
        bounds = np.ceil(self.cutoff * np.linalg.norm(inv, axis=0)).astype(int) + 1
        ranges = [range(-int(b), int(b) + 1) if periodic else [0] for b, periodic in zip(bounds, c.pbc)]
        count = np.prod([len(r) for r in ranges]) * len(c.species) ** 2
        if count > 5_000_000:
            raise ValueError("Cell too large/skewed for the reference evaluator; use smaller fitting cells.")
        ii, jj, vectors, distances = [], [], [], []
        n = len(c.species)
        for image in product(*ranges):
            delta = positions[None, :, :] - positions[:, None, :] + np.asarray(image) @ cell
            r = np.linalg.norm(delta, axis=2)
            mask = r < self.cutoff
            if image == (0, 0, 0):
                np.fill_diagonal(mask, False)
            if np.any(r[mask] < 1e-7):
                raise ValueError(f"{c.name}: overlapping atoms or periodic images.")
            i, j = np.nonzero(mask)
            ii.extend(i); jj.extend(j)
            vectors.extend(delta[i, j]); distances.extend(r[i, j])
        return (np.asarray(ii, int), np.asarray(jj, int), np.asarray(vectors, float).reshape(-1, 3), np.asarray(distances))

    def evaluate(self, c: Configuration, neighbors=None) -> Prediction:
        self.validate()
        i, j, vector, r = self.neighbors(c) if neighbors is None else neighbors
        f, density, pair = self._splines()
        species = np.asarray(c.species)
        rho = np.zeros(len(species))
        drho_j = np.zeros(len(r)); drho_i = np.zeros(len(r)); dphi = np.zeros(len(r))
        pair_energy = 0.0
        for a in self.symbols:
            for b in self.symbols:
                mask = (species[i] == a) & (species[j] == b)
                key_j = b if self.style == "eam/alloy" else density_key(a, b)
                key_i = a if self.style == "eam/alloy" else density_key(b, a)
                vals = density[key_j](r[mask])
                if self.style != "eam/he" and np.any(vals < -1e-10):
                    raise ValueError("Density interpolation became negative; refine the radial grid.")
                np.add.at(rho, i[mask], vals)
                drho_j[mask] = density[key_j](r[mask], 1)
                drho_i[mask] = density[key_i](r[mask], 1)
        if np.any(rho > self.density_grid[-1]) or np.any(rho < self.density_grid[0]):
            raise ValueError(f"{c.name}: density outside embedding grid; extend the domain explicitly.")
        fp = np.zeros(len(species)); energy = 0.0
        for s in self.symbols:
            mask = species == s
            energy += float(np.sum(f[s](rho[mask])))
            fp[mask] = f[s](rho[mask], 1)
        for a_index, a in enumerate(self.symbols):
            for b in self.symbols[:a_index + 1]:
                mask = ((species[i] == a) & (species[j] == b)) | ((species[i] == b) & (species[j] == a))
                spline = pair[pair_key(a, b)]
                radii = r[mask]
                rp = spline(radii)
                pair_energy += 0.5 * float(np.sum(rp / radii))
                dphi[mask] = (spline(radii, 1) * radii - rp) / radii ** 2
        coeff = fp[i] * drho_j + fp[j] * drho_i + dphi
        force_pairs = coeff[:, None] * vector / r[:, None]
        forces = np.zeros((len(species), 3))
        np.add.at(forces, i, force_pairs)
        stress = 0.5 * np.einsum("ni,nj->ij", vector, force_pairs) / abs(np.linalg.det(c.cell))
        return Prediction(energy + pair_energy, forces, stress, rho)

    def write_setfl(self, path, comment="GZ-toolkit eamGen candidate; validate before use"):
        self.validate()
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        header = f"{len(self.density_grid)} {self.density_grid[1] - self.density_grid[0]:.16e} {len(self.radial_grid)} {self.radial_grid[1]:.16e} {self.cutoff:.16e}"
        if self.style == "eam/he":
            header += f" {self.density_grid[-1]:.16e}"
        comments = [comment.replace("\n", " "), "Units: metal (eV, Angstrom)", self.style]
        if self.metadata.get("source_comments"):
            comments = [line.replace("\n", " ") for line in self.metadata["source_comments"][:3]]
            comments[0] += " | " + comment.replace("\n", " ")
        lines = comments + [
                 f"{len(self.elements)} " + " ".join(self.symbols),
                 header]
        def table(values):
            lines.extend(" ".join(f"{x:.16e}" for x in values[k:k + 5]) for k in range(0, len(values), 5))
        for e in self.elements:
            lines.append(f"{e.atomic_number} {e.mass:.16e} {e.lattice_constant:.16e} {e.lattice}")
            table(self.embedding[e.symbol])
            if self.style == "eam/alloy":
                table(self.density[e.symbol])
            else:
                # setfl sections are grouped by donor beta, receivers alpha inside.
                for receiver in self.symbols:
                    table(self.density[density_key(receiver, e.symbol)])
        for i, a in enumerate(self.symbols):
            for b in self.symbols[:i + 1]:
                key = pair_key(a, b)
                rp = np.asarray(self.radial_grid) * self.pair[key]
                rp[0] = self.metadata.get("rphi_at_zero", {}).get(key, 0.0)
                table(rp)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return path

    @classmethod
    def read_setfl(cls, path, style="eam/alloy"):
        """Import setfl with an explicit style; never guess FS/HE/funcfl formats.

        Internal phi(0) is extrapolated from the first two nonzero radii;
        setfl stores r*phi and cannot specify phi(0). No fitting point should
        approach that unresolved core. The global cutoff may lie inside the
        radial table or at most one interval beyond its final sample.
        """
        path = Path(path)
        lines = path.read_text(encoding="utf-8").splitlines()
        tokens = iter(" ".join(lines[3:]).split())
        try:
            n = int(next(tokens)); symbols = [next(tokens) for _ in range(n)]
            nrho, drho, nr, dr, cutoff = int(next(tokens)), float(next(tokens)), int(next(tokens)), float(next(tokens)), float(next(tokens))
            rho_max = float(next(tokens)) if style == "eam/he" else (nrho - 1) * drho
            if n < 1 or nrho < 4 or nr < 4 or style not in {"eam/alloy", "eam/fs", "eam/he"}:
                raise ValueError("Invalid setfl header/style.")
            elements, embedding, density, pairs, rp_zero = [], {}, {}, {}, {}
            def table(size):
                return [float(next(tokens)) for _ in range(size)]
            for s in symbols:
                elements.append(Element(s, int(next(tokens)), float(next(tokens)), float(next(tokens)), next(tokens)))
                embedding[s] = table(nrho)
                if style == "eam/alloy":
                    density[s] = table(nr)
                else:
                    for receiver in symbols:
                        density[density_key(receiver, s)] = table(nr)
            r = np.arange(nr) * dr
            for i, a in enumerate(symbols):
                for b in symbols[:i + 1]:
                    rp = np.asarray(table(nr))
                    phi = np.zeros(nr); phi[1:] = rp[1:] / r[1:]
                    phi[0] = 2 * phi[1] - phi[2]
                    pairs[pair_key(a, b)] = phi.tolist()
                    rp_zero[pair_key(a, b)] = float(rp[0])
            if next(tokens, None) is not None:
                raise ValueError("Extra tables found: wrong style or unsupported CD-EAM.")
        except (StopIteration, IndexError) as exc:
            raise ValueError("Incomplete setfl file.") from exc
        rho_grid = rho_max - (nrho - 1) * drho + np.arange(nrho) * drho
        model = cls(elements, r.tolist(), rho_grid.tolist(), embedding, density, pairs,
                    {"origin": "setfl", "source": str(path.resolve()), "validated": False,
                     "rphi_at_zero": rp_zero, "source_comments": lines[:3],
                     "source_sha256": sha256(path.read_bytes()).hexdigest()}, style, cutoff)
        model.validate()
        return model


def seed_model(elements: list[Element], cutoff=6.0, rho_max=100.0, nr=301, nrho=301, style="eam/alloy", rho_min=0.0) -> EAMModel:
    """Generic smooth starting guess, NOT a physical parameterization.

    Density gauge/shape stays fixed in the first fitter. All elemental and
    cross-pair functions are present, so adding an element never relies on
    an undocumented mixing rule. User/reference data determines corrections.
    """
    if not np.isfinite(cutoff) or not np.isfinite(rho_max) or cutoff <= 0 or rho_max <= 0:
        raise ValueError("cutoff and rho_max must be positive.")
    r = np.linspace(0, cutoff, nr); rho = np.linspace(rho_min, rho_max, nrho)
    taper = (1 - r / cutoff) ** 4
    density = {e.symbol: (np.exp(-r) * taper).tolist() for e in elements}
    embedding = {e.symbol: (-np.sqrt(np.maximum(rho, 0) + 0.01) + 0.1).tolist() for e in elements}
    if style != "eam/alloy":
        density = {density_key(a.symbol, b.symbol): density[b.symbol][:] for a in elements for b in elements}
    pairs = {pair_key(a.symbol, b.symbol): (0.1 * np.exp(-2 * r) * taper).tolist() for a in elements for b in elements}
    model = EAMModel(elements, r.tolist(), rho.tolist(), embedding, density, pairs,
                     {"origin": "generic_seed", "validated": False, "density_policy": "fixed"}, style)
    model.validate()
    return model
