"""Portable reference configurations. Units: eV, Angstrom, eV/Angstrom^3.

Cells contain lattice vectors as rows; stress is positive in tension.
Energy offsets are explicit per-species reference alignments, not fitted away.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
import numpy as np


@dataclass
class Configuration:
    name: str
    species: list[str]
    positions: list[list[float]]
    cell: list[list[float]]
    energy: float | None = None
    forces: list[list[float]] | None = None
    stress: list[list[float]] | None = None
    pbc: list[bool] = field(default_factory=lambda: [True, True, True])
    split: str = "train"
    group: str = ""
    weight: float = 1.0
    provenance: dict = field(default_factory=dict)

    def validate(self) -> None:
        n = len(self.species)
        if not self.name or n == 0 or any(not s for s in self.species):
            raise ValueError("Configuration needs a name and nonempty species.")
        for key, shape in (("positions", (n, 3)), ("cell", (3, 3)),
                           ("forces", (n, 3)), ("stress", (3, 3))):
            value = getattr(self, key)
            if value is not None:
                arr = np.asarray(value, dtype=float)
                if arr.shape != shape or not np.isfinite(arr).all():
                    raise ValueError(f"{self.name}: invalid {key}; expected finite {shape}.")
        if abs(np.linalg.det(self.cell)) < 1e-10:
            raise ValueError(f"{self.name}: cell is singular.")
        if len(self.pbc) != 3 or any(type(p) is not bool for p in self.pbc):
            raise ValueError("pbc must contain three booleans.")
        if self.split not in {"train", "validation", "test"}:
            raise ValueError("split must be train, validation or test.")
        if not np.isfinite(self.weight) or self.weight <= 0:
            raise ValueError("Configuration weight must be positive.")
        if self.energy is not None and not np.isfinite(self.energy):
            raise ValueError("Energy must be finite.")
        if self.stress is not None and (not all(self.pbc) or not np.allclose(self.stress, np.asarray(self.stress).T)):
            raise ValueError("Stress targets require a fully periodic cell and symmetric tensor.")


@dataclass
class Dataset:
    configurations: list[Configuration]
    energy_reference: str
    energy_offsets: dict[str, float] = field(default_factory=dict)
    schema_version: int = 1

    def validate(self, require_targets: bool = True) -> None:
        if self.schema_version != 1 or not self.energy_reference.strip():
            raise ValueError("Dataset needs schema_version=1 and a documented energy_reference.")
        if not self.configurations:
            raise ValueError("Dataset is empty.")
        if any(not np.isfinite(v) for v in self.energy_offsets.values()):
            raise ValueError("Energy offsets must be finite.")
        names, groups = set(), {}
        for c in self.configurations:
            c.validate()
            if c.name in names:
                raise ValueError(f"Duplicate configuration name: {c.name}")
            names.add(c.name)
            if not c.provenance:
                raise ValueError(f"{c.name}: record reference provenance before fitting.")
            if require_targets and c.energy is None and c.forces is None and c.stress is None:
                raise ValueError(f"{c.name}: no fitting targets.")
            if c.group:
                previous = groups.setdefault(c.group, c.split)
                if previous != c.split:
                    raise ValueError(f"Group {c.group!r} crosses splits; related configurations would leak.")

    def aligned_energy(self, c: Configuration) -> float | None:
        if c.energy is None:
            return None
        return float(c.energy - sum(self.energy_offsets.get(s, 0.0) for s in c.species))

    def save(self, path: str | Path) -> Path:
        self.validate(require_targets=False)
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2, allow_nan=False), encoding="utf-8")
        return path

    @classmethod
    def load(cls, path: str | Path) -> "Dataset":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        payload["configurations"] = [Configuration(**c) for c in payload["configurations"]]
        dataset = cls(**payload)
        dataset.validate(require_targets=False)
        return dataset


def compare_metrics(actual: dict[str, float], targets: list[dict]) -> list[dict]:
    """Compare testing metrics with sourced experimental/DFT targets.

    Targets: metric, value, tolerance, units, source (URL/DOI), conditions.
    Missing results fail review rather than silently passing.
    """
    rows = []
    for t in targets:
        for key in ("metric", "value", "tolerance", "units", "source", "conditions"):
            if key not in t or (key in {"metric", "units", "source", "conditions"} and not t[key]):
                raise ValueError(f"Comparison target needs {key}.")
        if not np.isfinite(t["value"]) or not np.isfinite(t["tolerance"]) or t["tolerance"] <= 0:
            raise ValueError("Target value/tolerance must be finite; tolerance must be positive.")
        value = actual.get(t["metric"])
        error = None if value is None else float(value - t["value"])
        rows.append({**t, "actual": value, "error": error,
                     "passed": error is not None and np.isfinite(error) and abs(error) <= t["tolerance"]})
    return rows
