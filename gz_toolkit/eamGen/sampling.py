"""Reproducible small FCC/BCC snapshots to send to reference calculators.

These are geometry proposals, not labelled data or validated defect structures.
Related volume/strain snapshots share a group and must stay in one split.
"""
from itertools import product
import numpy as np

from .data import Configuration


def crystal_snapshots(name, host, lattice_constant, lattice="fcc", repetitions=2,
                      composition=None, scales=(0.95, 1.0, 1.05), displacement=0.0,
                      seed=0, split="train", gas=None, gas_site="tetra", vacancy=False):
    """Propose EOS/alloy/gas snapshots with deterministic chemistry/displacements.

    composition uses fractions summing to 1; gas is inserted explicitly after
    substitution. FCC gas sites: tetra (1/4,1/4,1/4), octa (1/2,0,0).
    BCC sites: tetra (1/4,1/2,0), octa (1/2,1/2,0), in host-cell units.
    A vacancy-gas case replaces a selected host atom with gas at that site.
    """
    if lattice.lower() not in {"fcc", "bcc"} or not isinstance(repetitions, int) or repetitions < 1:
        raise ValueError("Snapshot sampler supports FCC/BCC and positive integer repetitions.")
    if not np.isfinite(lattice_constant) or lattice_constant <= 0 or not np.isfinite(displacement) or displacement < 0:
        raise ValueError("Positive lattice constant and nonnegative displacement required.")
    if not scales or any(not np.isfinite(s) or s <= 0 for s in scales):
        raise ValueError("Supply positive finite volume scales.")
    if not isinstance(seed, int):
        raise ValueError("Record an integer seed for reproducibility.")
    composition = {host: 1.0} if composition is None else composition
    if host not in composition or any(not np.isfinite(v) or v < 0 for v in composition.values()) or not np.isclose(sum(composition.values()), 1):
        raise ValueError("Composition fractions must be nonnegative, sum to 1 and include the host.")
    basis = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]]) if lattice.lower() == "fcc" else np.array([[0, 0, 0], [.5, .5, .5]])
    coordinates = np.vstack([basis + np.array(index) for index in product(range(repetitions), repeat=3)])
    rng = np.random.default_rng(seed)
    n = len(coordinates)
    # Largest remainder allocation gives the closest possible finite-cell composition.
    fractions = np.asarray(list(composition.values()))
    counts = np.floor(fractions * n).astype(int)
    for k in np.argsort(-(fractions*n - counts))[:n - counts.sum()]:
        counts[k] += 1
    chemistry = np.concatenate([np.repeat(s, count) for s, count in zip(composition, counts)])
    rng.shuffle(chemistry)
    chemistry = chemistry.tolist()
    if vacancy:
        if n <= 1:
            raise ValueError("Vacancy would empty the cell.")
        chemistry.pop(0); coordinates = coordinates[1:]
    if gas:
        if vacancy:
            site = np.zeros(3)
        elif gas_site == "tetra":
            site = np.array([.25, .25, .25] if lattice.lower() == "fcc" else [.25, .5, 0])
        elif gas_site == "octa":
            site = np.array([.5, 0, 0] if lattice.lower() == "fcc" else [.5, .5, 0])
        else:
            raise ValueError("Gas site must be tetra or octa.")
        coordinates = np.vstack((coordinates, site)); chemistry.append(gas)
    noise = rng.normal(scale=displacement, size=coordinates.shape)
    configurations = []
    for index, scale in enumerate(scales):
        a = lattice_constant * scale
        positions = coordinates * a + noise
        cell = np.eye(3) * repetitions * a
        c = Configuration(f"{name}_{index:03d}", chemistry[:], positions.tolist(), cell.tolist(), split=split,
                          group=name, provenance={"source": "geometry proposal; reference calculation required", "seed": seed,
                                                  "lattice": lattice, "composition": composition, "scale": scale,
                                                  "displacement_A": displacement, "gas_site": gas_site if gas else None,
                                                  "vacancy": vacancy})
        c.validate(); configurations.append(c)
    return configurations
