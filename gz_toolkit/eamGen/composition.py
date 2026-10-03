"""Explicit provisional composition of baseline functions, with provenance.

Composition is not a validated alloy potential or an automatic gauge alignment.
The primary model wins every overlap. Missing cross functions are deliberate
seed guesses, named in the model metadata for subsequent reference fitting.
"""
from copy import deepcopy
import numpy as np
from scipy.interpolate import CubicSpline

from .model import EAMModel, density_key, pair_key, seed_model


def compose_model(primary: EAMModel, secondary: EAMModel, symbols: list[str]) -> EAMModel:
    primary.validate(); secondary.validate()
    if primary.style != "eam/he":
        raise ValueError("Composition currently requires an eam/he primary baseline.")
    if len(set(symbols)) != len(symbols) or not symbols:
        raise ValueError("Choose unique output elements.")
    metadata = {e.symbol: e for e in secondary.elements}
    metadata.update({e.symbol: e for e in primary.elements})
    if not set(symbols).issubset(metadata):
        raise ValueError("Requested element absent from both baselines.")
    result = seed_model([deepcopy(metadata[s]) for s in symbols], cutoff=primary.cutoff,
                        rho_min=primary.density_grid[0], rho_max=primary.density_grid[-1],
                        nr=len(primary.radial_grid), nrho=len(primary.density_grid), style="eam/he")
    # Keep the primary grids exactly, including files whose cutoff is not their final sample.
    result.radial_grid = primary.radial_grid[:]
    result.density_grid = primary.density_grid[:]
    result.cutoff_radius = primary.cutoff
    sources, provisional = {}, []
    for component in ("embedding", "density", "pair"):
        tables = getattr(result, component)
        for key in tables:
            chosen, source_key = None, key
            for source in (primary, secondary):
                candidate_key = key
                if component == "density" and source.style == "eam/alloy":
                    receiver, donor = key.split("<-")
                    if receiver not in source.symbols or donor not in source.symbols:
                        continue
                    candidate_key = donor
                if candidate_key in getattr(source, component):
                    chosen, source_key = source, candidate_key
                    break
            if chosen is None:
                # Regenerate guess on the actual primary grids (seed grid may differ).
                grid = np.asarray(result.radial_grid)
                tables[key] = (np.exp(-grid) * np.maximum(1 - grid / result.cutoff, 0) ** 4).tolist()
                if component == "pair":
                    tables[key] = (0.1 * np.exp(-2 * grid) * np.maximum(1 - grid / result.cutoff, 0) ** 4).tolist()
                provisional.append(f"{component}:{key}")
                continue
            old_grid = np.asarray(chosen.density_grid if component == "embedding" else chosen.radial_grid)
            new_grid = np.asarray(result.density_grid if component == "embedding" else result.radial_grid)
            if chosen is primary:
                tables[key] = getattr(chosen, component)[source_key][:]
            else:
                # Secondary embedding has no defined extension onto negative HE density.
                # Reject extrapolation: keep a documented provisional seed until user supplies a function.
                if component == "embedding" and (new_grid[0] < old_grid[0] or new_grid[-1] > old_grid[-1]):
                    provisional.append(f"embedding:{key} (secondary domain incompatible; generic seed)")
                    continue
                spline = CubicSpline(old_grid, getattr(chosen, component)[source_key], bc_type="natural")
                if component == "embedding":
                    tables[key] = spline(new_grid).tolist()
                else:
                    values = np.zeros_like(new_grid)
                    mask = new_grid < chosen.cutoff
                    if np.any(new_grid[mask] > old_grid[-1] + np.diff(old_grid)[0] + 1e-8):
                        raise ValueError("Secondary radial domain cannot cover composition.")
                    values[mask] = spline(new_grid[mask])
                    tables[key] = values.tolist()
                provisional.append(f"{component}:{key} (secondary gauge; refit/review)")
            sources[f"{component}:{key}"] = chosen.metadata.get("source", chosen.metadata.get("origin", "unspecified"))
    result.metadata = {"origin": "provisional_composition", "validated": False, "sources": sources,
                       "provisional_functions": provisional,
                       "rphi_at_zero": {key: value for key, value in primary.metadata.get("rphi_at_zero", {}).items() if key in result.pair},
                       "source_comments": primary.metadata.get("source_comments", []),
                       "baseline_sha256": [primary.metadata.get("source_sha256"), secondary.metadata.get("source_sha256")],
                       "warning": "Baselines are not gauge aligned. Cross interactions require reference fitting and regression tests."}
    result.validate()
    return result
