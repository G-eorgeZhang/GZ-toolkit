"""Selective, regularized spline corrections to an existing EAM model.

With fixed density functions the fit is linear. Opt-in density corrections
use nonlinear least squares, including signed receiver/donor densities for HE.
The fitter never changes validation/test targets or silently aligns energies.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass, field
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import lsq_linear, least_squares

from .data import Dataset
from .model import EAMModel


@dataclass
class FitOptions:
    embedding_keys: list[str] | None = None  # None = all; [] = frozen
    pair_keys: list[str] | None = None
    density_keys: list[str] = field(default_factory=list)  # opt-in nonlinear fit
    knots: int = 8
    radial_min: float = 0.5  # preserve the baseline short-range core
    energy_scale: float = 0.02  # eV/atom; inverse weighting scales
    force_scale: float = 0.1  # eV/Angstrom
    stress_scale: float = 0.01  # eV/Angstrom^3
    regularization: float = 1e-3
    coefficient_bound: float = 10.0
    max_evaluations: int = 200
    max_parameters: int = 500

    def validate(self, model):
        if self.knots < 4 or self.max_evaluations < 1 or self.max_parameters < 1:
            raise ValueError("Need >=4 knots and positive evaluation/parameter limits.")
        values = [self.energy_scale, self.force_scale, self.stress_scale, self.coefficient_bound]
        if any(not np.isfinite(v) or v <= 0 for v in values):
            raise ValueError("Scales and coefficient_bound must be positive and finite.")
        if not np.isfinite(self.regularization) or self.regularization <= 0:
            raise ValueError("Positive regularization is required to control EAM gauge/nonidentifiability.")
        if not 0 <= self.radial_min < model.cutoff:
            raise ValueError("radial_min must be inside the cutoff.")


@dataclass
class FitResult:
    model: EAMModel
    report: dict
    coefficients: list[float]


def _correction_basis(model, options):
    basis, labels = [], []
    for component, requested in (("embedding", options.embedding_keys), ("pair", options.pair_keys),
                                 ("density", options.density_keys)):
        tables = getattr(model, component)
        keys = list(tables) if requested is None else requested
        if len(set(keys)) != len(keys) or not set(keys).issubset(tables):
            raise ValueError(f"Invalid/duplicate {component} keys: {keys}")
        embedding = component == "embedding"
        grid = np.asarray(model.density_grid if embedding else model.radial_grid)
        knots = np.linspace(grid[0], grid[-1], options.knots) if embedding else np.linspace(options.radial_min, model.cutoff, options.knots)
        # Preserve F(0) and radial core/cutoff. Densities may be signed in HE.
        excluded = int(np.argmin(abs(knots))) if embedding else None
        if embedding and not grid[0] <= 0 <= grid[-1]:
            raise ValueError("Embedding correction requires a density grid containing zero.")
        for key in keys:
            for k in range(options.knots):
                if embedding and k == excluded:
                    continue
                if not embedding and k in (0, options.knots - 1):
                    continue
                y = np.zeros(options.knots); y[k] = 1
                spline = CubicSpline(knots, y, bc_type="natural" if embedding else ((1, 0.0), (1, 0.0)))
                if embedding:
                    values = spline(grid) - spline(0.0)
                else:
                    values = np.zeros_like(grid)
                    mask = (grid >= options.radial_min) & (grid <= model.cutoff)
                    values[mask] = spline(grid[mask])
                basis.append((component, key, values))
                labels.append(f"{component}:{key}:k{k}")
    if not basis or len(basis) > options.max_parameters:
        raise ValueError(f"Fit has {len(basis)} parameters; select functions/knots within max_parameters.")
    return basis, labels


def _apply(model, basis, coefficients):
    candidate = deepcopy(model)
    for (component, key, values), coefficient in zip(basis, coefficients):
        tables = getattr(candidate, component)
        tables[key] = (np.asarray(tables[key]) + coefficient * values).tolist()
    return candidate


def evaluate_dataset(model: EAMModel, dataset: Dataset, include_test=False) -> dict:
    """Raw per-split RMSE plus per-case errors. Test split is opt-in.

    Train RMSE is a fitting diagnostic, not a claim of physical accuracy.
    """
    dataset.validate()
    residuals = {}; cases = []
    for c in dataset.configurations:
        if c.split == "test" and not include_test:
            continue
        p = model.evaluate(c)
        values = residuals.setdefault(c.split, {"energy_per_atom": [], "force": [], "stress": []})
        row = {"name": c.name, "split": c.split, "group": c.group,
               "density_range": [float(p.densities.min()), float(p.densities.max())]}
        if c.energy is not None:
            error = (p.energy - dataset.aligned_energy(c)) / len(c.species)
            values["energy_per_atom"].append(error)
            row["energy_error_per_atom"] = float(error)
        if c.forces is not None:
            error = p.forces - c.forces
            values["force"].extend(error.ravel())
            row["force_rmse"] = float(np.sqrt(np.mean(error ** 2)))
        if c.stress is not None:
            error = p.stress - c.stress
            values["stress"].extend(error[np.triu_indices(3)])
            row["stress_rmse"] = float(np.sqrt(np.mean(error ** 2)))
        cases.append(row)
    metrics = {split: {key + "_rmse": float(np.sqrt(np.mean(np.asarray(v) ** 2))) if v else None
                       for key, v in values.items()} for split, values in residuals.items()}
    return {"metrics": metrics, "cases": cases, "test_evaluated": include_test,
            "validation_present": "validation" in metrics,
            "units": {"energy_per_atom": "eV/atom", "force": "eV/Angstrom", "stress": "eV/Angstrom^3"}}


def fit_model(model: EAMModel, dataset: Dataset, options: FitOptions | None = None) -> FitResult:
    model.validate(); dataset.validate()
    options = options or FitOptions()
    options.validate(model)
    training = [c for c in dataset.configurations if c.split == "train"]
    if not training:
        raise ValueError("No training configurations.")
    basis, labels = _correction_basis(model, options)
    geometry = [model.neighbors(c) for c in training]
    def vector(candidate):
        out = []
        for c, neighbors in zip(training, geometry):
            p = candidate.evaluate(c, neighbors)
            w = np.sqrt(c.weight)
            if c.energy is not None:
                out.append(w * (p.energy - dataset.aligned_energy(c)) / len(c.species) / options.energy_scale)
            if c.forces is not None:
                out.extend((w * (p.forces - c.forces) / options.force_scale / np.sqrt(p.forces.size)).ravel())
            if c.stress is not None:
                out.extend(w * (p.stress - c.stress)[np.triu_indices(3)] / options.stress_scale / np.sqrt(6))
        return np.asarray(out)
    initial = vector(model)  # fail immediately for invalid starting densities
    zero = np.zeros(len(basis))
    if not options.density_keys:
        columns = []
        for k in range(len(basis)):
            x = zero.copy(); x[k] = 1
            columns.append(vector(_apply(model, basis, x)) - initial)
        design = np.column_stack(columns)
        ridge = np.sqrt(options.regularization) * np.eye(len(basis))
        solution = lsq_linear(np.vstack((design, ridge)), np.concatenate((-initial, zero)),
                              bounds=(-options.coefficient_bound, options.coefficient_bound),
                              max_iter=options.max_evaluations)
        rank = int(np.linalg.matrix_rank(design))
        evaluations = solution.nit
    else:
        def residual(x):
            try:
                values = vector(_apply(model, basis, x))
            except ValueError:
                # Invalid-density trials are rejected; no extrapolated embedding targets.
                values = np.full(initial.shape, 1e6)
            return np.concatenate((values, np.sqrt(options.regularization) * x))
        solution = least_squares(residual, zero, bounds=(-options.coefficient_bound, options.coefficient_bound),
                                 max_nfev=options.max_evaluations)
        rank = None
        evaluations = solution.nfev
    fitted = _apply(model, basis, solution.x)
    fitted.validate()
    fitted.metadata.update({"validated": False, "origin": "fit", "fit_options": asdict(options)})
    report = evaluate_dataset(fitted, dataset)
    report.update({"optimizer_success": bool(solution.success), "optimizer_message": solution.message,
                   "evaluations": evaluations, "parameters": len(basis), "data_rank": rank,
                   "parameter_labels": labels, "options": asdict(options), "energy_reference": dataset.energy_reference,
                   "energy_offsets": dataset.energy_offsets, "status": "candidate_requires_review",
                   "warnings": ["Internal interpolation differs from LAMMPS; verify exported file with LAMMPS.",
                                "Density functions are frozen." if not options.density_keys else "Density changes can alter host interactions; run regression cases.",
                                "No held-out validation configurations." if not report["validation_present"] else "Final test split has not been inspected."]})
    if rank is not None and rank < len(basis):
        report["warnings"].append("Training data does not identify all coefficients; regularization selects a solution.")
    return FitResult(fitted, report, solution.x.tolist())
