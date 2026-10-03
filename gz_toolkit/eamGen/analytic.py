"""Named analytic functions from Zhou et al., PRB 103, 014108 (2021).

Equations 1, 2, 4, 6 and 7, for extending the signed-density HE formulation.
Parameter signs follow the actual equation/table (Table II has beta_HeHe < 0).
The original Pd/H functions are retained from the supplied baseline tables.
"""
from __future__ import annotations
from copy import deepcopy
from dataclasses import asdict, dataclass, field
import numpy as np
from scipy.special import erfc, erfcinv
from scipy.optimize import least_squares

from .data import Dataset
from .model import EAMModel
from .fitting import FitOptions, FitResult, evaluate_dataset


def zhou_cutoff(r, rs, rc):
    """Eq. 2: fc(rs)=0.9, fc(rc-)=1e-5, then zero at/above rc."""
    if not np.isfinite(rs) or not np.isfinite(rc) or not 0 <= rs < rc:
        raise ValueError("Cutoff needs 0 <= rs < rc.")
    r = np.asarray(r, float)
    mu, nu = erfcinv(2e-5), erfcinv(1.8)
    return np.where(r < rc, 0.5 * erfc((mu * (r-rs) + nu * (rc-r)) / (rc-rs)), 0.0)


def analytic_values(form, grid, parameters):
    x = np.asarray(grid, float)
    p = parameters
    if not all(np.isfinite(v) for v in p.values()):
        raise ValueError("Analytic parameters must be finite.")
    if form == "zhou_repulsive_pair":
        if p["r0"] <= 0 or p["alpha"] < 0 or p["E0"] < 0:
            raise ValueError("Repulsive pair requires r0>0, alpha>=0, E0>=0.")
        values = p["E0"] * np.exp(-p["alpha"] * (x-p["r0"]) / p["r0"]) * zhou_cutoff(x, p["rs"], p["rc"])
    elif form == "zhou_density":
        if p["gamma"] < 0:
            raise ValueError("Density decay gamma must be nonnegative; amplitude may be signed.")
        values = p["amplitude"] * np.exp(-p["gamma"] * x) * zhou_cutoff(x, p["rs"], p["rc"])
    elif form == "zhou_he_core_pair":
        if p["r0"] <= 0:
            raise ValueError("Core radius must be positive.")
        delta = x - p["r0"]
        values = np.where(x < p["r0"], p["alpha"] * delta**2 + p["beta"] * delta**3, 0.0)
    elif form == "zhou_he_embedding":
        if p["rho0"] >= 0 or any(p[k] < 0 for k in ("F0", "F2", "F3")):
            raise ValueError("HE embedding requires negative rho0 and nonnegative F0/F2/F3.")
        values = np.empty_like(x)
        low = x <= p["rho0"]; middle = (x > p["rho0"]) & (x < 0); high = x >= 0
        values[low] = p["F0"]
        values[middle] = p["F0"] * (0.5 - 0.5 * np.cos(x[middle] / p["rho0"] * np.pi))
        values[high] = p["F2"] * x[high]**2 + p["F3"] * x[high]**3
    else:
        raise ValueError(f"Unknown analytic form: {form}")
    if not np.isfinite(values).all():
        raise ValueError("Analytic function overflowed.")
    return values


@dataclass
class AnalyticTerm:
    component: str
    key: str
    form: str
    parameters: dict[str, float]
    bounds: dict[str, list[float]] = field(default_factory=dict)  # only these parameters are free
    source: str = ""


def apply_analytic_terms(model: EAMModel, terms: list[AnalyticTerm]) -> EAMModel:
    result = deepcopy(model); seen = set()
    expected = {"zhou_repulsive_pair": "pair", "zhou_density": "density",
                "zhou_he_core_pair": "pair", "zhou_he_embedding": "embedding"}
    for term in terms:
        if expected.get(term.form) != term.component or not term.source:
            raise ValueError("Analytic term needs a matching component/form and documented source/hypothesis.")
        if (term.component, term.key) in seen:
            raise ValueError("Duplicate analytic function.")
        seen.add((term.component, term.key))
        if term.component == "embedding" and model.style != "eam/he":
            raise ValueError("Signed HE embedding requires eam/he.")
        tables = getattr(result, term.component)
        if term.key not in tables:
            raise ValueError(f"Unknown model function: {term.component}:{term.key}")
        if term.component != "embedding" and term.parameters.get("rc", term.parameters.get("r0", 0)) > model.cutoff:
            raise ValueError("Analytic radial cutoff cannot exceed the model cutoff.")
        grid = model.density_grid if term.component == "embedding" else model.radial_grid
        tables[term.key] = analytic_values(term.form, grid, term.parameters).tolist()
        if term.component == "pair":
            result.metadata.setdefault("rphi_at_zero", {})[term.key] = 0.0
    result.metadata.update({"analytic_terms": [asdict(t) for t in terms], "validated": False})
    result.validate()
    return result


def fit_analytic_model(model: EAMModel, dataset: Dataset, terms: list[AnalyticTerm], options=None) -> FitResult:
    """Fit bounded named parameters; all unlisted model functions remain frozen.

    Select spline fitting or analytic fitting per revision. A later revision can
    use the previous candidate as its baseline, but cannot fit both parameter
    families simultaneously with this first implementation.
    """
    model.validate(); dataset.validate()
    options = options or FitOptions(); options.validate(model)
    training = [c for c in dataset.configurations if c.split == "train"]
    if not training:
        raise ValueError("No training configurations.")
    initial_terms = deepcopy(terms)
    labels, locations, initial, lower, upper = [], [], [], [], []
    for i, term in enumerate(terms):
        for parameter, bounds in term.bounds.items():
            if parameter not in term.parameters or len(bounds) != 2:
                raise ValueError("Bounds must name an existing parameter and contain [lower, upper].")
            value = term.parameters[parameter]
            if not all(np.isfinite(v) for v in bounds) or not bounds[0] < bounds[1] or not bounds[0] <= value <= bounds[1]:
                raise ValueError("Finite strict bounds must contain the initial parameter.")
            labels.append(f"{term.component}:{term.key}:{parameter}"); locations.append((i, parameter))
            initial.append(value); lower.append(bounds[0]); upper.append(bounds[1])
    if not labels or len(labels) > options.max_parameters:
        raise ValueError("Supply bounded free parameters within max_parameters.")
    initial = np.asarray(initial); scales = np.maximum(abs(initial), 1e-3)
    geometry = [model.neighbors(c) for c in training]
    def candidate(x):
        selected = deepcopy(initial_terms)
        for (i, parameter), value in zip(locations, x):
            selected[i].parameters[parameter] = float(value)
        return apply_analytic_terms(model, selected)
    def vector(m):
        values = []
        for c, neighbors in zip(training, geometry):
            p = m.evaluate(c, neighbors); w = np.sqrt(c.weight)
            if c.energy is not None:
                values.append(w * (p.energy - dataset.aligned_energy(c)) / len(c.species) / options.energy_scale)
            if c.forces is not None:
                values.extend((w * (p.forces-c.forces) / options.force_scale / np.sqrt(p.forces.size)).ravel())
            if c.stress is not None:
                values.extend(w * (p.stress-c.stress)[np.triu_indices(3)] / options.stress_scale / np.sqrt(6))
        return np.asarray(values)
    initial_values = vector(candidate(initial))
    def residual(x):
        try:
            values = vector(candidate(x))
        except (ValueError, FloatingPointError):
            values = np.full(initial_values.shape, 1e6)
        return np.concatenate((values, np.sqrt(options.regularization) * (x-initial) / scales))
    solution = least_squares(residual, initial, bounds=(lower, upper), x_scale=scales, max_nfev=options.max_evaluations)
    fitted = candidate(solution.x)
    report = evaluate_dataset(fitted, dataset)
    report.update({"mode": "analytic", "optimizer_success": bool(solution.success), "optimizer_message": solution.message,
                   "evaluations": solution.nfev, "parameters": len(labels), "parameter_labels": labels,
                   "options": asdict(options), "analytic_terms": fitted.metadata["analytic_terms"],
                   "energy_reference": dataset.energy_reference, "energy_offsets": dataset.energy_offsets,
                   "status": "candidate_requires_review", "warnings": ["Verify exported file against LAMMPS and independently validate physical behavior."]})
    return FitResult(fitted, report, solution.x.tolist())
