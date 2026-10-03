"""Reproducible candidate revisions and potential-testing handoff."""
from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re
import shutil

from .data import Dataset
from .fitting import FitOptions, fit_model
from .model import EAMModel, Element, seed_model


def _json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False), encoding="utf-8")


def _safe_name(value):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value) or value in {".", ".."}:
        raise ValueError("Use a simple project/revision name without path separators.")
    return value


def init_project(root, elements: list[Element], style="eam/he", baseline: EAMModel | None = None) -> Path:
    root = Path(root)
    if root.exists() and any(root.iterdir()):
        raise FileExistsError("Fitting project must be new or empty.")
    root.mkdir(parents=True, exist_ok=True)
    model = baseline or seed_model(elements, style=style, rho_min=-10.0 if style == "eam/he" else 0.0)
    model.save(root / "baseline.json")
    (root / "revisions").mkdir()
    (root / "references").mkdir()
    _json(root / "fit_options.json", asdict(FitOptions()))
    _json(root / "dataset.json", {"schema_version": 1, "energy_reference": "REPLACE: reference method, energy zero and settings",
                                   "energy_offsets": {}, "configurations": []})
    _json(root / "targets.json", [])
    _json(root / "project.json", {"schema_version": 1, "elements": model.symbols, "style": model.style,
                                   "status": "reference_data_required", "intended_use": "Describe system, temperature, compositions and target phenomena.",
                                   "sources": [], "decisions": []})
    return root


def fit_revision(root, name, baseline=None) -> Path:
    root = Path(root)
    _safe_name(name)
    data_path = root / "dataset.json"
    model_path = Path(baseline) if baseline is not None else root / "baseline.json"
    options_path = root / "fit_options.json"
    dataset = Dataset.load(data_path)
    if dataset.energy_reference.startswith("REPLACE:"):
        raise ValueError("Document the reference energy convention before fitting.")
    options = FitOptions(**json.loads(options_path.read_text(encoding="utf-8")))
    destination = root / "revisions" / name
    if destination.exists():
        raise FileExistsError("Revision already exists; use a new name.")
    terms_path = root / "analytic_terms.json"
    if terms_path.exists():
        from .analytic import AnalyticTerm, fit_analytic_model
        terms = [AnalyticTerm(**t) for t in json.loads(terms_path.read_text(encoding="utf-8"))]
        result = fit_analytic_model(EAMModel.load(model_path), dataset, terms, options)
    else:
        result = fit_model(EAMModel.load(model_path), dataset, options)
    destination.mkdir(parents=True)
    result.model.save(destination / "model.json")
    suffix = result.model.style.replace("/", ".")
    result.model.write_setfl(destination / f"candidate.{suffix}")
    hashes = {}
    for path, target in ((data_path, "dataset.json"), (model_path, "baseline.json"), (options_path, "fit_options.json")):
        shutil.copy2(path, destination / target)
        hashes[target] = sha256(path.read_bytes()).hexdigest()
    if terms_path.exists():
        shutil.copy2(terms_path, destination / "analytic_terms.json")
        hashes["analytic_terms.json"] = sha256(terms_path.read_bytes()).hexdigest()
    result.report.update({"revision": name, "created": datetime.now(timezone.utc).isoformat(), "input_sha256": hashes})
    _json(destination / "report.json", result.report)
    _json(destination / "coefficients.json", result.coefficients)
    _json(destination / "review.json", {"status": "pending", "reviewer": None, "reason": None})
    return destination


def record_review(revision, status, reviewer, reason) -> Path:
    """Record a real user's decision; approval is for testing, not deployment."""
    if status not in {"approved_for_testing", "rejected"} or not reviewer.strip() or not reason.strip():
        raise ValueError("Supply status approved_for_testing/rejected, reviewer and reason.")
    revision = Path(revision)
    report = json.loads((revision / "report.json").read_text(encoding="utf-8"))
    if status == "approved_for_testing" and not report["optimizer_success"]:
        raise ValueError("Optimizer did not converge; revise the fit before approving for testing.")
    path = revision / "review.json"
    previous = json.loads(path.read_text(encoding="utf-8")) if path.exists() else None
    _json(path, {"status": status, "reviewer": reviewer, "reason": reason,
                 "recorded": datetime.now(timezone.utc).isoformat(), "previous": previous})
    return path


def export_testing_project(model: EAMModel, root, pot_name="eam_candidate", gases=None) -> Path:
    """Create a new potential_testing project; never submit or overwrite runs.

    Defaults identify He as gas, Pd/Ni as FCC from baseline metadata. Review
    all crystal choices, sizes and HPC settings in the generated config.
    """
    from gz_toolkit.potential_testing.config import PotentialConfig, PotentialMetadata, WorkflowOptions, HPCOptions, save_potential_config
    from gz_toolkit.potential_testing.project import write_driver_scripts
    model.validate(); _safe_name(pot_name)
    gases = ([s for s in model.symbols if s == "He"] if gases is None else gases)
    if not set(gases).issubset(model.symbols):
        raise ValueError("Gas species must be in the model.")
    metals = [s for s in model.symbols if s not in gases]
    if not metals:
        raise ValueError("Testing needs at least one host metal.")
    structures = {e.symbol: {"structure": e.lattice.lower()} for e in model.elements if e.symbol in metals}
    if any(d["structure"] not in {"bcc", "fcc", "hcp"} for d in structures.values()):
        raise ValueError("Testing supports bcc/fcc/hcp hosts; specify supported baseline metadata.")
    if any(e.lattice_constant <= 0 for e in model.elements if e.symbol in metals):
        raise ValueError("Host lattice constants must be positive for testing.")
    root = Path(root)
    if root.exists() and any(root.iterdir()):
        raise FileExistsError("Testing destination must be new or empty.")
    root.mkdir(parents=True, exist_ok=True)
    files = root / "potentials" / pot_name
    files.mkdir(parents=True)
    filename = f"{pot_name}.{model.style.replace('/', '.')}"
    model.write_setfl(files / filename)
    suites = [{"A": a, "B": b, "fractions_atpct": [25, 50, 75], "ordering": "random"}
              for i, a in enumerate(metals) for b in metals[i + 1:]]
    cfg = PotentialConfig(
        potential=PotentialMetadata(pot_name=pot_name,
                                    pot_lines=f"pair_style {model.style}\npair_coeff * * {filename} " + " ".join(model.symbols),
                                    type_map={s: i + 1 for i, s in enumerate(model.symbols)},
                                    masses={i + 1: e.mass for i, e in enumerate(model.elements)}, single_elements=metals,
                                    gases=gases, crystal_structures=structures, two_element_suites=suites,
                                    potential_files=[filename], lc_initial=model.elements[model.symbols.index(metals[0])].lattice_constant,
                                    gas_max_m=2, gas_max_n=2, gas_with_interstitial=False),
        workflow=WorkflowOptions(include_elastic=True, include_alloy_suite=bool(suites), include_gas_complexes=bool(gases)),
        hpc=HPCOptions())
    save_potential_config(cfg, root / "pot_inputs" / f"{pot_name}.json")
    write_driver_scripts(root)
    _json(root / "eamGen_provenance.json", {"model_metadata": model.metadata,
                                           "potential_sha256": sha256((files / filename).read_bytes()).hexdigest(),
                                           "status": "candidate_for_testing", "required": "LAMMPS parity, convergence, holdout and physical regression checks"})
    return root
