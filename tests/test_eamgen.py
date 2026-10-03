"""Numerical and workflow checks; synthetic fits do not validate a physical potential."""
from copy import deepcopy
from dataclasses import asdict
import json
import sys

import numpy as np
import pytest

from gz_toolkit.eamGen import (
    Configuration, Dataset, Element, EAMModel, FitOptions, seed_model, fit_model,
    compose_model, density_key, compare_metrics, init_project, fit_revision,
    record_review, export_testing_project, prepare_lammps, prepare_vasp,
    import_reference, run_reference,
    AnalyticTerm, analytic_values, zhou_cutoff, apply_analytic_terms, fit_analytic_model, crystal_snapshots,
)
from gz_toolkit.eamGen.cli import main


def model(style="eam/he"):
    m = seed_model([Element("Pd", 46, 106.42, 3.89, "FCC"), Element("He", 2, 4.0026, 0, "none")],
                   cutoff=4.0, rho_max=10, rho_min=-2 if style == "eam/he" else 0,
                   nr=81, nrho=121, style=style)
    if style == "eam/he":
        m.density["Pd<-He"] = (-0.1 * np.asarray(m.density["Pd<-He"])).tolist()
        m.density["He<-Pd"] = (0.8 * np.asarray(m.density["He<-Pd"])).tolist()
    return m


def configuration(distance=1.8, periodic=False):
    return Configuration("pair", ["Pd", "He"], [[0.6, 0.4, 0.2], [distance + 0.6, 0.5, 0.3]],
                         [[9, 0, 0], [0.6, 8, 0], [0.3, 0.4, 7]], pbc=[periodic] * 3,
                         provenance={"source": "synthetic regression fixture"})


@pytest.mark.parametrize("style", ["eam/alloy", "eam/fs", "eam/he"])
def test_setfl_roundtrip(tmp_path, style):
    m = model(style)
    m.metadata["rphi_at_zero"] = {"He-Pd": 2.3}
    path = m.write_setfl(tmp_path / "potential.eam")
    loaded = EAMModel.read_setfl(path, style)
    assert loaded.style == style
    assert loaded.symbols == m.symbols
    assert np.allclose(loaded.density_grid, m.density_grid)
    for key in m.density:
        assert np.allclose(loaded.density[key], m.density[key])
    for key in m.pair:
        assert np.allclose(loaded.pair[key][1:], m.pair[key][1:])
    assert loaded.metadata["rphi_at_zero"]["He-Pd"] == 2.3
    c = configuration()
    assert loaded.evaluate(c).energy == pytest.approx(m.evaluate(c).energy, abs=1e-10)
    assert np.allclose(loaded.evaluate(c).forces, m.evaluate(c).forces, atol=1e-10)
    with pytest.raises(ValueError):
        EAMModel.read_setfl(path, "eam/fs" if style != "eam/fs" else "eam/alloy")


@pytest.mark.parametrize("style", ["eam/alloy", "eam/fs", "eam/he"])
def test_forces_are_energy_gradient(style):
    m = model(style); c = configuration(); p = m.evaluate(c)
    h = 1e-6
    for atom in range(2):
        for axis in range(3):
            plus, minus = deepcopy(c), deepcopy(c)
            plus.positions[atom][axis] += h; minus.positions[atom][axis] -= h
            force = -(m.evaluate(plus).energy - m.evaluate(minus).energy) / (2 * h)
            assert p.forces[atom, axis] == pytest.approx(force, abs=2e-7)
    assert np.allclose(p.forces.sum(axis=0), 0, atol=1e-12)
    if style == "eam/he":
        assert p.densities[0] < 0 < p.densities[1]


def test_periodic_images_translation_and_stress():
    m = model()
    c = Configuration("small", ["Pd", "He"], [[0.3, 0.2, 0.1], [1.4, 1.1, 0.9]],
                      [[2.8, 0, 0], [0.7, 2.9, 0], [0.3, 0.4, 3.1]], provenance={"source": "synthetic"})
    p = m.evaluate(c)
    shifted = deepcopy(c); shifted.positions[1] = (np.asarray(c.positions[1]) + np.asarray(c.cell)[0]).tolist()
    assert m.evaluate(shifted).energy == pytest.approx(p.energy, abs=1e-10)
    assert len(m.neighbors(c)[0]) > 2  # multiple images, not minimum-image only
    h = 1e-6
    for axis in range(3):
        strains = []
        for sign in (1, -1):
            transform = np.eye(3); transform[axis, axis] += sign * h
            trial = deepcopy(c)
            trial.cell = (np.asarray(c.cell) @ transform).tolist()
            trial.positions = (np.asarray(c.positions) @ transform).tolist()
            strains.append(m.evaluate(trial).energy)
        derivative = (strains[0] - strains[1]) / (2 * h) / abs(np.linalg.det(c.cell))
        assert p.stress[axis, axis] == pytest.approx(derivative, abs=2e-7)
    assert np.allclose(p.forces.sum(axis=0), 0, atol=1e-10)
    energies = []
    for sign in (1, -1):
        transform = np.eye(3); transform[0, 1] = sign * h
        trial = deepcopy(c)
        trial.cell = (np.asarray(c.cell) @ transform).tolist()
        trial.positions = (np.asarray(c.positions) @ transform).tolist()
        energies.append(m.evaluate(trial).energy)
    assert p.stress[0, 1] == pytest.approx((energies[0]-energies[1])/(2*h)/abs(np.linalg.det(c.cell)), abs=2e-7)


def training_dataset(truth):
    configurations = []
    for k, distance in enumerate(np.linspace(1.1, 3.5, 12)):
        c = configuration(distance)
        c.name = f"pair_{k}"
        c.split = "validation" if k % 4 == 0 else "train"
        c.group = c.name
        p = truth.evaluate(c)
        c.energy = p.energy; c.forces = p.forces.tolist()
        configurations.append(c)
    return Dataset(configurations, "synthetic model energy zero")


def test_selective_fit_and_holdout():
    base = model(); truth = deepcopy(base)
    r = np.asarray(base.radial_grid)
    truth.pair["He-Pd"] = (np.asarray(base.pair["He-Pd"]) + 0.15 * np.maximum(r - 0.5, 0)**2 * (4-r)**2).tolist()
    data = training_dataset(truth)
    options = FitOptions(embedding_keys=[], pair_keys=["He-Pd"], knots=8, regularization=1e-8)
    result = fit_model(base, data, options)
    assert result.report["optimizer_success"]
    assert result.model.embedding == base.embedding
    assert result.model.density == base.density
    assert result.model.pair["Pd-Pd"] == base.pair["Pd-Pd"]
    assert result.report["metrics"]["validation"]["energy_per_atom_rmse"] < 0.015
    assert result.report["metrics"]["validation"]["force_rmse"] < 0.03
    assert result.model.metadata["validated"] is False
    changed = deepcopy(data)
    for c in changed.configurations:
        if c.split == "validation":
            c.energy += 100
    other = fit_model(base, changed, options)
    assert np.allclose(result.coefficients, other.coefficients)  # holdout never enters fit


def test_density_fit_is_opt_in():
    base = model(); truth = deepcopy(base)
    r = np.asarray(base.radial_grid)
    truth.density["Pd<-He"] = (np.asarray(base.density["Pd<-He"]) - 0.003 * np.maximum(r-0.5, 0)**2 * (4-r)**2).tolist()
    data = training_dataset(truth)
    options = FitOptions(embedding_keys=[], pair_keys=[], density_keys=["Pd<-He"], knots=5, regularization=1e-8,
                         coefficient_bound=0.1, max_evaluations=80)
    result = fit_model(base, data, options)
    assert result.report["optimizer_success"]
    assert result.model.density["He<-Pd"] == base.density["He<-Pd"]
    assert result.report["metrics"]["validation"]["energy_per_atom_rmse"] < 0.01


def test_validation_and_missing_metrics():
    data = training_dataset(model())
    data.configurations[1].group = data.configurations[0].group
    with pytest.raises(ValueError, match="crosses splits"):
        data.validate()
    target = {"metric": "lc_Pd", "value": 3.89, "tolerance": 0.02, "units": "Angstrom",
              "source": "synthetic fixture", "conditions": "0 K"}
    assert compare_metrics({}, [target])[0]["passed"] is False
    assert compare_metrics({"lc_Pd": 3.9}, [target])[0]["passed"] is True


def test_compose_keeps_primary_and_identifies_missing_cross_terms():
    primary = model()
    secondary = seed_model([Element("Pd", 46, 106.42, 3.89, "FCC"), Element("Ni", 28, 58.6934, 3.52, "FCC")],
                           cutoff=4, rho_max=10, nr=81, nrho=121)
    result = compose_model(primary, secondary, ["Pd", "He", "Ni"])
    assert result.style == "eam/he"
    assert result.embedding["Pd"] == primary.embedding["Pd"]
    assert result.density["Pd<-He"] == primary.density["Pd<-He"]
    assert result.pair["He-Pd"] == primary.pair["He-Pd"]
    assert "pair:He-Ni" in result.metadata["provisional_functions"]
    assert "density:Ni<-He" in result.metadata["provisional_functions"]
    assert len(result.density) == 9


def test_revision_review_testing_handoff(tmp_path):
    base = model()
    root = init_project(tmp_path / "fit", base.elements, baseline=base)
    training_dataset(base).save(root / "dataset.json")
    (root / "fit_options.json").write_text(json.dumps(asdict(FitOptions(embedding_keys=[], pair_keys=["He-Pd"], knots=4))))
    revision = fit_revision(root, "r001")
    assert (revision / "candidate.eam.he").exists()
    with pytest.raises(FileExistsError):
        fit_revision(root, "r001")
    with pytest.raises(ValueError, match="human review"):
        main(["handoff", "--revision", str(revision), "--root", str(tmp_path / "testing")])
    record_review(revision, "approved_for_testing", "test reviewer", "Synthetic regression only")
    assert main(["handoff", "--revision", str(revision), "--root", str(tmp_path / "testing")]) == 0
    cfg = json.loads((tmp_path / "testing/pot_inputs/eam_candidate.json").read_text())
    assert cfg["potential"]["pot_lines"].startswith("pair_style eam/he")
    assert cfg["potential"]["single_elements"] == ["Pd"]
    assert cfg["potential"]["gases"] == ["He"]
    with pytest.raises(FileExistsError):
        export_testing_project(base, tmp_path / "testing")


def test_reference_preparation_and_lammps_import(tmp_path):
    c = configuration(periodic=True)
    root = prepare_lammps(c, tmp_path / "lmp", {"Pd": 1, "He": 2}, {"Pd": 106.42, "He": 4}, "pair_style zero 4\npair_coeff * *")
    (root / "energy.txt").write_text("-3.5\n")
    (root / "stress.txt").write_text("1 2 3 0.1 0.2 0.3\n")
    (root / "forces.dump").write_text("ITEM: ATOMS id fx fy fz\n2 -1 -2 -3\n1 1 2 3\n")
    result = import_reference(root)
    assert result.energy == -3.5
    assert np.allclose(result.forces, [[1, 2, 3], [-1, -2, -3]])
    assert "output_sha256" in result.provenance
    templates = tmp_path / "templates"; templates.mkdir()
    for name in ("INCAR", "KPOINTS", "POTCAR"):
        (templates / name).write_text("user input\n")
    vasp = prepare_vasp(c, tmp_path / "vasp", templates, ["He", "Pd"])
    metadata = json.loads((vasp / "reference.json").read_text())
    assert metadata["permutation"] == [1, 0]
    assert "He Pd" in (vasp / "POSCAR").read_text()


def test_external_execution_records_failure(tmp_path):
    root = prepare_lammps(configuration(), tmp_path / "run", {"Pd": 1, "He": 2}, {"Pd": 106.42, "He": 4}, "pair_style zero 4\npair_coeff * *")
    with pytest.raises(RuntimeError, match="failed"):
        run_reference(root, [sys.executable, "-c", "raise SystemExit(3)"], timeout=10)
    assert json.loads((root / "execution.json").read_text())["returncode"] == 3
    with pytest.raises(ValueError, match="did not complete"):
        import_reference(root)


def test_invalid_density_and_overlap_are_rejected():
    m = model(); c = configuration()
    c.positions[1] = c.positions[0][:]
    with pytest.raises(ValueError, match="overlapping"):
        m.evaluate(c)
    m.density_grid = np.linspace(-0.001, 0.001, 121).tolist()
    with pytest.raises(ValueError, match="outside embedding"):
        m.evaluate(configuration())


def test_zhou_functions_and_parameter_fit():
    assert zhou_cutoff(np.array([2.0, 4.0]), 2.0, 4.0)[0] == pytest.approx(0.9)
    assert zhou_cutoff(np.array([2.0, 4.0]), 2.0, 4.0)[1] == 0
    p = {"F0": 1.4, "F2": 0.02, "F3": 0.001, "rho0": -1.0}
    v = analytic_values("zhou_he_embedding", [-2, -1, -0.5, 0, 1], p)
    assert np.allclose(v, [1.4, 1.4, .7, 0, .021])
    base = model()
    term = AnalyticTerm("pair", "He-Pd", "zhou_repulsive_pair",
                        {"E0": .05, "alpha": 2., "r0": 2., "rs": 3., "rc": 4.},
                        {"E0": [.001, .5]}, "synthetic analytic fixture")
    truth_term = deepcopy(term); truth_term.parameters["E0"] = .15
    truth = apply_analytic_terms(base, [truth_term])
    result = fit_analytic_model(base, training_dataset(truth), [term], FitOptions(regularization=1e-10))
    assert result.report["optimizer_success"]
    assert result.coefficients[0] == pytest.approx(.15, abs=1e-5)
    assert result.model.embedding == base.embedding
    assert result.model.density == base.density


def test_snapshot_groups_reproducibility_and_unlabelled_data():
    args = dict(name="PdNiHe", host="Pd", lattice_constant=3.8, composition={"Pd": .75, "Ni": .25},
                gas="He", vacancy=True, repetitions=2, displacement=.01, seed=24, split="validation")
    snapshots = crystal_snapshots(**args)
    assert snapshots == crystal_snapshots(**args)
    assert len(snapshots) == 3 and len(snapshots[0].species) == 32
    assert all(c.group == "PdNiHe" and c.split == "validation" for c in snapshots)
    Dataset(snapshots, "geometry only").validate(require_targets=False)
    with pytest.raises(ValueError, match="no fitting targets"):
        Dataset(snapshots, "geometry only").validate()


def test_testing_config_is_valid_for_ternary(tmp_path):
    from gz_toolkit.potential_testing.config import load_potential_config
    from gz_toolkit.potential_testing.validate import validate_config
    m = seed_model([Element("Pd", 46, 106.42, 3.89, "FCC"), Element("He", 2, 4, 0, "none"),
                    Element("Ni", 28, 58.69, 3.52, "FCC")], style="eam/he", rho_min=-10)
    root = export_testing_project(m, tmp_path / "ternary")
    cfg = load_potential_config(root / "pot_inputs/eam_candidate.json")
    assert validate_config(cfg) == []
    assert cfg.potential.single_elements == ["Pd", "Ni"]
    assert cfg.potential.type_map == {"Pd": 1, "He": 2, "Ni": 3}
    assert cfg.potential.gas_with_interstitial is False
