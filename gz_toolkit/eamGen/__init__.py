"""eamGen: human-guided EAM/FS/HE candidate fitting, references and validation."""
from .data import Configuration, Dataset, compare_metrics
from .model import Element, EAMModel, Prediction, seed_model, pair_key, density_key
from .composition import compose_model
from .fitting import FitOptions, FitResult, fit_model, evaluate_dataset
from .project import init_project, fit_revision, record_review, export_testing_project
from .references import prepare_lammps, prepare_vasp, run_reference, import_reference
from .sampling import crystal_snapshots
from .analytic import AnalyticTerm, analytic_values, zhou_cutoff, apply_analytic_terms, fit_analytic_model

__all__ = ["Configuration", "Dataset", "Element", "EAMModel", "Prediction", "seed_model", "pair_key", "density_key",
           "compose_model", "FitOptions", "FitResult", "fit_model", "evaluate_dataset", "compare_metrics",
           "init_project", "fit_revision", "record_review", "export_testing_project", "prepare_lammps",
           "prepare_vasp", "run_reference", "import_reference", "crystal_snapshots", "AnalyticTerm",
           "analytic_values", "zhou_cutoff", "apply_analytic_terms", "fit_analytic_model"]
