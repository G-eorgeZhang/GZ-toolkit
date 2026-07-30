"""
GZ-toolkit (import name: gz_toolkit) — LAMMPS structure manipulation toolkit.

Provides tools for building crystal structures, creating defects
(voids, dislocation loops, Frenkel pairs), and merging atomistic
configurations for LAMMPS molecular dynamics simulations.

Submodules are imported lazily (PEP 562): ``import gz_toolkit`` is cheap,
and heavy dependencies (pymatgen, pandas) are only loaded when the
corresponding attribute is first accessed.
"""

import importlib

__version__ = "0.1.0"

# Map of public attribute -> submodule that provides it.
_LAZY_ATTRS = {
    # Core
    "Modlmp_LmpData": "gz_toolkit.core.modlmp",
    # Crystal builder
    "Gen_crystal": "gz_toolkit.buildmtx.buildstr",
    # Defect builders
    "BCCDefect": "gz_toolkit.defect.bcc_defect",
    "BccDefect": "gz_toolkit.defect.bcc_defect",
    "FCCDefect": "gz_toolkit.defect.fcc_defect",
    "FccDefect": "gz_toolkit.defect.fcc_defect",
    # Jobs module
    "DirectoryManager": "gz_toolkit.jobs",
    "JobTemplate": "gz_toolkit.jobs",
    "ScatterSubmitter": "gz_toolkit.jobs",
    "PatchSubmitter": "gz_toolkit.jobs",
    "generate_folder_tree": "gz_toolkit.jobs",
    # Potential testing module
    "HPCOptions": "gz_toolkit.potential_testing",
    "PotentialConfig": "gz_toolkit.potential_testing",
    "PotentialMetadata": "gz_toolkit.potential_testing",
    "WorkflowOptions": "gz_toolkit.potential_testing",
    "load_potential_config": "gz_toolkit.potential_testing",
    "load_potential_configs": "gz_toolkit.potential_testing",
    "save_potential_config": "gz_toolkit.potential_testing",
    "init_potential_project": "gz_toolkit.potential_testing",
    "run_workflow": "gz_toolkit.potential_testing",
    "run_all": "gz_toolkit.potential_testing",
    "summarize_all": "gz_toolkit.potential_testing",
}

__all__ = list(_LAZY_ATTRS)


def __getattr__(name):
    module_path = _LAZY_ATTRS.get(name)
    if module_path is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = importlib.import_module(module_path)
    value = getattr(module, name)
    globals()[name] = value  # cache for subsequent lookups
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__))
