"""gz_toolkit.potential_testing — LAMMPS potential testing workflow scaffolding.

Top-level entry points the user typically calls from their project's
``main.py`` and ``summarize.py``::

    from gz_toolkit.potential_testing import run_all, summarize_all

    run_all(pot_inputs_dir="pot_inputs", run_dir=".")
    summarize_all(pot_inputs_dir="pot_inputs", run_dir=".")
"""

from gz_toolkit.potential_testing.config import (
    DEFAULT_POINT_DEFECTS,
    HPCOptions,
    PotentialConfig,
    PotentialMetadata,
    WorkflowOptions,
    expand_gas_pairs,
    load_potential_config,
    load_potential_configs,
    save_potential_config,
)
from gz_toolkit.potential_testing.project import init_potential_project, default_config, write_driver_scripts
from gz_toolkit.potential_testing.pipeline import (
    build_defects_post_reference,
    prepare_defect_stage,
    prepare_reference_stage,
    summarize_stage,
)
from gz_toolkit.potential_testing.parallel import emit_jobs, scatter_cases
from gz_toolkit.potential_testing.seakmc_runner import emit_seakmc_inputs
from gz_toolkit.potential_testing.elastic import build_elastic_cases, render_elastic_inputs
from gz_toolkit.potential_testing.summary import (
    build_pot_tag_json,
    write_summary,
    write_wide_summary,
    write_grouped_summaries,
)
from gz_toolkit.potential_testing.validate import validate_config, validate_project
from gz_toolkit.potential_testing.check import check_potential_files, check_case_dirs, check_project
from gz_toolkit.potential_testing.status import classify_case, collect_status, format_status_table
from gz_toolkit.potential_testing.wizard import run_init_wizard
from gz_toolkit.potential_testing.workflow import run_all, run_workflow, summarize_all

__all__ = [
    # Config
    "DEFAULT_POINT_DEFECTS",
    "PotentialConfig",
    "PotentialMetadata",
    "WorkflowOptions",
    "HPCOptions",
    "load_potential_config",
    "load_potential_configs",
    "save_potential_config",
    "expand_gas_pairs",
    # Project bootstrap
    "init_potential_project",
    "default_config",
    "write_driver_scripts",
    # Pipeline stages
    "prepare_reference_stage",
    "prepare_defect_stage",
    "build_defects_post_reference",
    "summarize_stage",
    # Parallel job generation
    "emit_jobs",
    "scatter_cases",
    # SEAKMC
    "emit_seakmc_inputs",
    # Elastic constants
    "build_elastic_cases",
    "render_elastic_inputs",
    # Summary
    "build_pot_tag_json",
    "write_summary",
    "write_wide_summary",
    "write_grouped_summaries",
    # Validation / status / wizard
    "validate_config",
    "validate_project",
    "check_potential_files",
    "check_case_dirs",
    "check_project",
    "classify_case",
    "collect_status",
    "format_status_table",
    "run_init_wizard",
    # High-level driver
    "run_all",
    "summarize_all",
    "run_workflow",
]
