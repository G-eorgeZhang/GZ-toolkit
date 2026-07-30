"""High-level driver: prepare every potential's reference dirs and emit a launcher.

For backwards compatibility, ``run_workflow(config, run_dir)`` is preserved as
the single-potential variant used by older tests.
"""

from __future__ import annotations

from pathlib import Path

from gz_toolkit.potential_testing.config import (
    PotentialConfig,
    load_potential_configs,
)
from gz_toolkit.potential_testing.parallel import emit_jobs
from gz_toolkit.potential_testing.pipeline import (
    build_defects_post_reference,
    prepare_defect_stage,
    prepare_reference_stage,
    summarize_stage,
)
from gz_toolkit.potential_testing.summary import write_wide_summary


def run_all(pot_inputs_dir: str | Path = "pot_inputs", run_dir: str | Path = ".") -> Path:
    """Set up every potential's reference dirs and emit ``submit_all.sh``.

    Returns the path of ``submit_all.sh``. The user runs it (e.g.
    ``bash submit_all.sh``) to actually push jobs to the scheduler.
    """
    project_root = Path(run_dir).resolve()
    pot_inputs = Path(pot_inputs_dir).resolve() if not Path(pot_inputs_dir).is_absolute() \
        else Path(pot_inputs_dir)
    configs = load_potential_configs(pot_inputs)
    for cfg in configs:
        prepare_reference_stage(cfg, run_dir=project_root)
    return emit_jobs(configs, project_root, pot_inputs)


def summarize_all(pot_inputs_dir: str | Path = "pot_inputs", run_dir: str | Path = ".") -> Path:
    """Aggregate every potential's results into one wide CSV."""
    return write_wide_summary(pot_inputs_dir=pot_inputs_dir, run_dir=run_dir)


def run_workflow(config: PotentialConfig, run_dir: str | Path = ".") -> Path:
    """Single-potential helper kept for older tests."""
    prepare_reference_stage(config, run_dir=run_dir)
    prepare_defect_stage(config, run_dir=run_dir)
    return summarize_stage(config, run_dir=run_dir)
