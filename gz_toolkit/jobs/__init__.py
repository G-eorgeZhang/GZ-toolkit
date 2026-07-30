"""
gz_toolkit.jobs - HPC job script generation for simulation workflows.

Linux-only:
  This module targets Linux HPC environments (SLURM-style job scripts).

Provides tools for:
  - Creating SLURM job scripts from templates.
  - Generating scatter or patch job scripts (no auto-submission).
"""

from gz_toolkit.jobs.directory import DirectoryManager
from gz_toolkit.jobs.template_engine import (
    JobTemplate,
    load_cluster_info,
    load_cluster_meta,
    list_available_machines,
)
from gz_toolkit.jobs.submission import ScatterSubmitter, PatchSubmitter
from gz_toolkit.jobs.tree import generate_folder_tree

__all__ = [
    "DirectoryManager",
    "JobTemplate",
    "load_cluster_info",
    "load_cluster_meta",
    "list_available_machines",
    "ScatterSubmitter",
    "PatchSubmitter",
    "generate_folder_tree",
]
