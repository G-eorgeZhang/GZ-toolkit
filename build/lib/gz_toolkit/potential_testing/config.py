"""Configuration schema for potential testing workflows.

A project contains a `pot_inputs/` directory with one JSON per potential.
Each file is loaded into a `PotentialConfig` (potential + workflow + hpc).
Defaults assume:
  * classical CPU potential
  * point defects only (loops / dislocations / alloys / gas complexes opt-in)
  * scatter-style submission
  * SEAKMC disabled
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Any


# Default catalog of point defects that are always built when run_energies is on.
DEFAULT_POINT_DEFECTS: list[str] = [
    "1vac",
    "2vac1nn",
    "2vac2nn",
    "2vac3nn",
    "2vac4nn",
    "1int_tetra",
    "1int_octa",
    "dumbbell_100",
    "dumbbell_111",
    "dumbbell_110",
    "2int",
    "3int",
]


@dataclass
class PotentialMetadata:
    pot_name: str
    pot_lines: str
    type_map: dict[str, int]
    masses: dict[int, float]
    single_elements: list[str] = field(default_factory=list)
    gases: list[str] = field(default_factory=list)
    # gas-defect complex enumeration ranges.
    gas_max_m: int = 4   # max number of vacancies (or interstitials) in the complex
    gas_max_n: int = 4   # max number of gas atoms in the complex
    # User-selected (metal, gas) pairs for which gas-defect complexes should be built.
    # Empty list -> use the cartesian product of single_elements x gases.
    gas_complex_pairs: list[list[str]] = field(default_factory=list)
    # Whether (gas)n(V)m and/or (gas)n(I)m complexes are wanted.
    gas_with_vacancy: bool = True
    gas_with_interstitial: bool = True
    two_element_suites: list[dict[str, Any]] = field(default_factory=list)
    # Single-point compositions for 3+ (or 2) element random-substitution alloys.
    # Each entry: {"composition": {el: at.% or null, ...}, "structure": "bcc",
    # "ordering": "random"}. See `resolve_composition`.
    multi_element_suites: list[dict[str, Any]] = field(default_factory=list)
    crystal_structures: dict[str, dict[str, Any]] = field(default_factory=dict)
    lc_initial: float = 2.85
    size_single: int = 20
    size_alloy: int = 20
    # Elastic-constant supercell replications. Kept small: the ELASTIC
    # workflow runs 12 deformations per case, and Cij converges fast with size.
    size_elastic: int = 4
    # Alloy elastic cells need to be a bit larger so the random composition
    # rounds reasonably; still one random realization -> treat Cij as estimate.
    size_elastic_alloy: int = 6
    potential_files: list[str] = field(default_factory=list)
    # Type of potential — selects HPC defaults and SEAKMC compatibility.
    # "classical" (eam/meam/hybrid), "ml" (mtp/snap), "universal" (chgnet/mace).
    kind: str = "classical"
    # Whether/where to copy this pot's <pot_name>.json into gz_toolkit/pot_infobank/
    # once its tests finish (see gz_toolkit.pot_infobank.promote_tag_json):
    #   False (default) -> never copy.
    #   True             -> copy to the cluster's path2gz_toolkit/pot_infobank/
    #                        (hpc.cluster must be set; see cluster_info/<CLUSTER>/meta.json).
    #   a path string    -> copy there verbatim, no validation.
    promote_to_infobank: bool | str = False


@dataclass
class HPCOptions:
    mode: str = "scatter"  # "scatter" or "batch"
    cluster: str | None = None
    partition_key: str | None = None
    walltime: str = "24:00:00"
    nodes: int = 1
    ntasks: int = 40
    output_name: str = "run.job"
    batch_max_per_file: int = 50
    # Hard upper limit on how many sims one batch .job can drive (site-imposed).
    # When a work group exceeds this, the group splits into _1, _2, ... files.
    max_jobs_per_batch: int = 200
    # Default LAMMPS run command; {{NTASKS}}, {{INPUT}}, {{LAMMPS_BIN}} are placeholders.
    run_command: str = "srun -n {{NTASKS}} {{LAMMPS_BIN}} -in {{INPUT}}"
    # CPU vs GPU pots (universal/ML are usually GPU).
    device: str = "cpu"          # "cpu" or "gpu"
    gpus_per_node: int = 0
    # Extra `module load` lines for the job script.
    modules: list[str] = field(default_factory=list)
    # Optional conda env activation (ignored if None).
    conda_env: str | None = None
    # Full path to a custom lmp_mpi binary (e.g. a custom-compiled MTP build).
    # Falls back to plain "lmp_mpi" on PATH if not set.
    lammps_binary_path: str | None = None
    # Directory prepended to LD_LIBRARY_PATH before the run command — needed
    # when lmp_mpi was built without RPATH and depends on shared libs (e.g.
    # the LAMMPS src/ dir of a custom MTP-enabled build).
    lammps_src_path: str | None = None


@dataclass
class WorkflowOptions:
    # Stage gates — kept for backwards compat with existing CLI subcommands.
    run_find_lc: bool = True
    run_reference_supercell: bool = True
    run_defect_generation: bool = True
    run_energies: bool = True
    run_summary: bool = True

    # SEAKMC — migration energy barriers for vacancy / single-interstitial only.
    seakmc_enabled: bool = False
    seakmc_temp_K: float = 600.0
    seakmc_template_path: str | None = None  # if None, gz_toolkit ships a default template
    # Which point defects to send to SEAKMC. By default, single point defects.
    seakmc_defects: list[str] = field(default_factory=lambda: ["1vac", "1int_tetra", "1int_octa"])

    # Defect catalog — point defects always run; everything else is opt-in.
    defect_catalog: list[str] = field(default_factory=lambda: list(DEFAULT_POINT_DEFECTS))
    include_dislocation_loops: bool = False
    include_dislocation_lines: bool = False
    include_alloy_suite: bool = False
    include_gas_complexes: bool = False

    # Alloy point-defect formation energies (two_element_suites and
    # multi_element_suites compositions). Separate opt-in from
    # include_alloy_suite -- materially more expensive (replicas x
    # defect_catalog x composition), since every alloy build is randomized
    # (see structure_ops.py) and needs averaging to be meaningful.
    include_alloy_defects: bool = False
    # Independent random realizations averaged per (composition, case) to
    # smooth out configurational noise from the random solute placement.
    alloy_defect_replicas: int = 3

    # Elastic constants (full 6x6 Cij via the LAMMPS examples/ELASTIC scheme).
    # Core test -> on by default for every element in single_elements.
    include_elastic: bool = True
    # Also compute Cij for each alloy composition in two_element_suites
    # (single random realization on a small cell — an estimate, off by default).
    elastic_for_alloys: bool = False
    # Finite-deformation strain magnitude (the ELASTIC `up` variable).
    elastic_strain: float = 1.0e-6

    # Geometry knobs previously hardcoded in structure_ops.py.
    dumbbell_sep_factor: float = 0.35   # dumbbell separation = factor * lc
    loop_radius_factor: float = 1.5     # SIL loop radius = factor * lc
    # Energy/force tolerance used in every generated minimize command.
    minimize_tol: float = 1.0e-18

    # Parallelization tier:
    #   1 = single job for all potentials
    #   2 = one job per potential
    #   3 = one job per (potential, work-group)
    #   4 = one job per simulation (scatter; auto-split by max_jobs_per_batch)
    paral_degree: int = 4


@dataclass
class PotentialConfig:
    potential: PotentialMetadata
    workflow: WorkflowOptions = field(default_factory=WorkflowOptions)
    hpc: HPCOptions = field(default_factory=HPCOptions)


# ---------------------------------------------------------------------------
# Loaders / dumpers
# ---------------------------------------------------------------------------


def _coerce_mass_keys(payload: dict[str, Any]) -> dict[str, Any]:
    """Convert masses {str: number} to {int: float} after JSON load."""
    potential = payload.get("potential", {})
    masses = potential.get("masses", {})
    if isinstance(masses, dict):
        potential["masses"] = {int(k): float(v) for k, v in masses.items()}
    payload["potential"] = potential
    return payload


def _filter_kwargs(cls, payload: dict[str, Any]) -> dict[str, Any]:
    """Drop keys that aren't dataclass fields so old/new JSONs both load."""
    valid = {f for f in cls.__dataclass_fields__}  # type: ignore[attr-defined]
    return {k: v for k, v in payload.items() if k in valid}


def load_potential_config(path: str | Path) -> PotentialConfig:
    """Load a single per-potential JSON file."""
    cfg_path = Path(path)
    text = cfg_path.read_text(encoding="utf-8")
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Invalid JSON in '{cfg_path}': {exc.msg} at line {exc.lineno}, "
            f"column {exc.colno}. Check for a trailing comma, a Python-style "
            f"comment (# ...), or an unquoted value near that spot."
        ) from exc
    payload = _coerce_mass_keys(payload)
    potential = PotentialMetadata(**_filter_kwargs(PotentialMetadata, payload["potential"]))
    workflow = WorkflowOptions(**_filter_kwargs(WorkflowOptions, payload.get("workflow", {})))
    hpc = HPCOptions(**_filter_kwargs(HPCOptions, payload.get("hpc", {})))
    return PotentialConfig(potential=potential, workflow=workflow, hpc=hpc)


def load_potential_configs(pot_inputs_dir: str | Path = "pot_inputs") -> list[PotentialConfig]:
    """Load every *.json / *.yaml in `pot_inputs_dir` as a separate potential config.

    Files are sorted alphabetically so the order is deterministic.
    YAML support is best-effort — only enabled if PyYAML is installed.
    """
    root = Path(pot_inputs_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"pot_inputs directory not found: {root}")

    configs: list[PotentialConfig] = []
    for f in sorted(root.iterdir()):
        if f.suffix.lower() == ".json":
            configs.append(load_potential_config(f))
        elif f.suffix.lower() in (".yaml", ".yml"):
            try:
                import yaml  # type: ignore
            except ImportError as exc:  # pragma: no cover - optional dep
                raise RuntimeError(
                    f"YAML config '{f}' requires PyYAML; install with `pip install pyyaml`."
                ) from exc
            payload = yaml.safe_load(f.read_text(encoding="utf-8"))
            payload = _coerce_mass_keys(payload)
            potential = PotentialMetadata(**_filter_kwargs(PotentialMetadata, payload["potential"]))
            workflow = WorkflowOptions(**_filter_kwargs(WorkflowOptions, payload.get("workflow", {})))
            hpc = HPCOptions(**_filter_kwargs(HPCOptions, payload.get("hpc", {})))
            configs.append(PotentialConfig(potential=potential, workflow=workflow, hpc=hpc))
    if not configs:
        raise RuntimeError(f"No potential config files found in {root}")
    return configs


def save_potential_config(config: PotentialConfig, path: str | Path) -> None:
    """Serialize a PotentialConfig to JSON."""
    cfg_path = Path(path)
    cfg_path.parent.mkdir(parents=True, exist_ok=True)
    cfg_path.write_text(json.dumps(asdict(config), indent=2), encoding="utf-8")


# ---------------------------------------------------------------------------
# Convenience helpers used by the pipeline
# ---------------------------------------------------------------------------


def expand_gas_pairs(meta: PotentialMetadata) -> list[tuple[str, str]]:
    """Return the (metal, gas) pairs requested for gas-defect complexes.

    If the user provided ``gas_complex_pairs`` in the JSON, only those are used.
    Otherwise the cartesian product (single_elements x gases) is returned.
    """
    if meta.gas_complex_pairs:
        return [(p[0], p[1]) for p in meta.gas_complex_pairs if len(p) == 2]
    return [(m, g) for m in meta.single_elements for g in meta.gases]


def resolve_composition(composition: dict[str, Any], tol: float = 1e-6) -> dict[str, float]:
    """Resolve a (possibly partial) at.% composition into concrete numbers.

    ``composition`` maps element -> at.% (a number) or ``None``/``null`` for
    "whatever is left". Any number of elements may be ``None``: the remaining
    at.% (100 minus the explicit entries) is split evenly across them, e.g.
    ``{"Cr": 3, "Fe": None}`` -> ``{"Cr": 3.0, "Fe": 97.0}`` and
    ``{"Fe": None, "Cr": None, "Ni": None}`` -> equal thirds.

    With no ``None`` entries, the explicit values must already sum to 100
    (within ``tol``). Key order is preserved — callers use the first key as
    the "host" lattice species for cell seeding.
    """
    if len(composition) < 2:
        raise ValueError(f"Composition {composition} needs at least 2 elements for an alloy.")

    explicit = {el: float(v) for el, v in composition.items() if v is not None}
    null_elements = [el for el, v in composition.items() if v is None]

    explicit_sum = sum(explicit.values())
    if null_elements:
        remainder = 100.0 - explicit_sum
        if remainder <= 0:
            raise ValueError(
                f"Composition {composition}: explicit fractions already sum to "
                f"{explicit_sum:g}, leaving nothing for {null_elements}."
            )
        share = remainder / len(null_elements)
        resolved = dict(explicit)
        for el in null_elements:
            resolved[el] = share
    else:
        if abs(explicit_sum - 100.0) > tol:
            raise ValueError(
                f"Composition {composition} sums to {explicit_sum:g}, not 100 "
                f"(no null entry to absorb the difference)."
            )
        resolved = explicit

    return {el: resolved[el] for el in composition}
