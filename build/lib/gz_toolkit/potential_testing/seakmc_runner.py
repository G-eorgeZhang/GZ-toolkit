"""SEAKMC input file generator.

For each point-defect case where SEAKMC is requested, emit:

    <case_dir>/meb/
    ├── input.yaml          # populated from gz_toolkit-shipped (or user-supplied) template
    ├── run_seakmc_p.py     # standard 1-line MPI driver
    └── relaxed.data        # symlink/copy of the relaxed defect data file

Only the fields that are *known to vary per (potential, element)* are filled
in by gz_toolkit:

  * potential.species      — list ordered by gz_toolkit's type_map
  * potential.pair_style
  * potential.FileName     — stripped from `pot_lines`
  * potential.bondlengths  — heuristic from lc + crystal structure
  * potential.coordnums    — heuristic from crystal structure
  * data.FileName          — relaxed defect file
  * kinetic_MC.Temp        — workflow.seakmc_temp_K

All other fields come from a baseline template (default ships with gz_toolkit;
user can override by setting ``workflow.seakmc_template_path`` in the JSON).
"""

from __future__ import annotations

import shutil
from pathlib import Path

from gz_toolkit.potential_testing.config import PotentialConfig


# ---------------------------------------------------------------------------
# Default input.yaml template (BCC defaults; FCC overrides applied programmatically)
# ---------------------------------------------------------------------------

_DEFAULT_TEMPLATE = """\
# Auto-generated SEAKMC input.yaml
system:
  Restart:
    WriteRestart: True
    LoadRestart: True
    AVStep4Restart: 1000
    KMCStep4Restart: 1
  Tolerance: 0.1

kinetic_MC:
  NSteps: 100
  Temp: {temp_K}
  AccStyle: NoAcc

potential:
  species:
{species_block}
  pair_style: {pair_style}
  FileName: {pot_filename}
  bondlengths:
{bondlengths_block}
  coordnums:
{coordnums_block}

data:
  FileName: {data_filename}
  atom_style: atomic
  Relaxed: True
  BoxRelax: False
  boundary: p p p

active_volume:
  Style: defects
  DActive: {dactive}
  DBuffer: 1.5
  DFixed: {dactive}
  FindDefects:
    Method: BLCN
    DCut4Def: 0.1
  PDReduction: True
  DCut4PDR: 3.5
  Overlapping: True
  NMin4AV: 40

spsearch:
  Method: dimer
  NSearch: 10
  NMax4Trans: 1000
  FConv: 1.0e-6
  TrialStepsize: 0.015
  MaxStepsize: 0.05
  DimerSep: 0.005
  TransHorizon: True
  LocalRelax:
    LocalRelax: True
  HandleVN:
    CenterVN: True
    RescaleVN: True
    RescaleStyle4LOGV: LOGNRAS

saddle_point:
  BarrierCut: 5.0
  BarrierMin: 0.0
  BackBarrierMin: 0.0
  ValidSPs:
    CheckConnectivity: True
    RealtimeDelete: True

dynamic_matrix:
  SNC: False
  CalPrefactor: False

defect_bank:
  Recycle: True
  SaveDB: True
  SavePath: DefectBank
  LoadDB: False

force_evaluator:
  Bin: pylammps
  NSteps4Relax: 10000
  Relaxation:
    BoxRelax: False

visual:
  Screen: True
  Log: True
  Write_SP_Summary: True
  Write_Data_SPs:
    Write_KMC_Data: True
"""


_RUN_DRIVER = """\
#!/usr/bin/env python
\"\"\"Wrapper that invokes the SEAKMC parallel runner.\"\"\"

if __name__ == "__main__":
    # The parallel build of SEAKMC_py exposes a function-style entry point.
    # If your installation is named differently, edit this single line.
    from seakmc_p.SEAKMC import run as seakmc_run
    seakmc_run()
"""


# ---------------------------------------------------------------------------
# Heuristics
# ---------------------------------------------------------------------------


def _bondlength_for_structure(structure: str, lc: float) -> float:
    """First-nearest-neighbour distance for a cubic structure."""
    s = structure.lower()
    if s == "bcc":
        return 0.866 * lc            # sqrt(3)/2
    if s == "fcc":
        return 0.7071 * lc           # sqrt(2)/2
    if s == "hcp":
        return lc                    # a (in-plane) — caller may override
    return lc


def _coord_number_for_structure(structure: str) -> int:
    s = structure.lower()
    if s == "bcc":
        return 8
    if s == "fcc":
        return 12
    if s == "hcp":
        return 12
    return 8


def _extract_pair_style_and_filename(pot_lines: str, potential_files: list[str]) -> tuple[str, str]:
    """Best-effort parse of ``pair_style ...`` and ``pair_coeff ...`` lines."""
    pair_style = "eam"
    pair_filename = ""
    for line in pot_lines.splitlines():
        line = line.strip()
        if line.lower().startswith("pair_style"):
            parts = line.split(None, 1)
            if len(parts) > 1:
                pair_style = parts[1].strip()
        elif line.lower().startswith("pair_coeff"):
            tokens = line.split()
            for tok in tokens[1:]:
                if "*" in tok or tok in ("pair_coeff",):
                    continue
                if tok in potential_files or "." in tok and "/" not in tok:
                    pair_filename = tok
                    break
    if not pair_filename and potential_files:
        pair_filename = potential_files[0]
    return pair_style, pair_filename


# ---------------------------------------------------------------------------
# Public emitter
# ---------------------------------------------------------------------------


def _bond_block(types: list[int], bondlength: float) -> str:
    """One ``- t1 t1 d`` line per type (homogeneous bond)."""
    return "\n".join(f"    - {t} {t} {bondlength:.4f}" for t in types)


def _coord_block(types: list[int], coordnum: int) -> str:
    return "\n".join(f"    - {t} {coordnum}" for t in types)


def _species_block(species_in_order: list[str]) -> str:
    return "\n".join(f"    - {sp}" for sp in species_in_order)


def emit_seakmc_inputs(
    case_dir: Path,
    config: PotentialConfig,
    element: str,
    structure_name: str,
    relaxed_lc: float,
    relaxed_data_path: Path,
    dactive: float | None = None,
) -> Path:
    """Write ``meb/input.yaml`` and ``meb/run_seakmc_p.py`` next to the relaxed defect.

    Returns the meb directory path.
    """
    meb_dir = case_dir / "meb"
    meb_dir.mkdir(parents=True, exist_ok=True)

    # Copy the relaxed data file into meb/ so SEAKMC and the rest of the
    # workflow are decoupled.
    target_data = meb_dir / relaxed_data_path.name
    if relaxed_data_path.exists():
        shutil.copy2(relaxed_data_path, target_data)

    # Choose template — user override or shipped default.
    template_text = _DEFAULT_TEMPLATE
    user_tpl = config.workflow.seakmc_template_path
    if user_tpl:
        p = Path(user_tpl)
        if p.is_file():
            template_text = p.read_text(encoding="utf-8")

    # Resolve fields.
    species_in_order = sorted(config.potential.type_map, key=lambda e: config.potential.type_map[e])
    types_in_order = [config.potential.type_map[e] for e in species_in_order]
    pair_style, pot_filename = _extract_pair_style_and_filename(
        config.potential.pot_lines, config.potential.potential_files
    )
    bond = _bondlength_for_structure(structure_name, relaxed_lc)
    cn = _coord_number_for_structure(structure_name)
    if dactive is None:
        # ~3 lattice constants is a reasonable default for point defects.
        dactive = max(7.0, 2.7 * relaxed_lc)

    # Render. ``str.format`` is sufficient since all replacement keys are
    # alphanumeric, and the default template uses no other curly braces.
    rendered = template_text.format(
        temp_K=config.workflow.seakmc_temp_K,
        species_block=_species_block(species_in_order),
        pair_style=pair_style,
        pot_filename=pot_filename,
        bondlengths_block=_bond_block(types_in_order, bond),
        coordnums_block=_coord_block(types_in_order, cn),
        data_filename=target_data.name,
        dactive=f"{dactive:.2f}",
    )
    (meb_dir / "input.yaml").write_text(rendered, encoding="utf-8")
    (meb_dir / "run_seakmc_p.py").write_text(_RUN_DRIVER, encoding="utf-8")

    return meb_dir


# Backward-compat helper used by the pipeline to decide whether SEAKMC should
# run for a given defect.
def is_seakmc_eligible(case_name: str, workflow_seakmc_defects: list[str]) -> bool:
    return case_name in set(workflow_seakmc_defects)
