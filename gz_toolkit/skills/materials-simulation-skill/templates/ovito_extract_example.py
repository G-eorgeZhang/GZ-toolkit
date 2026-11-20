"""
OVITO-based extractor functions for summarize.py (field type `ovito` or `python`).
Copy into reffiles/ (or workflow/) and point a summary field at it:

  summary:
    fields:
      - {name: n_vac, type: ovito, file: "dump.final", script: "reffiles/ovito_extract.py:ws_vacancies",
         reference: "reffiles/ref0_perfect.data"}
      - {name: n_dislo_len, type: ovito, file: "dump.final", script: "reffiles/ovito_extract.py:dxa_length"}

Each function gets (run_dir: pathlib.Path, field: dict) and returns one value (number/str).
Requires the `ovito` Python module (pip install ovito) in the workflow environment.
Verify modifier names/arguments against https://docs.ovito.org/python/ for the installed version.
"""
import glob
from pathlib import Path

from ovito.io import import_file
from ovito.modifiers import (CommonNeighborAnalysisModifier, DislocationAnalysisModifier,
                             WignerSeitzAnalysisModifier)


def _file(run_dir: Path, pattern: str) -> str:
    hits = sorted(glob.glob(str(run_dir / pattern)))
    if not hits:
        raise FileNotFoundError(pattern)
    return hits[-1]


def ws_vacancies(run_dir: Path, field: dict):
    """Number of vacancies from Wigner-Seitz analysis against a perfect reference."""
    pipe = import_file(_file(run_dir, field["file"]))
    ws = WignerSeitzAnalysisModifier()
    root = run_dir
    while not (root / "sim_plan.yaml").exists() and root != root.parent:
        root = root.parent
    ws.reference.load(str(root / field["reference"]))
    pipe.modifiers.append(ws)
    return int(pipe.compute().attributes["WignerSeitz.vacancy_count"])


def cna_fraction(run_dir: Path, field: dict):
    """Fraction of atoms in the structure type named in field['structure'] (default BCC)."""
    pipe = import_file(_file(run_dir, field["file"]))
    pipe.modifiers.append(CommonNeighborAnalysisModifier())
    data = pipe.compute()
    name = field.get("structure", "BCC").upper()
    return data.attributes[f"CommonNeighborAnalysis.counts.{name}"] / data.particles.count


def dxa_length(run_dir: Path, field: dict):
    """Total dislocation line length (Angstrom) from DXA."""
    pipe = import_file(_file(run_dir, field["file"]))
    pipe.modifiers.append(DislocationAnalysisModifier(
        input_crystal_structure=getattr(DislocationAnalysisModifier.Lattice, field.get("lattice", "BCC"))))
    return float(pipe.compute().attributes["DislocationAnalysis.total_line_length"])
