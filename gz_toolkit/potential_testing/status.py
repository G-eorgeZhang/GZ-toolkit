"""One-command progress check for a potential-testing project.

Walks every case directory of every work group and classifies it:

  pending — no log.lammps yet (job not started)
  running — log exists but LAMMPS hasn't reached a terminal marker
  done    — log contains RESULT lines (our summary hooks)
  failed  — log contains ERROR, or finished cleanly without RESULT lines
"""

from __future__ import annotations

from pathlib import Path

from gz_toolkit.potential_testing.config import PotentialConfig
from gz_toolkit.potential_testing.parallel import _list_case_dirs

STATES = ("done", "failed", "running", "pending")

# Work groups walked per potential: (group name, lammps input filename).
_GROUPS: list[tuple[str, str]] = [
    ("point_defects", "in.relax.lammps"),
    ("elastic", "in.elastic"),
    ("loops", "in.relax.lammps"),
    ("alloy_lc", "in.alloy.lammps"),
    ("gas_complexes", "in.relax.lammps"),
    ("seakmc", ""),
]


def classify_case(case_dir: Path) -> str:
    """Classify one case directory by inspecting its log.lammps."""
    log = Path(case_dir) / "log.lammps"
    if not log.is_file():
        return "pending"
    txt = log.read_text(encoding="utf-8", errors="ignore")
    if "ERROR" in txt:
        return "failed"
    if "RESULT " in txt:
        return "done"
    # LAMMPS prints "Total wall time" on every clean exit; reaching it
    # without any RESULT line means our prints never ran -> failed.
    if "Total wall time" in txt:
        return "failed"
    return "running"


def _classify_seakmc(meb_dir: Path) -> str:
    """SEAKMC cases have no log.lammps; look for its own outputs."""
    meb = Path(meb_dir)
    if (meb / "Seakmc_summary.csv").is_file() or (meb / "SPOut").is_dir():
        return "done"
    # any output beyond the inputs we wrote -> probably running
    started = any(p.name not in ("input.yaml", "run_seakmc_p.py")
                  for p in meb.iterdir()) if meb.is_dir() else False
    return "running" if started else "pending"


def collect_status(
    cfg: PotentialConfig,
    run_dir: str | Path = ".",
) -> dict[str, list[tuple[str, str]]]:
    """Return {group: [(relative_case_path, state), ...]} for one potential."""
    root = Path(run_dir).resolve()
    pot_dir = root / cfg.potential.pot_name
    out: dict[str, list[tuple[str, str]]] = {}

    # reference dirs first — they gate everything else
    ref_root = pot_dir / "reference"
    refs: list[tuple[str, str]] = []
    if ref_root.is_dir():
        for el_dir in sorted(ref_root.iterdir()):
            if el_dir.is_dir():
                refs.append((str(el_dir.relative_to(root)), classify_case(el_dir)))
    out["reference"] = refs

    for group, input_name in _GROUPS:
        cases = _list_case_dirs(pot_dir, group, input_name)
        entries: list[tuple[str, str]] = []
        for case_dir, _ in cases:
            state = _classify_seakmc(case_dir) if group == "seakmc" else classify_case(case_dir)
            entries.append((str(Path(case_dir).relative_to(root)), state))
        if entries:
            out[group] = entries
    return out


def format_status_table(
    status_by_pot: dict[str, dict[str, list[tuple[str, str]]]],
    verbose: bool = False,
) -> str:
    """Render the status as a compact text table.

    One row per (potential, group) with done/failed/running/pending counts;
    with ``verbose`` the failed and running case paths are listed underneath.
    """
    lines: list[str] = []
    header = f"{'potential':<20} {'group':<16} {'done':>5} {'fail':>5} {'run':>5} {'wait':>5}"
    lines.append(header)
    lines.append("-" * len(header))
    for pot, groups in status_by_pot.items():
        for group, entries in groups.items():
            counts = {s: 0 for s in STATES}
            for _, state in entries:
                counts[state] = counts.get(state, 0) + 1
            lines.append(
                f"{pot:<20} {group:<16} {counts['done']:>5} {counts['failed']:>5} "
                f"{counts['running']:>5} {counts['pending']:>5}"
            )
            if verbose:
                for path, state in entries:
                    if state in ("failed", "running"):
                        lines.append(f"    [{state}] {path}")
    return "\n".join(lines)
