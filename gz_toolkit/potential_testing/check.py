"""Pre-flight checks for the manual, error-prone step between writing a
pot_inputs/*.json and submitting jobs: dropping the actual potential files
into ``potentials/<pot_name>/`` and getting ``pot_lines``/``potential_files``
to agree with them.

Two kinds of checks:

  check_potential_files(cfg, root)
      Static check of pot_inputs/<pot>.json against potentials/<pot_name>/,
      before anything has been built. Catches:
        * potential_files entries with no matching file on disk (or the
          reverse: files on disk that no config lists, so they're silently
          never copied anywhere)
        * filenames mentioned in pot_lines (pair_style/pair_coeff) that
          aren't in potential_files - the #1 cause of "works nowhere"
        * for MTP potentials, prints mlip.ini's contents and flags it if the
          .mtp filename it points to doesn't match potential_files

  check_case_dirs(cfg, run_dir)
      After `run` / `build-defects` have built case directories, verify each
      one actually received every file in potential_files, and that the copy
      matches the current source file (catches stale copies after editing a
      potential file post-build).

Both return a list of (level, message) pairs; level is "error" or "info".
"""

from __future__ import annotations

import re
from pathlib import Path

from gz_toolkit.potential_testing.config import PotentialConfig
from gz_toolkit.potential_testing.parallel import _list_case_dirs

Finding = tuple[str, str]  # (level, message) - level is "error" or "info"

_FILENAME_TOKEN_RE = re.compile(r"[\w.\-]+\.[A-Za-z0-9]+$")
_MTP_REF_RE = re.compile(r"[\w.\-]+\.mtp")

# Same work groups status.py walks, minus seakmc (no potential.inc there).
_GROUPS: list[tuple[str, str]] = [
    ("point_defects", "in.relax.lammps"),
    ("elastic", "in.elastic"),
    ("loops", "in.relax.lammps"),
    ("alloy_lc", "in.alloy.lammps"),
    ("gas_complexes", "in.relax.lammps"),
]


def _tokens_in_pot_lines(pot_lines: str) -> set[str]:
    """Filename-looking tokens (has a dot, no path separators) in pot_lines."""
    tokens: set[str] = set()
    for line in pot_lines.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        for tok in stripped.split():
            if "/" in tok or "\\" in tok or tok == "*":
                continue
            if _FILENAME_TOKEN_RE.match(tok):
                tokens.add(tok)
    return tokens


def check_potential_files(cfg: PotentialConfig, root: str | Path = ".") -> list[Finding]:
    """Cross-check pot_lines / potential_files against potentials/<pot_name>/."""
    findings: list[Finding] = []
    pot = cfg.potential
    src_dir = Path(root) / "potentials" / pot.pot_name

    if not src_dir.is_dir():
        findings.append((
            "error",
            f"potentials/{pot.pot_name}/ does not exist - create it and drop the "
            f"potential file(s) there ({', '.join(pot.potential_files) or 'none listed'})."
        ))
        return findings

    on_disk = {p.name for p in src_dir.iterdir() if p.is_file() and not p.name.startswith(".")}
    listed = set(pot.potential_files)

    for fname in sorted(listed - on_disk):
        findings.append((
            "error",
            f"potential_files lists '{fname}' but potentials/{pot.pot_name}/{fname} "
            f"does not exist."
        ))
    for fname in sorted(on_disk - listed):
        findings.append((
            "error",
            f"potentials/{pot.pot_name}/{fname} exists on disk but is not listed in "
            f"potential_files - it will never be copied into a case directory."
        ))

    # pot_lines cross-check - the classic "filename typo'd in pair_coeff" bug.
    mentioned = _tokens_in_pot_lines(pot.pot_lines)
    for tok in sorted(mentioned - listed):
        findings.append((
            "error",
            f"pot_lines references '{tok}' but it is not in potential_files - "
            f"either add it there or fix the typo in pot_lines."
        ))

    # MTP: mlip.ini carries its own pointer to the .mtp file. We can't fully
    # validate it (format is potential-specific) so surface its contents and
    # flag an obvious name mismatch for the user to eyeball.
    mtp_files = [f for f in pot.potential_files if f.lower().endswith(".mtp")]
    ini_files = [f for f in pot.potential_files if f.lower().endswith(".ini")]
    for ini_name in ini_files:
        ini_path = src_dir / ini_name
        if not ini_path.is_file():
            continue
        ini_text = ini_path.read_text(encoding="utf-8", errors="ignore")
        findings.append((
            "info",
            f"contents of {ini_name} (check the MTP filename it points to by hand):\n"
            + ini_text.rstrip()
        ))
        referenced = set(_MTP_REF_RE.findall(ini_text))
        if mtp_files and referenced and not referenced & set(mtp_files):
            findings.append((
                "error",
                f"{ini_name} points to {sorted(referenced)}, but potential_files "
                f"only lists {mtp_files} - {ini_name} will load the wrong potential."
            ))

    return findings


def _all_case_dirs(cfg: PotentialConfig, run_dir: str | Path) -> list[Path]:
    root = Path(run_dir).resolve()
    pot_dir = root / cfg.potential.pot_name
    dirs: list[Path] = []
    ref_root = pot_dir / "reference"
    if ref_root.is_dir():
        dirs.extend(p for p in ref_root.iterdir() if p.is_dir())
    for group, input_name in _GROUPS:
        dirs.extend(d for d, _ in _list_case_dirs(pot_dir, group, input_name))
    return dirs


def check_case_dirs(cfg: PotentialConfig, run_dir: str | Path = ".") -> list[Finding]:
    """Verify built case directories actually hold the files potential_files promises."""
    findings: list[Finding] = []
    pot = cfg.potential
    root = Path(run_dir).resolve()
    src_dir = root / "potentials" / pot.pot_name

    case_dirs = _all_case_dirs(cfg, run_dir)
    if not case_dirs:
        return findings  # nothing built yet - nothing to check

    for fname in pot.potential_files:
        src = src_dir / fname
        src_size = src.stat().st_size if src.is_file() else None
        for case_dir in case_dirs:
            dst = case_dir / fname
            rel = dst.relative_to(root)
            if not dst.is_file():
                findings.append(("error", f"{rel} is missing - potential_files lists '{fname}'."))
            elif src_size is not None and dst.stat().st_size != src_size:
                findings.append((
                    "error",
                    f"{rel} differs in size from potentials/{pot.pot_name}/{fname} - "
                    f"looks like a stale copy (source was edited after this dir was built)."
                ))

    return findings


def check_project(pot_inputs_dir: str | Path = "pot_inputs", run_dir: str | Path = ".") -> dict[str, list[Finding]]:
    """Run both checks for every config in pot_inputs_dir. {pot_name: [(level, msg), ...]}."""
    from gz_toolkit.potential_testing.config import load_potential_configs

    root = Path(run_dir)
    report: dict[str, list[Finding]] = {}
    for cfg in load_potential_configs(pot_inputs_dir):
        findings = check_potential_files(cfg, root=root)
        findings.extend(check_case_dirs(cfg, run_dir=run_dir))
        report[cfg.potential.pot_name] = findings
    return report
