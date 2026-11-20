"""
checklib.py — evaluate the completion checks defined in sim_plan `checks:` for one run dir.

Check types (all paths relative to the run dir; globs allowed in `file`):
  {type: exists,       file: log.lammps}
  {type: missing,      file: core.*}                         # file must NOT exist
  {type: min_lines,    file: dump.final, n: 2010}
  {type: min_size,     file: out.data, bytes: 1000}
  {type: contains,     file: log.lammps, pattern: "Total wall time"}      # regex
  {type: not_contains, file: log.lammps, pattern: "ERROR|Lost atoms"}
  {type: last_line,    file: Seakmc.log, pattern: "finished"}              # regex on last non-empty line
  {type: count,        file: "SPOut/KMC_*_SPs.csv", min: 10}               # number of files matching glob
A run is FINISHED if every check passes, NOT_STARTED if none of the files referenced
exist, otherwise RUNNING/FAILED (UNFINISHED).
"""
from __future__ import annotations

import glob
import re
from pathlib import Path


def _files(run: Path, pattern: str):
    return sorted(Path(p) for p in glob.glob(str(run / pattern)))


def _tail_line(f: Path) -> str:
    with open(f, "rb") as fh:
        fh.seek(0, 2)
        size = fh.tell()
        fh.seek(max(0, size - 8192))
        lines = [l for l in fh.read().decode(errors="replace").splitlines() if l.strip()]
    return lines[-1] if lines else ""


def run_check(run: Path, c: dict) -> tuple[bool, str]:
    t = c["type"]
    fs = _files(run, c["file"])
    if t == "missing":
        return (not fs), f"{c['file']} should not exist"
    if t == "count":
        return len(fs) >= int(c.get("min", 1)), f"{len(fs)} files match {c['file']} (need {c.get('min', 1)})"
    if not fs:
        return False, f"{c['file']} not found"
    f = fs[-1]
    if t == "exists":
        return True, ""
    if t == "min_size":
        return f.stat().st_size >= int(c["bytes"]), f"{f.name} smaller than {c['bytes']} B"
    if t == "min_lines":
        with open(f, "rb") as fh:
            n = sum(1 for _ in fh)
        return n >= int(c["n"]), f"{f.name} has {n} lines (need {c['n']})"
    if t in ("contains", "not_contains"):
        rx = re.compile(c["pattern"])
        with open(f, errors="replace") as fh:
            hit = any(rx.search(line) for line in fh)
        ok = hit if t == "contains" else not hit
        return ok, f"{f.name} {'lacks' if t == 'contains' else 'contains'} /{c['pattern']}/"
    if t == "last_line":
        return bool(re.search(c["pattern"], _tail_line(f))), f"{f.name} last line lacks /{c['pattern']}/"
    raise ValueError(f"unknown check type {t}")


def status(run: Path, checks: list[dict]) -> tuple[str, list[str]]:
    if not checks:
        return "UNKNOWN", ["no checks defined in sim_plan"]
    fails = []
    for c in checks:
        ok, why = run_check(run, c)
        if not ok:
            fails.append(why)
    if not fails:
        return "FINISHED", []
    touched = any(_files(run, c["file"]) for c in checks if c["type"] not in ("missing",))
    return ("UNFINISHED" if touched else "NOT_STARTED"), fails
