"""
simflow_common.py — shared helpers for the materials-simulation workflow scripts.

Everything is driven by two files in the simulation ROOT directory:
  sim_plan.yaml (or sim_plan.json)  — the plan: variables, reffiles, job, checks, summary
  manifest.csv                      — one row per run dir, written by make_tree.py

Only the Python standard library is required. PyYAML is used if installed;
otherwise write the plan as sim_plan.json.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import sys
from pathlib import Path

PLAN_NAMES = ("sim_plan.yaml", "sim_plan.yml", "sim_plan.json")
MANIFEST = "manifest.csv"
PLACEHOLDER = re.compile(r"\{\{\s*([A-Za-z_][A-Za-z0-9_]*)\s*\}\}")


# --------------------------------------------------------------------------- plan
def find_root(start: str | os.PathLike | None = None) -> Path:
    """Walk upward from `start` (default cwd) until a sim_plan file is found."""
    p = Path(start or os.getcwd()).resolve()
    for d in [p, *p.parents]:
        if any((d / n).is_file() for n in PLAN_NAMES):
            return d
    sys.exit(f"[simflow] no sim_plan.yaml/json found in {p} or any parent directory")


def load_plan(root: Path) -> dict:
    for n in PLAN_NAMES:
        f = root / n
        if not f.is_file():
            continue
        text = f.read_text(encoding="utf-8")
        if f.suffix == ".json":
            return json.loads(text)
        try:
            import yaml  # type: ignore
        except ImportError:
            sys.exit("[simflow] PyYAML not installed: `pip install pyyaml`, or convert the plan to sim_plan.json")
        return yaml.safe_load(text)
    sys.exit(f"[simflow] no plan file in {root}")


def levels(plan: dict) -> list[dict]:
    return plan.get("tree", {}).get("levels", [])


def run_level(plan: dict) -> int:
    """Depth of the run dirs. root = 0, first variable = 1, ..., runs = n_vars + 1."""
    return len(levels(plan)) + 1


def fmt_dirname(spec: dict, value) -> str:
    pattern = spec.get("dirname", spec["name"] + "{value}")
    return pattern.format(value=value, name=spec["name"])


def run_names(plan: dict) -> list[str]:
    runs = plan.get("tree", {}).get("runs", 1)
    prefix = plan.get("tree", {}).get("run_prefix", "run")
    if isinstance(runs, int):
        return [f"{prefix}{i}" for i in range(1, runs + 1)]
    return [str(r) for r in runs]  # explicit list of run dir names


# ----------------------------------------------------------------------- manifest
def read_manifest(root: Path) -> list[dict]:
    f = root / MANIFEST
    if not f.is_file():
        sys.exit(f"[simflow] {f} missing — run make_tree.py first")
    with open(f, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_manifest(root: Path, rows: list[dict]) -> None:
    if not rows:
        sys.exit("[simflow] manifest would be empty")
    cols = list(rows[0].keys())
    with open(root / MANIFEST, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


def rel_parts(row: dict) -> list[str]:
    return Path(row["relpath"]).parts


def ancestor(row: dict, level: int) -> str:
    """Relative path of this run's directory at `level` (0 = root -> '.')."""
    parts = rel_parts(row)
    return "." if level == 0 else Path(*parts[:level]).as_posix()


# ----------------------------------------------------------------------- template
def stable_seed(text: str, base: int = 0) -> int:
    """Deterministic positive 31-bit seed derived from a string (e.g. the relpath)."""
    h = int(hashlib.sha256(text.encode()).hexdigest()[:8], 16)
    return (h + base) % 2_000_000_000 + 1


def context_for(row: dict, plan: dict) -> dict:
    """All {{KEY}} values available for one run dir."""
    ctx = {k: v for k, v in row.items()}
    ctx["RELPATH"] = row["relpath"]
    ctx["RUN"] = row["run"]
    ctx["RUN_INDEX"] = row["run_index"]
    ctx["SEED"] = str(stable_seed(row["relpath"], int(plan.get("seed_base", 0))))
    ctx["JOBNAME"] = job_name(plan, row["relpath"])
    for k, v in (plan.get("constants") or {}).items():
        ctx.setdefault(k, str(v))
    return ctx


def render(text: str, ctx: dict, where: str = "") -> str:
    missing = sorted({m for m in PLACEHOLDER.findall(text) if m not in ctx})
    if missing:
        raise KeyError(f"unfilled placeholders {missing} in {where}")
    return PLACEHOLDER.sub(lambda m: str(ctx[m.group(1)]), text)


def job_name(plan: dict, relpath: str) -> str:
    prefix = (plan.get("job") or {}).get("name_prefix", plan.get("project", "sim"))
    tail = re.sub(r"[^A-Za-z0-9._-]+", "_", relpath.replace(os.sep, "_").replace("/", "_")).strip("_.")
    return f"{prefix}_{tail}" if tail else prefix


# ------------------------------------------------------------------------- misc
def hms_to_s(t: str) -> int:
    """SLURM time formats: MM, MM:SS, HH:MM:SS, D-HH, D-HH:MM, D-HH:MM:SS."""
    days = 0
    if "-" in t:
        d, t = t.split("-", 1)
        days = int(d)
        parts = [int(x) for x in t.split(":")]
        while len(parts) < 3:
            parts.append(0)
        h, m, s = parts
    else:
        parts = [int(x) for x in t.split(":")]
        if len(parts) == 1:
            h, m, s = 0, parts[0], 0
        elif len(parts) == 2:
            h, (m, s) = 0, parts
        else:
            h, m, s = parts
    return ((days * 24 + h) * 60 + m) * 60 + s


def s_to_slurm(sec: int) -> str:
    d, rem = divmod(int(sec), 86400)
    h, rem = divmod(rem, 3600)
    m, s = divmod(rem, 60)
    return f"{d}-{h:02d}:{m:02d}:{s:02d}" if d else f"{h:02d}:{m:02d}:{s:02d}"


def log(msg: str) -> None:
    print(f"[simflow] {msg}", flush=True)
