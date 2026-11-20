#!/usr/bin/env python3
"""
make_structures.py — step 1: build the reference chain and write final structures into run dirs.

COPY this template to ROOT/workflow/make_structures.py and fill in the TODO(gz) blocks with
real GZ_toolkit calls (read the toolkit's API first — never invent function names).

Rules this script follows (from the materials-simulation skill):
  * ONE change per stage. Stage 0 = perfect crystal. Each later stage reads the previous
    stage's file and applies exactly one modification.
  * Reference stages are written to reffiles/ (ref0_*, ref1_*, ...). The LAST stage writes
    directly into every run dir as `structures.file` from sim_plan.
  * Each stage is run separately (--stage N) so the user can inspect it before the next.

    python workflow/make_structures.py --stage 0      # perfect crystal  -> reffiles/ref0_perfect.data
    python workflow/make_structures.py --stage 1      # + one defect     -> reffiles/ref1_*.data
    python workflow/make_structures.py --stage final  # per-run change   -> <run dir>/structure.data
    python workflow/make_structures.py --stage final --workers 16
"""
import argparse
import csv
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load_plan():
    for n in ("sim_plan.yaml", "sim_plan.yml", "sim_plan.json"):
        f = ROOT / n
        if f.is_file():
            if f.suffix == ".json":
                return json.loads(f.read_text())
            import yaml
            return yaml.safe_load(f.read_text())
    sys.exit("no sim_plan in " + str(ROOT))


PLAN = load_plan()
for _name, _path in (PLAN.get("environment", {}).get("packages") or {}).items():
    if _path and _path != "pip":
        sys.path.insert(0, _path)

# TODO(gz): import the toolkit, e.g.  import gz_toolkit as gz   (check the real package name)
REF = ROOT / PLAN.get("reffiles", {}).get("dir", "reffiles")
OUTNAME = PLAN.get("structures", {}).get("file", "structure.data")


# ------------------------------------------------------------------------ stages
def stage0_perfect():
    """Perfect crystal. Parameters come from the user, not from defaults."""
    out = REF / "ref0_perfect.data"
    # TODO(gz): build the lattice: element, lattice constant (relaxed with THIS potential),
    #           orientation (e.g. x=[11-2] y=[111] z=[1-10]), repeats, write LAMMPS data.
    raise NotImplementedError("fill in stage 0 with GZ_toolkit calls")
    return out


def stage1_reference():
    """Exactly one modification of stage 0, e.g. insert a screw dislocation."""
    src, out = REF / "ref0_perfect.data", REF / "ref1_defect.data"
    # TODO(gz): read src, apply ONE change, write out.
    raise NotImplementedError("fill in stage 1")
    return out


def final_for_run(row: dict) -> str:
    """Per-run modification of the last reference, written into the run dir.
    `row` has relpath, every tree variable (as strings), run, run_index."""
    src = REF / "ref1_defect.data"
    dst = ROOT / row["relpath"] / OUTNAME
    # Example of mapping variables -> structure parameters:
    #   cluster = row["cluster"]            # e.g. "He2V1"
    #   site    = SITES[int(row["run_index"]) - 1]
    # TODO(gz): read src, apply ONE change for this run, write dst.
    raise NotImplementedError("fill in the final stage")
    return str(dst)


# ------------------------------------------------------------------------ checks
def report(path: Path):
    """Print what the user needs to verify a structure: atoms, types, box. Works on LAMMPS data."""
    natoms = ntypes = None
    box = {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            s = line.split()
            if "atoms" in s and len(s) == 2:
                natoms = int(s[0])
            elif "atom" in s and "types" in s:
                ntypes = int(s[0])
            elif len(s) >= 4 and s[2].endswith("lo") and s[3].endswith("hi"):
                box[s[2][0]] = float(s[1]) - float(s[0])
            elif len(s) >= 3 and s[-3:] == ["xy", "xz", "yz"]:
                box["tilt"] = tuple(map(float, s[:3]))
            if i > 40:
                break
    print(f"  {path.relative_to(ROOT)}: atoms={natoms} types={ntypes} box={box}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, help="0, 1, ..., or final")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--only", default=None, help="relpath prefix filter for the final stage")
    a = ap.parse_args()

    if a.stage == "0":
        report(stage0_perfect())
    elif a.stage == "1":
        report(stage1_reference())
    elif a.stage == "final":
        rows = list(csv.DictReader(open(ROOT / "manifest.csv", newline="")))
        if a.only:
            rows = [r for r in rows if r["relpath"].startswith(a.only)]
        if a.workers > 1:
            with ProcessPoolExecutor(a.workers) as ex:
                outs = list(ex.map(final_for_run, rows))
        else:
            outs = [final_for_run(r) for r in rows]
        print(f"wrote {len(outs)} structures; sample:")
        for o in outs[:3]:
            report(Path(o))
    else:
        sys.exit("unknown stage")


if __name__ == "__main__":
    main()
