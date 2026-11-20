#!/usr/bin/env python3
"""
make_tree.py — build the variable directory tree and manifest.csv from sim_plan.

    python workflow/make_tree.py            # create dirs + manifest.csv
    python workflow/make_tree.py --dry-run  # print the tree only, touch nothing

Layout produced (example with 3 variables):
    ROOT/                         level 0
      reffiles/
      T300K/                      level 1  (variable 1)
        rho1e14/                  level 2  (variable 2)
          He1V1/                  level 3  (variable 3)
            run1/ run2/ ...       level 4  (run dirs: simulations happen here)

Every level's values can be a plain list, or a dict keyed by the PARENT dir name when
the values of a level depend on its parent (non-Cartesian trees):
    values: {T300K: [1, 2], T600K: [1, 2, 3]}
Existing directories are never deleted. manifest.csv is (re)written each time.
"""
import argparse
from pathlib import Path

from simflow_common import (find_root, fmt_dirname, levels, load_plan, log,
                            run_names, write_manifest)


def expand(plan):
    specs = levels(plan)
    rows = []

    def recurse(depth, parts, values):
        if depth == len(specs):
            for i, r in enumerate(run_names(plan), start=1):
                row = {"relpath": Path(*parts, r).as_posix() if parts else r}
                row.update(values)
                row["run"] = r
                row["run_index"] = i
                rows.append(row)
            return
        spec = specs[depth]
        vals = spec["values"]
        if isinstance(vals, dict):  # depends on parent dir name
            parent = parts[-1] if parts else "."
            vals = vals.get(parent, [])
        for v in vals:
            d = fmt_dirname(spec, v)
            recurse(depth + 1, parts + [d], {**values, spec["name"]: v})

    recurse(0, [], {})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=None, help="simulation root (default: search upward from cwd)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    root = find_root(a.root)
    plan = load_plan(root)
    rows = expand(plan)
    if not rows:
        raise SystemExit("[simflow] plan expands to zero run dirs — check tree.levels")

    # print a compact tree: variable dirs in full, run dirs collapsed to one line
    seen = set()
    runs_under = {}
    for row in rows:
        runs_under.setdefault(Path(row["relpath"]).parts[:-1], []).append(row["run"])
    print(f"{root.name}/   (level 0)")
    print("    reffiles/")
    for parent, runs in runs_under.items():
        for depth in range(1, len(parent) + 1):
            key = parent[:depth]
            if key not in seen:
                seen.add(key)
                print("    " * depth + key[-1] + f"/   (level {depth})")
        shown = runs if len(runs) <= 4 else runs[:2] + ["..."] + runs[-1:]
        print("    " * (len(parent) + 1) + " ".join(f"{r}/" if r != "..." else r for r in shown)
              + f"   ({len(runs)} runs, level {len(parent) + 1})")
    n_leaf_parents = len({str(Path(r['relpath']).parent) for r in rows})
    log(f"{len(rows)} run dirs under {n_leaf_parents} variable-combination dirs; run level = {len(levels(plan)) + 1}")

    if a.dry_run:
        return
    (root / plan.get("reffiles", {}).get("dir", "reffiles")).mkdir(exist_ok=True)
    for row in rows:
        (root / row["relpath"]).mkdir(parents=True, exist_ok=True)
    write_manifest(root, rows)
    log(f"created dirs and wrote {root / 'manifest.csv'}")


if __name__ == "__main__":
    main()
