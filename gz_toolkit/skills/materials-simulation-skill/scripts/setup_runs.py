#!/usr/bin/env python3
"""
setup_runs.py — populate every run dir from reffiles/ (step 3, "set up simulations").

    python workflow/setup_runs.py                 # all run dirs, parallel over local cores
    python workflow/setup_runs.py --workers 1     # serial
    python workflow/setup_runs.py --only T300K    # only run dirs whose relpath starts with T300K
    python workflow/setup_runs.py --dry-run       # report what would happen

What it does in each run dir (all driven by sim_plan `reffiles:`):
  copy:       files copied verbatim (globs allowed), e.g. potential files
  link:       files symlinked instead of copied (large, never-changing files)
  templates:  files copied with {{KEY}} placeholders filled from manifest.csv
              (each level variable, RUN, RUN_INDEX, SEED, RELPATH, JOBNAME, plan constants)
              an entry may be "in.lammps" or {src: in.lammps.tmpl, dst: in.lammps}
  require:    files that must already be in the run dir (the structure from step 1)

It never deletes anything. A run dir that has been submitted (.simflow_submitted marker)
is skipped unless --force, so finished or running simulations are not overwritten.
Run this on a compute node for very large trees: the work is split over --workers processes.
"""
import argparse
import glob
import os
import shutil
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from simflow_common import context_for, find_root, load_plan, log, read_manifest, render


def plan_entries(plan):
    rf = plan.get("reffiles") or {}
    refdir = rf.get("dir", "reffiles")
    tmpl = []
    for t in rf.get("templates") or []:
        tmpl.append((t, t) if isinstance(t, str) else (t["src"], t.get("dst", t["src"])))
    return refdir, rf.get("copy") or [], rf.get("link") or [], tmpl, rf.get("require") or []


def setup_one(args):
    root, plan, row, force, dry = args
    root = Path(root)
    run = root / row["relpath"]
    refdir, copies, links, templates, require = plan_entries(plan)
    ref = root / refdir
    msgs = []
    if (run / ".simflow_submitted").exists() and not force:
        return row["relpath"], "SKIP (already submitted; use --force)", []
    if not dry:
        run.mkdir(parents=True, exist_ok=True)
    ctx = context_for(row, plan)

    for pattern in copies:
        hits = sorted(glob.glob(str(ref / pattern)))
        if not hits:
            msgs.append(f"missing reffile {pattern}")
        for h in hits:
            if not dry:
                shutil.copy2(h, run / Path(h).name)
    for pattern in links:
        hits = sorted(glob.glob(str(ref / pattern)))
        if not hits:
            msgs.append(f"missing reffile {pattern}")
        for h in hits:
            dst = run / Path(h).name
            if dry:
                continue
            if dst.is_symlink() or dst.exists():
                dst.unlink()
            try:
                dst.symlink_to(os.path.relpath(h, run))
            except OSError:  # Windows without symlink privilege -> copy
                shutil.copy2(h, dst)
    for src, dst in templates:
        f = ref / src
        if not f.is_file():
            msgs.append(f"missing template {src}")
            continue
        try:
            text = render(f.read_text(encoding="utf-8"), ctx, where=f"{src} for {row['relpath']}")
        except KeyError as e:
            msgs.append(str(e))
            continue
        if not dry:
            (run / dst).write_text(text, encoding="utf-8")
    for req in require:
        if not dry and not glob.glob(str(run / req)):
            msgs.append(f"required file not present: {req}")
    return row["relpath"], ("OK" if not msgs else "PROBLEM"), msgs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root")
    ap.add_argument("--workers", type=int, default=min(32, os.cpu_count() or 1))
    ap.add_argument("--only", default=None, help="prefix filter on relpath")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    root = find_root(a.root)
    plan = load_plan(root)
    rows = read_manifest(root)
    if a.only:
        rows = [r for r in rows if r["relpath"].replace("\\", "/").startswith(a.only.rstrip("/"))]
    jobs = [(str(root), plan, r, a.force, a.dry_run) for r in rows]
    log(f"setting up {len(jobs)} run dirs with {a.workers} worker(s){' [dry-run]' if a.dry_run else ''}")

    if a.workers > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            results = list(ex.map(setup_one, jobs, chunksize=max(1, len(jobs) // (a.workers * 4))))
    else:
        results = [setup_one(j) for j in jobs]

    bad = [r for r in results if r[1] != "OK"]
    for rel, status, msgs in bad[:20]:
        print(f"  {status:8s} {rel}: {'; '.join(msgs)}")
    if len(bad) > 20:
        print(f"  ... and {len(bad) - 20} more")
    log(f"{len(results) - len(bad)} OK, {len(bad)} with problems/skipped")
    raise SystemExit(1 if any(s == "PROBLEM" for _, s, _ in bad) else 0)


if __name__ == "__main__":
    main()
