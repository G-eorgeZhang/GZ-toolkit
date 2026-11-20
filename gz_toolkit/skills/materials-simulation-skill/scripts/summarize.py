#!/usr/bin/env python3
"""
summarize.py — extract values from every run dir into one CSV (step 4, "summarize").

    python workflow/summarize.py                     # -> ROOT/<summary.output> (default summary.csv)
    python workflow/summarize.py --only-finished     # skip runs whose checks fail
    python workflow/summarize.py --workers 16        # parallel (useful with OVITO extractors)
    python workflow/summarize.py --average           # also write summary_avg.csv (mean/std over runs)

Columns: relpath, every tree variable, run, then one column per `summary.fields` entry.
Field types (paths relative to run dir; `which`: last | first | min | max | mean | all | <int index>):
  regex          {name: Eb, type: regex, file: log.lammps, pattern: "Eb = (\\S+)", which: last}
  lammps_thermo  {name: PE, type: lammps_thermo, file: log.lammps, column: PotEng, which: last, run: -1}
  csv_column     {name: Emin, type: csv_column, file: Seakmc_summary.csv, column: barrier, which: min}
  vasp_energy    {name: E0, type: vasp_energy, file: OSZICAR}                # last E0= value
  file_lines     {name: nlines, type: file_lines, file: dump.final}
  count_files    {name: nsteps, type: count_files, file: "DataOut/KMC_*.dat"}
  ovito          {name: nvac, type: ovito, file: final.dump, script: "reffiles/ovito_ws.py:count_vacancies"}
  python         {name: x, type: python, script: "reffiles/my_extract.py:func"}   # func(run_dir: Path, field) -> value
Missing files / failed extractions give an empty cell and a warning — the run is not dropped.
"""
import argparse
import csv
import glob
import importlib.util
import math
import re
import statistics
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from checklib import status
from simflow_common import find_root, levels, load_plan, log, read_manifest

NUM = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def pick(vals, which):
    vals = [v for v in vals if v is not None]
    if not vals:
        return None
    if which in (None, "last"):
        return vals[-1]
    if which == "first":
        return vals[0]
    if which == "all":
        return ";".join(map(str, vals))
    nums = [float(v) for v in vals]
    if which == "min":
        return min(nums)
    if which == "max":
        return max(nums)
    if which == "mean":
        return sum(nums) / len(nums)
    return vals[int(which)]


def one_file(run, pattern):
    hits = sorted(glob.glob(str(run / pattern)))
    if not hits:
        raise FileNotFoundError(pattern)
    return Path(hits[-1])


def lammps_thermo(path, column, run_index=-1):
    """Return the values of `column` from thermo block number `run_index` of a LAMMPS log."""
    blocks, cur, hdr = [], None, None
    with open(path, errors="replace") as fh:
        for line in fh:
            s = line.split()
            if s and s[0] == "Step" and len(s) > 1:
                hdr, cur = s, []
                blocks.append((hdr, cur))
                continue
            if cur is not None:
                if line.startswith("Loop time") or not s:
                    cur = None
                    continue
                if len(s) == len(hdr):
                    try:
                        cur.append([float(x) for x in s])
                    except ValueError:
                        cur = None
    blocks = [b for b in blocks if b[1]]
    if not blocks:
        raise ValueError("no thermo data")
    hdr, data = blocks[int(run_index)]
    if column not in hdr:
        raise KeyError(f"column {column} not in thermo header {hdr}")
    i = hdr.index(column)
    return [row[i] for row in data]


def load_func(root, spec):
    f, fn = spec.rsplit(":", 1)
    path = (root / f).resolve()
    sp = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(mod)
    return getattr(mod, fn)


def extract(root, run, fld):
    t = fld["type"]
    which = fld.get("which")
    if t == "regex":
        rx = re.compile(fld["pattern"])
        vals = []
        with open(one_file(run, fld["file"]), errors="replace") as fh:
            for line in fh:
                m = rx.search(line)
                if m:
                    vals.append(m.group(1) if m.groups() else m.group(0))
        return pick(vals, which)
    if t == "lammps_thermo":
        return pick(lammps_thermo(one_file(run, fld["file"]), fld["column"], fld.get("run", -1)), which)
    if t == "csv_column":
        with open(one_file(run, fld["file"]), newline="") as fh:
            rd = csv.DictReader(fh, delimiter=fld.get("delimiter", ","))
            vals = [row[fld["column"]].strip() for row in rd if row.get(fld["column"], "").strip()]
        return pick(vals, which)
    if t == "vasp_energy":
        f = one_file(run, fld.get("file", "OSZICAR"))
        key = fld.get("key", "E0=")
        vals = []
        with open(f, errors="replace") as fh:
            for line in fh:
                if key in line:
                    m = NUM.search(line.split(key, 1)[1])
                    if m:
                        vals.append(float(m.group(0)))
        return pick(vals, which)
    if t == "file_lines":
        with open(one_file(run, fld["file"]), "rb") as fh:
            return sum(1 for _ in fh)
    if t == "count_files":
        return len(glob.glob(str(run / fld["file"])))
    if t in ("python", "ovito"):
        func = load_func(root, fld["script"])
        return func(run, fld)
    raise ValueError(f"unknown field type {t}")


def do_row(args):
    root, plan, row, only_finished = args
    root = Path(root)
    run = root / row["relpath"]
    out = {"relpath": row["relpath"]}
    out.update({s["name"]: row[s["name"]] for s in levels(plan)})
    out["run"] = row["run"]
    if only_finished:
        st, _ = status(run, plan.get("checks") or [])
        if st != "FINISHED":
            return None, []
    warns = []
    for fld in (plan.get("summary") or {}).get("fields", []):
        try:
            out[fld["name"]] = extract(root, run, fld)
        except Exception as e:  # keep the row, flag the cell
            out[fld["name"]] = ""
            warns.append(f"{row['relpath']}: {fld['name']}: {type(e).__name__}: {e}")
    return out, warns


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root")
    ap.add_argument("--only-finished", action="store_true")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--average", action="store_true", help="also write mean/std over runs per variable combination")
    ap.add_argument("-o", "--output", default=None)
    a = ap.parse_args()

    root = find_root(a.root)
    plan = load_plan(root)
    fields = (plan.get("summary") or {}).get("fields", [])
    if not fields:
        raise SystemExit("[simflow] sim_plan has no summary.fields — ask the user which values to extract")
    rows = read_manifest(root)
    jobs = [(str(root), plan, r, a.only_finished) for r in rows]
    if a.workers > 1:
        with ProcessPoolExecutor(a.workers) as ex:
            res = list(ex.map(do_row, jobs))
    else:
        res = [do_row(j) for j in jobs]

    out_rows = [r for r, _ in res if r]
    warns = [w for _, ws in res for w in ws]
    outfile = root / (a.output or (plan.get("summary") or {}).get("output", "summary.csv"))
    if not out_rows:
        raise SystemExit("[simflow] nothing to write (no finished runs?)")
    with open(outfile, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    # group warnings: "field: ErrorType: msg" -> list of runs
    grouped = {}
    for wmsg in warns:
        rel, rest = wmsg.split(": ", 1)
        grouped.setdefault(rest, []).append(rel)
    for rest, rels in grouped.items():
        tail = f" ({rels[0]}{', ...' if len(rels) > 1 else ''})"
        print(f"  WARN {len(rels):4d} run(s): {rest}{tail}")
    log(f"wrote {len(out_rows)} rows x {len(fields)} fields -> {outfile.name}")

    if a.average:
        varnames = [s["name"] for s in levels(plan)]
        groups = {}
        for r in out_rows:
            groups.setdefault(tuple(r[v] for v in varnames), []).append(r)
        avg_rows = []
        for key, grp in groups.items():
            ar = dict(zip(varnames, key))
            ar["n_runs"] = len(grp)
            for fld in fields:
                vals = []
                for r in grp:
                    try:
                        v = float(r[fld["name"]])
                        if not math.isnan(v):
                            vals.append(v)
                    except (TypeError, ValueError):
                        pass
                ar[fld["name"] + "_mean"] = statistics.mean(vals) if vals else ""
                ar[fld["name"] + "_std"] = statistics.stdev(vals) if len(vals) > 1 else ""
                ar[fld["name"] + "_n"] = len(vals)
            avg_rows.append(ar)
        af = outfile.with_name(outfile.stem + "_avg.csv")
        with open(af, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(avg_rows[0].keys()))
            w.writeheader()
            w.writerows(avg_rows)
        log(f"wrote {len(avg_rows)} averaged rows -> {af.name}")


if __name__ == "__main__":
    main()
