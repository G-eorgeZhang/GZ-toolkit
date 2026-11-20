#!/usr/bin/env python3
"""
make_jobs.py — write SLURM job files at a chosen PARALLELIZATION LEVEL.

    python workflow/make_jobs.py --level 4            # scatter: one job per run dir
    python workflow/make_jobs.py --level 3            # one job per level-3 dir, its runs in serial
    python workflow/make_jobs.py --level 3 --pack 4   # ... 4 runs at a time inside the job
    python workflow/make_jobs.py --level 0            # one job for the whole tree
    python workflow/make_jobs.py --level 3 --only-unfinished   # resubmission: skip FINISHED runs

Levels: root = 0, variable 1 = 1, ..., run dirs = (number of variables + 1) = scatter.
  * level == run level  -> the reference job (reffiles/<job>) is rendered into every run dir.
  * level <  run level  -> one job per directory at that level. Its header (#SBATCH, module,
    source, export ...) comes from the reference job; its body is the reference job's run
    section repeated for every run dir below it:  cd run1 -> srun ... -> cd back -> cd run2 ...

Splitting header/body of the reference job: put a line `# @RUN` in it — everything after is
the per-run body. Without the marker, the header is the leading block of comments/#SBATCH/
module/source/export/conda/cd $SLURM_SUBMIT_DIR lines and the body is everything after.

Walltime: set job.time_per_run in sim_plan (e.g. "02:00:00") and the aggregated job's
--time becomes time_per_run x ceil(n_runs / pack) (+ job.time_margin, default 10 %).
Writes jobs_lv<L>.txt in the root (list of job files) for submit.py.
"""
import argparse
import math
import re
from pathlib import Path

from checklib import status
from simflow_common import (ancestor, context_for, find_root, hms_to_s, job_name, levels,
                            load_plan, log, read_manifest, render, run_level, s_to_slurm)

HEADER_CMDS = {"module", "ml", "source", "export", "conda", "set", "ulimit", "umask",
               "unset", "alias", "spack", ".", "mamba", "micromamba"}
SBATCH_NAME = re.compile(r"^(#SBATCH\s+)(--job-name[= ]\S+|-J\s*\S+)(.*)$")
SBATCH_TIME = re.compile(r"^(#SBATCH\s+)(--time[= ]\S+|-t\s*\S+)(.*)$")


def split_job(text: str):
    lines = text.splitlines()
    for i, l in enumerate(lines):
        if re.match(r"^\s*#\s*@RUN\b", l):
            return lines[:i], lines[i + 1:]
    i = 0
    for i, l in enumerate(lines):
        s = l.strip()
        if not s or s.startswith("#"):
            continue
        tok = s.split()[0]
        if tok in HEADER_CMDS or re.match(r"cd\s+\"?\$\{?SLURM_SUBMIT_DIR", s):
            continue
        break
    else:
        i = len(lines)
    return lines[:i], lines[i:]


def set_directive(header, rx, flag, value):
    out, done = [], False
    for l in header:
        if rx.match(l):
            out.append(f"#SBATCH {flag}={value}")
            done = True
        else:
            out.append(l)
    if not done:  # insert after the last #SBATCH line (or after shebang)
        idx = max([i for i, l in enumerate(out) if l.startswith("#SBATCH")] or [0])
        out.insert(idx + 1, f"#SBATCH {flag}={value}")
    return out


def fail_fast(header):
    # Keep all scheduler directives before the first executable command.
    index = max([i for i, line in enumerate(header) if line.startswith("#SBATCH")] or [0])
    return header[:index + 1] + ["set -e", "set -o pipefail"] + header[index + 1:]


def wait_for_pack():
    # A bare `wait` reports success even if one of several background runs fails.
    return ["simflow_failed=0", 'for simflow_pid in "${simflow_pids[@]}"; do',
            '  wait "$simflow_pid" || simflow_failed=1', "done",
            '[ "$simflow_failed" -eq 0 ] || exit 1', "simflow_pids=()"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root")
    ap.add_argument("--level", type=int, default=None)
    ap.add_argument("--pack", type=int, default=None, help="runs executed concurrently inside one job")
    ap.add_argument("--only-unfinished", action="store_true")
    ap.add_argument("--only", default=None, help="prefix filter on relpath")
    a = ap.parse_args()

    root = find_root(a.root)
    plan = load_plan(root)
    jobcfg = plan.get("job") or {}
    rl = run_level(plan)
    level = a.level if a.level is not None else int(jobcfg.get("level", rl))
    pack = a.pack if a.pack is not None else int(jobcfg.get("pack", 1))
    if pack < 1:
        raise SystemExit("[simflow] pack must be positive")
    if not 0 <= level <= rl:
        raise SystemExit(f"[simflow] level must be 0..{rl} (run level = {rl}, scatter)")

    refdir = root / (plan.get("reffiles") or {}).get("dir", "reffiles")
    refjob_name = (plan.get("reffiles") or {}).get("job")
    if not refjob_name or not (refdir / refjob_name).is_file():
        raise SystemExit(f"[simflow] reference job file reffiles/{refjob_name} not found — ask the user for it")
    ref_text = (refdir / refjob_name).read_text(encoding="utf-8")

    rows = read_manifest(root)
    if a.only:
        rows = [r for r in rows if r["relpath"].replace("\\", "/").startswith(a.only.rstrip("/"))]
    if a.only_unfinished:
        checks = plan.get("checks") or []
        rows = [r for r in rows if status(root / r["relpath"], checks)[0] != "FINISHED"]
        log(f"--only-unfinished: {len(rows)} runs left")

    written = []
    if level == rl:  # ---------------------------------------------------- scatter
        for r in rows:
            ctx = context_for(r, plan)
            text = render(ref_text, ctx, where=f"{refjob_name} for {r['relpath']}")
            hdr, body = split_job(text)
            hdr = set_directive(hdr, SBATCH_NAME, "--job-name", ctx["JOBNAME"])
            if jobcfg.get("time_per_run"):
                hdr = set_directive(hdr, SBATCH_TIME, "--time", jobcfg["time_per_run"])
            out = root / r["relpath"] / refjob_name
            with out.open("w", encoding="utf-8", newline="\n") as fh:
                fh.write("\n".join(fail_fast(hdr) + body) + "\n")
            written.append(out)
    else:  # ---------------------------------------------------------- aggregated
        groups = {}
        for r in rows:
            groups.setdefault(ancestor(r, level), []).append(r)
        stem = Path(refjob_name).stem
        varnames = [s["name"] for s in levels(plan)][:level]
        for gdir, grows in groups.items():
            gctx = {k: grows[0][k] for k in varnames}
            gctx.update({"RELPATH": gdir, "JOBNAME": job_name(plan, "" if gdir == "." else gdir),
                         "RUN": "multi", "RUN_INDEX": "0", "SEED": "0"})
            gctx.update({k: str(v) for k, v in (plan.get("constants") or {}).items()})
            hdr_raw, body_raw = split_job(ref_text)
            hdr = render("\n".join(hdr_raw), gctx, where=f"{refjob_name} header for {gdir}").splitlines()
            hdr = set_directive(hdr, SBATCH_NAME, "--job-name", gctx["JOBNAME"])
            if jobcfg.get("time_per_run"):
                per = hms_to_s(str(jobcfg["time_per_run"]))
                tot = per * math.ceil(len(grows) / pack) * (1 + float(jobcfg.get("time_margin", 0.1)))
                hdr = set_directive(hdr, SBATCH_TIME, "--time", s_to_slurm(int(tot)))
            lines = fail_fast(hdr) + ["", f"# ===== simflow: level-{level} job, {len(grows)} run(s), pack={pack} =====",
                           'JOBDIR="${SLURM_SUBMIT_DIR:-$(pwd)}"', ""]
            if pack > 1:
                lines.append("simflow_pids=()")
            for i, r in enumerate(grows, 1):
                sub = Path(r["relpath"]).relative_to(gdir).as_posix() if gdir != "." else r["relpath"]
                sub = sub.replace("\\", "/")
                body = render("\n".join(body_raw), context_for(r, plan), where=f"body for {r['relpath']}")
                lines += [f"# ---- {i}/{len(grows)}: {r['relpath']}",
                          "(",
                          f'  cd "$JOBDIR/{sub}" || {{ echo "[simflow] cannot cd {sub}"; exit 1; }}',
                          f'  echo "[simflow] start {sub} $(date)"']
                lines += ["  " + l if l.strip() else l for l in body.splitlines()]
                lines += [f'  echo "[simflow] end   {sub} $(date) exit=$?"',
                          ")" + (" &" if pack > 1 else "")]
                if pack > 1:
                    lines.append('simflow_pids+=("$!")')
                if pack > 1 and i % pack == 0:
                    lines.extend(wait_for_pack())
                lines.append("")
            if pack > 1:
                lines.extend(wait_for_pack())
            out = (root / gdir) / f"{stem}_lv{level}.job"
            with out.open("w", encoding="utf-8", newline="\n") as fh:
                fh.write("\n".join(lines) + "\n")
            written.append(out)

    listfile = root / f"jobs_lv{level}.txt"
    listfile.write_text("\n".join(str(p.relative_to(root)).replace("\\", "/") for p in written) + "\n")
    log(f"level {level} (run level {rl}): wrote {len(written)} job file(s); list in {listfile.name}")
    if level < rl and pack > 1:
        log("pack>1: make sure each srun in the body asks for its share of the node "
            "(e.g. srun -n $((SLURM_NTASKS/PACK)) --exact) or runs will oversubscribe cores")


if __name__ == "__main__":
    try:
        main()
    except KeyError as e:
        raise SystemExit(f"[simflow] {e.args[0]}\n  -> a placeholder in the reference job has no value at this level "
                         "(level-L headers only know variables of levels 1..L, JOBNAME, RELPATH and constants). "
                         "Remove it from the header or move it below the '# @RUN' line.")
