#!/usr/bin/env python3
"""
submit.py — sbatch the job files written by make_jobs.py.

    python workflow/submit.py --level 3 --dry-run    # show what would be submitted
    python workflow/submit.py --level 3              # submit
    python workflow/submit.py --level 4 --max 50     # submit at most 50 (queue limits)

Each job is submitted from its own directory (so SLURM_SUBMIT_DIR = that dir).
Job IDs are appended to jobs_submitted.csv in the root, and every run dir covered by a
submitted job gets a .simflow_submitted marker (setup_runs.py will not overwrite it).
"""
import argparse
import csv
import datetime as dt
import re
import shlex
import shutil
import subprocess
import time
from pathlib import Path

from simflow_common import find_root, load_plan, log, read_manifest, run_level


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root")
    ap.add_argument("--level", type=int, required=True)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--max", type=int, default=None)
    ap.add_argument("--delay", type=float, default=0.2, help="seconds between sbatch calls")
    ap.add_argument("--sbatch-args", default="", help="extra args, e.g. '--dependency=afterok:123'")
    a = ap.parse_args()

    root = find_root(a.root)
    plan = load_plan(root)
    lst = root / f"jobs_lv{a.level}.txt"
    if not lst.is_file():
        raise SystemExit(f"[simflow] {lst.name} missing — run make_jobs.py --level {a.level} first")
    jobs = [l.strip() for l in lst.read_text().splitlines() if l.strip()]
    if a.max:
        jobs = jobs[: a.max]
    if not a.dry_run and not shutil.which("sbatch"):
        raise SystemExit("[simflow] sbatch not found on this machine — are you on the cluster login node?")

    rows = read_manifest(root)
    rl = run_level(plan)
    if a.dry_run:
        for j in jobs:
            jf = root / j
            cmd = ["sbatch", "--parsable", *shlex.split(a.sbatch_args), jf.name]
            print(f"  (cd {jf.parent.relative_to(root)} && {shlex.join(cmd)})")
        log(f"would submit {len(jobs)} job(s)")
        return
    logf = root / "jobs_submitted.csv"
    new = not logf.exists()
    with open(logf, "a", newline="") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(["time", "level", "job_file", "job_id"])
        for j in jobs:
            jf = root / j
            cmd = ["sbatch", "--parsable", *shlex.split(a.sbatch_args), jf.name]
            res = subprocess.run(cmd, cwd=jf.parent, capture_output=True, text=True)
            m = re.fullmatch(r"(\d+)(?:;[^\s;]+)?", res.stdout.strip())
            jid = m.group(1) if (res.returncode == 0 and m) else "FAILED:" + (res.stderr.strip() or res.stdout.strip())
            w.writerow([dt.datetime.now().isoformat(timespec="seconds"), a.level, j, jid])
            fh.flush()
            print(f"  {j} -> {jid}")
            if jid.startswith("FAILED"):
                raise SystemExit(f"[simflow] submission stopped; failure recorded in {logf}")
            if not jid.startswith("FAILED"):
                jobdir = jf.parent.relative_to(root)
                for r in rows:
                    rp = Path(r["relpath"])
                    if a.level == rl and rp == jobdir or a.level < rl and (str(jobdir) == "." or jobdir in rp.parents):
                        (root / rp / ".simflow_submitted").write_text(jid + "\n")
            time.sleep(a.delay)
    log(f"submitted {len(jobs)} job(s)")


if __name__ == "__main__":
    main()
