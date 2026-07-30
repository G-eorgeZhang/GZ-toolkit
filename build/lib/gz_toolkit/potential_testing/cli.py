"""Command-line entry points for the potential-testing pipeline.

Usage::

    # Initialize a fresh project (creates pot_inputs/, main.py, summarize.py).
    python -m gz_toolkit.potential_testing.cli init --root . --pot-name eam_fs

    # Plan + emit reference dirs and submit_all.sh (top-level driver, called by main.py).
    python -m gz_toolkit.potential_testing.cli run --pot-inputs pot_inputs --run-dir .

    # Build defect / alloy / gas dirs from already-relaxed reference data.
    # Called from inside a SLURM job after reference jobs complete.
    python -m gz_toolkit.potential_testing.cli build-defects \\
        --pot-config pot_inputs/eam_fs.json --run-dir .

    # Used in paral_degree=4 only — write per-case .job files and sbatch them.
    python -m gz_toolkit.potential_testing.cli scatter-cases \\
        --pot-config pot_inputs/eam_fs.json --run-dir .

    # Aggregate results into a wide CSV (called by summarize.py).
    python -m gz_toolkit.potential_testing.cli summarize --pot-inputs pot_inputs --run-dir .
"""

from __future__ import annotations

import argparse
from pathlib import Path

from gz_toolkit.potential_testing.config import (
    load_potential_config,
    load_potential_configs,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="gz_toolkit.potential_testing")
    sub = parser.add_subparsers(dest="command", required=True)

    p_init = sub.add_parser("init", help="Initialize a project skeleton")
    p_init.add_argument("--root", default=".")
    p_init.add_argument("--pot-name", default="tao_2.31")
    p_init.add_argument("--interactive", action="store_true",
                        help="Answer plain-language questions instead of editing JSON by hand")

    p_val = sub.add_parser("validate", help="Check every config for mistakes before submitting")
    p_val.add_argument("--pot-inputs", default="pot_inputs")

    p_st = sub.add_parser("status", help="Show done/failed/running/pending counts per work group")
    p_st.add_argument("--pot-inputs", default="pot_inputs")
    p_st.add_argument("--run-dir", default=".")
    p_st.add_argument("-v", "--verbose", action="store_true",
                      help="Also list the failed and running case directories")

    p_run = sub.add_parser("run", help="Prepare reference dirs and emit submit_all.sh")
    p_run.add_argument("--pot-inputs", default="pot_inputs")
    p_run.add_argument("--run-dir", default=".")

    p_bd = sub.add_parser("build-defects", help="Build defect/alloy/gas dirs from relaxed reference")
    p_bd.add_argument("--pot-config", required=True)
    p_bd.add_argument("--run-dir", default=".")

    p_sc = sub.add_parser("scatter-cases", help="Write per-case .job files and submit them (paral=4)")
    p_sc.add_argument("--pot-config", required=True)
    p_sc.add_argument("--run-dir", default=".")
    p_sc.add_argument("--pot-inputs", default="pot_inputs")

    p_sum = sub.add_parser("summarize", help="Write wide-format summary.csv")
    p_sum.add_argument("--pot-inputs", default="pot_inputs")
    p_sum.add_argument("--run-dir", default=".")

    # Backwards-compat single-stage hooks (older tests)
    p_ref = sub.add_parser("reference", help="[compat] Single-pot reference stage")
    p_ref.add_argument("--config", required=True)
    p_ref.add_argument("--run-dir", default=".")
    p_def = sub.add_parser("defects", help="[compat] Single-pot defect stage")
    p_def.add_argument("--config", required=True)
    p_def.add_argument("--run-dir", default=".")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.command == "init":
        if args.interactive:
            from gz_toolkit.potential_testing.wizard import run_init_wizard
            run_init_wizard(root=args.root)
            return 0
        from gz_toolkit.potential_testing.project import init_potential_project
        root = init_potential_project(root_dir=args.root, pot_name=args.pot_name)
        print(f"Initialized potential testing project at: {root.resolve()}")
        return 0

    if args.command == "validate":
        from gz_toolkit.potential_testing.validate import validate_project
        report = validate_project(args.pot_inputs)
        n_problems = 0
        for fname, problems in report.items():
            if problems:
                n_problems += len(problems)
                print(f"{fname}:")
                for p in problems:
                    print(f"  - {p}")
            else:
                print(f"{fname}: OK")
        if n_problems:
            print(f"\n{n_problems} problem(s) found — fix them before submitting.")
            return 1
        print("\nAll configs valid.")
        return 0

    if args.command == "status":
        from gz_toolkit.potential_testing.status import collect_status, format_status_table
        configs = load_potential_configs(args.pot_inputs)
        status_by_pot = {
            cfg.potential.pot_name: collect_status(cfg, run_dir=args.run_dir)
            for cfg in configs
        }
        print(format_status_table(status_by_pot, verbose=args.verbose))
        return 0

    if args.command == "run":
        from gz_toolkit.potential_testing.workflow import run_all
        out = run_all(pot_inputs_dir=args.pot_inputs, run_dir=args.run_dir)
        print(f"Wrote launcher: {out}")
        return 0

    if args.command == "build-defects":
        from gz_toolkit.potential_testing.pipeline import build_defects_post_reference
        cfg = load_potential_config(args.pot_config)
        manifest = build_defects_post_reference(cfg, run_dir=args.run_dir)
        print(f"Built {len(manifest['point_defects'])} point defects, "
              f"{len(manifest.get('elastic', []))} elastic cases, "
              f"{len(manifest['loops'])} loops, {len(manifest['alloys'])} alloy cases, "
              f"{len(manifest['gas_complexes'])} gas-complex cases.")
        return 0

    if args.command == "scatter-cases":
        from gz_toolkit.potential_testing.parallel import scatter_cases
        cfg = load_potential_config(args.pot_config)
        written = scatter_cases(cfg, Path(args.run_dir).resolve(), Path(args.pot_inputs).resolve())
        print(f"Scatter wrote {len(written)} per-case .job files (and attempted sbatch).")
        return 0

    if args.command == "summarize":
        from gz_toolkit.potential_testing.summary import write_wide_summary
        out = write_wide_summary(pot_inputs_dir=args.pot_inputs, run_dir=args.run_dir)
        print(f"Summary written to: {out}")
        return 0

    if args.command == "reference":
        from gz_toolkit.potential_testing.pipeline import prepare_reference_stage
        cfg = load_potential_config(args.config)
        prepare_reference_stage(cfg, run_dir=args.run_dir)
        return 0

    if args.command == "defects":
        from gz_toolkit.potential_testing.pipeline import prepare_defect_stage
        cfg = load_potential_config(args.config)
        prepare_defect_stage(cfg, run_dir=args.run_dir)
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
