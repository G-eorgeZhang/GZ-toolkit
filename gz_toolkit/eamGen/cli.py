"""Agent-friendly tools: JSON inputs, JSON reports, explicit external execution."""
from __future__ import annotations
import argparse
from dataclasses import asdict
import json
from pathlib import Path

from .composition import compose_model
from .data import Configuration, Dataset, compare_metrics
from .fitting import evaluate_dataset
from .model import EAMModel, Element
from .project import init_project, fit_revision, record_review, export_testing_project
from .references import prepare_lammps, prepare_vasp, run_reference, import_reference
from .sampling import crystal_snapshots


def _load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _write(path, value):
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"Output exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False), encoding="utf-8")


def build_parser():
    parser = argparse.ArgumentParser(prog="gz-toolkit-eamgen")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("init"); p.add_argument("--root", required=True); p.add_argument("--elements", required=True, help="JSON list of Element metadata")
    p.add_argument("--style", choices=["eam/alloy", "eam/fs", "eam/he"], default="eam/he"); p.add_argument("--baseline", help="Model JSON")
    p = sub.add_parser("import-potential"); p.add_argument("--file", required=True); p.add_argument("--style", required=True, choices=["eam/alloy", "eam/fs", "eam/he"]); p.add_argument("--output", required=True)
    p = sub.add_parser("compose"); p.add_argument("--primary", required=True); p.add_argument("--secondary", required=True); p.add_argument("--elements", nargs="+", required=True); p.add_argument("--output", required=True)
    p = sub.add_parser("fit"); p.add_argument("--root", required=True); p.add_argument("--revision", required=True); p.add_argument("--baseline")
    p = sub.add_parser("evaluate"); p.add_argument("--model", required=True); p.add_argument("--dataset", required=True); p.add_argument("--include-test", action="store_true"); p.add_argument("--output", required=True)
    p = sub.add_parser("review"); p.add_argument("--revision", required=True); p.add_argument("--status", required=True, choices=["approved_for_testing", "rejected"]); p.add_argument("--reviewer", required=True); p.add_argument("--reason", required=True)
    p = sub.add_parser("handoff"); p.add_argument("--revision", required=True); p.add_argument("--root", required=True); p.add_argument("--pot-name", default="eam_candidate")
    p = sub.add_parser("prepare-reference"); p.add_argument("--spec", required=True, help="JSON: backend, configuration, and backend-specific options"); p.add_argument("--directory", required=True)
    p = sub.add_parser("run-reference"); p.add_argument("--directory", required=True); p.add_argument("--argv-file", required=True, help="JSON list: synchronous command argv"); p.add_argument("--timeout", type=float, default=3600)
    p = sub.add_parser("import-reference"); p.add_argument("--directory", required=True); p.add_argument("--output", required=True)
    p = sub.add_parser("compare"); p.add_argument("--metrics", required=True, help="JSON mapping metric -> number in target units"); p.add_argument("--targets", required=True); p.add_argument("--output", required=True)
    p = sub.add_parser("sample"); p.add_argument("--spec", required=True, help="JSON list of crystal_snapshots keyword objects"); p.add_argument("--output", required=True)
    p = sub.add_parser("collect-references"); p.add_argument("--directories", nargs="+", required=True); p.add_argument("--energy-reference", required=True); p.add_argument("--offsets", help="JSON species -> energy offset in eV"); p.add_argument("--output", required=True)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command == "init":
        baseline = EAMModel.load(args.baseline) if args.baseline else None
        result = init_project(args.root, [Element(**e) for e in _load(args.elements)], args.style, baseline)
    elif args.command == "import-potential":
        if Path(args.output).exists():
            raise FileExistsError(args.output)
        result = EAMModel.read_setfl(args.file, args.style).save(args.output)
    elif args.command == "compose":
        if Path(args.output).exists():
            raise FileExistsError(args.output)
        result = compose_model(EAMModel.load(args.primary), EAMModel.load(args.secondary), args.elements).save(args.output)
    elif args.command == "fit":
        result = fit_revision(args.root, args.revision, args.baseline)
    elif args.command == "evaluate":
        _write(args.output, evaluate_dataset(EAMModel.load(args.model), Dataset.load(args.dataset), args.include_test)); result = args.output
    elif args.command == "review":
        result = record_review(args.revision, args.status, args.reviewer, args.reason)
    elif args.command == "handoff":
        revision = Path(args.revision)
        if _load(revision / "review.json")["status"] != "approved_for_testing":
            raise ValueError("Record a human review approved_for_testing before handoff.")
        result = export_testing_project(EAMModel.load(revision / "model.json"), args.root, args.pot_name)
    elif args.command == "prepare-reference":
        spec = _load(args.spec); backend = spec.pop("backend"); c = Configuration(**spec.pop("configuration"))
        if backend == "lammps":
            result = prepare_lammps(c, args.directory, **spec)
        elif backend == "vasp":
            result = prepare_vasp(c, args.directory, **spec)
        else:
            raise ValueError("Backend must be lammps or vasp.")
    elif args.command == "run-reference":
        result = run_reference(args.directory, _load(args.argv_file), args.timeout)
    elif args.command == "import-reference":
        _write(args.output, asdict(import_reference(args.directory))); result = args.output
    elif args.command == "sample":
        configurations = [c for spec in _load(args.spec) for c in crystal_snapshots(**spec)]
        dataset = Dataset(configurations, "REPLACE: reference calculations and energy convention required")
        dataset.validate(require_targets=False)
        _write(args.output, asdict(dataset)); result = args.output
    elif args.command == "collect-references":
        dataset = Dataset([import_reference(d) for d in args.directories], args.energy_reference,
                          _load(args.offsets) if args.offsets else {})
        dataset.validate()
        _write(args.output, asdict(dataset)); result = args.output
    else:
        _write(args.output, compare_metrics(_load(args.metrics), _load(args.targets))); result = args.output
    print(json.dumps({"command": args.command, "output": str(result)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
