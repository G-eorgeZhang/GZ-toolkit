"""Project bootstrap helpers for potential-testing workflows.

Layout produced by ``init_potential_project``::

    <root>/
    ├── main.py            # one-liner driver
    ├── summarize.py       # one-liner summary
    ├── check.py           # user-initiated pre-flight check, run before main.py
    ├── pot_inputs/
    │   └── <pot_name>.json
    └── potentials/
        └── <pot_name>/
"""

from __future__ import annotations

from pathlib import Path

from gz_toolkit.potential_testing.config import (
    PotentialConfig,
    PotentialMetadata,
    save_potential_config,
)


DEFAULT_MAIN = """\"\"\"Entry point — submits all potential-testing jobs for every config in pot_inputs/.\"\"\"

from gz_toolkit.potential_testing import run_all


def main() -> None:
    run_all(pot_inputs_dir="pot_inputs", run_dir=".")


if __name__ == "__main__":
    main()
"""


DEFAULT_SUMMARY = """\"\"\"Aggregate results from every <pot_name>/ workspace into CSV(s).

Parses each potential's log.lammps files directly — safe to run any time,
even while some potentials are still running (unfinished cases just show up
as blank cells). Independent of <pot_name>.json / `summarize-pot` (that's a
separate, automatic step for gz_toolkit.pot_infobank — see the manual).

SUMMARY_GROUPS controls how potentials are split across output files:
  {}                                          -> single summary.csv, everyone
  {"summary1": ["pot1", "pot2"]}              -> summary1.csv with only those pots
  {"summary1": [...], "summary2": "rest"}     -> summary1.csv + summary2.csv (everyone else)
\"\"\"

from gz_toolkit.potential_testing import write_grouped_summaries

SUMMARY_GROUPS: dict[str, list[str] | str] = {}


def main() -> None:
    out_paths = write_grouped_summaries(pot_inputs_dir="pot_inputs", run_dir=".", groups=SUMMARY_GROUPS)
    for p in out_paths:
        print(f"Summary written to: {p}")


if __name__ == "__main__":
    main()
"""


DEFAULT_CHECK = """\"\"\"Sanity-check potential files before running main.py.

User-initiated only — nothing else in the pipeline calls this. Run it by hand
after editing pot_inputs/*.json and dropping files into potentials/<pot_name>/,
and again any time you touch either, to catch the usual mistakes: a filename
in pot_lines that doesn't match potential_files, a potential file listed but
never copied into potentials/<pot_name>/, or (for MTP) an mlip.ini pointing at
the wrong .mtp file.
\"\"\"

from gz_toolkit.potential_testing import check_project


def main() -> None:
    report = check_project(pot_inputs_dir="pot_inputs", run_dir=".")
    n_errors = 0
    for pot_name, findings in report.items():
        print(f"{pot_name}:")
        if not findings:
            print("  OK")
            continue
        for level, msg in findings:
            if level == "error":
                n_errors += 1
                prefix = "  - "
            else:
                prefix = "  [info] "
            print(prefix + msg.replace("\\n", "\\n    "))
    if n_errors:
        print(f"\\n{n_errors} problem(s) found - fix them before running main.py.")
    else:
        print("\\nAll potential files check out.")


if __name__ == "__main__":
    main()
"""


def default_config(pot_name: str = "tao_2.31") -> PotentialConfig:
    """Generate a default per-potential config (BCC Fe + He + Cr example)."""
    return PotentialConfig(
        potential=PotentialMetadata(
            pot_name=pot_name,
            pot_lines="pair_style    mlip mlip.ini\npair_coeff    *  *",
            type_map={"Fe": 1, "He": 2, "Cr": 3},
            masses={1: 55.845, 2: 4.0026, 3: 51.9961},
            single_elements=["Fe", "Cr"],
            gases=["He"],
            gas_max_m=4,
            gas_max_n=4,
            gas_complex_pairs=[],   # empty -> all (metal, gas) combos
            gas_with_vacancy=True,
            gas_with_interstitial=True,
            two_element_suites=[
                {"A": "Fe", "B": "Cr", "fractions_atpct": [3, 5, 8, 10, 50], "ordering": "random"}
            ],
            crystal_structures={"Fe": {"structure": "bcc"}, "Cr": {"structure": "bcc"}},
            lc_initial=2.85,
            size_single=20,
            size_alloy=20,
            potential_files=["mlip.ini"],
            kind="classical",
        )
    )


def write_driver_scripts(
    root_dir: str | Path = ".",
    force: bool = False,
    only: list[str] | None = None,
) -> list[Path]:
    """Drop ``main.py`` / ``summarize.py`` / ``check.py`` into ``root_dir``.

    Independent of ``pot_inputs/`` — safe to run against a project whose
    ``pot_inputs/*.json`` were written by hand (never went through ``init``),
    to pull in the same driver scripts ``init`` would have created.

    ``force=False`` (default): existing files are left untouched — idempotent,
    only fills in whatever's missing.
    ``force=True``: overwrites existing files with the current template too —
    use this to pull in template updates (e.g. a project scaffolded before
    ``SUMMARY_GROUPS`` existed in ``summarize.py``). Any hand edits to an
    overwritten file are lost.
    ``only``: restrict to a subset of ``{"main.py", "summarize.py", "check.py"}``
    (default: all three) — e.g. ``only=["summarize.py"]`` to upgrade just that
    one file without touching a hand-edited ``main.py``/``check.py``.

    Returns the paths actually written.
    """
    root = Path(root_dir)
    root.mkdir(parents=True, exist_ok=True)
    templates = {
        "main.py": DEFAULT_MAIN,
        "summarize.py": DEFAULT_SUMMARY,
        "check.py": DEFAULT_CHECK,
    }
    names = only if only is not None else list(templates)
    written: list[Path] = []
    for name in names:
        path = root / name
        if force or not path.exists():
            path.write_text(templates[name], encoding="utf-8")
            written.append(path)
    return written


def init_potential_project(root_dir: str | Path = ".", pot_name: str = "tao_2.31") -> Path:
    """Create the project skeleton: ``main.py``, ``summarize.py``, ``pot_inputs/``, ``potentials/``.

    A starter JSON for ``pot_name`` is dropped into ``pot_inputs/`` so the user
    can copy it for additional potentials.
    """
    root = Path(root_dir)
    root.mkdir(parents=True, exist_ok=True)

    # potentials/<pot_name>/ — where the user drops their pair-style files.
    pot_files_dir = root / "potentials" / pot_name
    pot_files_dir.mkdir(parents=True, exist_ok=True)
    (pot_files_dir / ".gitkeep").write_text("", encoding="utf-8")

    # pot_inputs/ — one JSON per potential.
    pot_inputs_dir = root / "pot_inputs"
    pot_inputs_dir.mkdir(parents=True, exist_ok=True)

    cfg_path = pot_inputs_dir / f"{pot_name}.json"
    if not cfg_path.exists():
        save_potential_config(default_config(pot_name=pot_name), cfg_path)

    write_driver_scripts(root)

    return root
