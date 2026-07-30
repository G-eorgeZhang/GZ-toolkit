"""Project bootstrap helpers for potential-testing workflows.

Layout produced by ``init_potential_project``::

    <root>/
    ├── main.py            # one-liner driver
    ├── summarize.py       # one-liner summary
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


DEFAULT_SUMMARY = """\"\"\"Aggregate results from every <pot_name>/ workspace into one CSV.\"\"\"

from gz_toolkit.potential_testing import summarize_all


def main() -> None:
    out = summarize_all(pot_inputs_dir="pot_inputs", run_dir=".")
    print(f"Summary written to: {out}")


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

    # Drivers.
    main_path = root / "main.py"
    if not main_path.exists():
        main_path.write_text(DEFAULT_MAIN, encoding="utf-8")
    summary_path = root / "summarize.py"
    if not summary_path.exists():
        summary_path.write_text(DEFAULT_SUMMARY, encoding="utf-8")

    return root
