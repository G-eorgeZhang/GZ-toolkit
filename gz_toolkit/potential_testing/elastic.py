"""Elastic-constant (full 6x6 Cij) case generation.

Uses the LAMMPS examples/ELASTIC finite-deformation scheme: four scripts
(in.elastic, init.mod, potential.mod, displace.mod) run 6 box deformations
(each +/-) and print the complete stiffness tensor. The templates shipped in
``templates/`` are adapted so the cell comes from ``structure.data`` (built
here at the reference-relaxed lattice constant) and the potential from the
per-case ``potential.inc`` — the same path therefore serves pure elements,
random alloys, and hcp cells.

Each elastic case is self-contained: ``in.elastic`` starts with a
``fix box/relax aniso 0.0`` minimize, so it only depends on the reference
stage (for the relaxed lc), never on other work groups.
"""

from __future__ import annotations

from importlib import resources
from pathlib import Path
from typing import Any

from gz_toolkit.potential_testing.config import PotentialConfig, resolve_composition
from gz_toolkit.potential_testing.structure_ops import (
    build_alloy_structure,
    build_multi_alloy_structure,
    build_reference_structure,
)

ELASTIC_INPUT_NAME = "in.elastic"

# The 9 main stiffness components surfaced in summary.csv; the remaining
# off-diagonal couplings (C14..C56) go to summary.json only.
MAIN_CIJ = ("C11", "C22", "C33", "C12", "C13", "C23", "C44", "C55", "C66")
ALL_CIJ = MAIN_CIJ + (
    "C14", "C15", "C16", "C24", "C25", "C26",
    "C34", "C35", "C36", "C45", "C46", "C56",
)


def _read_template(name: str) -> str:
    return (
        resources.files("gz_toolkit.potential_testing")
        .joinpath("templates", name)
        .read_text(encoding="utf-8")
    )


def render_elastic_inputs(config: PotentialConfig) -> dict[str, str]:
    """Render the four ELASTIC scripts for one case directory.

    Returns {filename: text} for in.elastic / init.mod / potential.mod /
    displace.mod. Only init.mod carries substitutions (strain magnitude and
    minimize settings); the other three are static.
    """
    tol = f"{config.workflow.minimize_tol:g}"
    init = (
        _read_template("init.mod.tpl")
        .replace("__UP__", f"{config.workflow.elastic_strain:g}")
        .replace("__ETOL__", tol)
        .replace("__FTOL__", tol)
        .replace("__MAXITER__", "100000")
        .replace("__MAXEVAL__", "100000")
    )
    return {
        "in.elastic": _read_template("in.elastic"),
        "init.mod": init,
        "potential.mod": _read_template("potential.mod"),
        "displace.mod": _read_template("displace.mod"),
    }


def _write_elastic_scripts(case_dir: Path, rendered: dict[str, str]) -> None:
    for name, text in rendered.items():
        (case_dir / name).write_text(text, encoding="utf-8")


def build_elastic_cases(
    config: PotentialConfig,
    pot_dir: Path,
    project_root: Path,
    per_element_lc: dict[str, float],
) -> list[dict[str, Any]]:
    """Create <pot>/elastic/<case>/ directories.

    Cases:
      * one per element in ``single_elements`` (cell at the relaxed lc,
        ``size_elastic`` replications);
      * one per alloy composition in ``two_element_suites`` and
        ``multi_element_suites`` when both ``include_elastic`` and
        ``elastic_for_alloys`` are set (single random realization at
        ``size_elastic_alloy`` — an estimate).

    Returns manifest entries mirroring the other work groups.
    """
    # Imported lazily: pipeline imports this module at load time.
    from gz_toolkit.potential_testing.pipeline import (
        _copy_potential_files,
        _write_potential_include,
    )

    manifest: list[dict[str, Any]] = []
    rendered = render_elastic_inputs(config)
    elements = config.potential.single_elements or [next(iter(config.potential.type_map))]

    def _setup(case_dir: Path) -> None:
        case_dir.mkdir(parents=True, exist_ok=True)
        _write_potential_include(case_dir, config.potential.pot_lines, config.potential.masses)
        _copy_potential_files(project_root, config.potential.pot_name, case_dir,
                              config.potential.potential_files)
        _write_elastic_scripts(case_dir, rendered)

    # ---- pure elements ----
    for element in elements:
        struct = config.potential.crystal_structures.get(element, {}).get("structure", "bcc")
        lc = per_element_lc.get(element, config.potential.lc_initial)
        case_dir = pot_dir / "elastic" / element
        _setup(case_dir)
        build_reference_structure(
            case_dir, config, element, struct,
            lc_value=lc, n_replicate=config.potential.size_elastic,
        )
        manifest.append({
            "element": element, "case": "elastic",
            "dir": str(case_dir.relative_to(project_root)),
        })

    # ---- alloys (opt-in) ----
    if config.workflow.elastic_for_alloys:
        for suite in config.potential.two_element_suites:
            A = suite.get("A")
            B = suite.get("B")
            if not (A and B):
                continue
            ordering = suite.get("ordering", "random")
            struct = config.potential.crystal_structures.get(A, {}).get("structure", "bcc")
            lc_A = per_element_lc.get(A, config.potential.lc_initial)
            for frac in suite.get("fractions_atpct", []):
                f_int = int(round(float(frac)))
                tag = f"{A}{B}_{A}{100 - f_int}_{B}{f_int}"
                case_dir = pot_dir / "elastic" / tag
                _setup(case_dir)
                build_alloy_structure(
                    case_dir, config, A, B, float(frac), struct,
                    ordering=ordering, relaxed_lc_A=lc_A,
                    n_replicate=config.potential.size_elastic_alloy,
                )
                manifest.append({
                    "A": A, "B": B, "fraction_B": frac, "case": "elastic_alloy",
                    "dir": str(case_dir.relative_to(project_root)),
                })

        for suite in config.potential.multi_element_suites:
            composition_raw = suite.get("composition", {})
            if not isinstance(composition_raw, dict) or len(composition_raw) < 2:
                continue
            try:
                composition = resolve_composition(composition_raw)
            except ValueError:
                continue  # already reported by validate_config
            suite_elements = list(composition.keys())
            host = suite_elements[0]
            ordering = suite.get("ordering", "random")
            struct = suite.get("structure") or config.potential.crystal_structures.get(host, {}).get("structure", "bcc")
            lc_host = per_element_lc.get(host, config.potential.lc_initial)
            pair_tag = "".join(suite_elements)
            comp_tag = "_".join(f"{el}{int(round(pct))}" for el, pct in composition.items())
            case_dir = pot_dir / "elastic" / f"{pair_tag}_{comp_tag}"
            _setup(case_dir)
            build_multi_alloy_structure(
                case_dir, config, composition, struct,
                ordering=ordering, relaxed_lc_host=lc_host,
                n_replicate=config.potential.size_elastic_alloy,
            )
            manifest.append({
                "composition": composition, "case": "elastic_alloy_multi",
                "dir": str(case_dir.relative_to(project_root)),
            })

    return manifest
