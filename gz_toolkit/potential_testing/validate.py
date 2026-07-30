"""Plain-language config validation.

Run BEFORE submitting anything to the cluster — every check produces a
sentence a materials scientist can act on without reading the source.
"""

from __future__ import annotations

import re
from pathlib import Path

from gz_toolkit.potential_testing.config import (
    PotentialConfig,
    load_potential_config,
    resolve_composition,
)

_WALLTIME_RE = re.compile(r"^\d+:\d{2}:\d{2}$|^\d+-\d+:\d{2}:\d{2}$")


def validate_config(cfg: PotentialConfig) -> list[str]:
    """Return a list of plain-language problems; an empty list means valid."""
    problems: list[str] = []
    pot = cfg.potential
    wf = cfg.workflow
    hpc = cfg.hpc

    # --- type_map / masses consistency -----------------------------------
    type_ids = set(pot.type_map.values())
    mass_ids = set(pot.masses.keys())
    for tid in sorted(type_ids - mass_ids):
        els = [e for e, i in pot.type_map.items() if i == tid]
        problems.append(
            f"Atom type {tid} (element {', '.join(els)}) is in type_map but has "
            f"no mass — add \"{tid}\": <mass in amu> to \"masses\"."
        )
    for tid in sorted(mass_ids - type_ids):
        problems.append(
            f"\"masses\" defines atom type {tid}, but no element in type_map "
            f"uses that type — remove it or fix type_map."
        )

    # --- elements ----------------------------------------------------------
    for el in pot.single_elements:
        if el not in pot.type_map:
            problems.append(
                f"Element '{el}' is listed in single_elements but missing from "
                f"type_map — add \"{el}\": <type id> to \"type_map\"."
            )
        if el not in pot.crystal_structures:
            problems.append(
                f"Element '{el}' is listed in single_elements but has no entry in "
                f"crystal_structures — add e.g. \"{el}\": {{\"structure\": \"bcc\"}}."
            )

    for gas in pot.gases:
        if gas not in pot.type_map:
            problems.append(
                f"Gas '{gas}' is listed in gases but missing from type_map — "
                f"add \"{gas}\": <type id> to \"type_map\"."
            )

    for pair in pot.gas_complex_pairs:
        if len(pair) != 2:
            problems.append(
                f"gas_complex_pairs entry {pair} must be a [metal, gas] pair."
            )
            continue
        for name in pair:
            if name not in pot.type_map:
                problems.append(
                    f"gas_complex_pairs mentions '{name}', which is not in type_map."
                )

    # --- alloys --------------------------------------------------------------
    for suite in pot.two_element_suites:
        A, B = suite.get("A"), suite.get("B")
        if not (A and B):
            problems.append(
                f"two_element_suites entry {suite} needs both \"A\" and \"B\" elements."
            )
            continue
        for name in (A, B):
            if name not in pot.type_map:
                problems.append(
                    f"Alloy element '{name}' (suite {A}-{B}) is not in type_map."
                )
        for frac in suite.get("fractions_atpct", []):
            try:
                f = float(frac)
            except (TypeError, ValueError):
                problems.append(
                    f"Alloy fraction '{frac}' in suite {A}-{B} is not a number."
                )
                continue
            if not (0.0 < f < 100.0):
                problems.append(
                    f"Alloy fraction {f} at.% in suite {A}-{B} must be between "
                    f"0 and 100 (exclusive)."
                )

    # --- multi-element (3+) alloys --------------------------------------------
    for suite in pot.multi_element_suites:
        composition = suite.get("composition")
        if not isinstance(composition, dict) or len(composition) < 2:
            problems.append(
                f"multi_element_suites entry {suite} needs a \"composition\" "
                f"dict with at least 2 elements, e.g. {{\"Fe\": null, \"Cr\": 3}}."
            )
            continue
        for name in composition:
            if name not in pot.type_map:
                problems.append(
                    f"multi_element_suites composition mentions '{name}', which "
                    f"is not in type_map."
                )
        try:
            resolve_composition(composition)
        except ValueError as exc:
            problems.append(str(exc))
        ordering = suite.get("ordering", "random")
        requested = [o.strip().lower() for o in str(ordering).split("|") if o.strip()]
        if requested and "random" not in requested:
            problems.append(
                f"multi_element_suites entry {suite} requests ordering "
                f"'{ordering}', but only 'random' is supported for 3+ element "
                f"alloys."
            )

    # --- SEAKMC ---------------------------------------------------------------
    if wf.seakmc_enabled:
        unknown = [d for d in wf.seakmc_defects if d not in wf.defect_catalog]
        for d in unknown:
            problems.append(
                f"seakmc_defects contains '{d}', which is not in defect_catalog — "
                f"SEAKMC can only run on defects that are actually built."
            )

    # --- numbers / formats -----------------------------------------------------
    if wf.paral_degree not in (1, 2, 3, 4):
        problems.append(
            f"paral_degree is {wf.paral_degree}, but must be 1, 2, 3 or 4 "
            f"(1 = everything in one job ... 4 = one job per simulation)."
        )
    if pot.lc_initial <= 0:
        problems.append(f"lc_initial must be positive (got {pot.lc_initial}).")
    for name in ("size_single", "size_alloy", "size_elastic", "size_elastic_alloy"):
        v = getattr(pot, name)
        if int(v) < 1:
            problems.append(f"{name} must be a positive integer (got {v}).")
    if wf.elastic_strain <= 0:
        problems.append(
            f"elastic_strain must be positive (got {wf.elastic_strain}); "
            f"the LAMMPS ELASTIC default is 1.0e-6."
        )
    if not _WALLTIME_RE.match(hpc.walltime):
        problems.append(
            f"walltime '{hpc.walltime}' does not look like a SLURM time "
            f"(expected HH:MM:SS or D-HH:MM:SS, e.g. 24:00:00)."
        )
    if hpc.ntasks < 1:
        problems.append(f"ntasks must be at least 1 (got {hpc.ntasks}).")

    if wf.elastic_for_alloys and not pot.two_element_suites and not pot.multi_element_suites:
        problems.append(
            "elastic_for_alloys is on, but two_element_suites and "
            "multi_element_suites are both empty — there are no alloy "
            "compositions to compute Cij for."
        )

    if wf.include_alloy_defects:
        if not pot.two_element_suites and not pot.multi_element_suites:
            problems.append(
                "include_alloy_defects is on, but two_element_suites and "
                "multi_element_suites are both empty — there are no alloy "
                "compositions to build point defects for."
            )
        if wf.alloy_defect_replicas < 1:
            problems.append(
                f"alloy_defect_replicas must be a positive integer (got "
                f"{wf.alloy_defect_replicas}) — each alloy build is randomized, "
                f"so at least one replica is needed."
            )

    return problems


def validate_project(pot_inputs_dir: str | Path = "pot_inputs") -> dict[str, list[str]]:
    """Validate every config file in a pot_inputs directory.

    Returns {filename: [problems]} — files that fail to load report the load
    error as their single problem. An all-empty dict of lists means the
    project is good to go.
    """
    root = Path(pot_inputs_dir)
    if not root.is_dir():
        return {str(root): [f"pot_inputs directory not found: {root}"]}

    report: dict[str, list[str]] = {}
    found = False
    for f in sorted(root.iterdir()):
        if f.suffix.lower() not in (".json", ".yaml", ".yml"):
            continue
        found = True
        try:
            cfg = load_potential_config(f)
        except Exception as exc:  # noqa: BLE001 - report any load failure
            report[f.name] = [f"Could not load this config file: {exc}"]
            continue
        report[f.name] = validate_config(cfg)
    if not found:
        report[str(root)] = ["No .json/.yaml config files found here."]
    return report
