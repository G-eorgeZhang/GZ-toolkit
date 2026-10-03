"""Reference -> defects -> energies -> SEAKMC pipeline.

The flow per potential::

    prepare_reference(cfg)
        # Builds <pot>/reference/<element>/ for each metal in single_elements
        # and writes box-relax LAMMPS inputs.

    build_defects_post_reference(cfg)
        # Reads each <pot>/reference/<element>/log.lammps to extract lc/Ecoh
        # (used only as a scale target -- every case below is built fresh
        # from a lc=1 seed, see structure_ops.py's module docstring), then builds:
        #   <pot>/point_defects/<element>/<defect>/
        #   <pot>/loops/<element>/<loop>/                  (opt-in, BCC only)
        #   <pot>/disloc_lines/<element>/<line>/           (opt-in, NotImplementedError)
        #   <pot>/alloy_lc/<A><B>/<A_X_B_Y>/               (opt-in; binary scan, two_element_suites)
        #   <pot>/alloy_lc/<A><B><C>.../<comp_tag>/        (opt-in; fixed N-element point, multi_element_suites)
        #   <pot>/point_defects_alloy/<pair>/<comp>/rep<k>/<case>/  (opt-in; alloy Ef, replica-averaged)
        #   <pot>/gas_complexes/<metal>-<gas>/<...>/       (opt-in)
        # Each case dir gets structure.data + in.relax.lammps + potential.inc.

    summarize(cfg)
        # Long-form CSV (kept for back-compat; see summary.py for wide CSV).

The job-script generation lives in ``parallel.py``. The end-to-end driver
(``run_all`` / ``summarize_all`` in workflow.py) wires these stages together.
"""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import re
import shutil
from typing import Any

from gz_toolkit.potential_testing.config import PotentialConfig, expand_gas_pairs, resolve_composition
from gz_toolkit.potential_testing.elastic import build_elastic_cases
from gz_toolkit.potential_testing.seakmc_runner import emit_seakmc_inputs, is_seakmc_eligible
from gz_toolkit.potential_testing.structure_ops import (
    build_alloy_structure,
    build_case_structure,
    build_gas_complex_structure,
    build_multi_alloy_structure,
    build_reference_structure,
    build_single_atom_structure,
    build_single_gas_in_bulk,
)
from gz_toolkit.potential_testing.summary import write_summary


# ---------------------------------------------------------------------------
# LAMMPS input + potential.inc writers
# ---------------------------------------------------------------------------


def _write_potential_include(case_dir: Path, pot_lines: str, masses: dict[int, float]) -> None:
    lines = ["# auto-generated potential include", ""]
    for tid in sorted(masses):
        lines.append(f"mass {tid} {masses[tid]}")
    lines.extend(["", pot_lines.strip(), ""])
    (case_dir / "potential.inc").write_text("\n".join(lines), encoding="utf-8")


def _copy_potential_files(base_root: Path, pot_name: str, case_dir: Path, potential_files: list[str]) -> None:
    src_root = base_root / "potentials" / pot_name
    for fname in potential_files:
        src = src_root / fname
        if src.exists():
            shutil.copy2(src, case_dir / fname)


def render_min_input(
    *,
    box_relax: bool,
    n_replicate: int | None = None,
    emit_ecoh: bool = False,
    minimize_tol: float = 1.0e-18,
    run_zero: bool = False,
) -> str:
    """Render a LAMMPS minimize input (shared by reference/relax/alloy writers).

    box_relax     : add ``fix box/relax iso 0.0`` around the minimize.
    n_replicate   : emit ``variable nrep`` + ``RESULT LC`` (lx / nrep).
    emit_ecoh     : emit ``RESULT ECOH`` (pe / natoms, i.e. energy per atom —
                    NOT true cohesive energy; combine with the isolated
                    single-atom run to get E_coh) — reference runs only.
    minimize_tol  : etol and ftol for the minimize command.
    run_zero      : insert ``run 0`` before the minimize (defect relax runs).

    Every rendered input reads ``structure.data``, includes ``potential.inc``,
    writes ``relaxed.data`` and prints ``RESULT NATOMS`` / ``RESULT PE_TOTAL``.
    """
    tol = f"{minimize_tol:g}"
    lines = [
        "units           metal",
        "dimension       3",
        "boundary        p p p",
        "atom_style      atomic",
        "atom_modify     map array sort 0 0.0",
        "",
        "read_data       structure.data",
        "include         potential.inc",
        "",
    ]
    if n_replicate is not None:
        lines += [f"variable nrep equal {int(n_replicate)}", ""]
    lines += [
        "thermo 100",
        "thermo_style    custom step atoms pe ke etotal temp press vol pxx pyy pzz lx ly lz",
        "thermo_modify   lost ignore format 8 %15.15g",
        "",
        "neighbor        1.0 bin",
        "neigh_modify    every 1 delay 0 check yes",
    ]
    if box_relax:
        lines.append("fix             1 all box/relax iso 0.0 vmax 0.001")
    lines.append("min_style       cg")
    if run_zero:
        lines.append("run             0")
    lines.append(f"minimize        {tol} {tol} 100000 100000")
    if box_relax:
        lines.append("unfix           1")
    lines += [
        "",
        "# write a relaxed data file for downstream construction",
        "write_data      relaxed.data",
        "",
        "variable        natoms equal atoms",
        "variable        pe_tot equal pe",
    ]
    if emit_ecoh:
        lines.append("variable        ecoh   equal pe/v_natoms")
    if n_replicate is not None:
        lines.append("variable        lc_rel equal lx/v_nrep")
        lines.append('print           "RESULT LC ${lc_rel}"')
    if emit_ecoh:
        lines.append('print           "RESULT ECOH ${ecoh}"')
    lines += [
        'print           "RESULT NATOMS ${natoms}"',
        'print           "RESULT PE_TOTAL ${pe_tot}"',
        "",
    ]
    return "\n".join(lines)


def _write_reference_input(case_dir: Path, n_replicate: int, tol: float = 1.0e-18) -> None:
    """Reference input: box/relax + write_data so the supercell is reusable."""
    text = render_min_input(box_relax=True, n_replicate=n_replicate, emit_ecoh=True, minimize_tol=tol)
    (case_dir / "in.reference.lammps").write_text(text, encoding="utf-8")


def _write_relax_input(case_dir: Path, tol: float = 1.0e-18) -> None:
    """Relax input for defect cases: minimize at fixed box, write energies."""
    text = render_min_input(box_relax=False, run_zero=True, minimize_tol=tol)
    (case_dir / "in.relax.lammps").write_text(text, encoding="utf-8")


def _write_alloy_lc_input(case_dir: Path, n_replicate: int, tol: float = 1.0e-18) -> None:
    """Alloy box-relax input — same as reference but for the alloy supercell."""
    text = render_min_input(box_relax=True, n_replicate=n_replicate, minimize_tol=tol)
    (case_dir / "in.alloy.lammps").write_text(text, encoding="utf-8")


def render_single_atom_input() -> str:
    """One isolated atom, no minimization needed (nothing to relax) — just
    read the box, compute the potential, and print its energy."""
    lines = [
        "units           metal",
        "dimension       3",
        "boundary        p p p",
        "atom_style      atomic",
        "atom_modify     map array sort 0 0.0",
        "",
        "read_data       structure.data",
        "include         potential.inc",
        "",
        "thermo 1",
        "thermo_style    custom step atoms pe",
        "",
        "neighbor        1.0 bin",
        "neigh_modify    every 1 delay 0 check yes",
        "",
        "run             0",
        "",
        "variable        pe_single equal pe",
        'print           "RESULT PE_SINGLE ${pe_single}"',
        "",
    ]
    return "\n".join(lines)


def _write_single_atom_input(case_dir: Path) -> None:
    text = render_single_atom_input()
    (case_dir / "in.single_atom.lammps").write_text(text, encoding="utf-8")


# ---------------------------------------------------------------------------
# Reference reading
# ---------------------------------------------------------------------------


def read_reference_lc(reference_dir: Path, default_lc: float) -> float:
    """Pull `RESULT LC` from log.lammps; fall back to default_lc."""
    logf = reference_dir / "log.lammps"
    if not logf.exists():
        return default_lc
    txt = logf.read_text(encoding="utf-8", errors="ignore")
    m = re.findall(r"RESULT\s+LC\s+([0-9Ee+\-\.]+)", txt)
    if not m:
        return default_lc
    try:
        return float(m[-1])
    except ValueError:
        return default_lc


def read_reference_e_per_atom(reference_dir: Path) -> float | None:
    """Read the ``RESULT ECOH`` line, which is really just ``pe/natoms``
    (energy per atom of the relaxed bulk cell) — not the true cohesive
    energy. Combine with :func:`read_reference_single_atom` to get the
    actual E_coh (see ``read_reference_ecoh``)."""
    logf = reference_dir / "log.lammps"
    if not logf.exists():
        return None
    txt = logf.read_text(encoding="utf-8", errors="ignore")
    m = re.findall(r"RESULT\s+ECOH\s+([0-9Ee+\-\.]+)", txt)
    if not m:
        return None
    try:
        return float(m[-1])
    except ValueError:
        return None


def read_reference_single_atom(reference_dir: Path) -> float | None:
    """Read the isolated single-atom energy from ``<reference_dir>/single_atom/log.lammps``."""
    logf = reference_dir / "single_atom" / "log.lammps"
    if not logf.exists():
        return None
    txt = logf.read_text(encoding="utf-8", errors="ignore")
    m = re.findall(r"RESULT\s+PE_SINGLE\s+([0-9Ee+\-\.]+)", txt)
    if not m:
        return None
    try:
        return float(m[-1])
    except ValueError:
        return None


def read_reference_ecoh(reference_dir: Path) -> float | None:
    """True cohesive energy: E_per_atom - E_single_atom.

    Returns ``None`` unless both the bulk reference run and the isolated
    single-atom run have results available.
    """
    e_per_atom = read_reference_e_per_atom(reference_dir)
    e_single = read_reference_single_atom(reference_dir)
    if e_per_atom is None or e_single is None:
        return None
    return e_per_atom - e_single


# ---------------------------------------------------------------------------
# Stage 1: reference per element
# ---------------------------------------------------------------------------


def prepare_reference_stage(config: PotentialConfig, run_dir: str | Path = ".") -> Path:
    """Build one reference cell per element listed in ``single_elements``."""
    root = Path(run_dir).resolve()
    pot_dir = root / config.potential.pot_name
    pot_dir.mkdir(parents=True, exist_ok=True)

    elements = config.potential.single_elements or [next(iter(config.potential.type_map))]
    for element in elements:
        struct = config.potential.crystal_structures.get(element, {}).get("structure", "bcc")
        ref_dir = pot_dir / "reference" / element
        ref_dir.mkdir(parents=True, exist_ok=True)
        build_reference_structure(ref_dir, config, element, struct)
        _write_potential_include(ref_dir, config.potential.pot_lines, config.potential.masses)
        _copy_potential_files(root, config.potential.pot_name, ref_dir, config.potential.potential_files)
        _write_reference_input(ref_dir, n_replicate=config.potential.size_single,
                               tol=config.workflow.minimize_tol)

        # Isolated-atom sanity check: same element, single atom, huge box.
        single_dir = ref_dir / "single_atom"
        single_dir.mkdir(parents=True, exist_ok=True)
        build_single_atom_structure(single_dir, config, element)
        _write_potential_include(single_dir, config.potential.pot_lines, config.potential.masses)
        _copy_potential_files(root, config.potential.pot_name, single_dir, config.potential.potential_files)
        _write_single_atom_input(single_dir)

    # Persist a snapshot of the config so later stages can verify what was run.
    snapshot = pot_dir / "config_snapshot.json"
    snapshot.write_text(json.dumps(asdict(config), indent=2), encoding="utf-8")
    return pot_dir


# Backwards-compat single-element entry point.
def prepare_reference_stage_single_element(config: PotentialConfig, run_dir: str | Path = ".") -> Path:
    pot_dir = prepare_reference_stage(config, run_dir=run_dir)
    elements = config.potential.single_elements or [next(iter(config.potential.type_map))]
    return pot_dir / "reference" / elements[0]


# ---------------------------------------------------------------------------
# Stage 2: build all defect / alloy / gas case directories
# ---------------------------------------------------------------------------


def _setup_case_dir(
    case_dir: Path,
    config: PotentialConfig,
    project_root: Path,
    input_writer=None,
) -> None:
    case_dir.mkdir(parents=True, exist_ok=True)
    _write_potential_include(case_dir, config.potential.pot_lines, config.potential.masses)
    _copy_potential_files(project_root, config.potential.pot_name, case_dir, config.potential.potential_files)
    if input_writer is None:
        _write_relax_input(case_dir, tol=config.workflow.minimize_tol)
    else:
        input_writer(case_dir)


def _iter_alloy_compositions(config: PotentialConfig):
    """Yield ``(pair_tag, comp_tag, composition, structure_name)`` for every
    composition named in ``two_element_suites`` (each scanned fraction) and
    ``multi_element_suites`` — the set of alloy points that
    ``include_alloy_defects`` builds point-defect/Ef cases for.
    """
    for suite in config.potential.two_element_suites:
        A = suite.get("A")
        B = suite.get("B")
        if not (A and B):
            continue
        struct = config.potential.crystal_structures.get(A, {}).get("structure", "bcc")
        for frac in suite.get("fractions_atpct", []):
            f = float(frac)
            composition = {A: 100.0 - f, B: f}
            pair_tag = f"{A}{B}"
            comp_tag = f"{A}{int(round(100.0 - f))}_{B}{int(round(f))}"
            yield pair_tag, comp_tag, composition, struct

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
        struct = suite.get("structure") or config.potential.crystal_structures.get(host, {}).get("structure", "bcc")
        pair_tag = "".join(suite_elements)
        comp_tag = "_".join(f"{el}{int(round(pct))}" for el, pct in composition.items())
        yield pair_tag, comp_tag, composition, struct


def build_defects_post_reference(config: PotentialConfig, run_dir: str | Path = ".") -> dict[str, Any]:
    """After reference has run, build every requested defect / alloy / gas case.

    Returns a manifest dict describing what was built — used by parallel.py to
    decide how to batch the resulting jobs.
    """
    root = Path(run_dir).resolve()
    pot_dir = root / config.potential.pot_name
    elements = config.potential.single_elements or [next(iter(config.potential.type_map))]

    manifest: dict[str, Any] = {
        "pot_name": config.potential.pot_name,
        "reference": {},
        "point_defects": [],
        "elastic": [],
        "loops": [],
        "disloc_lines": [],
        "alloys": [],
        "point_defects_alloy": [],
        "gas_complexes": [],
        "seakmc": [],
    }

    # ---- collect per-element relaxed lc / ecoh ----
    per_el_lc: dict[str, float] = {}
    per_el_ecoh: dict[str, float] = {}
    for element in elements:
        ref_dir = pot_dir / "reference" / element
        lc = read_reference_lc(ref_dir, config.potential.lc_initial)
        e_per_atom = read_reference_e_per_atom(ref_dir)
        e_single = read_reference_single_atom(ref_dir)
        ecoh = e_per_atom - e_single if e_per_atom is not None and e_single is not None else None
        per_el_lc[element] = lc
        if ecoh is not None:
            per_el_ecoh[element] = ecoh
        manifest["reference"][element] = {
            "lc": lc, "ecoh": ecoh, "e_per_atom": e_per_atom, "e_single_atom": e_single,
        }

    (pot_dir / "reference_values.json").write_text(json.dumps(manifest["reference"], indent=2), encoding="utf-8")

    # ---- point defects per element ----
    for element in elements:
        struct = config.potential.crystal_structures.get(element, {}).get("structure", "bcc")
        lc = per_el_lc[element]
        n_single = config.potential.size_single
        for case in config.workflow.defect_catalog:
            case_dir = pot_dir / "point_defects" / element / case
            _setup_case_dir(case_dir, config, root)
            try:
                build_case_structure(
                    case_dir, config, {element: 100.0}, case, struct,
                    target_lc=lc, n_replicate=n_single,
                )
                manifest["point_defects"].append({
                    "element": element, "case": case, "dir": str(case_dir.relative_to(root)),
                })
            except NotImplementedError as exc:
                manifest["point_defects"].append({
                    "element": element, "case": case, "skipped": str(exc),
                })

    # ---- elastic constants (6x6 Cij; elements always, alloys opt-in) ----
    if config.workflow.include_elastic:
        manifest["elastic"] = build_elastic_cases(config, pot_dir, root, per_el_lc)

    # ---- loops (BCC only, opt-in) ----
    if config.workflow.include_dislocation_loops:
        for element in elements:
            struct = config.potential.crystal_structures.get(element, {}).get("structure", "bcc")
            if struct.lower() != "bcc":
                continue
            lc = per_el_lc[element]
            for loop_case in ("SIL111", "SIL100"):
                case_dir = pot_dir / "loops" / element / loop_case
                _setup_case_dir(case_dir, config, root)
                build_case_structure(
                    case_dir, config, {element: 100.0}, loop_case, struct,
                    target_lc=lc, n_replicate=config.potential.size_single,
                )
                manifest["loops"].append({
                    "element": element, "case": loop_case, "dir": str(case_dir.relative_to(root)),
                })

    # ---- dislocation lines (placeholder; opt-in) ----
    if config.workflow.include_dislocation_lines:
        for element in elements:
            for line_case in ("edgedislo111", "edgedislo100", "screw111"):
                case_dir = pot_dir / "disloc_lines" / element / line_case
                case_dir.mkdir(parents=True, exist_ok=True)
                (case_dir / "TODO.txt").write_text(
                    "Dislocation-line construction is not implemented yet.\n"
                    "See gz_toolkit/potential_testing/structure_ops.py::_apply_dislocation_case.\n",
                    encoding="utf-8",
                )
                manifest["disloc_lines"].append({
                    "element": element, "case": line_case, "skipped": "not_implemented",
                })

    # ---- alloys (LC only) ----
    if config.workflow.include_alloy_suite:
        for suite in config.potential.two_element_suites:
            A = suite.get("A")
            B = suite.get("B")
            ordering = suite.get("ordering", "random")
            fractions = suite.get("fractions_atpct", [])
            if not (A and B):
                continue
            struct = config.potential.crystal_structures.get(A, {}).get("structure", "bcc")
            lc_A = per_el_lc.get(A, config.potential.lc_initial)
            for frac in fractions:
                f_int = int(round(float(frac)))
                tag = f"{A}{100 - f_int}_{B}{f_int}"
                case_dir = pot_dir / "alloy_lc" / f"{A}{B}" / tag
                _setup_case_dir(
                    case_dir,
                    config,
                    root,
                    input_writer=lambda d, n=config.potential.size_alloy,
                    t=config.workflow.minimize_tol: _write_alloy_lc_input(d, n, tol=t),
                )
                status = build_alloy_structure(
                    case_dir, config, A, B, float(frac), struct,
                    ordering=ordering, relaxed_lc_A=lc_A,
                )
                manifest["alloys"].append({
                    "A": A, "B": B, "fraction_B": frac, "ordering": status,
                    "dir": str(case_dir.relative_to(root)),
                })

        # ---- multi-element (3+) alloys: one fixed composition per entry ----
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
            lc_host = per_el_lc.get(host, config.potential.lc_initial)
            pair_tag = "".join(suite_elements)
            comp_tag = "_".join(f"{el}{int(round(pct))}" for el, pct in composition.items())
            case_dir = pot_dir / "alloy_lc" / pair_tag / comp_tag
            _setup_case_dir(
                case_dir,
                config,
                root,
                input_writer=lambda d, n=config.potential.size_alloy,
                t=config.workflow.minimize_tol: _write_alloy_lc_input(d, n, tol=t),
            )
            status = build_multi_alloy_structure(
                case_dir, config, composition, struct,
                ordering=ordering, relaxed_lc_host=lc_host,
            )
            manifest["alloys"].append({
                "composition": composition, "ordering": status,
                "dir": str(case_dir.relative_to(root)),
            })

    # ---- alloy point-defect formation energies (opt-in, replica-averaged) ----
    # Separate gate from include_alloy_suite: every alloy build is an
    # independent random draw (see structure_ops.py), so a meaningful Ef
    # needs several replicas averaged together -- materially more expensive
    # than the plain lc scan above.
    if config.workflow.include_alloy_defects:
        n_alloy = config.potential.size_alloy
        for pair_tag, comp_tag, composition, struct in _iter_alloy_compositions(config):
            host = next(iter(composition))
            lc_host = per_el_lc.get(host, config.potential.lc_initial)
            base = pot_dir / "point_defects_alloy" / pair_tag / comp_tag
            for rep in range(config.workflow.alloy_defect_replicas):
                for case in ["bulk"] + config.workflow.defect_catalog:
                    case_dir = base / f"rep{rep}" / case
                    _setup_case_dir(case_dir, config, root)
                    try:
                        build_case_structure(
                            case_dir, config, composition, case, struct,
                            target_lc=lc_host, n_replicate=n_alloy,
                        )
                        manifest["point_defects_alloy"].append({
                            "pair": pair_tag, "composition": comp_tag, "replica": rep,
                            "case": case, "dir": str(case_dir.relative_to(root)),
                        })
                    except NotImplementedError as exc:
                        manifest["point_defects_alloy"].append({
                            "pair": pair_tag, "composition": comp_tag, "replica": rep,
                            "case": case, "skipped": str(exc),
                        })

    # ---- gas complexes ----
    if config.workflow.include_gas_complexes and config.potential.gases:
        pairs = expand_gas_pairs(config.potential)
        for metal, gas in pairs:
            struct = config.potential.crystal_structures.get(metal, {}).get("structure", "bcc")
            lc = per_el_lc.get(metal, config.potential.lc_initial)
            base = pot_dir / "gas_complexes" / f"{metal}-{gas}"

            # Single-gas references (m=0) needed for binding-energy formulas.
            for site in ("tetra", "octa"):
                case_dir = base / f"{gas}1_{site}"
                _setup_case_dir(case_dir, config, root)
                build_single_gas_in_bulk(case_dir, metal, gas, site, config, struct, relaxed_lc=lc)
                manifest["gas_complexes"].append({
                    "metal": metal, "gas": gas, "n": 1, "m": 0, "kind": f"single_{site}",
                    "dir": str(case_dir.relative_to(root)),
                })

            kinds: list[str] = []
            if config.potential.gas_with_vacancy:
                kinds.append("vacancy")
            if config.potential.gas_with_interstitial:
                kinds.append("interstitial")

            for kind in kinds:
                vi = "V" if kind == "vacancy" else "I"
                for n in range(1, config.potential.gas_max_n + 1):
                    for m in range(1, config.potential.gas_max_m + 1):
                        case_dir = base / f"{gas}{n}_{vi}{m}"
                        _setup_case_dir(case_dir, config, root)
                        build_gas_complex_structure(
                            case_dir, metal, gas, n, m, kind,
                            config, struct, relaxed_lc=lc,
                        )
                        manifest["gas_complexes"].append({
                            "metal": metal, "gas": gas, "n": n, "m": m, "kind": kind,
                            "dir": str(case_dir.relative_to(root)),
                        })

    # ---- SEAKMC inputs for selected point defects ----
    if config.workflow.seakmc_enabled:
        for element in elements:
            struct = config.potential.crystal_structures.get(element, {}).get("structure", "bcc")
            lc = per_el_lc[element]
            for case in config.workflow.defect_catalog:
                if not is_seakmc_eligible(case, config.workflow.seakmc_defects):
                    continue
                case_dir = pot_dir / "point_defects" / element / case
                # Use case_dir/relaxed.data after the relax stage; if missing,
                # fall back to the un-relaxed structure.data so the file at
                # least exists. Job ordering ensures relaxed.data appears
                # before SEAKMC runs.
                relaxed = case_dir / "relaxed.data"
                if not relaxed.exists():
                    relaxed = case_dir / "structure.data"
                meb_dir = emit_seakmc_inputs(
                    case_dir, config, element, struct, lc, relaxed,
                )
                manifest["seakmc"].append({
                    "element": element, "case": case,
                    "dir": str(meb_dir.relative_to(root)),
                })

    (pot_dir / "build_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


# ---------------------------------------------------------------------------
# Backwards-compatible single-pot entry points (used by older tests/CLI)
# ---------------------------------------------------------------------------


def prepare_defect_stage(config: PotentialConfig, run_dir: str | Path = ".") -> Path:
    """Run stage 2 for a single potential (assumes stage 1 already completed)."""
    root = Path(run_dir).resolve()
    pot_dir = root / config.potential.pot_name
    build_defects_post_reference(config, run_dir=run_dir)
    return pot_dir


def summarize_stage(config: PotentialConfig, run_dir: str | Path = ".") -> Path:
    """Long-form CSV one-pot summary (kept for older tests)."""
    root = Path(run_dir).resolve()
    pot_dir = root / config.potential.pot_name
    rows: list[dict[str, Any]] = []
    if not pot_dir.is_dir():
        out = root / "summary.csv"
        write_summary(rows, out)
        return out

    for element_dir in sorted((pot_dir / "reference").glob("*")) if (pot_dir / "reference").is_dir() else []:
        if not element_dir.is_dir():
            continue
        lc = read_reference_lc(element_dir, config.potential.lc_initial)
        e_per_atom = read_reference_e_per_atom(element_dir)
        e_single = read_reference_single_atom(element_dir)
        ecoh = e_per_atom - e_single if e_per_atom is not None and e_single is not None else None
        case = f"reference/{element_dir.name}"
        rows.append({"pot_name": config.potential.pot_name, "case": case, "metric": "lc", "value": lc})
        if e_per_atom is not None:
            rows.append({"pot_name": config.potential.pot_name, "case": case,
                         "metric": "E_per_atom", "value": e_per_atom})
        if e_single is not None:
            rows.append({"pot_name": config.potential.pot_name, "case": case,
                         "metric": "E_single_atom", "value": e_single})
        if ecoh is not None:
            rows.append({"pot_name": config.potential.pot_name, "case": case,
                         "metric": "E_coh", "value": ecoh})

    out = root / "summary.csv"
    write_summary(rows, out)
    return out
