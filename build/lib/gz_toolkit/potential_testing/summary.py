"""Summary generation for potential-testing workflows.

Two independent consumers of ``_collect_pot_metrics`` (the log parser), kept
deliberately decoupled so one can't go stale/missing and break the other:

* ``build_pot_tag_json(cfg, run_dir)`` — writes the standardized per-pot
  record to ``<run_dir>/<pot_name>/<pot_name>.json`` once that pot's own
  tests finish (see ``gz_toolkit.potential_testing.cli`` ``summarize-pot``,
  wired to run automatically at the end of each pot's job chain). This is
  "JSON for code" — what ``gz_toolkit.pot_infobank.promote_tag_json`` copies
  elsewhere (same filename) for ``buildmtx``/``defect``/``analyze`` to read a
  potential's numbers back out by tag.
* ``write_wide_summary`` / ``write_grouped_summaries(pot_inputs_dir, run_dir,
  groups)`` — parses each pot's logs directly (same as ``build_pot_tag_json``,
  independently) and pivots into CSV(s) for humans. Doesn't touch or require
  that per-pot JSON — works even if that pot's automatic summarize-pot step
  hasn't run yet (or ever), as long as the raw logs exist.

Wide-format layout — rows are grouped by *system* (pure element, then gas,
then alloy pair, then gas-defect-complex pair — each in the order first seen
across ``pot_inputs/``), with a ``# <system>`` section marker row ahead of
each group::

    metric,                        pot1,    pot2,    pot3
    # Fe
    lc_Fe,                         2.8557,  2.8553,  2.8612
    Ecoh_Fe,                       -4.122,  ...
    E_per_atom_Fe,                 -4.122,  ...
    Esingle_Fe,                    0.0,     ...
    Ef_1vac_Fe,                    1.85,    ...
    Emig_1vac_Fe,                  0.62,    ...
    # Cr
    lc_Cr,                         2.8841,  ...
    ...
    # He
    ...
    # FeCr
    alloy_lc_FeCr_Fe95Cr5,         2.859,   ...
    # FeHe
    Ef_He_in_Fe_tetra,             ...
    Eb_V_He1_V1_Fe,                ...
"""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
import re
from pathlib import Path
from typing import Any

from gz_toolkit.analyze.energetics import (
    MAIN_CIJ,
    binding_energy,
    cohesive_energy,
    formation_energy,
    voigt_moduli,
)
from gz_toolkit.potential_testing.config import (
    PotentialConfig,
    expand_gas_pairs,
    load_potential_configs,
)


# ---------------------------------------------------------------------------
# Long-form (legacy)
# ---------------------------------------------------------------------------


def write_summary(rows: list[dict[str, Any]], output_csv: str | Path) -> None:
    """Long-form summary writer (kept for back-compat with old tests)."""
    out_path = Path(output_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        rows = [{"metric_group": "workflow", "metric_name": "status", "value": "no-data"}]
    fieldnames = sorted({key for row in rows for key in row.keys()})
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# Log parsing helpers
# ---------------------------------------------------------------------------


_RES_PATTERNS = {
    "natoms":     re.compile(r"RESULT\s+NATOMS\s+([0-9Ee+\-\.]+)"),
    "pe":         re.compile(r"RESULT\s+PE_TOTAL\s+([0-9Ee+\-\.]+)"),
    "lc":         re.compile(r"RESULT\s+LC\s+([0-9Ee+\-\.]+)"),
    # LAMMPS still tags this line "ECOH" but it is really just pe/natoms
    # (energy per atom); true cohesive energy is this minus the isolated
    # single-atom energy ("pe_single" below) — computed in _collect_pot_metrics.
    "e_per_atom": re.compile(r"RESULT\s+ECOH\s+([0-9Ee+\-\.]+)"),
    "pe_single":  re.compile(r"RESULT\s+PE_SINGLE\s+([0-9Ee+\-\.]+)"),
}


def _extract_results(log_path: Path) -> dict[str, float]:
    if not log_path.is_file():
        return {}
    txt = log_path.read_text(encoding="utf-8", errors="ignore")
    out: dict[str, float] = {}
    for key, pat in _RES_PATTERNS.items():
        m = pat.findall(txt)
        if m:
            try:
                out[key] = float(m[-1])
            except ValueError:
                pass
    return out


_ELASTIC_PATTERN = re.compile(r"RESULT\s+(C\d\d)\s+([0-9Ee+\-\.]+)")


def _extract_elastic(log_path: Path) -> dict[str, float]:
    """Parse ``RESULT Cij value`` lines from an in.elastic log (last wins)."""
    if not log_path.is_file():
        return {}
    txt = log_path.read_text(encoding="utf-8", errors="ignore")
    out: dict[str, float] = {}
    for comp, value in _ELASTIC_PATTERN.findall(txt):
        try:
            out[comp] = float(value)
        except ValueError:
            pass
    return out


# Row order for elastic_* keys within one tag: the 9 Cij first, then the
# Voigt moduli as Young's (E), Shear (G), Bulk (B), Poisson's ratio (nu).
# Voigt-modulus math itself lives in gz_toolkit.analyze.energetics.voigt_moduli.
_ELASTIC_COMPONENT_ORDER = {c: i for i, c in enumerate(MAIN_CIJ)}
_ELASTIC_COMPONENT_ORDER.update({
    "E": len(MAIN_CIJ),
    "G": len(MAIN_CIJ) + 1,
    "B": len(MAIN_CIJ) + 2,
    "nu": len(MAIN_CIJ) + 3,
})


def _seakmc_min_barrier(meb_dir: Path) -> float | None:
    """Return the minimum forward barrier from Seakmc_summary.csv."""
    summary = meb_dir / "Seakmc_summary.csv"
    if not summary.is_file():
        # Some SEAKMC builds write to SPOut/KMC_*_SPs.csv per step.
        return _seakmc_min_from_spout(meb_dir)
    barriers: list[float] = []
    with summary.open("r", encoding="utf-8", errors="ignore") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            for key in ("barrier", "Barrier", "forward_barrier"):
                if key in row and row[key]:
                    try:
                        barriers.append(float(row[key]))
                    except ValueError:
                        pass
                    break
    return min(barriers) if barriers else None


def _seakmc_min_from_spout(meb_dir: Path) -> float | None:
    spout = meb_dir / "SPOut"
    if not spout.is_dir():
        return None
    barriers: list[float] = []
    for csvfile in spout.glob("KMC_*_SPs.csv"):
        with csvfile.open("r", encoding="utf-8", errors="ignore") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                for key in ("barrier", "Barrier"):
                    if key in row and row[key]:
                        try:
                            barriers.append(float(row[key]))
                        except ValueError:
                            pass
                        break
    return min(barriers) if barriers else None


# ---------------------------------------------------------------------------
# Per-potential metric extraction
# ---------------------------------------------------------------------------


def _natoms_change_for_case(case_name: str) -> int:
    """Return the change in atom count vs the perfect cell, signed.

    Vacancies are negative; interstitials positive. Used in formation-energy
    formula:  E_f = E_def - (N_def / N_ref) * E_ref_total.
    """
    n = 0
    if case_name == "1vac":
        n = -1
    elif case_name.startswith("2vac"):
        n = -2
    elif case_name == "1int_tetra" or case_name == "1int_octa":
        n = +1
    elif case_name.startswith("dumbbell"):
        n = +1
    elif case_name == "2int":
        n = +2
    elif case_name == "3int":
        n = +3
    return n


def _gas_complex_atom_count(metal: str, gas: str, n_gas: int, m_def: int, kind: str,
                            type_map: dict[str, int]) -> tuple[int, int]:
    """Return (delta_metal, delta_gas) vs the perfect-metal reference cell."""
    if kind == "vacancy":
        return (-m_def, n_gas)
    return (m_def, n_gas)


def _collect_pot_metrics(
    cfg: PotentialConfig, run_dir: Path
) -> tuple[dict[str, float | str], dict[str, str]]:
    """Walk a single potential's workspace.

    Returns ``(metrics, system_of)``: ``metrics`` maps metric label -> value
    (unchanged shape, still dumped to summary.json); ``system_of`` maps that
    same label -> the "system" it belongs to (a pure element, a gas, an alloy
    pair tag like "FeCr", or a gas-complex pair tag like "FeHe"), recorded
    exactly where each metric is built rather than re-parsed from the label
    string afterwards (alloy/gas labels aren't unambiguous to reverse-parse).
    """
    pot_dir = run_dir / cfg.potential.pot_name
    out: dict[str, float | str] = {}
    system_of: dict[str, str] = {}
    if not pot_dir.is_dir():
        return out, system_of

    def put(label: str, value: float | str, system: str) -> None:
        out[label] = value
        system_of[label] = system

    elements = cfg.potential.single_elements or [next(iter(cfg.potential.type_map))]

    # Reference: lc, E/atom, isolated single-atom energy, true Ecoh, total PE, N_atoms.
    #   E_per_atom = pe / natoms of the relaxed bulk reference cell
    #   Esingle    = energy of one isolated atom (same potential, huge box)
    #   Ecoh       = E_per_atom - Esingle  (the actual cohesive energy)
    ref_total_pe: dict[str, float] = {}
    ref_natoms: dict[str, int] = {}
    for el in elements:
        log = pot_dir / "reference" / el / "log.lammps"
        res = _extract_results(log)
        if "lc" in res:
            put(f"lc_{el}", res["lc"], el)
        if "pe" in res:
            ref_total_pe[el] = res["pe"]
        if "natoms" in res:
            ref_natoms[el] = int(res["natoms"])

        single_log = pot_dir / "reference" / el / "single_atom" / "log.lammps"
        single_res = _extract_results(single_log)
        e_single = single_res.get("pe_single")
        if e_single is not None:
            put(f"Esingle_{el}", e_single, el)

        if "e_per_atom" in res:
            e_per_atom = res["e_per_atom"]
            put(f"E_per_atom_{el}", e_per_atom, el)
            if e_single is not None:
                put(f"Ecoh_{el}", cohesive_energy(e_per_atom, e_single), el)

    # Point defects: formation energies.
    pd_root = pot_dir / "point_defects"
    if pd_root.is_dir():
        for el_dir in sorted(pd_root.iterdir()):
            if not el_dir.is_dir():
                continue
            el = el_dir.name
            n_ref = ref_natoms.get(el)
            e_ref = ref_total_pe.get(el)
            if n_ref is None or e_ref is None:
                continue
            for case_dir in sorted(el_dir.iterdir()):
                if not case_dir.is_dir():
                    continue
                res = _extract_results(case_dir / "log.lammps")
                if "pe" not in res or "natoms" not in res:
                    continue
                n_def = int(res["natoms"])
                e_def = res["pe"]
                ef = formation_energy(e_def, n_def, e_ref, n_ref)
                put(f"Ef_{case_dir.name}_{el}", ef, el)

                # SEAKMC barrier.
                meb = case_dir / "meb"
                if meb.is_dir():
                    bar = _seakmc_min_barrier(meb)
                    if bar is not None:
                        put(f"Emig_{case_dir.name}_{el}", bar, el)

    # Elastic constants: 6x6 Cij (+ Voigt B/G/nu) per element / alloy case.
    elastic_root = pot_dir / "elastic"
    if elastic_root.is_dir():
        for case_dir in sorted(elastic_root.iterdir()):
            if not case_dir.is_dir():
                continue
            cij = _extract_elastic(case_dir / "log.lammps")
            if not cij:
                continue
            tag = case_dir.name  # element symbol or "<A><B>_<comp>" alloy tag
            # Pure-element case dirs are named exactly the element symbol;
            # alloy case dirs are "<A><B>_<comp>" — the pair tag is the part
            # before the first underscore.
            system = tag if tag in elements else tag.split("_", 1)[0]
            for comp in MAIN_CIJ:
                if comp in cij:
                    put(f"elastic_{comp}_{tag}", cij[comp], system)
            for comp, val in cij.items():
                if comp not in MAIN_CIJ:
                    # off-diagonal couplings: kept in summary.json only
                    put(f"elastic_raw_{comp}_{tag}", val, system)
            for name, val in voigt_moduli(cij).items():
                put(f"elastic_{name}_{tag}", val, system)

    # Loops: formation energy.
    loops_root = pot_dir / "loops"
    if loops_root.is_dir():
        for el_dir in sorted(loops_root.iterdir()):
            el = el_dir.name
            n_ref = ref_natoms.get(el)
            e_ref = ref_total_pe.get(el)
            if n_ref is None or e_ref is None:
                continue
            for case_dir in sorted(el_dir.iterdir()):
                if not case_dir.is_dir():
                    continue
                res = _extract_results(case_dir / "log.lammps")
                if "pe" not in res or "natoms" not in res:
                    continue
                ef = formation_energy(res["pe"], int(res["natoms"]), e_ref, n_ref)
                put(f"Ef_loop_{case_dir.name}_{el}", ef, el)

    # Alloys: lattice constant.
    alloy_root = pot_dir / "alloy_lc"
    if alloy_root.is_dir():
        for pair_dir in sorted(alloy_root.iterdir()):
            for case_dir in sorted(pair_dir.iterdir()):
                res = _extract_results(case_dir / "log.lammps")
                if "lc" in res:
                    put(f"alloy_lc_{pair_dir.name}_{case_dir.name}", res["lc"], pair_dir.name)

    # Alloy point defects: formation energy, averaged over independent random
    # replicas (each alloy build is a fresh solute placement — see
    # structure_ops.py — so a single realization is too noisy on its own).
    # Each replica pairs its own "bulk" (no-defect) case with each requested
    # defect case; Ef is per-atom (not per-species, unlike the pure-element
    # block above), since the defect cell's matrix is itself multi-species.
    alloy_pd_root = pot_dir / "point_defects_alloy"
    if alloy_pd_root.is_dir():
        for pair_dir in sorted(alloy_pd_root.iterdir()):
            if not pair_dir.is_dir():
                continue
            for comp_dir in sorted(pair_dir.iterdir()):
                if not comp_dir.is_dir():
                    continue
                per_case_efs: dict[str, list[float]] = {}
                for rep_dir in sorted(comp_dir.iterdir()):
                    if not rep_dir.is_dir():
                        continue
                    bulk_res = _extract_results(rep_dir / "bulk" / "log.lammps")
                    if "pe" not in bulk_res or "natoms" not in bulk_res:
                        continue
                    e_bulk = bulk_res["pe"]
                    n_bulk = int(bulk_res["natoms"])
                    for case_dir in sorted(rep_dir.iterdir()):
                        if case_dir.name == "bulk" or not case_dir.is_dir():
                            continue
                        res = _extract_results(case_dir / "log.lammps")
                        if "pe" not in res or "natoms" not in res:
                            continue
                        ef = formation_energy(res["pe"], int(res["natoms"]), e_bulk, n_bulk)
                        per_case_efs.setdefault(case_dir.name, []).append(ef)
                for case_name, efs in per_case_efs.items():
                    if not efs:
                        continue
                    mean_ef = sum(efs) / len(efs)
                    label = f"Ef_{case_name}_{pair_dir.name}_{comp_dir.name}"
                    put(label, mean_ef, pair_dir.name)
                    if len(efs) > 1:
                        variance = sum((x - mean_ef) ** 2 for x in efs) / (len(efs) - 1)
                        put(f"{label}_std", variance ** 0.5, pair_dir.name)

    # Gas complexes: formation + binding energy.
    gas_root = pot_dir / "gas_complexes"
    gas_ef: dict[tuple[str, str, int, int, str], float] = {}
    if gas_root.is_dir():
        for pair_dir in sorted(gas_root.iterdir()):
            # pair_dir name is "<metal>-<gas>"
            try:
                metal, gas = pair_dir.name.split("-", 1)
            except ValueError:
                continue
            gas_system = f"{metal}{gas}"
            n_ref = ref_natoms.get(metal)
            e_ref = ref_total_pe.get(metal)
            if n_ref is None or e_ref is None:
                continue
            for case_dir in sorted(pair_dir.iterdir()):
                res = _extract_results(case_dir / "log.lammps")
                if "pe" not in res or "natoms" not in res:
                    continue
                # Parse case_dir name: "<gas>n_<V|I>m" or "<gas>1_<site>" for singles.
                name = case_dir.name
                m_match = re.match(rf"{re.escape(gas)}(\d+)_([VI])(\d+)$", name)
                single_match = re.match(rf"{re.escape(gas)}1_(tetra|octa)$", name)
                e_def = res["pe"]
                n_def = int(res["natoms"])
                if single_match:
                    # Single gas in bulk: formation energy w.r.t. metal ref.
                    # N_metal_in_cell = N_ref (only the gas atom was added),
                    # so this degenerates to e_def - e_ref. Gas chemical
                    # potential is taken as zero (classical pots) — user can
                    # post-process.
                    ef = formation_energy(e_def, n_ref, e_ref, n_ref)
                    put(f"Ef_{gas}_in_{metal}_{single_match.group(1)}", ef, gas_system)
                    gas_ef[(metal, gas, 1, 0, single_match.group(1))] = ef
                elif m_match:
                    n_g = int(m_match.group(1))
                    kind = "vacancy" if m_match.group(2) == "V" else "interstitial"
                    m_d = int(m_match.group(3))
                    dM, _ = _gas_complex_atom_count(metal, gas, n_g, m_d, kind,
                                                    cfg.potential.type_map)
                    n_metal_in_cell = n_def - n_g
                    ef = formation_energy(e_def, n_metal_in_cell, e_ref, n_ref)
                    label = f"Ef_{gas}{n_g}_{m_match.group(2)}{m_d}_{metal}"
                    put(label, ef, gas_system)
                    gas_ef[(metal, gas, n_g, m_d, kind)] = ef

        # Binding energies: Eb_He^V = E_f(HenVm) - E_f(Hen V_{m-1}) - E_f(V)
        for (metal, gas, n_g, m_d, kind), ef in gas_ef.items():
            if m_d == 0 or n_g == 0:
                continue
            gas_system = f"{metal}{gas}"
            vi = "V" if kind == "vacancy" else "I"
            # Vacancy/interstitial single-defect Ef (need to find the right key)
            single_def_key = f"Ef_1vac_{metal}" if kind == "vacancy" else f"Ef_1int_octa_{metal}"
            if single_def_key not in out:
                continue
            ef_single_defect = out[single_def_key]
            # Binding via removing one defect atom (m -> m-1 vacancies/interstitials).
            ef_minus = gas_ef.get((metal, gas, n_g, m_d - 1, kind))
            if ef_minus is not None and isinstance(ef_minus, float):
                eb_v = binding_energy(ef, ef_minus, ef_single_defect)
                put(f"Eb_{vi}_{gas}{n_g}_{vi}{m_d}_{metal}", eb_v, gas_system)
            # Binding via removing one gas atom (n -> n-1 gas).
            ef_minus_g = gas_ef.get((metal, gas, n_g - 1, m_d, kind))
            single_gas_tetra = gas_ef.get((metal, gas, 1, 0, "tetra"))
            single_gas_octa = gas_ef.get((metal, gas, 1, 0, "octa"))
            ef_single_gas = single_gas_tetra if single_gas_tetra is not None else single_gas_octa
            if ef_minus_g is not None and ef_single_gas is not None:
                eb_g = binding_energy(ef, ef_minus_g, ef_single_gas)
                put(f"Eb_{gas}_{gas}{n_g}_{vi}{m_d}_{metal}", eb_g, gas_system)

    return out, system_of


# ---------------------------------------------------------------------------
# Per-pot <pot_name>.json — the standardized record, built once when a pot finishes
# ---------------------------------------------------------------------------


def build_pot_tag_json(cfg: PotentialConfig, run_dir: str | Path = ".") -> Path:
    """Parse this pot's logs and write ``<run_dir>/<pot_name>/<pot_name>.json``.

    Named after the potential (not a fixed ``tag.json``) so the exact same
    filename is used once promoted into ``gz_toolkit/pot_infobank/`` (see
    ``gz_toolkit.pot_infobank.promote_tag_json``) — no rename needed there.

    This is a minimal record — just what structure-building code needs to
    go from ``pot="gao2011"`` to an actual cell: ``pot_name``, ``type_map``,
    ``masses``, and ``lc``. ``lc`` has pure elements as top-level keys
    (``lc["Fe"]``) and alloy compositions nested under ``lc["alloys"]``,
    keyed by pair then composition tag (``lc["alloys"]["FeCr"]["Fe95_Cr5"]``)
    — see ``gz_toolkit.pot_infobank.get_lc`` for a lookup helper that
    resolves either form, including "nearest available composition" for an
    alloy fraction that wasn't exactly tested. Everything else
    ``_collect_pot_metrics`` can produce (Ecoh, elastic constants,
    point-defect formation energies, SEAKMC barriers, ...) is deliberately
    left out of this file — it's all still in ``summary.csv``
    (``write_wide_summary``/``write_grouped_summaries``), which is where to
    look for potential-characterization data. This file is only for "what
    lattice constant do I build potential X at."

    This is the only function that ever re-parses raw LAMMPS logs for
    summary purposes — ``gz_toolkit.pot_infobank`` only ever reads the
    resulting file. Call this once a pot's own tests are done (see the
    ``summarize-pot`` CLI subcommand); re-running it just overwrites the
    file with fresh numbers.
    """
    root = Path(run_dir).resolve()
    metrics, _system_of = _collect_pot_metrics(cfg, root)
    lc: dict[str, float | dict] = {}
    alloy_lc: dict[str, dict[str, float]] = {}
    for key, value in metrics.items():
        if key.startswith("lc_"):
            lc[key[len("lc_"):]] = value
        elif key.startswith("alloy_lc_"):
            # "<pair><comp>" e.g. "FeCr_Fe95_Cr5" — pair never contains an
            # underscore (it's a bare element-symbol concatenation), so the
            # first "_" always splits pair from the composition tag, even
            # though the composition tag itself has an internal "_".
            pair, _, comp = key[len("alloy_lc_"):].partition("_")
            alloy_lc.setdefault(pair, {})[comp] = value
    if alloy_lc:
        lc["alloys"] = alloy_lc
    record = {
        "pot_name": cfg.potential.pot_name,
        "type_map": cfg.potential.type_map,
        "masses": cfg.potential.masses,
        "lc": lc,
        "updated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    out_dir = root / cfg.potential.pot_name
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{cfg.potential.pot_name}.json"
    out_path.write_text(json.dumps(record, indent=2), encoding="utf-8")
    return out_path


# ---------------------------------------------------------------------------
# Wide-format aggregator
# ---------------------------------------------------------------------------


def _ordered_systems(configs: list[PotentialConfig]) -> list[str]:
    """Systems in display order: elements, then gases, then alloy pairs
    ("FeCr"), then gas-defect-complex pairs ("FeHe") — each de-duplicated in
    first-seen order across every config in ``pot_inputs/``."""
    systems: list[str] = []
    seen: set[str] = set()

    def add(tag: str) -> None:
        if tag and tag not in seen:
            seen.add(tag)
            systems.append(tag)

    for cfg in configs:
        for el in cfg.potential.single_elements:
            add(el)
    for cfg in configs:
        for gas in cfg.potential.gases:
            add(gas)
    for cfg in configs:
        for suite in cfg.potential.two_element_suites:
            A, B = suite.get("A"), suite.get("B")
            if A and B:
                add(f"{A}{B}")
    for cfg in configs:
        for suite in cfg.potential.multi_element_suites:
            composition = suite.get("composition", {})
            if isinstance(composition, dict) and len(composition) >= 2:
                add("".join(composition.keys()))
    for cfg in configs:
        for metal, gas in expand_gas_pairs(cfg.potential):
            add(f"{metal}{gas}")
    return systems


def _write_one_summary(
    configs: list[PotentialConfig],
    project_root: Path,
    output_csv: str | Path,
) -> Path:
    """Pivot every config in ``configs`` into one wide CSV.

    Parses each pot's logs directly (independent of ``<pot_name>.json`` / the
    automatic ``summarize-pot`` step — see the module docstring) so this
    always works off whatever's actually on disk in ``<pot_name>/``.
    """
    pot_names = [cfg.potential.pot_name for cfg in configs]
    metrics_by_pot: dict[str, dict[str, float | str]] = {}
    system_of: dict[str, str] = {}
    for cfg in configs:
        metrics, sys_of = _collect_pot_metrics(cfg, project_root)
        metrics_by_pot[cfg.potential.pot_name] = metrics
        system_of.update(sys_of)

    all_keys: set[str] = set()
    for m in metrics_by_pot.values():
        all_keys.update(m.keys())
    # elastic_raw_* (off-diagonal couplings) are dropped from the CSV — and
    # from <pot_name>.json too (build_pot_tag_json curates just like this
    # does). Only _collect_pot_metrics(cfg, run_dir) called directly still
    # has them.
    all_keys = {k for k in all_keys if not k.startswith("elastic_raw_")}

    def sort_key(k: str) -> tuple:
        if k.startswith("lc_"):
            return (0, k)
        if k.startswith("Ecoh_"):
            return (1, k)
        if k.startswith("E_per_atom_"):
            return (2, k)
        if k.startswith("Esingle_"):
            return (3, k)
        if k.startswith("elastic_"):
            # "elastic_<COMP>_<tag>" — group by tag, then Cij (C11, C22, ...,
            # C66) in canonical order, then Young's/Shear/Bulk/Poisson's
            # right after the Cij block for that same tag.
            rest = k[len("elastic_"):]
            comp, _, tag = rest.partition("_")
            return (4, tag, _ELASTIC_COMPONENT_ORDER.get(comp, 999), comp)
        if k.startswith("alloy_lc_"):
            return (5, k)
        if k.startswith("Ef_"):
            return (6, k)
        if k.startswith("Emig_"):
            return (7, k)
        if k.startswith("Eb_"):
            return (8, k)
        return (9, k)

    systems = _ordered_systems(configs)
    keys_by_system: dict[str, list[str]] = {s: [] for s in systems}
    other_keys: list[str] = []
    for k in all_keys:
        s = system_of.get(k)
        (keys_by_system[s] if s in keys_by_system else other_keys).append(k)
    for s in systems:
        keys_by_system[s].sort(key=sort_key)
    other_keys.sort(key=sort_key)

    out_path = project_root / output_csv
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["metric"] + pot_names)
        blank_row = [""] * len(pot_names)
        for system in [*systems, "other"] if other_keys else systems:
            keys = other_keys if system == "other" else keys_by_system[system]
            if not keys:
                continue
            writer.writerow([f"# {system}"] + blank_row)
            for key in keys:
                row = [key]
                for pot in pot_names:
                    v = metrics_by_pot.get(pot, {}).get(key, "")
                    row.append(v if v == "" else f"{v}")
                writer.writerow(row)
    return out_path


def write_wide_summary(
    pot_inputs_dir: str | Path = "pot_inputs",
    run_dir: str | Path = ".",
    output_csv: str | Path = "summary.csv",
) -> Path:
    """Aggregate every potential's logs into a single wide CSV.

    Parses ``<pot_name>/`` directly — doesn't need ``<pot_name>.json`` to exist. The
    output has metric labels as the first column and one column per
    potential (in the same order as ``pot_inputs/`` files load). Rows are
    grouped by system — pure element, gas, alloy pair, gas-complex pair, in
    that order — with a ``# <system>`` marker row ahead of each group; within
    a system, rows follow the usual lc/Ecoh/Esingle/elastic/alloy/Ef/Emig/Eb
    ordering.
    """
    project_root = Path(run_dir).resolve()
    configs = load_potential_configs(pot_inputs_dir)
    return _write_one_summary(configs, project_root, output_csv)


_REST_SENTINELS = {"rest", "the rest"}


def write_grouped_summaries(
    pot_inputs_dir: str | Path = "pot_inputs",
    run_dir: str | Path = ".",
    groups: dict[str, list[str] | str] | None = None,
) -> list[Path]:
    """Write one wide-format summary CSV per group in ``groups``.

    Each value is either a list of potential names — that CSV contains only
    those pots — or the sentinel ``"rest"`` (case-insensitive) — that CSV
    contains every potential not named in any other group's list. Output
    files are named ``<group_key>.csv``.

    ``groups`` empty/``None`` writes a single ``summary.csv`` covering every
    potential in ``pot_inputs_dir`` — the default, unfiltered behaviour.
    """
    project_root = Path(run_dir).resolve()
    configs = load_potential_configs(pot_inputs_dir)
    if not groups:
        return [_write_one_summary(configs, project_root, "summary.csv")]

    by_name = {cfg.potential.pot_name: cfg for cfg in configs}
    claimed: set[str] = set()
    explicit: dict[str, list[PotentialConfig]] = {}
    rest_keys: list[str] = []
    for name, members in groups.items():
        if isinstance(members, str) and members.strip().lower() in _REST_SENTINELS:
            rest_keys.append(name)
            continue
        selected: list[PotentialConfig] = []
        for pot_name in members:
            cfg = by_name.get(pot_name)
            if cfg is None:
                raise ValueError(f"Group '{name}' references unknown potential '{pot_name}'")
            selected.append(cfg)
            claimed.add(pot_name)
        explicit[name] = selected

    out_paths = [_write_one_summary(selected, project_root, f"{name}.csv")
                 for name, selected in explicit.items()]
    if rest_keys:
        rest_configs = [cfg for cfg in configs if cfg.potential.pot_name not in claimed]
        out_paths.extend(_write_one_summary(rest_configs, project_root, f"{name}.csv")
                          for name in rest_keys)
    return out_paths
