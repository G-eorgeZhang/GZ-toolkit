"""Summary generation for potential-testing workflows.

Two functions are exposed:

* ``write_summary(rows, output_csv)``     — long-form CSV used by the legacy
  ``summarize_stage`` helper. One row per (pot, case, metric).
* ``write_wide_summary(pot_inputs_dir, run_dir)`` — pivoted CSV with metric
  labels as rows and potentials as columns. This is what users get when they
  run ``python summarize.py``.

Wide-format layout::

    metric,                        pot1,    pot2,    pot3
    lc_Fe,                         2.8557,  2.8553,  2.8612
    lc_Cr,                         2.8841,  ...
    Ecoh_Fe,                       -4.122,  ...
    Ef_1vac_Fe,                    1.85,    ...
    Ef_2vac1nn_Fe,                 ...
    Emig_1vac_Fe,                  0.62,    ...
    alloy_lc_FeCr_Fe95Cr5,         2.859,   ...
    Ef_He1V1_Fe,                   ...
    Eb_He1V1_Fe,                   ...
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Any

from gz_toolkit.potential_testing.config import (
    PotentialConfig,
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
    "natoms": re.compile(r"RESULT\s+NATOMS\s+([0-9Ee+\-\.]+)"),
    "pe":     re.compile(r"RESULT\s+PE_TOTAL\s+([0-9Ee+\-\.]+)"),
    "lc":     re.compile(r"RESULT\s+LC\s+([0-9Ee+\-\.]+)"),
    "ecoh":   re.compile(r"RESULT\s+ECOH\s+([0-9Ee+\-\.]+)"),
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

# Components surfaced in summary.csv; remaining couplings go to summary.json only.
_MAIN_CIJ = ("C11", "C22", "C33", "C12", "C13", "C23", "C44", "C55", "C66")


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


def _voigt_moduli(cij: dict[str, float]) -> dict[str, float]:
    """Voigt-average bulk/shear modulus + Poisson ratio from the 9 main Cij."""
    if not all(k in cij for k in _MAIN_CIJ):
        return {}
    diag = cij["C11"] + cij["C22"] + cij["C33"]
    off = cij["C12"] + cij["C13"] + cij["C23"]
    shear = cij["C44"] + cij["C55"] + cij["C66"]
    B = (diag + 2.0 * off) / 9.0
    G = (diag - off + 3.0 * shear) / 15.0
    out = {"B": B, "G": G}
    denom = 2.0 * (3.0 * B + G)
    if denom != 0:
        out["nu"] = (3.0 * B - 2.0 * G) / denom
    return out


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


def _collect_pot_metrics(cfg: PotentialConfig, run_dir: Path) -> dict[str, float | str]:
    """Walk a single potential's workspace and return {metric_label: value}."""
    pot_dir = run_dir / cfg.potential.pot_name
    out: dict[str, float | str] = {}
    if not pot_dir.is_dir():
        return out

    elements = cfg.potential.single_elements or [next(iter(cfg.potential.type_map))]

    # Reference: lc, Ecoh, total PE, N_atoms.
    ref_total_pe: dict[str, float] = {}
    ref_natoms: dict[str, int] = {}
    for el in elements:
        log = pot_dir / "reference" / el / "log.lammps"
        res = _extract_results(log)
        if "lc" in res:
            out[f"lc_{el}"] = res["lc"]
        if "ecoh" in res:
            out[f"Ecoh_{el}"] = res["ecoh"]
        if "pe" in res:
            ref_total_pe[el] = res["pe"]
        if "natoms" in res:
            ref_natoms[el] = int(res["natoms"])

    # Point defects: formation energies.
    pd_root = pot_dir / "point_defects"
    ef_table: dict[tuple[str, str], float] = {}
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
                ef = e_def - (n_def / n_ref) * e_ref
                label = f"Ef_{case_dir.name}_{el}"
                out[label] = ef
                ef_table[(case_dir.name, el)] = ef

                # SEAKMC barrier.
                meb = case_dir / "meb"
                if meb.is_dir():
                    bar = _seakmc_min_barrier(meb)
                    if bar is not None:
                        out[f"Emig_{case_dir.name}_{el}"] = bar

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
            for comp in _MAIN_CIJ:
                if comp in cij:
                    out[f"elastic_{comp}_{tag}"] = cij[comp]
            for comp, val in cij.items():
                if comp not in _MAIN_CIJ:
                    # off-diagonal couplings: kept in summary.json only
                    out[f"elastic_raw_{comp}_{tag}"] = val
            for name, val in _voigt_moduli(cij).items():
                out[f"elastic_{name}_{tag}"] = val

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
                ef = res["pe"] - (int(res["natoms"]) / n_ref) * e_ref
                out[f"Ef_loop_{case_dir.name}_{el}"] = ef

    # Alloys: lattice constant.
    alloy_root = pot_dir / "alloy_lc"
    if alloy_root.is_dir():
        for pair_dir in sorted(alloy_root.iterdir()):
            for case_dir in sorted(pair_dir.iterdir()):
                res = _extract_results(case_dir / "log.lammps")
                if "lc" in res:
                    out[f"alloy_lc_{pair_dir.name}_{case_dir.name}"] = res["lc"]

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
                    # Use cohesive energy of metal for the metal contribution
                    # since N_metal_in_cell = N_ref. Gas chemical potential is
                    # taken as zero (classical pots) — user can post-process.
                    ef = e_def - e_ref  # ΔE for adding the gas atom (ignoring μ_He)
                    out[f"Ef_{gas}_in_{metal}_{single_match.group(1)}"] = ef
                    gas_ef[(metal, gas, 1, 0, single_match.group(1))] = ef
                elif m_match:
                    n_g = int(m_match.group(1))
                    kind = "vacancy" if m_match.group(2) == "V" else "interstitial"
                    m_d = int(m_match.group(3))
                    dM, _ = _gas_complex_atom_count(metal, gas, n_g, m_d, kind,
                                                    cfg.potential.type_map)
                    n_metal_in_cell = n_def - n_g
                    ef = e_def - (n_metal_in_cell / n_ref) * e_ref
                    label = f"Ef_{gas}{n_g}_{m_match.group(2)}{m_d}_{metal}"
                    out[label] = ef
                    gas_ef[(metal, gas, n_g, m_d, kind)] = ef

        # Binding energies: Eb_He^V = E_f(HenVm) - E_f(Hen V_{m-1}) - E_f(V)
        for (metal, gas, n_g, m_d, kind), ef in gas_ef.items():
            if m_d == 0 or n_g == 0:
                continue
            vi = "V" if kind == "vacancy" else "I"
            # Vacancy/interstitial single-defect Ef (need to find the right key)
            single_def_key = f"Ef_1vac_{metal}" if kind == "vacancy" else f"Ef_1int_octa_{metal}"
            if single_def_key not in out:
                continue
            ef_single_defect = out[single_def_key]
            # Binding via removing one defect atom (m -> m-1 vacancies/interstitials).
            ef_minus = gas_ef.get((metal, gas, n_g, m_d - 1, kind))
            if ef_minus is not None and isinstance(ef_minus, float):
                eb_v = ef - ef_minus - ef_single_defect
                out[f"Eb_{vi}_{gas}{n_g}_{vi}{m_d}_{metal}"] = eb_v
            # Binding via removing one gas atom (n -> n-1 gas).
            ef_minus_g = gas_ef.get((metal, gas, n_g - 1, m_d, kind))
            single_gas_tetra = gas_ef.get((metal, gas, 1, 0, "tetra"))
            single_gas_octa = gas_ef.get((metal, gas, 1, 0, "octa"))
            ef_single_gas = single_gas_tetra if single_gas_tetra is not None else single_gas_octa
            if ef_minus_g is not None and ef_single_gas is not None:
                eb_g = ef - ef_minus_g - ef_single_gas
                out[f"Eb_{gas}_{gas}{n_g}_{vi}{m_d}_{metal}"] = eb_g

    return out


# ---------------------------------------------------------------------------
# Wide-format aggregator
# ---------------------------------------------------------------------------


def write_wide_summary(
    pot_inputs_dir: str | Path = "pot_inputs",
    run_dir: str | Path = ".",
    output_csv: str | Path = "summary.csv",
) -> Path:
    """Aggregate every potential's metrics into a single wide CSV.

    The output has metric labels as the first column and one column per
    potential (in the same order as ``pot_inputs/`` files load).
    """
    project_root = Path(run_dir).resolve()
    configs = load_potential_configs(pot_inputs_dir)

    pot_names = [cfg.potential.pot_name for cfg in configs]
    metrics_by_pot: dict[str, dict[str, float | str]] = {}
    for cfg in configs:
        metrics_by_pot[cfg.potential.pot_name] = _collect_pot_metrics(cfg, project_root)

    # Union of metric keys, ordered: lc_*, Ecoh_*, Ef_*, Emig_*, alloy_*, Eb_*, others.
    all_keys: set[str] = set()
    for m in metrics_by_pot.values():
        all_keys.update(m.keys())

    def sort_key(k: str) -> tuple[int, str]:
        if k.startswith("lc_"):
            return (0, k)
        if k.startswith("Ecoh_"):
            return (1, k)
        if k.startswith("elastic_"):
            return (2, k)
        if k.startswith("alloy_lc_"):
            return (3, k)
        if k.startswith("Ef_"):
            return (4, k)
        if k.startswith("Emig_"):
            return (5, k)
        if k.startswith("Eb_"):
            return (6, k)
        return (9, k)

    # elastic_raw_* (off-diagonal couplings) stay in summary.json only.
    ordered_keys = sorted(
        (k for k in all_keys if not k.startswith("elastic_raw_")),
        key=sort_key,
    )
    out_path = project_root / output_csv
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["metric"] + pot_names)
        for key in ordered_keys:
            row = [key]
            for pot in pot_names:
                v = metrics_by_pot.get(pot, {}).get(key, "")
                row.append(v if v == "" else f"{v}")
            writer.writerow(row)

    # Also dump raw per-pot metrics as JSON for downstream scripts.
    (project_root / "summary.json").write_text(
        json.dumps(metrics_by_pot, indent=2),
        encoding="utf-8",
    )
    return out_path
