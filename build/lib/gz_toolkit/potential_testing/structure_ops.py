"""Structure and defect generation helpers.

These helpers all operate on already-built supercells (or build them) and
write a `structure.data` LAMMPS file into the supplied case directory.

Two distinct entry points:

* ``build_reference_structure(...)``        — builds an N×N×N pure-element supercell.
* ``build_case_from_reference(...)``        — applies one defect to a copy of the
                                              relaxed reference and writes the result.
* ``build_alloy_structure(...)``            — random-substitution alloy supercell.
* ``build_gas_complex_structure(...)``      — (gas)n(V)m or (gas)n(I)m complex.

All functions accept a `relaxed_lc` value that should come from the reference
LAMMPS box-relax output. If unavailable, callers fall back to ``lc_initial``.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np

from gz_toolkit.buildmtx.buildstr import Gen_crystal
from gz_toolkit.core.modlmp import Modlmp_LmpData
from gz_toolkit.defect.bcc_defect import BCCDefect
from gz_toolkit.potential_testing.config import PotentialConfig


# ---------------------------------------------------------------------------
# Prototype helpers
# ---------------------------------------------------------------------------


def _prototype_from_structure(name: str) -> str:
    s = name.lower()
    if s == "bcc":
        return "A2"
    if s == "fcc":
        return "A1"
    if s == "hcp":
        return "A3"
    # Ceramics / non-elemental prototypes are not supported in this version.
    raise NotImplementedError(
        f"Crystal structure '{name}' is not supported by potential-testing yet. "
        "Currently only bcc/fcc/hcp pure-element cells are handled."
    )


def _basis_count(proto: str) -> int:
    return {"A1": 4, "A2": 2, "A3": 2}.get(proto, 2)


# ---------------------------------------------------------------------------
# Atom geometry helpers
# ---------------------------------------------------------------------------


def _nearest_index(data: Modlmp_LmpData, center: np.ndarray) -> int:
    xyz = data.atoms[["x", "y", "z"]].to_numpy(dtype=float)
    d = np.linalg.norm(xyz - center.reshape(1, 3), axis=1)
    return int(np.argmin(d))


def _delete_nearest_n(data: Modlmp_LmpData, center: np.ndarray, n: int) -> int:
    """Delete the n atoms closest to ``center``."""
    xyz = data.atoms[["x", "y", "z"]].to_numpy(dtype=float)
    d = np.linalg.norm(xyz - center.reshape(1, 3), axis=1)
    order = np.argsort(d)[: max(0, n)]
    ids = data.atoms.index.to_numpy()[order]
    data.atoms = data.atoms.drop(index=ids).copy()
    data.initialization(normalization=False, style=1)
    return int(len(ids))


def _delete_2vac_nnn(data: Modlmp_LmpData, center: np.ndarray, n_shell: int, tol: float = 1e-3) -> int:
    """Delete a divacancy: one atom near ``center`` plus its n-th-nearest neighbour."""
    xyz = data.atoms[["x", "y", "z"]].to_numpy(dtype=float)
    ids = data.atoms.index.to_numpy()
    if len(ids) < 2:
        return 0

    d0 = np.linalg.norm(xyz - center.reshape(1, 3), axis=1)
    i0 = int(np.argmin(d0))
    p0 = xyz[i0]
    id0 = ids[i0]

    keep = np.ones(len(ids), dtype=bool)
    keep[i0] = False
    xyz_r = xyz[keep]
    ids_r = ids[keep]
    d = np.linalg.norm(xyz_r - p0.reshape(1, 3), axis=1)
    if len(d) == 0:
        data.atoms = data.atoms.drop(index=[id0]).copy()
        data.initialization(normalization=False, style=1)
        return 1

    order = np.argsort(d)
    d_sorted = d[order]
    ids_sorted = ids_r[order]

    shells: list[list[int]] = []
    shell_d: list[float] = []
    for dist, aid in zip(d_sorted, ids_sorted):
        if not shell_d or abs(dist - shell_d[-1]) > tol:
            shell_d.append(float(dist))
            shells.append([int(aid)])
        else:
            shells[-1].append(int(aid))

    target = max(1, int(n_shell))
    if target <= len(shells):
        id1 = shells[target - 1][0]
    else:
        id1 = int(ids_sorted[0])

    data.atoms = data.atoms.drop(index=[int(id0), int(id1)]).copy()
    data.initialization(normalization=False, style=1)
    return 2


def _assert_full_force_field(data: Modlmp_LmpData, config: PotentialConfig) -> None:
    """Register every type in type_map on a freshly loaded data object.

    Reference/seed data files only declare the types present in them (usually
    one). Adding or re-typing atoms to a higher type id would then break
    mylammps' get_data_info() — registering the full force field first keeps
    the Masses table consistent with the potential.
    """
    ff_elements = list(config.potential.type_map.keys())
    atomic_masses = [config.potential.masses[config.potential.type_map[e]] for e in ff_elements]
    data.assert_force_field(ff_elements, atomic_masses=atomic_masses)


def _add_interstitial_cart(data: Modlmp_LmpData, center: np.ndarray, offset: np.ndarray, atom_type: int) -> None:
    data.add_atoms([center + offset], atom_type=atom_type, wrap=True, reinit=True)


def _make_dumbbell(data: Modlmp_LmpData, center: np.ndarray, direction: np.ndarray, atom_type: int, sep: float) -> None:
    """Replace the atom nearest ``center`` with two atoms separated along ``direction``."""
    idx = _nearest_index(data, center)
    atom_id = data.atoms.index.to_numpy()[idx]
    p = data.atoms.loc[atom_id, ["x", "y", "z"]].to_numpy(dtype=float)
    u = np.asarray(direction, dtype=float).reshape(3,)
    u = u / np.linalg.norm(u)
    p1 = p - 0.5 * sep * u
    p2 = p + 0.5 * sep * u
    data.atoms.loc[atom_id, ["x", "y", "z"]] = p1
    data.add_atoms([p2], atom_type=atom_type, wrap=True, reinit=False)
    data.initialization(normalization=False, style=1)


# ---------------------------------------------------------------------------
# Crystal-specific interstitial sites
# ---------------------------------------------------------------------------


def _interstitial_offset(structure: str, kind: str, lc: float) -> np.ndarray:
    """Return the cartesian offset (relative to a lattice atom) for an interstitial site.

    BCC tetrahedral: (a/4, a/2, 0)
    BCC octahedral : (a/2, a/2, 0)
    FCC tetrahedral: (a/4, a/4, a/4)
    FCC octahedral : (a/2, 0,   0)
    """
    s = structure.lower()
    k = kind.lower()
    if s == "bcc" and k == "tetra":
        return np.array([0.25, 0.5, 0.0]) * lc
    if s == "bcc" and k == "octa":
        return np.array([0.5, 0.5, 0.0]) * lc
    if s == "fcc" and k == "tetra":
        return np.array([0.25, 0.25, 0.25]) * lc
    if s == "fcc" and k == "octa":
        return np.array([0.5, 0.0, 0.0]) * lc
    # Fallback — unsupported combo, use a generic tetrahedral-ish offset.
    return np.array([0.25, 0.25, 0.0]) * lc


# ---------------------------------------------------------------------------
# Point-defect dispatcher
# ---------------------------------------------------------------------------


def _apply_point_defect_case(
    data: Modlmp_LmpData,
    case_name: str,
    structure_name: str,
    lc: float,
    atom_type: int,
    dumbbell_sep_factor: float = 0.35,
) -> str:
    c = np.asarray(data.get_center("cart"), dtype=float).reshape(3,)
    case = case_name.lower()
    dumbbell_sep = float(dumbbell_sep_factor) * lc

    if case == "1vac":
        _delete_nearest_n(data, c, 1)
        return "vacancy_1"
    if case == "2vac1nn":
        _delete_2vac_nnn(data, c, 1)
        return "2vac1nn"
    if case == "2vac2nn":
        _delete_2vac_nnn(data, c, 2)
        return "2vac2nn"
    if case == "2vac3nn":
        _delete_2vac_nnn(data, c, 3)
        return "2vac3nn"
    if case == "2vac4nn":
        _delete_2vac_nnn(data, c, 4)
        return "2vac4nn"
    if case == "1int_tetra":
        _add_interstitial_cart(data, c, _interstitial_offset(structure_name, "tetra", lc), atom_type)
        return "int_tetra"
    if case == "1int_octa":
        _add_interstitial_cart(data, c, _interstitial_offset(structure_name, "octa", lc), atom_type)
        return "int_octa"
    if case == "dumbbell_100":
        _make_dumbbell(data, c, np.array([1, 0, 0]), atom_type, dumbbell_sep)
        return "dumbbell_100"
    if case == "dumbbell_111":
        _make_dumbbell(data, c, np.array([1, 1, 1]), atom_type, dumbbell_sep)
        return "dumbbell_111"
    if case == "dumbbell_110":
        _make_dumbbell(data, c, np.array([1, 1, 0]), atom_type, dumbbell_sep)
        return "dumbbell_110"
    if case == "2int":
        _add_interstitial_cart(data, c, _interstitial_offset(structure_name, "tetra", lc), atom_type)
        _add_interstitial_cart(data, c + np.array([0.5, 0.0, 0.0]) * lc,
                               _interstitial_offset(structure_name, "octa", lc), atom_type)
        return "int_2"
    if case == "3int":
        _add_interstitial_cart(data, c, _interstitial_offset(structure_name, "tetra", lc), atom_type)
        _add_interstitial_cart(data, c + np.array([0.5, 0.0, 0.0]) * lc,
                               _interstitial_offset(structure_name, "octa", lc), atom_type)
        _add_interstitial_cart(data, c + np.array([0.0, 0.5, 0.0]) * lc,
                               _interstitial_offset(structure_name, "tetra", lc), atom_type)
        return "int_3"
    return "no_point_defect"


# ---------------------------------------------------------------------------
# Loops (BCC only) and dislocation-line placeholders
# ---------------------------------------------------------------------------


def _apply_loop_case(
    data: Modlmp_LmpData,
    case_name: str,
    structure_name: str,
    config: PotentialConfig,
    element: str,
    lc: float,
) -> str:
    if structure_name.lower() != "bcc":
        return "skipped_non_bcc"

    ff_elements = list(config.potential.type_map.keys())
    masses = [config.potential.masses[config.potential.type_map[e]] for e in ff_elements]
    c = np.asarray(data.get_center("cart"), dtype=float).reshape(3,)
    defects = BCCDefect(data)

    loop_radius = float(config.workflow.loop_radius_factor) * lc

    if case_name == "SIL111":
        defects.add_111_loop(
            loop_type="sil",
            radius=loop_radius,
            lattice_const=lc,
            ff_elements=ff_elements,
            atomic_masses=masses,
            center=c,
            habit_plane=(1, 1, 1),
            n_repeats=1,
        )
        return "loop_111_sia"

    if case_name == "SIL100":
        defects.add_100_loop(
            loop_type="sil",
            radius=loop_radius,
            lattice_const=lc,
            ff_elements=ff_elements,
            atomic_masses=masses,
            center=c,
            habit_plane=(1, 0, 0),
            n_repeats=1,
        )
        return "loop_100_sia"

    return "no_loop"


def _apply_dislocation_case(case_name: str) -> str:
    """Dislocation-line constructions are deferred — see README/TODO."""
    raise NotImplementedError(
        f"Dislocation-line case '{case_name}' is not implemented yet. "
        "Edge / screw lines need a multi-step build (cut & shift, anisotropic "
        "pre-displacement, then box-aware relax) that hasn't been wired in. "
        "Disable include_dislocation_lines until this lands."
    )


# ---------------------------------------------------------------------------
# Public builders
# ---------------------------------------------------------------------------


def build_reference_structure(
    reference_dir: Path,
    config: PotentialConfig,
    element: str,
    structure_name: str,
    lc_value: float | None = None,
    n_replicate: int | None = None,
) -> Path:
    """Build an N×N×N pure-element supercell at ``lc_value`` (or ``lc_initial``).

    ``n_replicate`` overrides ``size_single`` (used e.g. by the smaller
    elastic-constant cells).
    """
    proto = _prototype_from_structure(structure_name)
    seed_file = reference_dir / "seed.data"
    out_data = reference_dir / "structure.data"
    gen = Gen_crystal(atom_style="atomic")
    gen.seed_crystal(proto, [element] * _basis_count(proto), str(seed_file))
    n = int(n_replicate if n_replicate is not None else config.potential.size_single)
    lc = float(lc_value if lc_value is not None else config.potential.lc_initial)
    gen.replicate(str(seed_file), n, n, n, lc, str(out_data))
    return out_data


def build_case_from_reference(
    reference_data: Path,
    case_dir: Path,
    case_name: str,
    config: PotentialConfig,
    element: str,
    structure_name: str,
    relaxed_lc: float | None = None,
) -> str:
    """Apply a single defect to a copy of the relaxed reference cell."""
    out_data = case_dir / "structure.data"
    data = Modlmp_LmpData.from_file(str(reference_data), "atomic")
    _assert_full_force_field(data, config)
    atom_type = int(config.potential.type_map[element])
    lc = float(relaxed_lc if relaxed_lc is not None else config.potential.lc_initial)

    if case_name in {"SIL111", "SIL100"}:
        status = _apply_loop_case(data, case_name, structure_name, config, element, lc)
    elif case_name in {"edgedislo111", "edgedislo100", "screw111"}:
        status = _apply_dislocation_case(case_name)  # raises NotImplementedError
    else:
        status = _apply_point_defect_case(
            data, case_name, structure_name, lc, atom_type,
            dumbbell_sep_factor=config.workflow.dumbbell_sep_factor,
        )

    data.to_file(str(out_data))
    return status


def build_alloy_structure(
    case_dir: Path,
    config: PotentialConfig,
    element_A: str,
    element_B: str,
    fraction_B_atpct: float,
    structure_name: str,
    ordering: str = "random",
    relaxed_lc_A: float | None = None,
    seed: int | None = 0,
    n_replicate: int | None = None,
) -> str:
    """Build a binary-alloy supercell.

    ``ordering`` may include "random" and/or "B2" separated by ``|``; the first
    matching ordering is used. B2 only makes sense at 50 at.% — at other
    fractions B2 falls back to random with a warning baked into the status.
    ``n_replicate`` overrides ``size_alloy`` (used by the elastic cells).
    """
    proto = _prototype_from_structure(structure_name)
    seed_file = case_dir / "seed.data"
    out_data = case_dir / "structure.data"
    gen = Gen_crystal(atom_style="atomic")
    gen.seed_crystal(proto, [element_A] * _basis_count(proto), str(seed_file))
    n = int(n_replicate if n_replicate is not None else config.potential.size_alloy)
    lc = float(relaxed_lc_A if relaxed_lc_A is not None else config.potential.lc_initial)
    gen.replicate(str(seed_file), n, n, n, lc, str(out_data))

    data = Modlmp_LmpData.from_file(str(out_data), "atomic")
    _assert_full_force_field(data, config)
    type_A = int(config.potential.type_map[element_A])
    type_B = int(config.potential.type_map[element_B])

    requested = [o.strip().lower() for o in ordering.split("|") if o.strip()]
    use_b2 = "b2" in requested and abs(fraction_B_atpct - 50.0) < 1e-6 and proto == "A2"

    n_atoms = len(data.atoms)
    if use_b2:
        # B2: assign one sublattice (the (1/2,1/2,1/2) basis atoms) to type B.
        # Replicated seed has basis atoms interleaved: even index -> corner (A),
        # odd index -> body-centre (B).
        idx = np.arange(n_atoms)
        b_mask = (idx % 2) == 1
        data.atoms.loc[data.atoms.index[b_mask], "type"] = type_B
        status = f"alloy_B2_{element_A}{50}_{element_B}{50}"
    else:
        rng = np.random.default_rng(seed)
        n_B = int(round(n_atoms * fraction_B_atpct / 100.0))
        choice = rng.choice(n_atoms, size=n_B, replace=False)
        ids = data.atoms.index.to_numpy()
        data.atoms.loc[ids[choice], "type"] = type_B
        f_int = int(round(fraction_B_atpct))
        status = f"alloy_random_{element_A}{100 - f_int}_{element_B}{f_int}"

    # Force types A to be set explicitly so the data file's mass section is consistent.
    # (Atoms not picked above retain type_A from the seed crystal.)
    data.initialization(normalization=False, style=1)
    data.to_file(str(out_data))
    return status


def build_gas_complex_structure(
    reference_data: Path,
    case_dir: Path,
    metal: str,
    gas: str,
    n_gas: int,
    m_defect: int,
    defect_kind: str,
    config: PotentialConfig,
    structure_name: str,
    relaxed_lc: float | None = None,
) -> str:
    """Build (gas)_n(V)_m or (gas)_n(I)_m near the cell centre.

    ``defect_kind`` is "vacancy" or "interstitial". Vacancies are the m nearest
    atoms to the centre; interstitials are placed at successive offsets from
    a tetrahedral site around the centre. Gas atoms are then placed close to
    the defect cluster — the LAMMPS minimisation in the case dir resolves them
    into the lowest-energy configuration.
    """
    out_data = case_dir / "structure.data"
    data = Modlmp_LmpData.from_file(str(reference_data), "atomic")
    _assert_full_force_field(data, config)
    metal_type = int(config.potential.type_map[metal])
    gas_type = int(config.potential.type_map[gas])
    lc = float(relaxed_lc if relaxed_lc is not None else config.potential.lc_initial)
    c = np.asarray(data.get_center("cart"), dtype=float).reshape(3,)

    # --- defect side ---
    if defect_kind == "vacancy":
        if m_defect > 0:
            _delete_nearest_n(data, c, m_defect)
    elif defect_kind == "interstitial":
        # Place m_defect metal interstitials around the centre on a small lattice.
        offsets = [
            np.array([0.5, 0.0, 0.0]),
            np.array([0.0, 0.5, 0.0]),
            np.array([0.0, 0.0, 0.5]),
            np.array([0.5, 0.5, 0.0]),
            np.array([0.5, 0.0, 0.5]),
        ]
        for k in range(m_defect):
            off = offsets[k % len(offsets)] * lc
            _add_interstitial_cart(data, c, off, metal_type)
    else:
        raise ValueError(f"defect_kind must be 'vacancy' or 'interstitial', got {defect_kind!r}")

    # --- gas atoms: cluster near the defect centre ---
    # Generate a small spiral/grid of offsets so atoms don't overlap exactly
    # but are close enough that minimisation pulls them into the defect.
    gas_offsets = [
        np.array([0.0, 0.0, 0.0]),
        np.array([0.25, 0.0, 0.0]),
        np.array([0.0, 0.25, 0.0]),
        np.array([0.0, 0.0, 0.25]),
        np.array([0.25, 0.25, 0.0]),
        np.array([0.25, 0.0, 0.25]),
        np.array([0.0, 0.25, 0.25]),
        np.array([-0.25, 0.0, 0.0]),
        np.array([0.0, -0.25, 0.0]),
    ]
    for k in range(n_gas):
        off = gas_offsets[k % len(gas_offsets)] * lc
        # Tiny irrational shift to avoid exact overlap when k > len(gas_offsets).
        off = off + np.array([1e-3, 1e-3, 1e-3]) * (k // len(gas_offsets))
        _add_interstitial_cart(data, c, off, gas_type)

    data.initialization(normalization=False, style=1)
    data.to_file(str(out_data))
    return f"{gas}{n_gas}_{'V' if defect_kind == 'vacancy' else 'I'}{m_defect}"


def build_single_gas_in_bulk(
    reference_data: Path,
    case_dir: Path,
    metal: str,
    gas: str,
    site: str,
    config: PotentialConfig,
    structure_name: str,
    relaxed_lc: float | None = None,
) -> str:
    """One gas atom on a single tetrahedral or octahedral site — needed for E_He^f."""
    out_data = case_dir / "structure.data"
    data = Modlmp_LmpData.from_file(str(reference_data), "atomic")
    _assert_full_force_field(data, config)
    gas_type = int(config.potential.type_map[gas])
    lc = float(relaxed_lc if relaxed_lc is not None else config.potential.lc_initial)
    c = np.asarray(data.get_center("cart"), dtype=float).reshape(3,)
    off = _interstitial_offset(structure_name, site, lc)
    _add_interstitial_cart(data, c, off, gas_type)
    data.to_file(str(out_data))
    return f"{gas}_in_{metal}_{site}"
