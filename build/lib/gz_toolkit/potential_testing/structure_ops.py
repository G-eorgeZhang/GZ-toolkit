"""Structure and defect generation helpers.

Every structure below is built the same way, regardless of whether it's a
pure element or an N-element alloy:

  1. Seed + replicate at lattice parameter **1** (dimensionless) -- every atom
     starts out typed as the "host" (first) element of the composition. All
     the geometry helpers below (interstitial offsets, dumbbell separation,
     loop radius) are already expressed as ``fraction * lc``, so building at
     ``lc=1`` and passing ``lc=1.0`` to them just works.
  2. Apply the case-specific edit (vacancy deletion, interstitial/dumbbell/
     loop insertion) in that dimensionless frame. Anything inserted is typed
     as the placeholder (host) type -- deciding its real chemistry is
     deferred to step 4.
  3. Scale the whole cell (box + atoms, affinely) up to the real lattice
     constant.
  4. Only now decide chemistry: for an alloy, every placeholder-typed atom
     (i.e. everyone except any deliberately-typed gas atom) is randomly
     reassigned according to the target composition. For a pure element this
     is a no-op -- every atom is already the only element in play.

This way, an inserted atom (dumbbell partner, interstitial, loop atom) is
drawn from the same random distribution as every other atom -- there's never
a point where the code has to decide "what type should this new atom be"
ahead of the composition being applied.

Entry points:

* ``build_reference_structure(...)``   -- N x N x N pure-element supercell.
* ``build_single_atom_structure(...)`` -- one isolated atom (sanity check).
* ``build_case_structure(...)``        -- one point-defect/loop/bulk case,
                                          pure element or alloy.
* ``build_alloy_structure(...)``       -- binary-alloy bulk supercell.
* ``build_multi_alloy_structure(...)`` -- N-element (N>=2) bulk supercell.
* ``build_gas_complex_structure(...)`` -- (gas)n(V)m or (gas)n(I)m complex.

All functions accept a `relaxed_lc`/`target_lc` value that should come from
the reference LAMMPS box-relax output. If unavailable, callers fall back to
``lc_initial``.
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
        # A3_ORTHO, not A3: the primitive hexagonal cell has gamma = 120 deg,
        # which cannot be replicated into an untilted LAMMPS box. The
        # orthogonal 4-atom cell is the same crystal in a rectangular box.
        return "A3_ORTHO"
    # Ceramics / non-elemental prototypes are not supported in this version.
    raise NotImplementedError(
        f"Crystal structure '{name}' is not supported by potential-testing yet. "
        "Currently only bcc/fcc/hcp pure-element cells are handled."
    )


def _basis_count(proto: str) -> int:
    return {"A1": 4, "A2": 2, "A3": 2, "A3_ORTHO": 4}.get(proto, 2)


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


def _write_structure(data: Modlmp_LmpData, out_data: Path) -> None:
    """Reset atom IDs and write a `structure.data`.

    Defect/alloy edits (deletions, insertions, type reassignment) leave gaps
    or out-of-order atom IDs; LAMMPS data files expect a dense 1..N range.
    """
    data.reset_atom_ids()
    data.to_file(str(out_data))


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
        # loop_atom_type left as the default (None): the new atoms copy the
        # type of the slab atom they're duplicated from, which is still the
        # placeholder type at this point in the build -- see module docstring.
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
# Shared build recipe: seed@lc=1 -> edit -> scale -> decorate
# ---------------------------------------------------------------------------


def _build_unit_cell(
    case_dir: Path,
    config: PotentialConfig,
    structure_name: str,
    host_element: str,
    n_replicate: int,
) -> Modlmp_LmpData:
    """Seed + replicate at lattice parameter 1 (dimensionless unit cell).

    Every atom starts out typed as ``host_element`` -- the placeholder type
    that a later ``_decorate_alloy_types`` call reassigns. Building at a=1
    means every defect/loop offset (already expressed as ``fraction * lc``)
    can be computed with ``lc=1.0`` and stays correct once the cell is
    scaled up afterward.
    """
    proto = _prototype_from_structure(structure_name)
    seed_file = case_dir / "seed.data"
    unit_data = case_dir / "unit_cell.data"
    gen = Gen_crystal(atom_style="atomic")
    gen.seed_crystal(proto, [host_element] * _basis_count(proto), str(seed_file))
    n = int(n_replicate)
    gen.replicate(str(seed_file), n, n, n, 1.0, str(unit_data))

    data = Modlmp_LmpData.from_file(str(unit_data), "atomic")
    # Gen_crystal/pymatgen always numbers a single-species cell as type 1,
    # regardless of what atom-type id `host_element` actually has in
    # type_map — force every atom to the real type id before registering the
    # full force field.
    data.atoms["type"] = int(config.potential.type_map[host_element])
    _assert_full_force_field(data, config)
    return data


def _scale_cell_to_lc(data: Modlmp_LmpData, lc: float) -> None:
    """Scale a unit cell (built at a0 = 1) up to the physical lattice constant.

    ``_build_unit_cell`` replicates a seed whose *a* axis is 1, so each box
    edge is already the right multiple of a0 **for that prototype** -- ``n`` on
    every axis for cubic, but ``(n, sqrt(3) n, (c/a) n)`` for the orthogonal
    hcp cell. The rescale is therefore *multiplicative*: every axis is
    stretched by the same factor ``lc``, which reproduces ``n * lc`` for cubic
    while preserving the cell's shape (and hence c/a) for any prototype.

    Do NOT set the axes to an absolute ``n * lc``: that is only correct when
    a == b == c, and would squash an hcp cell into a cube.

    Atoms are carried affinely with the box.
    """
    factor = float(lc)
    if factor <= 0.0:
        raise ValueError(f"Lattice constant must be positive, got {lc}.")
    lengths = np.asarray(data.box.lengths, dtype=float).reshape(3)
    data.scale_box(
        axes=[0, 1, 2],
        values=(lengths * factor).tolist(),
        mode="lc",
        affine=True,
    )


def _decorate_alloy_types(
    data: Modlmp_LmpData,
    composition_atpct: dict[str, float],
    type_map: dict[str, int],
    placeholder_type: int,
    rng: np.random.Generator,
) -> None:
    """Randomly assign real chemistry to every placeholder-typed atom.

    Operates only on ``mask = (type == placeholder_type)`` -- any atom
    already given a final, real type at insert time (gas atoms) is
    untouched, since it was never placeholder-typed to begin with. Non-host
    elements are allocated largest-fraction-first so rounding remainders
    land on the majority species; the host (composition's first element,
    which is what ``placeholder_type`` already is) absorbs whatever is left.

    For a pure element (``composition_atpct`` has one entry), every eligible
    atom is already that element's real type, so this is a no-op.
    """
    elements = list(composition_atpct.keys())
    host = elements[0]
    mask = (data.atoms["type"] == placeholder_type).to_numpy()
    idx = data.atoms.index.to_numpy()[mask]
    n_eligible = len(idx)
    if n_eligible == 0:
        return
    pool = idx[rng.permutation(n_eligible)]

    others = [(el, composition_atpct[el]) for el in elements if el != host]
    others.sort(key=lambda kv: kv[1], reverse=True)

    consumed = 0
    for el, pct in others:
        n_el = min(int(round(n_eligible * pct / 100.0)), len(pool) - consumed)
        chosen = pool[consumed: consumed + n_el]
        consumed += n_el
        data.atoms.loc[chosen, "type"] = int(type_map[el])
    # host keeps whatever's left over -- already placeholder_type.


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
    out_data = reference_dir / "structure.data"
    n = int(n_replicate if n_replicate is not None else config.potential.size_single)
    lc = float(lc_value if lc_value is not None else config.potential.lc_initial)

    data = _build_unit_cell(reference_dir, config, structure_name, element, n)
    _scale_cell_to_lc(data, lc)
    _write_structure(data, out_data)
    return out_data


def build_single_atom_structure(
    case_dir: Path,
    config: PotentialConfig,
    element: str,
    box_length: float = 100.0,
) -> Path:
    """One isolated atom in a large periodic box.

    Sanity check on the potential itself: an atom with no neighbors within
    cutoff should read a well-defined (often ~0, but potential-dependent)
    per-atom energy. ``box_length`` (100 A) is far beyond any realistic
    interatomic cutoff, so the atom never sees its own periodic images.
    """
    from pymatgen.core import Lattice, Structure
    from pymatgen.io.lammps.data import LammpsData

    out_data = case_dir / "structure.data"
    struct = Structure(Lattice.cubic(box_length), [element], [[0.5, 0.5, 0.5]])
    LammpsData.from_structure(struct, atom_style="atomic").write_file(str(out_data))

    data = Modlmp_LmpData.from_file(str(out_data), "atomic")
    data.atoms["type"] = int(config.potential.type_map[element])
    _assert_full_force_field(data, config)
    _write_structure(data, out_data)
    return out_data


def build_case_structure(
    case_dir: Path,
    config: PotentialConfig,
    composition_atpct: dict[str, float],
    case_name: str,
    structure_name: str,
    target_lc: float,
    n_replicate: int,
    seed: int | None = None,
) -> str:
    """Build one point-defect/loop/bulk case, for a pure element or an alloy.

    ``case_name == "bulk"`` builds the undisturbed cell (used as the Ef
    reference for a replica); any other name is dispatched through
    ``_apply_point_defect_case``/``_apply_loop_case``. See the module
    docstring for the full seed@lc=1 -> edit -> scale -> decorate recipe.

    For a pure element, pass ``composition_atpct = {element: 100.0}`` --
    decoration is then a no-op. For an alloy, every call is an independent
    random draw (default ``seed=None``): alloy defect energies are meant to
    be averaged over several such replicas, not reproduced exactly.
    """
    elements = list(composition_atpct.keys())
    if not elements:
        raise ValueError("composition_atpct must have at least one element.")
    host = elements[0]
    placeholder_type = int(config.potential.type_map[host])
    out_data = case_dir / "structure.data"

    data = _build_unit_cell(case_dir, config, structure_name, host, n_replicate)

    if case_name == "bulk":
        status = "bulk"
    elif case_name in {"SIL111", "SIL100"}:
        status = _apply_loop_case(data, case_name, structure_name, config, host, 1.0)
    elif case_name in {"edgedislo111", "edgedislo100", "screw111"}:
        status = _apply_dislocation_case(case_name)  # raises NotImplementedError
    else:
        status = _apply_point_defect_case(
            data, case_name, structure_name, 1.0, placeholder_type,
            dumbbell_sep_factor=config.workflow.dumbbell_sep_factor,
        )

    _scale_cell_to_lc(data, target_lc)
    if len(elements) > 1:
        rng = np.random.default_rng(seed)
        _decorate_alloy_types(data, composition_atpct, config.potential.type_map, placeholder_type, rng)

    data.initialization(normalization=False, style=1)
    _write_structure(data, out_data)
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
    seed: int | None = None,
    n_replicate: int | None = None,
) -> str:
    """Build a binary-alloy supercell (bulk, no defect).

    ``ordering`` may include "random" and/or "B2" separated by ``|``; the first
    matching ordering is used. B2 only makes sense at 50 at.% — at other
    fractions B2 falls back to random with a warning baked into the status.
    ``n_replicate`` overrides ``size_alloy`` (used by the elastic cells).
    ``seed`` defaults to ``None`` (fresh random draw every call).
    """
    n = int(n_replicate if n_replicate is not None else config.potential.size_alloy)
    lc = float(relaxed_lc_A if relaxed_lc_A is not None else config.potential.lc_initial)
    out_data = case_dir / "structure.data"

    requested = [o.strip().lower() for o in ordering.split("|") if o.strip()]
    use_b2 = (
        "b2" in requested
        and abs(fraction_B_atpct - 50.0) < 1e-6
        and _prototype_from_structure(structure_name) == "A2"
    )

    data = _build_unit_cell(case_dir, config, structure_name, element_A, n)

    if use_b2:
        # B2: assign one sublattice (the (1/2,1/2,1/2) basis atoms) to type B.
        # Replicated seed has basis atoms interleaved: even index -> corner (A),
        # odd index -> body-centre (B).
        type_B = int(config.potential.type_map[element_B])
        n_atoms = len(data.atoms)
        idx = np.arange(n_atoms)
        b_mask = (idx % 2) == 1
        data.atoms.loc[data.atoms.index[b_mask], "type"] = type_B
        status = f"alloy_B2_{element_A}{50}_{element_B}{50}"
    else:
        rng = np.random.default_rng(seed)
        _decorate_alloy_types(
            data,
            {element_A: 100.0 - fraction_B_atpct, element_B: fraction_B_atpct},
            config.potential.type_map,
            int(config.potential.type_map[element_A]),
            rng,
        )
        f_int = int(round(fraction_B_atpct))
        status = f"alloy_random_{element_A}{100 - f_int}_{element_B}{f_int}"

    _scale_cell_to_lc(data, lc)
    data.initialization(normalization=False, style=1)
    _write_structure(data, out_data)
    return status


def build_multi_alloy_structure(
    case_dir: Path,
    config: PotentialConfig,
    composition_atpct: dict[str, float],
    structure_name: str,
    ordering: str = "random",
    relaxed_lc_host: float | None = None,
    seed: int | None = None,
    n_replicate: int | None = None,
) -> str:
    """Build an N-element (N >= 2) random-substitution alloy supercell.

    ``composition_atpct`` must already be fully resolved (no ``None`` values,
    see ``config.resolve_composition``) and sum to ~100. The first key is
    used as the host lattice species for seeding/replication. Non-host
    elements are allocated largest-fraction-first so rounding remainders
    land on the biggest species; the host absorbs whatever is left over, so
    only ``ordering="random"`` is supported (no B2-style sublattice ordering
    for N > 2 species). ``seed`` defaults to ``None`` (fresh random draw
    every call).
    """
    elements = list(composition_atpct.keys())
    if len(elements) < 2:
        raise ValueError("composition_atpct needs at least 2 elements for an alloy.")
    host = elements[0]

    requested = [o.strip().lower() for o in ordering.split("|") if o.strip()]
    if requested and "random" not in requested:
        raise ValueError(
            f"Only 'random' ordering is supported for {len(elements)}-element "
            f"alloys (got {ordering!r})."
        )

    n = int(n_replicate if n_replicate is not None else config.potential.size_alloy)
    lc = float(relaxed_lc_host if relaxed_lc_host is not None else config.potential.lc_initial)
    out_data = case_dir / "structure.data"

    data = _build_unit_cell(case_dir, config, structure_name, host, n)
    rng = np.random.default_rng(seed)
    _decorate_alloy_types(
        data, composition_atpct, config.potential.type_map,
        int(config.potential.type_map[host]), rng,
    )

    tag = "_".join(f"{el}{int(round(composition_atpct[el]))}" for el in elements)
    status = f"alloy_random_{tag}"

    _scale_cell_to_lc(data, lc)
    data.initialization(normalization=False, style=1)
    _write_structure(data, out_data)
    return status


def build_gas_complex_structure(
    case_dir: Path,
    metal: str,
    gas: str,
    n_gas: int,
    m_defect: int,
    defect_kind: str,
    config: PotentialConfig,
    structure_name: str,
    relaxed_lc: float | None = None,
    n_replicate: int | None = None,
) -> str:
    """Build (gas)_n(V)_m or (gas)_n(I)_m near the cell centre.

    ``defect_kind`` is "vacancy" or "interstitial". Vacancies are the m nearest
    atoms to the centre; interstitials are placed at successive offsets from
    a tetrahedral site around the centre. Gas atoms are then placed close to
    the defect cluster — the LAMMPS minimisation in the case dir resolves them
    into the lowest-energy configuration. Gas atoms are always inserted with
    their real, final ``gas`` type -- they are never part of the metal-side
    chemistry decoration.
    """
    out_data = case_dir / "structure.data"
    n = int(n_replicate if n_replicate is not None else config.potential.size_single)
    lc = float(relaxed_lc if relaxed_lc is not None else config.potential.lc_initial)
    metal_type = int(config.potential.type_map[metal])
    gas_type = int(config.potential.type_map[gas])

    data = _build_unit_cell(case_dir, config, structure_name, metal, n)
    c = np.asarray(data.get_center("cart"), dtype=float).reshape(3,)

    # --- defect side (built at lc=1; scaled below) ---
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
            off = offsets[k % len(offsets)]
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
        off = gas_offsets[k % len(gas_offsets)]
        # Tiny irrational shift to avoid exact overlap when k > len(gas_offsets).
        off = off + np.array([1e-3, 1e-3, 1e-3]) * (k // len(gas_offsets))
        _add_interstitial_cart(data, c, off, gas_type)

    _scale_cell_to_lc(data, lc)
    data.initialization(normalization=False, style=1)
    _write_structure(data, out_data)
    return f"{gas}{n_gas}_{'V' if defect_kind == 'vacancy' else 'I'}{m_defect}"


def build_single_gas_in_bulk(
    case_dir: Path,
    metal: str,
    gas: str,
    site: str,
    config: PotentialConfig,
    structure_name: str,
    relaxed_lc: float | None = None,
    n_replicate: int | None = None,
) -> str:
    """One gas atom on a single tetrahedral or octahedral site — needed for E_He^f."""
    out_data = case_dir / "structure.data"
    n = int(n_replicate if n_replicate is not None else config.potential.size_single)
    lc = float(relaxed_lc if relaxed_lc is not None else config.potential.lc_initial)
    gas_type = int(config.potential.type_map[gas])

    data = _build_unit_cell(case_dir, config, structure_name, metal, n)
    c = np.asarray(data.get_center("cart"), dtype=float).reshape(3,)
    off = _interstitial_offset(structure_name, site, 1.0)
    _add_interstitial_cart(data, c, off, gas_type)

    _scale_cell_to_lc(data, lc)
    _write_structure(data, out_data)
    return f"{gas}_in_{metal}_{site}"
