"""Crystal-agnostic defect machinery shared by the per-lattice builders.

``LatticeDefect`` holds everything that does not depend on the Bravais
lattice: point defects, voids, the dislocation orientation/dispatch logic,
the prismatic-loop pipeline (footprint selection -> SIA insert / vacancy
closure), Frenkel pairs and coordinate transforms.

A concrete subclass supplies the crystallography through three hooks:

    _default_burgerm(lattice_const, burgers) -> |b|
        Magnitude of the Burgers vector for that lattice's common
        dislocations.

    _family_loop_geometry(family, habit_plane, a) -> (n, |b|, k, d_layer)
        Habit-plane unit normal, Burgers magnitude, the number ``k`` of
        atomic layers spanned by one Burgers repeat, and the interplanar
        spacing. The invariant ``k * d_layer == |b|`` must hold -- it is
        what makes the inserted platelet continue the stacking.

    _interstitial_offsets(lattice_const) -> (m, 3) Cartesian offsets
        Candidate interstitial sites relative to a host atom.

See ``bcc_defect.BCCDefect`` and ``fcc_defect.FCCDefect`` for the two
implementations.
"""

from __future__ import annotations

import warnings
from typing import Optional

import numpy as np
from pymatgen.core.lattice import Lattice
from mylammps.inputs.data import find_symmop_lattices

from gz_toolkit.core.modlmp import Modlmp_LmpData


class LatticeDefect:
    """Base defect builder; subclass per Bravais lattice."""

    #: human-readable lattice name, set by subclasses (used in messages)
    lattice_name = "generic"

    def __init__(self, data: Modlmp_LmpData):
        self.data = data

    # ==================================================================
    # Crystallography hooks -- subclasses must override
    # ==================================================================

    @staticmethod
    def _default_burgerm(lattice_const: float, burgers) -> float:
        raise NotImplementedError("Subclass must provide _default_burgerm.")

    #: Prismatic-loop families available on THIS lattice. Loops are
    #: structure-dependent -- bcc grows <111> and <100> loops, fcc grows
    #: <111> Frank and <110> perfect loops -- so each subclass declares its
    #: own registry and asking for a family the lattice does not have is an
    #: error rather than silently producing a wrong platelet.
    #:
    #: Each entry:
    #:   allowed_abs : sorted |hkl| signature the habit plane must match
    #:   label       : family name for error messages ("100", not "001")
    #:   burgers     : a -> |b|
    #:   k           : atomic layers spanned by one Burgers repeat
    #:   d_layer     : a -> interplanar spacing
    #:   description : short human-readable summary
    #: The invariant k * d_layer == |b| is enforced at lookup time.
    LOOP_FAMILIES: dict[str, dict] = {}

    @classmethod
    def loop_families(cls) -> dict[str, str]:
        """Loop family name -> description, for this lattice."""
        return {k: v.get("description", "") for k, v in cls.LOOP_FAMILIES.items()}

    @classmethod
    def _family_loop_geometry(cls, family, habit_plane, a: float):
        """(unit cubic normal, |b|, k layers, d_layer) for a loop family.

        Resolved from this lattice's LOOP_FAMILIES registry.
        """
        fam = str(family).replace("<", "").replace(">", "").strip()
        spec = cls.LOOP_FAMILIES.get(fam)
        if spec is None:
            avail = ", ".join(
                f"<{k}> ({v.get('description', '')})" for k, v in cls.LOOP_FAMILIES.items()
            ) or "none"
            raise ValueError(
                f"{cls.lattice_name} has no <{fam}> prismatic-loop family. "
                f"Available for {cls.lattice_name}: {avail}."
            )

        a = float(a)
        _hkl, n = cls._parse_family(
            habit_plane, spec["allowed_abs"], "habit_plane", label=spec["label"],
        )
        b_mag = float(spec["burgers"](a))
        k = int(spec["k"])
        d_layer = float(spec["d_layer"](a))

        # b must be a whole number of atomic layers, or the inserted platelet
        # is not a lattice translation and the stacking will not join up.
        if abs(k * d_layer - b_mag) > 1e-9 * max(1.0, b_mag):
            raise ValueError(
                f"LOOP_FAMILIES['{fam}'] for {cls.lattice_name} is inconsistent: "
                f"k*d_layer = {k * d_layer:.6f} but |b| = {b_mag:.6f}."
            )
        return n, b_mag, k, d_layer

    @staticmethod
    def _interstitial_offsets_frac() -> np.ndarray:
        """Interstitial sites as fractions of the cubic lattice constant."""
        raise NotImplementedError("Subclass must provide _interstitial_offsets_frac.")

    @classmethod
    def _interstitial_offsets(cls, lattice_const: float) -> np.ndarray:
        """Cartesian interstitial offsets (Angstrom) from a host atom."""
        return cls._interstitial_offsets_frac() * float(lattice_const)

    # ==================================================================
    # Point defects
    # ==================================================================

    def vacancy(self, coords, *, tolerance: float = 1e-3, reinit: bool = True):
        return self.data.delete_atoms("coords", coords=coords, tolerance=tolerance, reinit=reinit)

    def interstitial(self, sites, *, atom_type: int, wrap: bool = False, reinit: bool = True):
        return len(self.data.add_atoms(sites, atom_type=atom_type, wrap=wrap, reinit=reinit))

    def interstitial_cluster(self, sites, *, atom_type: int, wrap: bool = False, reinit: bool = True):
        return self.interstitial(sites, atom_type=atom_type, wrap=wrap, reinit=reinit)

    def dumbbell(self, center, *, atom_type: int, direction=(1, 1, 1), separation: float = 1.0,
                 remove_center_atom: bool = True, tolerance: float = 1e-3, wrap: bool = True,
                 reinit: bool = True):
        direction = self._unit(direction)
        center = np.asarray(center, dtype=float).reshape(3,)
        if remove_center_atom:
            n_removed = self.data.delete_atoms(
                "coords", coords=center, tolerance=tolerance, reinit=False,
            )
            if not n_removed:
                # Silently skipping the removal turns a dumbbell (net +1 atom)
                # into an interstitial pair (net +2), which is easy to miss in
                # a formation-energy calculation. Usually means `center` is not
                # a lattice site -- e.g. get_center("cart") returns the atom
                # centroid, which need not coincide with an atom.
                warnings.warn(
                    f"dumbbell(remove_center_atom=True) found no atom within "
                    f"tolerance={tolerance} of center={center.tolist()}; nothing "
                    f"was removed, so this adds 2 atoms rather than replacing 1. "
                    f"Pass an actual lattice site, or increase tolerance.",
                    stacklevel=2,
                )
        sites = np.vstack([
            center - 0.5 * float(separation) * direction,
            center + 0.5 * float(separation) * direction,
        ])
        self.data.add_atoms(sites, atom_type=atom_type, wrap=wrap, reinit=False)
        if reinit:
            self.data.initialization(normalization=False, style=1)
        return 2

    def dumbbell_100(self, center, **kwargs):
        return self.dumbbell(center, direction=(1, 0, 0), **kwargs)

    def dumbbell_110(self, center, **kwargs):
        return self.dumbbell(center, direction=(1, 1, 0), **kwargs)

    def dumbbell_111(self, center, **kwargs):
        return self.dumbbell(center, direction=(1, 1, 1), **kwargs)

    def impurities(self, region_type: str, *, new_type: int, ff_elements, atomic_masses, **kwargs):
        return self.data.mod_atom_type(
            region_type, new_type=new_type, ff_elements=ff_elements,
            atomic_masses=atomic_masses, **kwargs,
        )

    def impurity_atom(self, site, *, atom_type: int, wrap: bool = False, reinit: bool = True):
        return self.interstitial(site, atom_type=atom_type, wrap=wrap, reinit=reinit)

    # ==================================================================
    # Voids
    # ==================================================================

    def void(self, center=None, radius: float = 0.0, *, region_type: str = "sphere",
             reinit: bool = True, **kwargs):
        return self.create_void(center, radius, region_type=region_type, reinit=reinit, **kwargs)

    def create_void(self, center, radius: float, *, region_type: str = "sphere",
                    reinit: bool = True, **kwargs):
        """Delete atoms to create a void.

        Parameters
        ----------
        center : (3,) Cartesian center (None -> box center)
        radius : float
            Radius for sphere/cylinder, half-side for cube.
        region_type : {"sphere", "cube"/"box", "cylinder"}
        kwargs :
            Forwarded to region_mask(), e.g. plane="xy", height=20 for cylinder.
        """
        data = self.data
        if data.atoms is None or len(data.atoms) == 0:
            return 0

        if center is None:
            center = data.get_center("cart")
        else:
            center = np.asarray(center, dtype=float).reshape(3,)
        rt = str(region_type).lower()

        if rt == "sphere":
            n_deleted = data.delete_atoms("sphere", center=center, radius=float(radius), reinit=reinit)
        elif rt in ("cube", "box"):
            cube_params = {"center": center, "side_length": 2.0 * float(radius)}
            n_deleted = data.delete_atoms("cube", cube_params=cube_params, reinit=reinit)
        elif rt == "cylinder":
            if "plane" not in kwargs:
                raise ValueError("cylinder requires plane='xy'/'xz'/'yz'")
            if "height" not in kwargs:
                raise ValueError("cylinder requires height=...")
            n_deleted = data.delete_atoms(
                "cylinder", plane=kwargs["plane"], center=center,
                radius=float(radius), height=float(kwargs["height"]), reinit=reinit,
            )
        else:
            raise ValueError(
                f"region_type must be one of: 'sphere', 'cube', 'box', "
                f"'cylinder'. Got '{region_type}'."
            )
        return int(n_deleted)

    # ==================================================================
    # Dislocations (myLAMMPS lmpData backend)
    #
    # myLAMMPS conventions the data must satisfy BEFORE the call:
    #   edge : y || Burgers vector, z || glide-plane normal, x || line
    #   screw: z || Burgers vector (= line direction for a screw)
    #
    # The methods below handle that automatically: give the crystallographic
    # directions currently on x/y/z (current_axes) plus the desired
    # burgers/glide spec and the needed swap_axes permutation is derived.
    # ==================================================================

    def dislocation(self, kind: str = "edge", **kwargs):
        """Dispatch to edge_dislocation / screw_dislocation by kind."""
        k = str(kind).lower()
        if k == "edge":
            return self.edge_dislocation(**kwargs)
        if k == "screw":
            return self.screw_dislocation(**kwargs)
        raise ValueError("kind must be 'edge' or 'screw'.")

    @staticmethod
    def _parse_hkl(spec) -> np.ndarray:
        """Parse a crystallographic direction: (1,1,1), [1,-1,0], or '1 -1 0'."""
        if isinstance(spec, str):
            s = spec.replace(",", " ").replace("(", " ").replace(")", " ")
            s = s.replace("[", " ").replace("]", " ").replace("<", " ").replace(">", " ")
            parts = [p for p in s.split() if p]
            if len(parts) != 3:
                raise ValueError(f"Direction must have 3 integers, got '{spec}'.")
            return np.array([int(p) for p in parts], dtype=float)
        return np.asarray(spec, dtype=float).reshape(3,)

    @classmethod
    def _find_parallel_axis(cls, target, axes_unit: list[np.ndarray], what: str) -> int:
        """Index of the (anti)parallel axis among axes_unit, or raise."""
        t = cls._unit(target)
        hits = [i for i, ax in enumerate(axes_unit)
                if abs(abs(float(np.dot(t, ax))) - 1.0) < 1e-6]
        if len(hits) != 1:
            axes_str = [a.round(4).tolist() for a in axes_unit]
            raise ValueError(
                f"{what} direction {np.asarray(target).tolist()} matches "
                f"{len(hits)} of the current axes {axes_str}. It must be parallel "
                f"(or antiparallel) to exactly one box vector — swap_axes can only "
                f"permute existing axes, so build the supercell with this "
                f"orientation first."
            )
        return hits[0]

    def _orient_for_dislocation(self, *, burgers, glide_plane, current_axes,
                                newaxis, burgers_target: int) -> list[int]:
        """Resolve the swap_axes permutation for a dislocation geometry.

        burgers_target: which axis index the Burgers vector must end up on
        (1 = y for edge, 2 = z for screw). For edge, the glide-plane normal
        goes to z and the remaining axis to x. For screw (glide_plane=None),
        the remaining two axes keep their original relative order.
        """
        if newaxis is not None:
            perm = [int(i) for i in newaxis]
            if sorted(perm) != [0, 1, 2]:
                raise ValueError(f"newaxis must be a permutation of [0,1,2], got {newaxis}.")
            return perm

        if current_axes is None:
            return [0, 1, 2]  # trust the caller: data already oriented

        axes = [self._unit(self._parse_hkl(ax)) for ax in current_axes]
        if len(axes) != 3:
            raise ValueError("current_axes must give the 3 directions on x, y, z.")

        ib = self._find_parallel_axis(self._parse_hkl(burgers), axes, "Burgers")

        if glide_plane is not None:
            if burgers_target != 1:
                raise ValueError("glide_plane is only meaningful for edge (burgers_target=1).")
            ig = self._find_parallel_axis(self._parse_hkl(glide_plane), axes, "Glide-plane normal")
            if ig == ib:
                raise ValueError("Burgers vector and glide-plane normal map to the same axis.")
            iline = ({0, 1, 2} - {ib, ig}).pop()
            return [iline, ib, ig]

        others = [i for i in (0, 1, 2) if i != ib]
        perm = [0, 0, 0]
        perm[burgers_target] = ib
        rest = [i for i in (0, 1, 2) if i != burgers_target]
        for slot, src in zip(rest, others):
            perm[slot] = src
        return perm

    def edge_dislocation(self, *, ff_elements: list[str], lattice_const: float | None = None,
                         burgers=None, glide_plane=None, current_axes: list | None = None,
                         newaxis: list[int] | None = None, burgerm: float | None = None,
                         nedges: int = 1, fy_start: float | None = None,
                         add_vacuum: bool = False, vacuum_direction: int = 2,
                         lvac: float = 20.0, restore_axes: bool = False,
                         reset_ids: bool = True) -> dict:
        """Insert edge dislocation(s) via lmpData.create_edge_dislocation.

        The data must already be a built supercell. Orientation is handled
        here: pass current_axes (+ burgers/glide_plane) for automatic axis
        permutation, pass newaxis=[...] directly, or pass neither if the data
        already has y || b, z || glide normal, x || line.
        """
        data = self.data
        if data.atoms is None or len(data.atoms) == 0:
            raise ValueError("No atoms loaded.")
        data.assert_force_field(ff_elements)

        if burgers is None:
            burgers = self.default_burgers
        if glide_plane is None:
            glide_plane = self.default_glide_plane

        if burgerm is None:
            if lattice_const is None:
                raise ValueError("Provide lattice_const (for auto |b|) or burgerm=...")
            burgerm = self._default_burgerm(lattice_const, burgers)
        burgerm = float(burgerm)

        perm = self._orient_for_dislocation(
            burgers=burgers, glide_plane=glide_plane, current_axes=current_axes,
            newaxis=newaxis, burgers_target=1,
        )
        n_before = int(len(data.atoms))

        if perm != [0, 1, 2]:
            data.swap_axes(perm)

        data.create_edge_dislocation(
            burgerm, nedges=int(nedges), fy_start=fy_start, add_vacuum=add_vacuum,
            direction=int(vacuum_direction), lvac=float(lvac), reset_ids=reset_ids,
        )

        restored = False
        if restore_axes and perm != [0, 1, 2]:
            data.swap_axes(np.argsort(perm).tolist())
            restored = True

        return {
            "kind": "edge", "burgerm": burgerm, "nedges": int(nedges), "newaxis": perm,
            "n_atoms_before": n_before, "n_atoms_after": int(len(data.atoms)),
            "n_atoms_removed": n_before - int(len(data.atoms)), "axes_restored": restored,
        }

    def screw_dislocation(self, *, ff_elements: list[str], lattice_const: float | None = None,
                          burgers=None, current_axes: list | None = None,
                          newaxis: list[int] | None = None, burgerm: float | None = None,
                          nscrews: int = 1, style: str | None = None,
                          handle_pbc: str = "tilt", orientation: bool = True,
                          add_vacuum: bool = False, vacuum_direction: int = 0,
                          lvac: float = 20.0, restore_axes: bool = False,
                          center_offset=(0.0, 0.0)) -> dict:
        """Insert screw dislocation(s) via lmpData.create_screw_dislocation.

        Convention: z || Burgers vector (= line direction). Orientation
        handling mirrors edge_dislocation (b is pinned to z; the other two
        axes keep their relative order).

        The core defaults to the atom nearest the box center (or the
        dipole/quadrupole positions for nscrews=2/4). center_offset nudges
        the core in the local x-y plane (Cartesian Angstrom) after that snap.
        """
        data = self.data
        if data.atoms is None or len(data.atoms) == 0:
            raise ValueError("No atoms loaded.")
        data.assert_force_field(ff_elements)

        if burgers is None:
            burgers = self.default_burgers
        if style is None:
            style = self.screw_style

        if burgerm is None:
            if lattice_const is None:
                raise ValueError("Provide lattice_const (for auto |b|) or burgerm=...")
            burgerm = self._default_burgerm(lattice_const, burgers)
        burgerm = float(burgerm)

        perm = self._orient_for_dislocation(
            burgers=burgers, glide_plane=None, current_axes=current_axes,
            newaxis=newaxis, burgers_target=2,
        )
        n_before = int(len(data.atoms))

        if perm != [0, 1, 2]:
            data.swap_axes(perm)

        data.create_screw_dislocation(
            burgerm, int(nscrews), style=str(style), handle_pbc=str(handle_pbc),
            orientation=orientation, add_vacuum=add_vacuum,
            direction=int(vacuum_direction), lvac=float(lvac), center_offset=center_offset,
        )

        restored = False
        if restore_axes and perm != [0, 1, 2]:
            data.swap_axes(np.argsort(perm).tolist())
            restored = True

        return {
            "kind": "screw", "burgerm": burgerm, "nscrews": int(nscrews), "newaxis": perm,
            "n_atoms_before": n_before, "n_atoms_after": int(len(data.atoms)),
            "axes_restored": restored,
        }

    # ==================================================================
    # Prismatic loops -- shared pipeline
    #
    # A loop is a coherent platelet of k atomic layers along the habit
    # normal. Both endings share the selection stage; only the "apply"
    # step differs:
    #   sil -> copy the slab and split it +/- b/2 along the normal
    #   vl  -> delete the slab and move the closure layers inward
    # ==================================================================

    def add_loop(self, *, family: str, radius: float, lattice_const: float,
                 ff_elements: list[str], atomic_masses: list[float], center,
                 loop_type: str = "sil", habit_plane=None, box_axes=None,
                 normal_box=None, shape: str = "round", n_repeats: int = 1,
                 triangle_rotation_deg: float = 0.0, in_plane_ref=None,
                 loop_atom_type: int | None = None, closure_layers: int = 2,
                 closure_disp_frac: float = 0.15, select_types=None,
                 wrap: bool = True, reinit: bool = True):
        """Add a prismatic loop of one of THIS lattice's families.

        ``family`` must be a key of this class's ``LOOP_FAMILIES`` (see
        :meth:`loop_families`), so a bcc builder accepts '111'/'100' and an
        fcc builder '111'/'110'. ``habit_plane`` defaults to the family's
        own signature (e.g. (1,1,1) for '111').

        The named wrappers on each subclass (``add_111_loop``,
        ``add_frank_loop``, ...) are thin shims over this.
        """
        fam = str(family).replace("<", "").replace(">", "").strip()
        if habit_plane is None:
            spec = self.LOOP_FAMILIES.get(fam)
            if spec is None:
                # let _family_loop_geometry raise the informative error
                self._family_loop_geometry(fam, (1, 1, 1), lattice_const)
            habit_plane = spec["default_habit"]
        return self._add_loop(
            family=fam, loop_type=loop_type, radius=radius,
            lattice_const=lattice_const, ff_elements=ff_elements,
            atomic_masses=atomic_masses, center=center, habit_plane=habit_plane,
            box_axes=box_axes, normal_box=normal_box, shape=shape,
            n_repeats=n_repeats, loop_atom_type=loop_atom_type,
            closure_layers=closure_layers, closure_disp_frac=closure_disp_frac,
            select_types=select_types, wrap=wrap, reinit=reinit,
            triangle_rotation_deg=triangle_rotation_deg, in_plane_ref=in_plane_ref,
        )

    def _add_loop(self, *, family, loop_type, radius, lattice_const, ff_elements,
                  atomic_masses, center, habit_plane, box_axes, normal_box, shape,
                  n_repeats, loop_atom_type, closure_layers, closure_disp_frac,
                  select_types, wrap, reinit, triangle_rotation_deg=0.0,
                  in_plane_ref=None):
        """Shared core: select the k-layer platelet, then dispatch sil/vl."""
        lt = str(loop_type).lower()
        is_sil = lt in ("sil", "sia", "interstitial")
        is_vl = lt in ("vl", "va", "vacancy")
        if not (is_sil or is_vl):
            raise ValueError("loop_type must be 'sil'/'interstitial' or 'vl'/'vacancy'.")

        ctx = self._select_loop_atoms(
            family=family, radius=radius, lattice_const=lattice_const,
            ff_elements=ff_elements, atomic_masses=atomic_masses, center=center,
            habit_plane=habit_plane, shape=shape, n_repeats=n_repeats,
            select_types=select_types, box_axes=box_axes, normal_box=normal_box,
            triangle_rotation_deg=triangle_rotation_deg, in_plane_ref=in_plane_ref,
        )
        if ctx is None:
            return 0

        if is_sil:
            return self._apply_sia_insert(
                ctx, loop_atom_type=loop_atom_type, wrap=wrap, reinit=reinit,
            )
        return self._apply_vacancy_closure(
            ctx, burgers_vector=None, closure_layers=closure_layers,
            closure_disp_frac=closure_disp_frac, wrap_move_atoms=wrap, reinit=reinit,
        )

    @classmethod
    def _resolve_loop_normal(cls, n_cubic, box_axes, normal_box):
        """Box-frame unit normal.

        normal_box (given directly) > rotate n_cubic via box_axes (rows =
        crystal directions on box x/y/z; orthonormal => pure rotation) >
        n_cubic (cubic-aligned cell).
        """
        if normal_box is not None:
            n = np.asarray(normal_box, dtype=float).reshape(3,)
            return n / np.linalg.norm(n)
        if box_axes is not None:
            A = np.asarray(box_axes, dtype=float).reshape(3, 3)
            M = A / np.linalg.norm(A, axis=1, keepdims=True)
            n = M @ np.asarray(n_cubic, dtype=float).reshape(3,)
            return n / np.linalg.norm(n)
        return np.asarray(n_cubic, dtype=float).reshape(3,)

    @staticmethod
    def _triangle_mask(cu, cv, r, rotation_deg=0.0):
        """Interior of an equilateral triangle of circumradius r.

        Three inward half-plane constraints p . m_k <= r/2, with the edge
        normals m_k 120 deg apart. The inradius is r/2, and the vertices sit
        exactly on the boundary. ``rotation_deg`` spins the triangle in the
        habit plane so its edges can be aligned with a chosen direction
        (for fcc Frank loops / SFT precursors the edges run along <110>).
        """
        rot = np.deg2rad(float(rotation_deg))
        mask = np.ones(cu.shape, dtype=bool)
        for k in range(3):
            phi = -0.5 * np.pi + k * (2.0 * np.pi / 3.0) + rot
            mask &= (cu * np.cos(phi) + cv * np.sin(phi)) <= 0.5 * float(r)
        return mask

    def _select_loop_atoms(self, *, family: str, radius: float, lattice_const: float,
                           ff_elements: list[str], atomic_masses: list[float], center,
                           habit_plane, shape: str, n_repeats: int, select_types,
                           box_axes=None, normal_box=None,
                           triangle_rotation_deg: float = 0.0,
                           in_plane_ref=None) -> dict | None:
        """Resolve the box-frame normal and select EXACTLY k*n_repeats layers.

        k and the layer spacing come from the lattice's _family_loop_geometry.
        The selection snaps onto the actual atomic planes (circular-mean
        phase) and takes nlay = k*n_repeats consecutive layers centered on
        the loop. plane_tol = 0.30 * d_layer. PBC-safe via minimum-image
        relative vectors (triclinic-safe).
        """
        data = self.data

        if ff_elements is None or atomic_masses is None:
            raise ValueError("Must provide ff_elements and atomic_masses.")
        if len(ff_elements) != len(atomic_masses):
            raise ValueError("ff_elements and atomic_masses must have the same length.")
        data.assert_force_field(ff_elements, atomic_masses=atomic_masses)

        if center is None:
            raise ValueError("center must be provided in Cartesian coordinates (shape (3,)).")
        c = np.asarray(center, dtype=float).reshape(3,)

        n_repeats = int(n_repeats)
        if n_repeats < 1:
            raise ValueError("n_repeats must be >= 1.")

        shape = str(shape).lower()
        if shape not in ("round", "square", "triangle"):
            raise ValueError("shape must be 'round', 'square' or 'triangle'.")

        a = float(lattice_const)
        r = float(radius)

        n_cubic, b_mag, k, d_layer = self._family_loop_geometry(family, habit_plane, a)
        n = self._resolve_loop_normal(n_cubic, box_axes, normal_box)
        nlay = k * n_repeats

        # A triangular footprint must be locked to the lattice rows, so default
        # its in-plane reference to a <110> direction in the habit plane and map
        # it into the box frame the same way the normal is mapped.
        if in_plane_ref is None and shape == "triangle":
            in_plane_ref = self._perpendicular_110(n_cubic)
        ref_box = None
        if in_plane_ref is not None:
            ref_box = self._resolve_loop_normal(
                self._unit(in_plane_ref), box_axes, normal_box,
            )

        u, v = self._plane_basis_from_normal(n, ref_box)
        plane_tol = 0.30 * d_layer

        pos_all = data.atoms[["x", "y", "z"]].to_numpy(dtype=float)
        rel_all = data.min_image_rel(pos_all, c)

        type_mask = None
        if select_types is not None:
            type_mask = data.select_type(select_types)

        dist = rel_all @ n                     # signed distance along the normal
        proj = rel_all - np.outer(dist, n)     # in-plane projection
        cu = proj @ u
        cv = proj @ v
        if shape == "round":
            footprint = (cu * cu + cv * cv) <= (r * r)
        elif shape == "square":
            footprint = (np.abs(cu) <= r) & (np.abs(cv) <= r)
        else:
            footprint = self._triangle_mask(cu, cv, r, triangle_rotation_deg)

        cand = footprint.copy()
        if type_mask is not None:
            cand &= type_mask
        if not np.any(cand):
            return None

        # Anchor on the ACTUAL atomic planes (circular mean of dist mod
        # d_layer), snap the center plane to the nearest layer, then take
        # exactly nlay consecutive layers. A plain width window is
        # phase-dependent and grabs k-1 or k layers; this reliably gives nlay.
        ang = 2.0 * np.pi * dist[cand] / d_layer
        phase = float(np.angle(np.mean(np.exp(1j * ang)))) / (2.0 * np.pi) * d_layer
        s_center = phase + round((0.0 - phase) / d_layer) * d_layer
        offsets = [s_center + (i - nlay // 2) * d_layer for i in range(nlay)]

        donor = np.zeros(len(pos_all), dtype=bool)
        for s in offsets:
            on_layer = footprint & (np.abs(dist - s) <= plane_tol)
            if type_mask is not None:
                on_layer &= type_mask
            donor |= on_layer

        sel = np.flatnonzero(donor)
        if sel.size == 0:
            return None

        return {
            "sel": sel, "footprint": footprint, "type_mask": type_mask,
            "rel_all": rel_all, "n": n, "u": u, "v": v, "a": a, "b_mag": b_mag,
            "k": k, "d_layer": d_layer, "plane_tol": plane_tol,
        }

    def _apply_sia_insert(self, ctx: dict, *, loop_atom_type: int | None,
                          wrap: bool, reinit: bool) -> int:
        """SIL ending: copy the selected k-layer slab and split it +/- b/2.

        Each selected atom -> r - b/2 (original moves) and a NEW atom at
        r + b/2, where b = |b| * n is the full Burgers vector along the habit
        normal (a lattice translation), so the inserted layers continue the
        stacking. Returns the number of atoms inserted.
        """
        data = self.data
        sel = ctx["sel"]
        half = 0.5 * ctx["b_mag"] * ctx["n"]

        sel_df = data.atoms.iloc[sel][["type", "x", "y", "z"]].copy()
        old_pos = sel_df[["x", "y", "z"]].to_numpy(dtype=float)
        sel_index = data.atoms.index.to_numpy()[sel]

        moved_old = old_pos - half.reshape(1, 3)
        new_pos = old_pos + half.reshape(1, 3)
        if wrap:
            moved_old = data.wrap_cart(moved_old)
            new_pos = data.wrap_cart(new_pos)
        data.atoms.loc[sel_index, ["x", "y", "z"]] = moved_old

        if loop_atom_type is None:
            new_types = sel_df["type"].to_numpy(dtype=int)
        else:
            valid_types = set(int(i) for i in data.masses.index.to_numpy(dtype=int))
            if int(loop_atom_type) not in valid_types:
                raise ValueError(
                    f"loop_atom_type={loop_atom_type} not in Masses after assert_force_field."
                )
            new_types = np.full(len(sel_df), int(loop_atom_type), dtype=int)

        for t_id in np.unique(new_types):
            data.add_atoms(new_pos[new_types == t_id], atom_type=int(t_id), wrap=False, reinit=False)

        if reinit:
            data.initialization(normalization=False, style=1)

        return int(sel.size)

    def _apply_vacancy_closure(self, ctx: dict, *, burgers_vector, closure_layers: int,
                               closure_disp_frac: float, wrap_move_atoms: bool,
                               reinit: bool) -> int:
        """Vacancy ending: delete selected atoms + move closure layers inward.

        Direction convention: signed distance is along b_hat; atoms ABOVE the
        deleted slab move by -delta*b_hat, atoms BELOW by +delta*b_hat.
        """
        data = self.data
        sel = ctx["sel"]
        footprint = ctx["footprint"]
        type_mask = ctx["type_mask"]
        rel_all = ctx["rel_all"]
        n = ctx["n"]
        a = ctx["a"]
        d_layer = ctx["d_layer"]
        plane_tol = ctx["plane_tol"]

        if burgers_vector is None:
            b_hat = np.asarray(n, dtype=float).reshape(3,)
        else:
            b_hat = np.asarray(burgers_vector, dtype=float).reshape(3,)
            nb = np.linalg.norm(b_hat)
            if nb == 0:
                raise ValueError("burgers_vector must be non-zero.")
            b_hat = b_hat / nb

        delta = float(closure_disp_frac) * a
        n_close = int(closure_layers)
        if n_close < 0:
            raise ValueError("closure_layers must be >= 0.")

        dist_b = rel_all @ b_hat
        s_bot_b = float(np.min(dist_b[sel]))
        s_top_b = float(np.max(dist_b[sel]))

        if n_close > 0 and delta != 0.0:
            N = len(data.atoms)
            del_mask = np.zeros(N, dtype=bool)
            del_mask[sel] = True

            move_above = np.zeros(N, dtype=bool)
            move_below = np.zeros(N, dtype=bool)

            for j in range(1, n_close + 1):
                sA = s_top_b + j * d_layer
                sB = s_bot_b - j * d_layer
                onA = (np.abs(dist_b - sA) <= plane_tol) & footprint
                onB = (np.abs(dist_b - sB) <= plane_tol) & footprint
                if type_mask is not None:
                    onA &= type_mask
                    onB &= type_mask
                move_above |= onA
                move_below |= onB

            move_above &= ~del_mask
            move_below &= ~del_mask

            if np.any(move_above):
                idxA = data.atoms.index.to_numpy()[move_above]
                posA = data.atoms.loc[idxA, ["x", "y", "z"]].to_numpy(dtype=float)
                posA = posA - delta * b_hat.reshape(1, 3)
                if wrap_move_atoms:
                    posA = data.wrap_cart(posA)
                data.atoms.loc[idxA, ["x", "y", "z"]] = posA

            if np.any(move_below):
                idxB = data.atoms.index.to_numpy()[move_below]
                posB = data.atoms.loc[idxB, ["x", "y", "z"]].to_numpy(dtype=float)
                posB = posB + delta * b_hat.reshape(1, 3)
                if wrap_move_atoms:
                    posB = data.wrap_cart(posB)
                data.atoms.loc[idxB, ["x", "y", "z"]] = posB

        mask = np.zeros(len(data.atoms), dtype=bool)
        mask[sel] = True
        data.atoms = data.atoms.loc[~mask].copy()

        if reinit:
            data.initialization(normalization=False, style=1)

        return int(sel.size)

    # ==================================================================
    # Random interstitials
    # ==================================================================

    def random_interstitials(self, *, fraction: float, atom_type: int,
                             lattice_const: float | None = None, min_dist: float = 1.2,
                             max_trials_per_atom: int = 20,
                             rng: np.random.Generator | None = None,
                             reinit: bool = True) -> int:
        """Add random interstitial atoms on this lattice's tetra/octa sites.

        Candidate sites are random lattice-specific offsets around random
        host atoms. A candidate is accepted only if it is at least
        ``min_dist`` from every existing atom and every previously placed
        interstitial (minimum-image). Existing-atom checks use a periodic
        KD-tree when the box is orthogonal; non-orthogonal boxes fall back to
        a vectorized minimum-image sweep.
        """
        data = self.data
        if data.atoms is None or len(data.atoms) == 0:
            return 0
        if fraction <= 0:
            return 0
        if rng is None:
            rng = np.random.default_rng()

        pos = data.atoms[["x", "y", "z"]].to_numpy(dtype=float)
        n_atoms = len(pos)
        n_add = int(fraction * n_atoms)
        if n_add < 1:
            return 0

        if lattice_const is None:
            # Legacy path: scale the fractional sites by the BOX matrix. That is
            # only right when the box IS one unit cell -- for an N x N x N
            # supercell every offset comes out N times too large, so the
            # "interstitials" land at random spots rather than tetra/octa sites.
            # Kept for backward compatibility; pass lattice_const to get it right.
            warnings.warn(
                "random_interstitials() without lattice_const scales the "
                "interstitial offsets by the box matrix, which is only correct "
                "if the box is a single unit cell. Pass lattice_const=<a> to "
                "place atoms on true tetrahedral/octahedral sites.",
                stacklevel=2,
            )
            offsets_cart = self._interstitial_offsets_frac() @ data.box.matrix.T
        else:
            offsets_cart = self._interstitial_offsets(lattice_const)

        tree = None
        lengths = None
        if bool(getattr(data.box, "is_orthogonal", False)):
            from scipy.spatial import cKDTree
            lengths = np.asarray(data.box.lengths, dtype=float)
            tree = cKDTree(np.mod(data.wrap_cart(pos), lengths), boxsize=lengths)

        def _clashes_existing(q: np.ndarray) -> bool:
            if tree is not None:
                return bool(tree.query_ball_point(np.mod(q, lengths), r=min_dist))
            return bool(np.any(data.min_image_dists(pos, q) < min_dist))

        new_positions: list[np.ndarray] = []
        for _ in range(n_add):
            for _ in range(max_trials_per_atom):
                r0 = pos[rng.integers(0, n_atoms)]
                candidate = data.wrap_cart(r0 + offsets_cart[rng.integers(len(offsets_cart))])
                if _clashes_existing(candidate):
                    continue
                if new_positions and np.any(
                        data.min_image_dists(np.vstack(new_positions), candidate) < min_dist):
                    continue
                new_positions.append(candidate)
                break

        if new_positions:
            data.add_atoms(np.vstack(new_positions), atom_type=atom_type, wrap=False, reinit=False)
        if reinit:
            data.initialization(normalization=False, style=1)

        return int(len(new_positions))

    # ==================================================================
    # Placeholders
    # ==================================================================

    def precipitate(self, *args, **kwargs):
        raise NotImplementedError(
            "Precipitate generation is a placeholder for a future implementation."
        )

    def grain_boundary(self, *args, **kwargs):
        raise NotImplementedError(
            "Grain-boundary/polycrystal generation is a placeholder for a future implementation."
        )

    def gas(self, *, mode: str = "interstitial", atom_type: int, sites=None,
            center=None, radius=None, **kwargs):
        mode = str(mode).lower()
        if mode in ("interstitial", "single"):
            if sites is None:
                if center is None:
                    raise ValueError("gas interstitial mode requires sites=... or center=...")
                sites = center
            return self.interstitial(sites, atom_type=atom_type, **kwargs)
        if mode in ("bubble", "void"):
            if center is None or radius is None:
                raise ValueError("gas bubble mode requires center=... and radius=...")
            removed = self.void(center=center, radius=radius, reinit=False, **kwargs)
            self.interstitial(center, atom_type=atom_type, reinit=True)
            return {"removed_atoms": removed, "gas_atoms": 1}
        raise ValueError("mode must be 'interstitial', 'single', 'bubble', or 'void'.")

    # ==================================================================
    # Geometry helpers
    # ==================================================================

    @staticmethod
    def _unit(v):
        v = np.asarray(v, dtype=float).reshape(3,)
        n = np.linalg.norm(v)
        if n == 0:
            raise ValueError("Zero vector.")
        return v / n

    @classmethod
    def _plane_basis_from_normal(cls, n, in_plane_ref=None):
        """Orthonormal basis (u, v) spanning the plane with unit normal n.

        With ``in_plane_ref=None`` the basis is *arbitrary* (some vector not
        parallel to n). That is fine for a circular or square footprint,
        which are rotationally symmetric enough not to care, but NOT for a
        triangular one: an arbitrary u leaves the triangle misaligned with
        the lattice rows, and the platelet then grows in jumps instead of
        row by row. Pass ``in_plane_ref`` (a direction that should become u,
        projected into the plane) whenever the footprint's orientation
        carries crystallographic meaning.
        """
        n = cls._unit(n)
        if in_plane_ref is not None:
            ref = np.asarray(in_plane_ref, dtype=float).reshape(3,)
            ref = ref - np.dot(ref, n) * n          # project into the plane
            if np.linalg.norm(ref) < 1e-8:
                raise ValueError(
                    "in_plane_ref is parallel to the plane normal; it must have "
                    "a component in the plane."
                )
            u = cls._unit(ref)
        else:
            a = np.array([1.0, 0.0, 0.0])
            if abs(np.dot(a, n)) > 0.9:
                a = np.array([0.0, 1.0, 0.0])
            u = cls._unit(np.cross(n, a))
        v = cls._unit(np.cross(n, u))
        return u, v

    @classmethod
    def _perpendicular_110(cls, n_cubic):
        """A <110> direction lying in the plane with cubic normal n_cubic.

        Used to align triangular loop footprints with the close-packed rows.
        For a {111} normal there are three such directions, all equivalent
        by symmetry, so the first is as good as any.
        """
        n = cls._unit(n_cubic)
        for v in ([1, -1, 0], [0, 1, -1], [-1, 0, 1],
                  [1, 1, 0], [0, 1, 1], [1, 0, 1]):
            vv = np.asarray(v, dtype=float)
            if abs(float(np.dot(cls._unit(vv), n))) < 1e-8:
                return vv
        raise ValueError(f"No <110> direction is perpendicular to {n_cubic}.")

    @classmethod
    def _parse_family(cls, spec, allowed_abs, what: str, label: str | None = None):
        """Parse and validate a direction against an allowed |hkl| signature.

        e.g. allowed_abs=[1, 1, 1] accepts (1,1,1), (-1,1,1), ...
        ``label`` is the family name used in the error message (defaults to
        the digits of allowed_abs); pass it so <100> is not reported as <001>.
        Returns (hkl_tuple, unit_normal).
        """
        if label is None:
            label = "".join(str(a) for a in allowed_abs)
        if isinstance(spec, str):
            s = spec.replace(",", " ").replace("(", " ").replace(")", " ")
            s = s.replace("[", " ").replace("]", " ")
            parts = [p for p in s.split() if p]
            if len(parts) != 3:
                raise ValueError(f"{what} must have 3 integers, got '{spec}'.")
            hkl = [int(p) for p in parts]
        else:
            arr = np.asarray(spec, dtype=int).reshape(3,)
            hkl = [int(arr[0]), int(arr[1]), int(arr[2])]

        if sorted(abs(x) for x in hkl) != sorted(allowed_abs):
            raise ValueError(f"{what} must be a <{label}> direction, got {hkl}.")
        return tuple(hkl), cls._unit(np.array(hkl, dtype=float))

    # ==================================================================
    # Frenkel pairs
    # ==================================================================

    def add_FrenkelPair(self, n: int, *, min_sep_from_vac: float = 1.5,
                        min_sep_from_atoms: float = 1.3, min_sep_between_new: float = 1.3,
                        offset_range: tuple[float, float] = (10.0, 20.0),
                        max_trials: int = 200, rng: Optional[np.random.Generator] = None,
                        eligible_region_type: str | None = None,
                        eligible_region_kwargs: dict | None = None,
                        exclude_region_type: str | None = None,
                        exclude_region_kwargs: dict | None = None,
                        atom_mask: Optional[np.ndarray] = None,
                        wrap_new: bool = True) -> dict:
        """Create n Frenkel pairs by moving random atoms to nearby interstitial sites.

        Region definitions use Modlmp_LmpData.region_mask(), so any of:
        "all", "sphere", "cube"/"box", "cylinder", "plane", "type", "atom_id".

        Distance checks against the unmoved atoms use a periodic KD-tree when
        the box is orthogonal; otherwise a vectorized minimum-image sweep.
        """
        data = self.data
        if data.atoms is None or len(data.atoms) == 0:
            return {"selected_indices": [], "moved_indices": [],
                    "skipped_indices": [], "new_positions": {}}

        df = data.atoms
        if not all(c in df.columns for c in ("x", "y", "z")):
            raise ValueError("data.atoms must contain Cartesian columns: x, y, z.")

        all_idx = df.index.to_numpy(dtype=int)
        n_atoms = len(all_idx)
        elig_mask = np.ones(n_atoms, dtype=bool)

        if atom_mask is not None:
            atom_mask = np.asarray(atom_mask, dtype=bool)
            if len(atom_mask) != n_atoms:
                raise ValueError("atom_mask length must equal len(data.atoms).")
            elig_mask &= atom_mask

        if eligible_region_type is not None:
            rk = {} if eligible_region_kwargs is None else dict(eligible_region_kwargs)
            m = data.region_mask(eligible_region_type, **rk)
            if m is None or len(m) != n_atoms:
                raise ValueError("eligible region_mask() must return a mask aligned with data.atoms.")
            elig_mask &= np.asarray(m, dtype=bool)

        if exclude_region_type is not None:
            rk = {} if exclude_region_kwargs is None else dict(exclude_region_kwargs)
            m = data.region_mask(exclude_region_type, **rk)
            if m is None or len(m) != n_atoms:
                raise ValueError("exclude region_mask() must return a mask aligned with data.atoms.")
            elig_mask &= ~np.asarray(m, dtype=bool)

        elig = all_idx[elig_mask]
        if len(elig) == 0:
            raise ValueError("No eligible atoms to select from (mask/region excluded everything).")

        if rng is None:
            rng = np.random.default_rng()

        n_pick = int(min(int(n), len(elig)))
        chosen = rng.choice(elig, size=n_pick, replace=False)

        pos = df[["x", "y", "z"]].to_numpy(dtype=float)
        index_to_row = {int(idx): i for i, idx in enumerate(all_idx)}

        fixed_mask = np.ones(n_atoms, dtype=bool)
        fixed_mask[[index_to_row[int(i)] for i in chosen]] = False
        fixed_positions = pos[fixed_mask]

        fixed_tree = None
        lengths = None
        if fixed_positions.size and bool(getattr(data.box, "is_orthogonal", False)):
            from scipy.spatial import cKDTree
            lengths = np.asarray(data.box.lengths, dtype=float)
            fixed_tree = cKDTree(np.mod(data.wrap_cart(fixed_positions), lengths), boxsize=lengths)

        def _clashes_fixed(q_cart: np.ndarray) -> bool:
            if not fixed_positions.size:
                return False
            qw = data.wrap_cart(q_cart)
            if fixed_tree is not None:
                return bool(fixed_tree.query_ball_point(np.mod(qw, lengths), r=min_sep_from_atoms))
            return bool(np.any(data.min_image_dists(fixed_positions, qw) < min_sep_from_atoms))

        placed_positions = []
        placed_map: dict[int, tuple] = {}
        moved, skipped = [], []

        r_lo, r_hi = float(offset_range[0]), float(offset_range[1])
        if not (r_lo > 0 and r_hi > r_lo):
            raise ValueError("offset_range must be positive with hi > lo.")

        for atom_id in chosen:
            atom_id = int(atom_id)
            r0 = pos[index_to_row[atom_id]].copy()
            ok = False
            for _ in range(int(max_trials)):
                v = rng.normal(size=3)
                nv = np.linalg.norm(v)
                if nv == 0:
                    continue
                v /= nv
                q = r0 + rng.uniform(r_lo, r_hi) * v
                if wrap_new:
                    q = data.wrap_cart(q)
                if np.linalg.norm(data.min_image_rel(q, r0)) < min_sep_from_vac:
                    continue
                if _clashes_fixed(q):
                    continue
                if placed_positions and np.any(
                        data.min_image_dists(np.vstack(placed_positions), q) < min_sep_between_new):
                    continue
                placed_positions.append(q)
                placed_map[atom_id] = (float(q[0]), float(q[1]), float(q[2]))
                moved.append(atom_id)
                ok = True
                break
            if not ok:
                skipped.append(atom_id)

        if moved:
            for atom_id, q in placed_map.items():
                df.loc[int(atom_id), ["x", "y", "z"]] = q
            data.atoms = df
            data.initialization(normalization=False, style=1)

        report = {
            "selected_indices": list(map(int, chosen.tolist())),
            "moved_indices": moved, "skipped_indices": skipped,
            "new_positions": placed_map, "eligible_count": int(len(elig)),
            "eligible_region_type": eligible_region_type,
            "exclude_region_type": exclude_region_type,
        }
        if skipped:
            print(f"[add_FrenkelPair] placed {len(moved)}/{len(chosen)}; "
                  f"skipped {len(skipped)} (max_trials reached).")
        else:
            print(f"[add_FrenkelPair] placed {len(moved)}/{len(chosen)} Frenkel pairs.")
        return report

    # ==================================================================
    # Coordinate transforms
    # ==================================================================

    def transform(self, old_sys, new_sys, by="supercell", tol=1e-8):
        """Apply a coordinate transformation to the structure.

        by="transmat"  : apply SymmOp to atom coordinates (no replication)
        by="supercell" : apply integer supercell matrix (returns NEW data)
        """
        data = self.data
        old_sys = np.asarray(old_sys, float).reshape(3, 3)
        new_sys = np.asarray(new_sys, float).reshape(3, 3)

        if by == "transmat":
            symmop = find_symmop_lattices(Lattice(old_sys), Lattice(new_sys))
            data.atoms = data.modify_by_symmetry(data.atoms, symmop, normalization=True)
            return data

        if by == "supercell":
            S = new_sys @ np.linalg.inv(old_sys)
            S_round = np.rint(S)
            if not np.allclose(S, S_round, atol=tol):
                raise ValueError(f"new_sys is not an integer combination of old_sys.\nS=\n{S}")
            return data.make_supercell(S_round.astype(int).tolist())

        raise ValueError("by must be 'transmat' or 'supercell'")
