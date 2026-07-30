"""BCC defect creation tools.

Method/class tree
-----------------
bcc_defect.py
`-- BCCDefect
    |-- point defects
    |   |-- vacancy(coords)
    |   |-- interstitial(sites)
    |   |-- interstitial_cluster(sites)
    |   |-- random_interstitials(...)
    |   |-- dumbbell(center, direction)
    |   |-- dumbbell_111(center)
    |   |-- dumbbell_100(center)
    |   |-- dumbbell_110(center)
    |   |-- impurities(region_type)
    |   `-- impurity_atom(site)
    |-- larger defects / sinks
    |   |-- void(center, radius)
    |   |-- dislocation(kind, ...)            # dispatch -> edge/screw
    |   |-- edge_dislocation(...)             # half-plane removal (myLAMMPS backend)
    |   |-- screw_dislocation(...)            # myLAMMPS backend (OVITO variant planned)
    |   |-- add_111_loop(loop_type="sil"|"vl", ...)   # <111> loop (3 layers)
    |   |-- add_100_loop(loop_type="sil"|"vl", ...)   # <100> loop (2 layers)
    |   |-- precipitate(...)                  # placeholder
    |   |-- gas(mode)
    |   `-- grain_boundary(...)               # placeholder
    `-- internal geometry/implementation methods
        |-- _family_loop_geometry(...)    family -> normal, |b|, k layers, d_layer
        |-- _resolve_loop_normal(...)     rotate the cubic normal into the box frame
        |-- _add_loop(...)                shared core (dispatch sil/vl)
        |-- _select_loop_atoms(...)       footprint + snap to exactly k*n_repeats layers
        |-- _apply_sia_insert(...)        SIL ending (copy slab, split +/- b/2)
        `-- _apply_vacancy_closure(...)   VL ending (delete + close gap)  [separate path]
"""

from gz_toolkit.core.modlmp import Modlmp_LmpData
from typing import Optional
import numpy as np
from pymatgen.core.lattice import Lattice
from mylammps.inputs.data import find_symmop_lattices

########################################################################################################################
class BCCDefect:
    """BCC defect builder using Modlmp_LmpData selection and edit methods."""

    def __init__(self, data: Modlmp_LmpData):
        self.data = data

    def vacancy(self, coords, *, tolerance: float = 1e-3, reinit: bool = True):
        return self.data.delete_atoms("coords", coords=coords, tolerance=tolerance, reinit=reinit)

    def interstitial(self, sites, *, atom_type: int, wrap: bool = False, reinit: bool = True):
        return len(self.data.add_atoms(sites, atom_type=atom_type, wrap=wrap, reinit=reinit))

    def interstitial_cluster(self, sites, *, atom_type: int, wrap: bool = False, reinit: bool = True):
        return self.interstitial(sites, atom_type=atom_type, wrap=wrap, reinit=reinit)

    def random_interstitials(self, **kwargs):
        return self.add_interstitial_bcc(**kwargs)

    def dumbbell(self,center,*,atom_type: int,direction=(1, 1, 1),separation: float = 1.0,remove_center_atom: bool = True,tolerance: float = 1e-3,wrap: bool = True,
        reinit: bool = True,):

        direction = self._unit(direction)
        center = np.asarray(center, dtype=float).reshape(3,)
        if remove_center_atom:
            self.data.delete_atoms("coords", coords=center, tolerance=tolerance, reinit=False)
        sites = np.vstack([
            center - 0.5 * float(separation) * direction,
            center + 0.5 * float(separation) * direction,
        ])
        self.data.add_atoms(sites, atom_type=atom_type, wrap=wrap, reinit=False)
        if reinit:
            self.data.initialization(normalization=False, style=1)
        return 2

    def dumbbell_111(self, center, **kwargs):
        return self.dumbbell(center, direction=(1, 1, 1), **kwargs)

    def dumbbell_100(self, center, **kwargs):
        return self.dumbbell(center, direction=(1, 0, 0), **kwargs)

    def dumbbell_110(self, center, **kwargs):
        return self.dumbbell(center, direction=(1, 1, 0), **kwargs)

    def impurities(self, region_type: str, *, new_type: int, ff_elements, atomic_masses, **kwargs):
        return self.data.mod_atom_type(
            region_type,
            new_type=new_type,
            ff_elements=ff_elements,
            atomic_masses=atomic_masses,
            **kwargs,
        )

    def impurity_atom(self, site, *, atom_type: int, wrap: bool = False, reinit: bool = True):
        return self.interstitial(site, atom_type=atom_type, wrap=wrap, reinit=reinit)

    def void(self, center=None, radius: float = 0.0, *, region_type: str = "sphere", reinit: bool = True, **kwargs):
        return self.create_void(center, radius, region_type=region_type, reinit=reinit, **kwargs)

    def dislocation(self, kind: str = "edge", **kwargs):
        """Dispatch to edge_dislocation / screw_dislocation by kind."""
        k = str(kind).lower()
        if k == "edge":
            return self.edge_dislocation(**kwargs)
        if k == "screw":
            return self.screw_dislocation(**kwargs)
        raise ValueError("kind must be 'edge' or 'screw'.")

    # ------------------------------------------------------------------
    # Dislocations (myLAMMPS lmpData backend)
    #
    # myLAMMPS conventions the data must satisfy BEFORE the call:
    #   edge : y || Burgers vector, z || glide-plane normal, x || line
    #          (create_edge_dislocation removes a half-plane slab of width
    #           ~b/2 along y in the upper half of z, then shrinks y)
    #   screw: z || Burgers vector (= line direction for a screw)
    #
    # The methods below handle that orientation automatically: give the
    # crystallographic directions currently on x/y/z (current_axes) plus
    # the desired burgers/glide spec, and the needed swap_axes permutation
    # is derived. Pass newaxis=[...] to override, or nothing if the data
    # is already in the required frame.
    # ------------------------------------------------------------------

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

    @staticmethod
    def _find_parallel_axis(target, axes_unit: list[np.ndarray], what: str) -> int:
        """Index of the (anti)parallel axis among axes_unit, or raise."""
        t = BCCDefect._unit(target)
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

    @staticmethod
    def _default_burgerm(lattice_const: float, burgers) -> float:
        """|b| from the lattice constant for common bcc Burgers vectors."""
        hkl = sorted(abs(int(round(x))) for x in BCCDefect._parse_hkl(burgers))
        a = float(lattice_const)
        if hkl == [1, 1, 1]:
            return a * np.sqrt(3.0) / 2.0   # ½<111>
        if hkl == [0, 0, 1]:
            return a                         # <100>
        raise ValueError(
            f"No default |b| for burgers {burgers}; pass burgerm=... explicitly."
        )

    def _orient_for_dislocation(
        self,
        *,
        burgers,
        glide_plane,
        current_axes,
        newaxis,
        burgers_target: int,
    ) -> list[int]:
        """Resolve the swap_axes permutation for a dislocation geometry.

        burgers_target: which axis index the Burgers vector must end up on
        (1 = y for edge, 2 = z for screw). For edge, the glide-plane normal
        goes to z and the remaining axis to x. For screw (glide_plane=None),
        the remaining two axes keep their original relative order.
        Returns the permutation (identity = [0, 1, 2]).
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
            # edge geometry: x = line, y = burgers, z = glide-plane normal
            if burgers_target != 1:
                raise ValueError("glide_plane is only meaningful for edge (burgers_target=1).")
            ig = self._find_parallel_axis(self._parse_hkl(glide_plane), axes, "Glide-plane normal")
            if ig == ib:
                raise ValueError("Burgers vector and glide-plane normal map to the same axis.")
            iline = ({0, 1, 2} - {ib, ig}).pop()
            return [iline, ib, ig]

        # screw: only the Burgers axis is pinned; others keep original order
        others = [i for i in (0, 1, 2) if i != ib]
        perm = [0, 0, 0]
        perm[burgers_target] = ib
        rest = [i for i in (0, 1, 2) if i != burgers_target]
        for slot, src in zip(rest, others):
            perm[slot] = src
        return perm

    def edge_dislocation(
        self,
        *,
        ff_elements: list[str],
        lattice_const: float | None = None,
        burgers=(1, 1, 1),                 # Burgers DIRECTION (½<111> assumed for |b|)
        glide_plane=(1, 1, 0),             # glide-plane normal, e.g. (1,1,0) or (1,1,2)
        current_axes: list | None = None,  # directions now on x,y,z e.g. [[1,1,0],[1,1,2],[1,1,1]]
        newaxis: list[int] | None = None,  # explicit swap_axes override
        burgerm: float | None = None,      # explicit |b| override
        nedges: int = 1,
        fy_start: float | None = None,
        add_vacuum: bool = False,
        vacuum_direction: int = 2,
        lvac: float = 20.0,
        restore_axes: bool = False,
        reset_ids: bool = True,
    ) -> dict:
        """Insert edge dislocation(s) via lmpData.create_edge_dislocation.

        The data must already be a built supercell (read/scale/replicate
        beforehand). Orientation is handled here: either pass current_axes
        (+ burgers/glide_plane) for automatic axis permutation, pass
        newaxis=[...] directly, or pass neither if the data already has
        y || b, z || glide normal, x || line.

        Example (⟨111⟩{110} edge in a x=[110], y=[112], z=[111] supercell)::

            defect.edge_dislocation(
                ff_elements=["Fe"],
                lattice_const=2.8553,
                burgers=(1, 1, 1),
                glide_plane=(1, 1, 0),
                current_axes=[[1, 1, 0], [1, 1, 2], [1, 1, 1]],
            )   # derives newaxis=[1, 2, 0] like the manual workflow

        Returns
        -------
        dict
            burgerm, newaxis used, atom counts before/after, restored flag.
        """
        data = self.data

        if data.atoms is None or len(data.atoms) == 0:
            raise ValueError("No atoms loaded.")
        data.assert_force_field(ff_elements)

        if burgerm is None:
            if lattice_const is None:
                raise ValueError("Provide lattice_const (for auto |b|) or burgerm=...")
            burgerm = self._default_burgerm(lattice_const, burgers)
        burgerm = float(burgerm)

        perm = self._orient_for_dislocation(
            burgers=burgers, glide_plane=glide_plane,
            current_axes=current_axes, newaxis=newaxis,
            burgers_target=1,  # edge: b -> y
        )
        n_before = int(len(data.atoms))

        if perm != [0, 1, 2]:
            data.swap_axes(perm)

        data.create_edge_dislocation(
            burgerm,
            nedges=int(nedges),
            fy_start=fy_start,
            add_vacuum=add_vacuum,
            direction=int(vacuum_direction),
            lvac=float(lvac),
            reset_ids=reset_ids,
        )

        restored = False
        if restore_axes and perm != [0, 1, 2]:
            inv = np.argsort(perm).tolist()
            data.swap_axes(inv)
            restored = True

        return {
            "kind": "edge",
            "burgerm": burgerm,
            "nedges": int(nedges),
            "newaxis": perm,
            "n_atoms_before": n_before,
            "n_atoms_after": int(len(data.atoms)),
            "n_atoms_removed": n_before - int(len(data.atoms)),
            "axes_restored": restored,
        }

    def screw_dislocation(
        self,
        *,
        ff_elements: list[str],
        lattice_const: float | None = None,
        burgers=(1, 1, 1),
        current_axes: list | None = None,  # directions now on x,y,z
        newaxis: list[int] | None = None,
        burgerm: float | None = None,
        nscrews: int = 1,
        style: str = "bcc",
        handle_pbc: str = "tilt",
        orientation: bool = True,
        add_vacuum: bool = False,
        vacuum_direction: int = 0,
        lvac: float = 20.0,
        restore_axes: bool = False,
    ) -> dict:
        """Insert screw dislocation(s) via lmpData.create_screw_dislocation.

        myLAMMPS convention: z || Burgers vector (= line direction); for bcc
        the recommended frame is x=[1,-1,0], y=[1,1,-2], z=[1,1,1]/2.
        Orientation handling mirrors edge_dislocation (b is pinned to z; the
        other two axes keep their relative order).

        Note: an OVITO-based screw builder is planned as an alternative
        backend; this method wraps the myLAMMPS implementation.
        """
        data = self.data

        if data.atoms is None or len(data.atoms) == 0:
            raise ValueError("No atoms loaded.")
        data.assert_force_field(ff_elements)

        if burgerm is None:
            if lattice_const is None:
                raise ValueError("Provide lattice_const (for auto |b|) or burgerm=...")
            burgerm = self._default_burgerm(lattice_const, burgers)
        burgerm = float(burgerm)

        perm = self._orient_for_dislocation(
            burgers=burgers, glide_plane=None,
            current_axes=current_axes, newaxis=newaxis,
            burgers_target=2,  # screw: b -> z
        )
        n_before = int(len(data.atoms))

        if perm != [0, 1, 2]:
            data.swap_axes(perm)

        data.create_screw_dislocation(
            burgerm,
            int(nscrews),
            style=str(style),
            handle_pbc=str(handle_pbc),
            orientation=orientation,
            add_vacuum=add_vacuum,
            direction=int(vacuum_direction),
            lvac=float(lvac),
        )

        restored = False
        if restore_axes and perm != [0, 1, 2]:
            inv = np.argsort(perm).tolist()
            data.swap_axes(inv)
            restored = True

        return {
            "kind": "screw",
            "burgerm": burgerm,
            "nscrews": int(nscrews),
            "newaxis": perm,
            "n_atoms_before": n_before,
            "n_atoms_after": int(len(data.atoms)),
            "axes_restored": restored,
        }

    # ------------------------------------------------------------------
    # Prismatic dislocation loops (the only two public entry points)
    #
    # A bcc loop is a coherent platelet of k atomic layers along the habit
    # normal:  <111> -> b = 1/2<111>, |b| = a*sqrt(3)/2, k = 3 layers
    #          <100> -> b = a<100>,   |b| = a,           k = 2 layers
    #
    # Both methods take loop_type = "sil" (self-interstitial, the focus here)
    # or "vl" (vacancy, a SEPARATE code path). For oriented cells, pass the
    # matrix coordinate system via box_axes (crystal directions on box x/y/z)
    # so the requested plane is rotated into the box frame before selection;
    # or pass normal_box to give the box-frame normal directly.
    # ------------------------------------------------------------------

    def add_111_loop(
        self,
        *,
        loop_type: str = "sil",
        radius: float,
        lattice_const: float,
        ff_elements: list[str],
        atomic_masses: list[float],
        center,
        habit_plane=(1, 1, 1),
        box_axes=None,
        normal_box=None,
        shape: str = "round",
        n_repeats: int = 1,
        loop_atom_type: int | None = None,       # SIL: type for the new atoms
        closure_layers: int = 2,                 # VL only
        closure_disp_frac: float = 0.15,         # VL only
        select_types=None,
        wrap: bool = True,
        reinit: bool = True,
    ):
        """Add a <111> prismatic loop (3 (111) layers; b = 1/2<111>)."""
        return self._add_loop(
            family="111", loop_type=loop_type, radius=radius,
            lattice_const=lattice_const, ff_elements=ff_elements,
            atomic_masses=atomic_masses, center=center, habit_plane=habit_plane,
            box_axes=box_axes, normal_box=normal_box, shape=shape,
            n_repeats=n_repeats, loop_atom_type=loop_atom_type,
            closure_layers=closure_layers, closure_disp_frac=closure_disp_frac,
            select_types=select_types, wrap=wrap, reinit=reinit,
        )

    def add_100_loop(
        self,
        *,
        loop_type: str = "sil",
        radius: float,
        lattice_const: float,
        ff_elements: list[str],
        atomic_masses: list[float],
        center,
        habit_plane=(1, 0, 0),
        box_axes=None,
        normal_box=None,
        shape: str = "round",
        n_repeats: int = 1,
        loop_atom_type: int | None = None,
        closure_layers: int = 2,
        closure_disp_frac: float = 0.15,
        select_types=None,
        wrap: bool = True,
        reinit: bool = True,
    ):
        """Add a <100> prismatic loop (2 (200) layers; b = a<100>)."""
        return self._add_loop(
            family="100", loop_type=loop_type, radius=radius,
            lattice_const=lattice_const, ff_elements=ff_elements,
            atomic_masses=atomic_masses, center=center, habit_plane=habit_plane,
            box_axes=box_axes, normal_box=normal_box, shape=shape,
            n_repeats=n_repeats, loop_atom_type=loop_atom_type,
            closure_layers=closure_layers, closure_disp_frac=closure_disp_frac,
            select_types=select_types, wrap=wrap, reinit=reinit,
        )

    def _add_loop(
        self, *, family, loop_type, radius, lattice_const, ff_elements,
        atomic_masses, center, habit_plane, box_axes, normal_box, shape,
        n_repeats, loop_atom_type, closure_layers, closure_disp_frac,
        select_types, wrap, reinit,
    ):
        """Shared core: select the k-layer disc, then dispatch sil/vl."""
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
        )
        if ctx is None:
            return 0

        if is_sil:
            return self._apply_sia_insert(
                ctx, loop_atom_type=loop_atom_type, wrap=wrap, reinit=reinit,
            )
        # vacancy: separate path (delete + close the gap), unchanged
        return self._apply_vacancy_closure(
            ctx, burgers_vector=None, closure_layers=closure_layers,
            closure_disp_frac=closure_disp_frac, wrap_move_atoms=wrap, reinit=reinit,
        )

    @staticmethod
    def _family_loop_geometry(family, habit_plane, a: float):
        """(unit cubic normal, |b|, k layers, d_layer) for the loop family.

        k = bcc atomic layers spanned by one Burgers repeat |b| (= k * d_layer).
        """
        fam = str(family).replace("<", "").replace(">", "").strip()
        if fam == "111":
            _hkl, n = BCCDefect._parse_111(habit_plane)
            return n, a * np.sqrt(3.0) / 2.0, 3, a / (2.0 * np.sqrt(3.0))
        if fam == "100":
            _hkl, n = BCCDefect._parse_100(habit_plane)
            return n, a, 2, a / 2.0
        raise ValueError("family must be '111' or '100'.")

    @staticmethod
    def _resolve_loop_normal(n_cubic, box_axes, normal_box):
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

    def precipitate(self, *args, **kwargs):
        raise NotImplementedError("Precipitate generation is a placeholder for a future implementation.")

    def gas(self, *, mode: str = "interstitial", atom_type: int, sites=None, center=None, radius=None, **kwargs):
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

    def grain_boundary(self, *args, **kwargs):
        raise NotImplementedError("Grain-boundary/polycrystal generation is a placeholder for a future implementation.")

    def add_interstitial_bcc(self,*,fraction: float,atom_type: int,min_dist: float = 1.2,max_trials_per_atom: int = 20,
        rng: np.random.Generator | None = None,reinit: bool = True,):
        """Add random BCC tetrahedral/octahedral interstitial atoms.

        Candidate sites are random tetra/octa offsets around random host atoms.
        A candidate is accepted only if it is at least ``min_dist`` away from
        every existing atom and every previously placed interstitial
        (minimum-image distances). Existing-atom checks use a periodic KD-tree
        when the box is orthogonal; non-orthogonal boxes fall back to a
        vectorized minimum-image sweep.
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

        matT = data.box.matrix.T

        tetra = []
        base = [
            (0.25, 0.5, 0.0),
            (0.25, -0.5, 0.0),
            (-0.25, 0.5, 0.0),
            (-0.25, -0.5, 0.0),
        ]
        perms = [
            (0, 1, 2),
            (0, 2, 1),
            (1, 0, 2),
            (1, 2, 0),
            (2, 0, 1),
            (2, 1, 0),
        ]
        for p in perms:
            for v in base:
                tetra.append([v[p[0]], v[p[1]], v[p[2]]])
        tetra = np.unique(np.array(tetra), axis=0)

        octa = np.array([
            [0.5, 0.0, 0.0],
            [-0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.0, -0.5, 0.0],
            [0.0, 0.0, 0.5],
            [0.0, 0.0, -0.5],
        ])

        offsets_cart = np.vstack([tetra, octa]) @ matT

        # Fast neighbor lookup against EXISTING atoms:
        # orthogonal box -> periodic cKDTree; otherwise min-image sweep.
        # np.mod guards against wrap_cart returning exactly L (float rounding),
        # which the periodic cKDTree rejects.
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

    def create_void(self, center,radius: float,*,region_type: str = "sphere",reinit: bool = True,**kwargs,):
        """
        Delete atoms to create a void.

        Parameters
        ----------
        center : (3,) Cartesian center
        radius : float
            Used for sphere/cylinder or half-length for cube.
        region_type : str
            Must be one of:
                "sphere"
                "cube" / "box"
                "cylinder"
        kwargs :
            Extra keywords passed directly to region_mask(), e.g.
                plane="xy", height=20 for cylinder
                cube_params={...} for box
        """
        data = self.data

        if data.atoms is None or len(data.atoms) == 0:
            return 0

        if center is None:
            center = data.get_center("cart")
        else:
            center = np.asarray(center, dtype=float).reshape(3,)
        rt = str(region_type).lower()

        # ---------------------------------
        # Sphere  → select_sphere(center, radius)
        # ---------------------------------
        if rt == "sphere":
            n_deleted = data.delete_atoms("sphere",center=center,radius=float(radius),reinit=reinit)

        # ---------------------------------
        # Cube/Box → select_cube(cube_params)
        # ---------------------------------
        elif rt in ("cube", "box"):
            cube_params = {
                "center": center,
                "side_length": 2.0 * float(radius)
            }
            n_deleted = data.delete_atoms("cube",cube_params=cube_params,reinit=reinit)

        # ---------------------------------
        # Cylinder → select_cylinder(...)
        # requires plane + height
        # ---------------------------------
        elif rt == "cylinder":
            if "plane" not in kwargs:
                raise ValueError("cylinder requires plane='xy'/'xz'/'yz'")
            if "height" not in kwargs:
                raise ValueError("cylinder requires height=...")

            n_deleted = data.delete_atoms("cylinder",plane=kwargs["plane"],center=center,radius=float(radius),height=float(kwargs["height"]),reinit=reinit)

        else:
            raise ValueError(
                f"region_type must be one of: "
                f"'sphere', 'cube', 'box', 'cylinder'. "
                f"Got '{region_type}'."
            )

        return int(n_deleted)

    ########################################################################################################################
    @staticmethod
    def _unit(v):
        v = np.asarray(v, dtype=float).reshape(3,)
        n = np.linalg.norm(v)
        if n == 0:
            raise ValueError("Zero vector.")
        return v / n


    @staticmethod
    def _plane_basis_from_normal(n):
        """Orthonormal basis (u,v) spanning the plane with unit normal n."""
        n = BCCDefect._unit(n)
        a = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(a, n)) > 0.9:
            a = np.array([0.0, 1.0, 0.0])
        u = BCCDefect._unit(np.cross(n, a))
        v = BCCDefect._unit(np.cross(n, u))
        return u, v


    @staticmethod
    def _parse_111(habit_plane):
        if isinstance(habit_plane, str):
            s = habit_plane.replace(",", " ").replace("(", " ").replace(")", " ")
            parts = [p for p in s.split() if p]
            if len(parts) != 3:
                raise ValueError("habit_plane must have 3 integers, e.g. '1 1 1' or '-1 1 1'.")
            hkl = [int(p) for p in parts]
        else:
            arr = np.asarray(habit_plane, dtype=int).reshape(3,)
            hkl = [int(arr[0]), int(arr[1]), int(arr[2])]

        if sorted([abs(x) for x in hkl]) != [1, 1, 1]:
            raise ValueError(f"habit_plane must be a <111> normal, got {hkl}.")

        n = BCCDefect._unit(np.array(hkl, dtype=float))
        return tuple(hkl), n


    @staticmethod
    def _parse_100(habit_plane):
        """
        Parse and validate a <100> normal: (±1,0,0), (0,±1,0), (0,0,±1).
        Returns (hkl_tuple, unit_normal).
        """
        if isinstance(habit_plane, str):
            s = habit_plane.replace(",", " ").replace("(", " ").replace(")", " ")
            parts = [p for p in s.split() if p]
            if len(parts) != 3:
                raise ValueError("habit_plane must have 3 integers, e.g. '1 0 0' or '0 -1 0'.")
            hkl = [int(p) for p in parts]
        else:
            arr = np.asarray(habit_plane, dtype=int).reshape(3,)
            hkl = [int(arr[0]), int(arr[1]), int(arr[2])]

        absvals = [abs(x) for x in hkl]
        if sorted(absvals) != [0, 0, 1]:
            raise ValueError(f"habit_plane must be a <100> normal, got {hkl}.")

        n = BCCDefect._unit(np.array(hkl, dtype=float))
        return tuple(hkl), n

    # ------------------------------------------------------------------
    # Shared prismatic-loop machinery
    #
    # Both public builders (add_111_loop / add_100_loop) share this pipeline:
    # resolve the box-frame normal from the matrix coord system, build the
    # footprint, and snap to EXACTLY k*n_repeats atomic layers. Only the final
    # "apply" step differs:
    #   sil -> copy the slab and split it +/- b/2 along the normal (insert)
    #   vl  -> delete the slab + move each column toward the gap (separate)
    # ------------------------------------------------------------------

    def _select_loop_atoms(
        self,
        *,
        family: str,                  # "111" or "100"
        radius: float,
        lattice_const: float,
        ff_elements: list[str],
        atomic_masses: list[float],
        center,
        habit_plane,
        shape: str,
        n_repeats: int,
        select_types,
        box_axes=None,
        normal_box=None,
    ) -> dict | None:
        """Resolve the box-frame normal and select EXACTLY k*n_repeats layers.

        k and the layer spacing come from the family (see
        _family_loop_geometry): <111> -> k=3, d_layer = a/(2*sqrt(3));
        <100> -> k=2, d_layer = a/2. The selection snaps onto the actual atomic
        planes (circular-mean phase) and takes nlay = k*n_repeats consecutive
        layers centered on the loop. plane_tol = 0.30 * d_layer. PBC-safe via
        minimum-image relative vectors (triclinic-safe).

        Returns
        -------
        dict or None
            Selection context consumed by _apply_sia_insert /
            _apply_vacancy_closure, or None if no atoms were selected. Keys:
            sel, footprint, type_mask, rel_all, n, u, v, a, b_mag, k, d_layer,
            plane_tol.
        """
        data = self.data

        # -------------------------
        # 0) force-field bookkeeping
        # -------------------------
        if ff_elements is None or atomic_masses is None:
            raise ValueError("Must provide ff_elements and atomic_masses.")
        if len(ff_elements) != len(atomic_masses):
            raise ValueError("ff_elements and atomic_masses must have the same length.")
        data.assert_force_field(ff_elements, atomic_masses=atomic_masses)

        # -------------------------
        # 1) validate inputs
        # -------------------------
        if center is None:
            raise ValueError("center must be provided in Cartesian coordinates (shape (3,)).")
        c = np.asarray(center, dtype=float).reshape(3,)

        n_repeats = int(n_repeats)
        if n_repeats < 1:
            raise ValueError("n_repeats must be >= 1.")

        shape = str(shape).lower()
        if shape not in ("round", "square"):
            raise ValueError("shape must be 'round' or 'square'.")

        a = float(lattice_const)
        r = float(radius)

        # Know the matrix coordinate system first: get the family geometry, then
        # rotate the requested cubic normal into the box frame (box_axes /
        # normal_box) so the correct plane is selected for this cell.
        n_cubic, b_mag, k, d_layer = BCCDefect._family_loop_geometry(family, habit_plane, a)
        n = BCCDefect._resolve_loop_normal(n_cubic, box_axes, normal_box)
        nlay = k * n_repeats

        u, v = BCCDefect._plane_basis_from_normal(n)
        plane_tol = 0.30 * d_layer

        # -------------------------
        # 2) footprint cylinder + snap to EXACTLY nlay (= k*n_repeats) layers
        # -------------------------
        pos_all = data.atoms[["x", "y", "z"]].to_numpy(dtype=float)
        rel_all = data.min_image_rel(pos_all, c)

        type_mask = None
        if select_types is not None:
            type_mask = data.select_type(select_types)

        dist = rel_all @ n                     # signed distance along the normal
        proj = rel_all - np.outer(dist, n)     # in-plane projection (cylinder along n)
        cu = proj @ u
        cv = proj @ v
        if shape == "round":
            footprint = (cu * cu + cv * cv) <= (r * r)
        else:
            footprint = (np.abs(cu) <= r) & (np.abs(cv) <= r)

        cand = footprint.copy()
        if type_mask is not None:
            cand &= type_mask
        if not np.any(cand):
            return None

        # Anchor on the ACTUAL atomic planes (circular mean of dist mod d_layer),
        # snap the center plane to the nearest layer, then take exactly nlay
        # consecutive layers centered there. A plain width window is phase-
        # dependent and grabs k-1 or k layers; this reliably gives nlay.
        ang = 2.0 * np.pi * dist[cand] / d_layer
        phase = float(np.angle(np.mean(np.exp(1j * ang)))) / (2.0 * np.pi) * d_layer
        s_center = phase + round((0.0 - phase) / d_layer) * d_layer
        # integer layer steps centered on s_center: works for even and odd nlay.
        # (half-integer steps would land BETWEEN planes and select nothing.)
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
            "sel": sel,
            "footprint": footprint,
            "type_mask": type_mask,
            "rel_all": rel_all,
            "n": n,
            "u": u,
            "v": v,
            "a": a,
            "b_mag": b_mag,
            "k": k,
            "d_layer": d_layer,
            "plane_tol": plane_tol,
        }

    def _apply_sia_insert(
        self,
        ctx: dict,
        *,
        loop_atom_type: int | None,
        wrap: bool,
        reinit: bool,
    ) -> int:
        """SIL ending: copy the selected k-layer slab and split it +/- b/2.

        Each selected atom -> r - b/2 (original moves) and a NEW atom at
        r + b/2, where b = |b| * n is the full Burgers vector along the habit
        normal (a lattice translation), so the inserted layers continue the
        stacking. Returns the number of atoms inserted (= selected slab size).
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
                raise ValueError(f"loop_atom_type={loop_atom_type} not in Masses after assert_force_field.")
            new_types = np.full(len(sel_df), int(loop_atom_type), dtype=int)

        for t_id in np.unique(new_types):
            data.add_atoms(new_pos[new_types == t_id], atom_type=int(t_id), wrap=False, reinit=False)

        if reinit:
            data.initialization(normalization=False, style=1)

        return int(sel.size)

    def _apply_vacancy_closure(
        self,
        ctx: dict,
        *,
        burgers_vector,
        closure_layers: int,
        closure_disp_frac: float,
        wrap_move_atoms: bool,
        reinit: bool,
    ) -> int:
        """Vacancy ending: delete selected atoms + move closure layers inward.

        Direction convention:
          - signed distance is computed along b_hat (unit burgers direction)
          - atoms ABOVE the deleted slab (larger dist) move by -delta*b_hat
          - atoms BELOW the deleted slab (smaller dist) move by +delta*b_hat
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

        # Burgers direction for "above/below" + closure motion
        # default: use plane normal (prismatic loop b || n)
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

        # signed distance along burgers direction; deleted slab extent
        dist_b = rel_all @ b_hat
        s_bot_b = float(np.min(dist_b[sel]))
        s_top_b = float(np.max(dist_b[sel]))

        # move closure layers toward the void
        if n_close > 0 and delta != 0.0:
            N = len(data.atoms)
            del_mask = np.zeros(N, dtype=bool)
            del_mask[sel] = True

            move_above = np.zeros(N, dtype=bool)
            move_below = np.zeros(N, dtype=bool)

            # d_layer is also the correct step along b when b || n (prismatic loops)
            for j in range(1, n_close + 1):
                sA = s_top_b + j * d_layer   # above deleted region along +b
                sB = s_bot_b - j * d_layer   # below deleted region along -b

                onA = (np.abs(dist_b - sA) <= plane_tol)
                onB = (np.abs(dist_b - sB) <= plane_tol)

                # apply same footprint (same loop radius/shape in the habit plane)
                onA &= footprint
                onB &= footprint

                if type_mask is not None:
                    onA &= type_mask
                    onB &= type_mask

                move_above |= onA
                move_below |= onB

            # don't move atoms that will be deleted
            move_above &= ~del_mask
            move_below &= ~del_mask

            if np.any(move_above):
                idxA = data.atoms.index.to_numpy()[move_above]
                posA = data.atoms.loc[idxA, ["x", "y", "z"]].to_numpy(dtype=float)
                # ABOVE -> move toward void: -b_hat
                posA = posA - delta * b_hat.reshape(1, 3)
                if wrap_move_atoms:
                    posA = data.wrap_cart(posA)
                data.atoms.loc[idxA, ["x", "y", "z"]] = posA

            if np.any(move_below):
                idxB = data.atoms.index.to_numpy()[move_below]
                posB = data.atoms.loc[idxB, ["x", "y", "z"]].to_numpy(dtype=float)
                # BELOW -> move toward void: +b_hat
                posB = posB + delta * b_hat.reshape(1, 3)
                if wrap_move_atoms:
                    posB = data.wrap_cart(posB)
                data.atoms.loc[idxB, ["x", "y", "z"]] = posB

        # Now actually delete the loop atoms
        mask = np.zeros(len(data.atoms), dtype=bool)
        mask[sel] = True
        data.atoms = data.atoms.loc[~mask].copy()

        if reinit:
            data.initialization(normalization=False, style=1)

        return int(sel.size)


    ########################################################################################################################

    def add_FrenkelPair(
        self,
        n: int,
        *,
        min_sep_from_vac: float = 1.5,
        min_sep_from_atoms: float = 1.3,
        min_sep_between_new: float = 1.3,
        offset_range: tuple[float, float] = (10.0, 20.0),
        max_trials: int = 200,
        rng: Optional[np.random.Generator] = None,

        # region-based eligibility/exclusion (built on data.region_mask/select_*)
        eligible_region_type: str | None = None,
        eligible_region_kwargs: dict | None = None,
        exclude_region_type: str | None = None,
        exclude_region_kwargs: dict | None = None,

        # existing mask still supported
        atom_mask: Optional[np.ndarray] = None,

        wrap_new: bool = True,
    ) -> dict:
        """
        Create n Frenkel (vacancy-interstitial) pairs by moving randomly chosen atoms
        to nearby interstitial positions under distance constraints.

        Parameters
        ----------
        eligible_region_type / eligible_region_kwargs :
            Restrict which atoms are allowed to be selected (e.g., avoid surfaces).
        exclude_region_type / exclude_region_kwargs :
            Forbid selection in a region (e.g., frozen slab).
        atom_mask :
            Combined with region masks if provided.

        Region definitions use Modlmp_LmpData.region_mask(), so any of:
          "all", "sphere", "cube"/"box", "cylinder", "plane", "type", "atom_id", ...

        Notes
        -----
        Distance checks against the unmoved (fixed) atoms use a periodic
        KD-tree when the box is orthogonal; non-orthogonal boxes fall back
        to a vectorized minimum-image sweep.
        """
        data = self.data

        if data.atoms is None or len(data.atoms) == 0:
            return {
                "selected_indices": [],
                "moved_indices": [],
                "skipped_indices": [],
                "new_positions": {},
            }

        df = data.atoms
        if not all(c in df.columns for c in ("x", "y", "z")):
            raise ValueError("data.atoms must contain Cartesian columns: x, y, z.")

        all_idx = df.index.to_numpy(dtype=int)
        n_atoms = len(all_idx)

        # -------------------------
        # Build final eligibility mask (row-aligned)
        # -------------------------
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
                raise ValueError("eligible region_mask() must return a boolean mask aligned with data.atoms.")
            elig_mask &= np.asarray(m, dtype=bool)

        if exclude_region_type is not None:
            rk = {} if exclude_region_kwargs is None else dict(exclude_region_kwargs)
            m = data.region_mask(exclude_region_type, **rk)
            if m is None or len(m) != n_atoms:
                raise ValueError("exclude region_mask() must return a boolean mask aligned with data.atoms.")
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

        # fixed atoms (not moved)
        fixed_mask = np.ones(n_atoms, dtype=bool)
        fixed_mask[[index_to_row[int(i)] for i in chosen]] = False
        fixed_positions = pos[fixed_mask]

        # Fast neighbor lookup against FIXED atoms.
        # np.mod guards against wrap_cart returning exactly L (float rounding),
        # which the periodic cKDTree rejects.
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
            row = index_to_row[atom_id]
            r0 = pos[row].copy()

            ok = False
            for _ in range(int(max_trials)):
                v = rng.normal(size=3)
                nv = np.linalg.norm(v)
                if nv == 0:
                    continue
                v /= nv

                rr = rng.uniform(r_lo, r_hi)
                q = r0 + rr * v
                if wrap_new:
                    q = data.wrap_cart(q)

                d_vac = np.linalg.norm(data.min_image_rel(q, r0))
                if d_vac < min_sep_from_vac:
                    continue

                if _clashes_fixed(q):
                    continue

                if placed_positions:
                    d_new = data.min_image_dists(np.vstack(placed_positions), q)
                    if np.any(d_new < min_sep_between_new):
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
            "moved_indices": moved,
            "skipped_indices": skipped,
            "new_positions": placed_map,
            "eligible_count": int(len(elig)),
            "eligible_region_type": eligible_region_type,
            "exclude_region_type": exclude_region_type,
        }

        if skipped:
            print(f"[add_FrenkelPair] placed {len(moved)}/{len(chosen)}; skipped {len(skipped)} (max_trials reached).")
        else:
            print(f"[add_FrenkelPair] placed {len(moved)}/{len(chosen)} Frenkel pairs.")
        return report


    def transform(self, old_sys, new_sys, by="supercell", tol=1e-8):
        """
        Apply a coordinate transformation to the structure.

        Parameters
        ----------
        old_sys, new_sys : array-like, shape (3,3)
            3x3 matrices with rows = basis vectors.
        by : str
            - "transmat": apply SymmOp to atom coordinates (no replication)
            - "supercell": apply integer supercell matrix (replicates; returns NEW data)
        tol : float
            Tolerance for integer check of supercell matrix.
        """
        data = self.data

        old_sys = np.asarray(old_sys, float).reshape(3, 3)
        new_sys = np.asarray(new_sys, float).reshape(3, 3)

        if by == "transmat":
            symmop = find_symmop_lattices(Lattice(old_sys), Lattice(new_sys))
            data.atoms = data.modify_by_symmetry(data.atoms, symmop, normalization=True)
            return data

        elif by == "supercell":
            # compute integer supercell matrix S such that new_sys = S * old_sys
            S = new_sys @ np.linalg.inv(old_sys)
            S_round = np.rint(S)
            if not np.allclose(S, S_round, atol=tol):
                raise ValueError(f"new_sys is not an integer combination of old_sys.\nS=\n{S}")
            S_int = S_round.astype(int)

            # IMPORTANT: make_supercell returns a NEW object
            return data.make_supercell(S_int.tolist())

        else:
            raise ValueError("by must be 'transmat' or 'supercell'")


BccDefect = BCCDefect
