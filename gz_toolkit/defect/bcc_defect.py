"""BCC defect creation tools.

Everything structure-agnostic (point defects, voids, dislocation
orientation, the prismatic-loop pipeline, Frenkel pairs, transforms) lives
in ``base.LatticeDefect``. This module supplies only the bcc
crystallography.

Method/class tree
-----------------
bcc_defect.py
`-- BCCDefect(LatticeDefect)
    |-- point defects                            [base]
    |   |-- vacancy(coords)
    |   |-- interstitial(sites) / interstitial_cluster(sites)
    |   |-- random_interstitials(...) / add_interstitial_bcc(...)
    |   |-- dumbbell(center, direction)
    |   |-- dumbbell_111 / dumbbell_100 / dumbbell_110(center)
    |   `-- impurities(region_type) / impurity_atom(site)
    |-- larger defects / sinks
    |   |-- void(center, radius)                 [base]
    |   |-- dislocation(kind, ...)               [base] -> edge/screw
    |   |-- edge_dislocation(...)                [base, bcc |b| + defaults]
    |   |-- screw_dislocation(...)               [base, style="bcc"]
    |   |-- add_loop(family=...)                      [base] validates family
    |   |-- add_111_loop(loop_type="sil"|"vl", n_layers=3|1, ...)  # <111> loop
    |   |-- add_100_loop(loop_type="sil"|"vl", ...)   # <100> loop (2 layers)
    |   |-- precipitate(...) / grain_boundary(...)    # placeholders
    |   `-- gas(mode)
    `-- bcc crystallography (this module)
        |-- LOOP_FAMILIES                 <111> and <100>; no <110> (that is fcc)
        |-- _default_burgerm(...)         1/2<111> -> a*sqrt(3)/2; <100> -> a
        `-- _interstitial_offsets_frac()  tetrahedral + octahedral sites
"""

from __future__ import annotations

import warnings
from importlib.resources import files, as_file
from contextlib import ExitStack, redirect_stdout
import io
from pathlib import Path

import numpy as np

from gz_toolkit.core.modlmp import Modlmp_LmpData      # re-exported for callers
from gz_toolkit.defect.base import LatticeDefect


########################################################################################################################
class BCCDefect(LatticeDefect):
    """BCC defect builder using Modlmp_LmpData selection and edit methods."""

    lattice_name = "bcc"

    def add_c15_cluster(self, nsias: int, *, lattice_const: float, center=None,
                        ff_elements=None, atomic_masses=None, template_0=None,
                        template_1=None, min_distance_factor: float = 0.3) -> dict:
        """Insert an unrelaxed C15 seed in an axis-aligned, pure BCC host.

        Wraps myLAMMPS.create_C15_BCC transactionally: reject an incorrect
        excess atom count, changed cell, invalid types or overlaps before
        updating self.data. Templates are Z12 coordinate seeds, not a claim
        that a potential stabilizes C15. Relax and characterize the result.
        nsias is the NET excess, not the number of atoms in the C15 region.
        """
        from scipy.spatial import cKDTree
        if isinstance(nsias, bool) or not isinstance(nsias, int) or nsias < 2:
            raise ValueError("nsias must be an integer >= 2.")
        a = float(lattice_const)
        if not np.isfinite(a) or a <= 0 or not 0 < min_distance_factor < 1:
            raise ValueError("Positive lattice_const and 0 < min_distance_factor < 1 are required.")
        data = self.data
        matrix = np.asarray(data.box.matrix)
        if not np.allclose(matrix, np.diag(np.diag(matrix))) or data.atoms['type'].nunique() != 1:
            raise ValueError("C15 insertion currently requires an axis-aligned, pure BCC host.")
        repeats = np.diag(matrix) / a
        if not np.allclose(repeats, np.rint(repeats)) or len(data.atoms) != int(round(2 * np.prod(repeats))):
            raise ValueError('C15 insertion requires a complete periodic BCC supercell without pre-existing defects.')
        # Verify actual BCC registry (a rectangular rotated cell is insufficient).
        xyz = data.atoms[['x', 'y', 'z']].to_numpy(dtype=float)
        registry = np.mod((xyz - xyz[0]) / a, 1.0)
        registry[np.isclose(registry, 1, atol=1e-5)] = 0
        if not np.all(np.any(np.all(np.isclose(registry[:, None, :], np.array([[0, 0, 0], [.5, .5, .5]])[None], atol=1e-5), axis=2), axis=1)):
            raise ValueError("Host coordinates are not on the expected cubic BCC registry.")
        if ff_elements is None:
            raise ValueError("Provide ff_elements in force-field atom-type order.")
        before = len(data.atoms)
        with ExitStack() as stack:
            seeds = []
            for supplied, name in ((template_0, 'c15_z12_0.POSCAR'), (template_1, 'c15_z12_1.POSCAR')):
                path = Path(supplied) if supplied else stack.enter_context(as_file(files('gz_toolkit.defect').joinpath('templates', name)))
                if not path.is_file():
                    raise FileNotFoundError(path)
                seeds.append(str(path))
            with redirect_stdout(io.StringIO()):
                candidate = data.create_C15_BCC(a, nsias, to_center=(center if center is not None else [.5, .5, .5]),
                    is_cartesian=center is not None, fC15_0=seeds[0], fC15_1=seeds[1],
                    ff_elements=list(ff_elements), atomic_masses=atomic_masses, reset_ids=True)
        excess = len(candidate.atoms) - before
        if excess != nsias or not np.allclose(candidate.box.matrix, matrix):
            raise ValueError(f"C15 backend returned excess={excess} (requested {nsias}) or changed the cell; host unchanged.")
        positions = candidate.atoms[['x', 'y', 'z']].to_numpy(dtype=float)
        if not np.isfinite(positions).all():
            raise ValueError("C15 backend returned nonfinite coordinates; host unchanged.")
        lower = np.asarray(candidate.box.bounds)[:, 0]
        lengths = np.diag(matrix)
        shifted = np.mod(positions - lower, lengths)
        distances = cKDTree(shifted, boxsize=lengths).query(shifted, k=2)[0][:, 1]
        if np.min(distances) < min_distance_factor * a:
            raise ValueError("C15 backend produced overlapping atoms; enlarge the cell or inspect templates.")
        host_type = int(data.atoms['type'].iloc[0])
        candidate.atoms['type'] = host_type
        data.__dict__.update(candidate.__dict__)
        return {"kind": "c15", "n_atoms_before": before, "n_atoms_after": len(data.atoms),
                "net_interstitials": excess, "min_distance": float(np.min(distances)), "relaxed": False}

    #: b = 1/2<111> on {110} is the dominant bcc slip system
    default_burgers = (1, 1, 1)
    default_glide_plane = (1, 1, 0)
    #: myLAMMPS create_screw_dislocation style
    screw_style = "bcc"

    # ------------------------------------------------------------------
    # bcc crystallography
    # ------------------------------------------------------------------

    @staticmethod
    def _default_burgerm(lattice_const: float, burgers) -> float:
        """|b| from the lattice constant for common bcc Burgers vectors."""
        hkl = sorted(abs(int(round(x))) for x in BCCDefect._parse_hkl(burgers))
        a = float(lattice_const)
        if hkl == [1, 1, 1]:
            return a * np.sqrt(3.0) / 2.0   # 1/2<111>
        if hkl == [0, 0, 1]:
            return a                         # <100>
        raise ValueError(
            f"No default |b| for burgers {burgers}; pass burgerm=... explicitly."
        )

    #: The two prismatic-loop families observed in bcc metals. <111> loops
    #: (b = 1/2<111>) dominate at low dose; <100> loops (b = a<100>) are the
    #: sessile ones that survive to higher dose in Fe. There is no <110>
    #: family here -- that is an fcc habit.
    LOOP_FAMILIES = {
        "111": {
            "allowed_abs": [1, 1, 1],
            "label": "111",
            "default_habit": (1, 1, 1),
            "burgers": lambda a: a * np.sqrt(3.0) / 2.0,   # 1/2<111>
            "k": 3,
            "d_layer": lambda a: a / (2.0 * np.sqrt(3.0)),
            "description": "b = 1/2<111>, glissile, 3 (111) layers",
            # 1 layer = one extra/missing (111) plane: a faulted platelet with
            # a/6<111> displacement (not a lattice vector). See _partial_insert_offset.
            "allowed_layers": (3, 1),
        },
        "100": {
            "allowed_abs": [0, 0, 1],
            "label": "100",
            "default_habit": (1, 0, 0),
            "burgers": lambda a: a,                        # a<100>
            "k": 2,
            "d_layer": lambda a: a / 2.0,
            "description": "b = a<100>, sessile, 2 (200) layers",
        },
    }

    def _partial_insert_offset(self, family, n_cubic, a: float, n_layers: int) -> np.ndarray:
        """Cubic-frame shift that puts a single inserted (111) plane on a free site.

        bcc (111) planes stack ABCABC with spacing d = a/(2 sqrt3). A copy of
        layer C moved only along the normal would sit directly over C, ~d
        away from it. Instead the copy takes the in-plane pattern of the layer
        below (B), which is offset from both neighbours (C below, A above), so
        every new atom lands ~a*sqrt3/2 from its nearest neighbours before
        relaxation. The returned vector is a nearest-neighbour vector a/2<111>
        pointing one layer DOWN; the caller keeps only its in-plane part.
        """
        fam = str(family).replace("<", "").replace(">", "").strip()
        if fam != "111" or int(n_layers) != 1:
            return super()._partial_insert_offset(family, n_cubic, a, n_layers)
        h = np.sign(np.asarray(n_cubic, dtype=float).reshape(3,))
        v = -h.copy()
        v[2] = h[2]              # v . h = -1  ->  one (111) layer below
        return 0.5 * float(a) * v

    @staticmethod
    def _parse_111(habit_plane):
        return BCCDefect._parse_family(habit_plane, [1, 1, 1], "habit_plane", label="111")

    @staticmethod
    def _parse_100(habit_plane):
        return BCCDefect._parse_family(habit_plane, [0, 0, 1], "habit_plane", label="100")

    @staticmethod
    def _interstitial_offsets_frac() -> np.ndarray:
        """Tetrahedral + octahedral offsets from a bcc host atom, in units of a.

        Tetrahedral sites sit at 1/4<210>-type positions (12 per cell),
        octahedral at 1/2<100> face/edge centres (6 per cell).
        """
        base = [(0.25, 0.5, 0.0), (0.25, -0.5, 0.0), (-0.25, 0.5, 0.0), (-0.25, -0.5, 0.0)]
        perms = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
        tetra = np.unique(
            np.array([[v[p[0]], v[p[1]], v[p[2]]] for p in perms for v in base]), axis=0
        )
        octa = np.array([
            [0.5, 0.0, 0.0], [-0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0], [0.0, -0.5, 0.0],
            [0.0, 0.0, 0.5], [0.0, 0.0, -0.5],
        ])
        return np.vstack([tetra, octa])

    # ------------------------------------------------------------------
    # Legacy alias
    # ------------------------------------------------------------------

    def add_interstitial_bcc(self, **kwargs) -> int:
        """Deprecated alias for :meth:`random_interstitials`."""
        return self.random_interstitials(**kwargs)

    # ------------------------------------------------------------------
    # Prismatic dislocation loops (public entry points)
    #
    # A bcc loop is a coherent platelet of k atomic layers along the habit
    # normal. Both take loop_type = "sil" (self-interstitial) or "vl"
    # (vacancy, a separate code path). For oriented cells, pass box_axes
    # (crystal directions on box x/y/z) so the requested plane is rotated
    # into the box frame, or normal_box to give the box-frame normal.
    # ------------------------------------------------------------------

    def add_111_loop(self, *, loop_type: str = "sil", radius: float, lattice_const: float,
                     ff_elements: list[str], atomic_masses: list[float], center,
                     habit_plane=(1, 1, 1), box_axes=None, normal_box=None,
                     shape: str = "round", n_repeats: int = 1,
                     loop_atom_type: int | None = None, closure_layers: int = 2,
                     closure_disp_frac: float = 0.15, select_types=None,
                     wrap: bool = True, reinit: bool = True, n_layers: int = 3,
                     plane_tol_factor: float = 0.30):
        """Add a <111> prismatic loop on (111) layers.

        n_layers=3 (default): one full Burgers repeat, b = 1/2<111>, a perfect
        glissile loop; ``n_repeats`` stacks several such repeats.
        n_layers=1: a single extra ("sil") or missing ("vl") (111) plane. Its
        displacement a/6<111> is not a lattice vector, so the loop is faulted;
        it needs n_repeats=1, and for "sil" the crystal orientation (a cubic
        cell or ``box_axes``, not ``normal_box``) to place the new plane.
        """
        return self._add_loop(
            family="111", loop_type=loop_type, radius=radius,
            lattice_const=lattice_const, ff_elements=ff_elements,
            atomic_masses=atomic_masses, center=center, habit_plane=habit_plane,
            box_axes=box_axes, normal_box=normal_box, shape=shape,
            n_repeats=n_repeats, loop_atom_type=loop_atom_type,
            closure_layers=closure_layers, closure_disp_frac=closure_disp_frac,
            select_types=select_types, wrap=wrap, reinit=reinit, n_layers=n_layers,
            plane_tol_factor=plane_tol_factor,
        )

    def add_100_loop(self, *, loop_type: str = "sil", radius: float, lattice_const: float,
                     ff_elements: list[str], atomic_masses: list[float], center,
                     habit_plane=(1, 0, 0), box_axes=None, normal_box=None,
                     shape: str = "round", n_repeats: int = 1,
                     loop_atom_type: int | None = None, closure_layers: int = 2,
                     closure_disp_frac: float = 0.15, select_types=None,
                     wrap: bool = True, reinit: bool = True,
                     plane_tol_factor: float = 0.30):
        """Add a <100> prismatic loop (2 (200) layers; b = a<100>).

        With shape="square", radius may be a scalar (square, side 2*radius)
        or a (half_u, half_v) pair for a rectangle -- see
        BCCDefectBase._plane_basis_from_normal for how u, v map onto box x/y/z.
        """
        return self._add_loop(
            family="100", loop_type=loop_type, radius=radius,
            lattice_const=lattice_const, ff_elements=ff_elements,
            atomic_masses=atomic_masses, center=center, habit_plane=habit_plane,
            box_axes=box_axes, normal_box=normal_box, shape=shape,
            n_repeats=n_repeats, loop_atom_type=loop_atom_type,
            closure_layers=closure_layers, closure_disp_frac=closure_disp_frac,
            select_types=select_types, wrap=wrap, reinit=reinit,
            plane_tol_factor=plane_tol_factor,
        )


BccDefect = BCCDefect
